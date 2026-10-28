"""Build editable, source-derived cytoskeleton animations from JSON templates.

Blender 5.x; requires prepared assets from prepare_filament_assets.py. All
motion is authored teaching motion. No stochastic kinetics or MD is implied.
"""
import argparse
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import collection, duplicate, load_chain, material, save_and_render, studio


def translate(x=0, y=0, z=0):
    return Matrix.Translation((x, y, z))


def smooth(t):
    t = min(1., max(0., t))
    return t*t*t*(10+t*(-15+6*t))


class SceneBuilder:
    def __init__(self, args, config, profile):
        self.args, self.config, self.profile = args, config, profile
        self.motion = json.loads(args.motion_profile.read_text())['paths']
        self.groups = {k: collection(v) for k, v in {
            'units': '01 Source-derived complete subunits',
            'rigs': '02 Editable assembly rigs',
            'context': '03 Context — illustrative positions',
            'labels': '04 Explanatory labels'}.items()}
        self.prototypes, self.units, self.events, self.asset_records = {}, [], [], {}
        pal = profile['palette']
        self.mats = {k: material('Polymer '+k, color, profile['roughness'], profile['bump']) for k, color in {
            'a': pal['alpha'], 'b': pal['beta'], 'new': pal['protein'],
            'depart': pal['dna1'], 'junction': pal['accent'], 'context': pal['heads']}.items()}

    def info(self, kind):
        if kind not in self.asset_records:
            self.asset_records[kind] = json.loads((self.args.asset_dir/kind/'manifest.json').read_text())
        return self.asset_records[kind]

    def unit(self, kind, chain, name, transform=None, color='a', parent=None, context=False):
        key = (kind, chain)
        group = self.groups['context' if context else 'units']
        if key not in self.prototypes:
            obj, _ = load_chain(self.args.asset_dir/kind, chain, self.mats[color], group, name)
            self.prototypes[key] = obj
        else:
            obj = duplicate(self.prototypes[key], name, group)
            # An animated source object must not donate its action to another unit.
            obj.animation_data_clear()
        obj.name = name
        obj.parent = parent
        obj.matrix_basis = transform if transform is not None else Matrix.Identity(4)
        obj.material_slots[0].link = 'OBJECT'
        obj.material_slots[0].material = self.mats[color]
        info = self.info(kind)
        obj['polymer_unit_id'] = name
        obj['asset_kind'], obj['source_chain'], obj['source_pdb'] = kind, chain, info['source']
        obj['role'] = 'context' if context else 'polymer_unit'
        obj['geometry_evidence'] = 'experimental-coordinate-derived display surface'
        obj['placement_evidence'] = 'schematic unless retained source assembly'
        if not context:
            self.units.append(obj)
        return obj

    def rig(self, name):
        obj = bpy.data.objects.new(name, None)
        self.groups['rigs'].objects.link(obj)
        obj['motion_evidence'] = 'schematic — not measured trajectory or molecular dynamics'
        return obj

    @staticmethod
    def interpolate(start, end, arc=0.):
        a, qa, _ = start.decompose()
        b, qb, _ = end.decompose()
        def pose(u):
            position = a.lerp(b, u) + Vector((0, 0, arc*math.sin(math.pi*u)**2))
            return Matrix.LocRotScale(position, qa.slerp(qb, u), Vector((1,1,1)))
        return pose

    def perform(self, obj, segments):
        prepared = []
        for seconds, pose, action in segments:
            first = 1 + round(seconds[0]*self.args.fps)
            last = min(self.args.frames, 1+round(seconds[1]*self.args.fps))
            if last <= first:
                raise ValueError('Motion interval must contain more than one frame')
            prepared.append((first,last,pose))
            self.events.append({'object': obj.name, 'action': action, 'first_frame': first, 'last_frame': last,
                                'start_local_matrix': [list(r) for r in pose(0)],
                                'end_local_matrix': [list(r) for r in pose(1)]})
        obj.rotation_mode = 'QUATERNION'
        for frame in range(1, self.args.frames+1):
            first,last,pose = next((s for s in prepared if frame<=s[1]),prepared[-1])
            matrix = pose(smooth((frame-first)/(last-first)))
            obj.location, obj.rotation_quaternion, obj.scale = matrix.decompose()
            obj.keyframe_insert('location', frame=frame)
            obj.keyframe_insert('rotation_quaternion', frame=frame)
        if obj.animation_data:
            for fcurve in obj.animation_data.action.layers[0].strips[0].channelbags[0].fcurves:
                for key in fcurve.keyframe_points:
                    key.interpolation = 'LINEAR'
    def move(self, obj, start, end, seconds, arc=0., action='approach'):
        self.perform(obj, [(seconds,self.interpolate(start,end,arc),action)])

    def clearance_pose(self, key, t):
        guide = self.motion[key]
        center, offset = Vector(guide['center']), Vector(guide['far_offset'])
        vector = Vector(guide['far_rotation_vector'])
        angle = vector.length * t**guide.get('rotation_exponent',1)
        rotation = Quaternion(vector.normalized(),angle).to_matrix().to_4x4()
        return translate(*(offset*t**guide.get('translation_exponent',1))) @ translate(*center) @ rotation @ translate(*(-center))

    def incoming(self, obj, target, seconds, index=0, distance=18, key='barbed', pool_row=0):
        near = target @ self.clearance_pose(key,1)
        start = translate(distance+pool_row*4, -18-index*11, 14+pool_row*12) @ target
        start = start @ Matrix.Rotation(math.radians(22), 4, 'Z')
        approach_fraction = .64 if self.config['asset_family']=='tubulin' else .58
        middle = seconds[0]+approach_fraction*(seconds[1]-seconds[0])
        self.perform(obj, [((seconds[0],middle),self.interpolate(start,near,2.5),'approach'),
                           ((middle,seconds[1]),lambda u: target@self.clearance_pose(key,1-u),'addition')])

    def actin(self):
        conf = self.config
        operation = conf['operation']
        step = np.array(self.info('actin')['step_matrix_nm'])
        count, cycles = conf['seed_subunits'], conf['cycles']
        offset = translate(-.5*(count-1)*self.info('actin')['rise_nm'])
        matrices = [offset @ Matrix(np.linalg.matrix_power(step, i)) for i in range(count+cycles)]
        rigs = None
        if operation in ('sever', 'anneal'):
            rigs = [self.rig('Left complete fragment'), self.rig('Right complete fragment')]
        units = []
        for i in range(count):
            parent = rigs[int(i >= count//2)] if rigs else None
            color = 'depart' if operation in ('depolymerize', 'treadmill') and i < cycles else ('a' if i%2 else 'b')
            units.append(self.unit('actin', 'D', f'Actin seed {i:03}', matrices[i], color, parent))
        if operation in ('polymerize', 'treadmill'):
            for j in range(cycles):
                obj = self.unit('actin', 'D', f'Actin incoming {j:03}', color='new')
                times = (1+j*3, 3.6+j*3)
                self.incoming(obj, matrices[count+j], times, j)
        if operation in ('depolymerize', 'treadmill'):
            for j in range(cycles):
                # Loss from the pointed side; internal units remain at fixed world positions.
                bound = matrices[j]
                near = bound @ self.clearance_pose('pointed',1)
                target = near.copy()
                # Continue out along this subunit's own cleared radial path;
                # a common exit lane makes successive released units cross.
                target.translation = bound.translation + (near.translation-bound.translation)*3.5
                self.perform(units[j], [((1+j*3,2+j*3),lambda u,m=bound: m@self.clearance_pose('pointed',u),'dissociation'),
                                        ((2+j*3,3.6+j*3),self.interpolate(near,target,-2),'departure')])
        if operation in ('sever', 'anneal'):
            anchor = matrices[count//2]
            near = self.clearance_pose('sever',1)
            far = near.copy()
            far.translation += 2*Vector(self.motion['sever']['far_offset'])
            local = self.interpolate(near,far)
            def world(matrix):
                return anchor@matrix@anchor.inverted()
            if operation=='sever':
                segments = [((2,5),lambda u: world(self.clearance_pose('sever',u)),'interface_separation'),
                            ((5,8),lambda u: world(local(u)),'fragment_departure')]
            else:
                segments = [((2,5),lambda u: world(local(1-u)),'fragment_approach'),
                            ((5,8),lambda u: world(self.clearance_pose('sever',1-u)),'interface_annealing')]
            self.perform(rigs[1],segments)
        if operation == 'bundle':
            for layer in (-1, 1):
                rig = self.rig(f'Complete filament {layer:+}')
                for i in range(count):
                    self.unit('actin', 'D', f'Bundle {layer:+} actin {i:03}', matrices[i], 'b' if layer<0 else 'new', rig)
                end = translate(0, 0, layer*9)
                start = translate(layer*8, 0, layer*14) @ Matrix.Rotation(math.radians(layer*10), 4, 'Y')
                self.move(rig, start, end, (1, 9), 0, 'parallel_bundle_alignment')
        return {'position': (2, -95, 33), 'target': (2, 0, 1), 'scale': 100 if operation=='bundle' else 84, 'light_scale': 2.0}

    def branch(self):
        info = self.info('branch')
        for chain in info['chains']:
            color = 'junction' if chain in info['arp_chains'] else ('a' if chain in info['mother_chains'] else 'b')
            self.unit('branch', chain, f'Branch source chain {chain}', color=color)
        step = np.array(info['step_matrix_nm'])
        for j in range(self.config['cycles']):
            target = Matrix(np.linalg.matrix_power(step, j+1))
            obj = self.unit('branch', 'j', f'Daughter new actin {j:03}', color='new')
            # Approach along the daughter direction, away from the mother filament.
            start = translate(7+j*3, -4, 17+j*4) @ target
            self.move(obj, start, target, (1+j*3, 3.6+j*3), 2, 'daughter_extension')
        return {'position': (4, -90, 23), 'target': (1, 0, 5), 'scale': 68, 'light_scale': 1.6}

    def cap(self):
        info = self.info('cap')
        for chain in info['actin_chains']:
            self.unit('cap', chain, f'Capped-reference actin {chain}', color='a' if ord(chain)%2 else 'b')
        rig = self.rig('CapZ heterodimer — both complete subunits')
        for chain in info['cap_chains']:
            self.unit('cap', chain, f'CapZ chain {chain}', color='new', parent=rig)
        self.incoming(rig, Matrix.Identity(4), (1,8), key='cap')
        return {'position': (6, -75, 24), 'target': (6, 0, 1), 'scale': 62, 'light_scale': 1.5}

    def tubulin(self):
        info = self.info('tubulin')
        repeat = Vector(info['repeat_nm'])
        count, cycles = self.config['seed_subunits'], self.config['cycles']
        operation = self.config['operation']
        rigs = {}

        def dimer(row, pf, label, color):
            rig = self.rig(label)
            rig.matrix_basis = translate(*((row-(count-1)/2)*repeat))
            for ci, chain in enumerate(info['dimers'][pf]):
                self.unit('tubulin', chain, label+' '+('alpha' if ci==0 else 'beta'), color=color if color else ('a' if ci==0 else 'b'), parent=rig)
            return rig

        for row in range(count):
            for pf in range(3):
                rigs[row, pf] = dimer(row, pf, f'Seed dimer row {row} pf {pf}', None)
        if operation in ('catastrophe', 'rescue'):
            for j in range(cycles):
                # Release the outer dimers before the central dimer, leaving
                # its local exit open. This is an authored teaching order.
                for rank, pf in enumerate((2, 1, 0)):
                    rig = rigs[count-1-j, pf]
                    start = rig.matrix_basis.copy()
                    near = start@self.clearance_pose(f'tub{pf}',1)
                    end = start.copy()
                    end.translation = Vector((32+j*14,24+pf*13,(pf-1)*14))
                    first = 1+j*4.8+rank*1.4
                    times = (first, first+3.6)
                    middle = first+1
                    self.perform(rig,[((times[0],middle),lambda u,m=start,k=f'tub{pf}': m@self.clearance_pose(k,u),'tip_dimer_release'),
                                      ((middle,times[1]),self.interpolate(near,end,1),'departure')])
        if operation in ('grow', 'rescue'):
            base = count if operation == 'grow' else count-cycles
            for j in range(cycles):
                for pf in range(3):
                    rig = dimer(base+j, pf, f'New dimer {j} pf {pf}', 'new')
                    target = rig.matrix_basis.copy()
                    restart = 1 if operation=='grow' else 2+cycles*4.8
                    first = restart+j*4.8+pf*1.4
                    times = (first, first+2.6)
                    self.incoming(rig, target, times, pf, 22, f'tub{pf}', j)
        return {'position': (5, -105, 29), 'target': (5, 0, 0), 'scale': 106, 'light_scale': 2.5}

    def context(self):
        # Repeat the same source and conformation: no unrelated species hidden in decorative crowding.
        if self.config['asset_family'] != 'actin' or not self.config.get('environment', True):
            return
        light_background = sum(int(self.profile['background'].lstrip('#')[i:i+2],16) for i in (0,2,4)) >= 400
        positions = [(-26,18,5),(28,20,3),(-31,28,-29),(31,28,-29)] if light_background else [(-26,18,5),(28,20,3),(-24,26,-22),(27,25,-21)]
        for index, position in enumerate(positions):
            transform = translate(*position) @ Matrix.Rotation(index*.9, 4, 'Y')
            obj = self.unit('actin', 'D', f'Context intact actin {index}', transform, 'context', context=True)
            obj['placement_evidence'] = 'authored background; count, location and free-monomer conformation are not measured'

    def labels(self, camera):
        args, config = self.args, self.config
        font = bpy.data.fonts.load(str(args.font)) if args.font else None
        bg = self.profile['background'].lstrip('#')
        dark = sum(int(bg[i:i+2], 16) for i in (0,2,4)) < 400
        mat = material('Legible annotations', 'DCE9ED' if dark else '223741', 1, 0)
        shader = mat.node_tree.nodes.get('Principled BSDF')
        shader.inputs['Emission Color'].default_value = mat.diffuse_color
        shader.inputs['Emission Strength'].default_value = .8
        if not dark:
            # Unlit dark type stays readable under the much brighter teaching lights.
            mat.node_tree.nodes.clear()
            emission = mat.node_tree.nodes.new('ShaderNodeEmission')
            emission.inputs['Color'].default_value = mat.diffuse_color
            output = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
            mat.node_tree.links.new(emission.outputs[0],output.inputs['Surface'])
        width = camera.data.ortho_scale
        height = width*args.height/args.width
        entries = [(config['title'], -.44*width, .405*height, .026*width),
                   (config['subtitle'], -.44*width, .34*height, .017*width),
                   (config['footer'], -.44*width, (-.43 if dark else -.48)*height, .014*width)]
        if config['asset_family'] in ('actin', 'tubulin'):
            left, right = ('尖端（−）', '倒钩端（+）') if config['asset_family']=='actin' else ('负端方向', '正端方向')
            label_y = -.35 if config['operation']=='bundle' else -.23
            entries += [(left, -.30*width, label_y*height, .017*width),
                        (right, .20*width, label_y*height, .017*width)]
        for i, (text, x, y, size) in enumerate(entries):
            data = bpy.data.curves.new(f'Annotation {i}', 'FONT')
            data.body, data.size = text, size
            if not dark:
                data.offset = size*.008
            if font:
                data.font = font
            data.materials.append(mat)
            obj = bpy.data.objects.new(f'Annotation {i}', data)
            self.groups['labels'].objects.link(obj)
            obj.parent = camera
            obj.location = (x, y, -55)
            obj['role'] = 'annotation'
        if font:
            font.pack()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--template', type=Path, required=True)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--profile', type=Path, default=Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json')
    p.add_argument('--font', type=Path, help='Font supporting the template label language; packed in project, not the skill')
    p.add_argument('--motion-profile', type=Path, default=Path(__file__).resolve().parents[1]/'assets/polymer-motion.json')
    p.add_argument('--engine', choices=['cycles','eevee'], default='cycles')
    p.add_argument('--compositor', choices=['clean','reference','monochrome'])
    p.add_argument('--width', type=int, default=1920)
    p.add_argument('--samples', type=int, default=48)
    p.add_argument('--frame', type=int)
    p.add_argument('--render', action='store_true')
    p.add_argument('--animate', action='store_true')
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists():
        p.error('Use a new scene directory')
    config = json.loads(args.template.read_text())
    if config.get('template_type') != 'executable_blender':
        p.error('This is a storyboard, not an executable polymer template; implement its documented project workflow')
    allowed = {'actin': {'polymerize','depolymerize','treadmill','sever','anneal','bundle'},
               'branch': {'branch_growth'}, 'cap': {'cap'}, 'tubulin': {'grow','catastrophe','rescue'}}
    family = config['asset_family']
    if family not in allowed or config['operation'] not in allowed[family]:
        p.error('This is not an executable cytoskeleton template; use its documented implementation workflow')
    if not 1 <= config['cycles'] <= 3 or not 6 <= config['seed_subunits'] <= 18 and family=='actin':
        p.error('This composition supports 1–3 cycles and 6–18 seed actin subunits')
    args.fps = config['fps']
    args.frames = round(config['duration_seconds']*args.fps)
    if args.frames < args.fps*11:
        p.error('Template event layout needs at least 11 seconds')
    if family=='tubulin':
        if config['operation']=='rescue':
            required = 3.6+9.6*config['cycles']
        else:
            required = (3.6 if config['operation']=='catastrophe' else 2)+4.8*config['cycles']
        if config['duration_seconds'] < required:
            p.error(f'This microtubule timeline needs duration_seconds >= {required:.1f}')
    args.frame = args.frame or args.frames
    args.height = round(args.width*5/8)
    if not 1 <= args.frame <= args.frames:
        p.error('Preview frame outside duration')
    profile = json.loads(args.profile.read_text())
    bpy.ops.wm.read_factory_settings(use_empty=True)
    builder = SceneBuilder(args, config, profile)
    settings = getattr(builder, {'actin':'actin','branch':'branch','cap':'cap','tubulin':'tubulin'}[family])()
    builder.context()
    camera = studio(profile, args, settings)
    builder.labels(camera)
    bpy.context.scene.frame_set(args.frames)
    bpy.context.view_layer.update()
    manifest = {'template': config, 'asset_records': builder.asset_records,
                'coordinate_units': '1 Blender unit = 1 nm', 'motion_evidence': 'schematic',
                'unit_count': len(builder.units), 'unit_ids': [o.name for o in builder.units],
                'events': builder.events,
                'endpoints_world': {o.name: [list(r) for r in o.matrix_world] for o in builder.units},
                'limitations': config['limitations'],
                'runtime': 'Linked molecular meshes; baked editable object/rig keys; no external drivers or handlers'}
    save_and_render(args, manifest)


if __name__ == '__main__':
    main()
