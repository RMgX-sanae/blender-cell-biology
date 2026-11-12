"""Display complete selected source-state molecules in separate static panels."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import collection, material, mesh, save_and_render, studio
from audit_source_states import assets


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font-zh', type=Path, required=True)
    p.add_argument('--font-en', type=Path, required=True)
    p.add_argument('--profile', type=Path, default=Path(__file__).resolve().parents[1]/'assets/style-structure-materials.json')
    p.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    p.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'], default='clean')
    p.add_argument('--width', type=int, default=1800)
    p.add_argument('--samples', type=int, default=32)
    p.add_argument('--render', action='store_true')
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():
        p.error('Choose a fresh output directory')
    a.frames, a.fps, a.frame, a.animate = 1, 24, 1, False
    manifest_path = a.asset_dir/'manifest.json'
    m = assets(a.asset_dir)
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile = json.loads(a.profile.read_text())
    models, labels = collection('01 Separate source states'), collection('02 Readable labels')
    fonts = [bpy.data.fonts.load(str(f)) for f in (a.font_zh, a.font_en)]
    for font in fonts: font.pack()
    ink = material('High contrast labels', 'EEF2EB', .7, 0)
    shader = ink.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Emission Color'].default_value = ink.diffuse_color
    shader.inputs['Emission Strength'].default_value = .8
    def text(body, x, z, size, color=ink):
        data = bpy.data.curves.new(body, 'FONT')
        data.body, data.font, data.font_bold = body, fonts[0], fonts[1]
        # Use the alternate font slot for Times Regular, not a bold treatment.
        for char, form in zip(body, data.body_format): form.use_bold = ord(char) < 128 or char in '·—'
        data.size, data.align_x = size, 'CENTER'
        obj = bpy.data.objects.new(body, data); labels.objects.link(obj)
        obj.location, obj.rotation_euler = (x, -10, z), (np.pi/2, 0, 0)
        data.materials.append(color)
    element_colors = {'C': 'B4C5CA', 'N': '7CA8F3', 'O': 'F48A86', 'P': 'E9B45F',
                      'S': 'E7D36C', 'Fe': 'D4895A', 'Mg': '95D6C2', 'Ca': 'B6D7A8'}
    spheres, element_materials = {}, {}
    for element, color in element_colors.items():
        mat = element_materials[element] = material('Element '+element, color, .45, 0,
            settings=profile.get('surface_roles', {}).get('cofactor'))
        mat['surface_role'] = 'cofactor'
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, radius=1)
        template = bpy.context.object
        template.data.materials.append(mat)
        spheres[element] = template.data
        bpy.data.objects.remove(template, do_unlink=True)
    bounds = np.concatenate([np.load(a.asset_dir/actor['file'], allow_pickle=False)['vertices']
                             for state in m['states'] for actor in state['actors']])
    low, high = bounds.min(0), bounds.max(0)
    center = (low+high)/2
    panel_width = float(high[0]-low[0])+2.5
    body_half_height = float(high[2]-low[2])/2
    width = len(m['states'])*panel_width
    height = 2*body_half_height+8.5
    placements = []
    for i, state in enumerate(m['states']):
        x = (i-(len(m['states'])-1)/2)*panel_width
        rig = bpy.data.objects.new(state['pdb']+' source panel', None)
        models.objects.link(rig); rig.location = (x-center[0], 0, -center[2])
        rig['source_state_panel'] = state['pdb']; placements.append(list(rig.location))
        legend = []
        for actor in state['actors']:
            data = np.load(a.asset_dir/actor['file'], allow_pickle=False)
            role = 'nucleic_acid' if actor.get('polymer_type', '').startswith(('polyribonucleotide', 'polydeoxyribonucleotide')) else 'protein'
            mat = material(state['pdb']+' '+actor['label'], actor['color'], .48, 0,
                           settings=profile.get('surface_roles', {}).get(role))
            mat['surface_role'] = role
            obj = mesh(state['pdb']+' '+actor['label'], data['vertices'].tolist(), data['faces'].tolist(), [mat], models)
            obj.parent = rig
            obj['source_state'], obj['source_file'] = state['pdb'], actor['file']
            obj['source_author_chain'] = actor['chain']
            for key in ('instance_id', 'label_chain', 'operator'):
                if key in actor: obj['source_'+key] = actor[key]
            # No smoothing modifier, deformation, section or per-protein recentering.
            label_mat = material(state['pdb']+' label '+actor['label'], actor['color'], .7, 0)
            shader = label_mat.node_tree.nodes.get('Principled BSDF')
            shader.inputs['Emission Color'].default_value = label_mat.diffuse_color
            shader.inputs['Emission Strength'].default_value = 1
            legend.append((actor['label'], label_mat))
        for j, ligand in enumerate(state['ligands']):
            for k, (point, element) in enumerate(zip(ligand['positions_nm'], ligand['elements'])):
                if element not in spheres: raise ValueError('Add a reviewed element display: '+element)
                obj = bpy.data.objects.new(f"{state['pdb']} {ligand['component']} {j}:{k}", spheres[element])
                models.objects.link(obj); obj.parent = rig
                obj.location = point; obj.scale = (.13, .13, .13)
                obj['source_ligand_state'], obj['ligand_index'], obj['atom_index'] = state['pdb'], j, k
                for key in ('instance_id', 'label_chain', 'operator', 'owner_instance'):
                    if key in ligand: obj['source_'+key] = ligand[key]
        text(state['title'], x, body_half_height+2.0, .9)
        text(state['state_label'], x, body_half_height+.95, .65)
        text(state['pdb'], x, -body_half_height-2.4, .75)
        for j, (label, mat) in enumerate(legend):
            text(label, x+(j-(len(legend)-1)/2)*3.7, -body_half_height-1.25, .62, mat)
    text('独立结构状态对照', 0, body_half_height+3.6, .95)
    studio(profile, a, {'position': (0, -62, 0), 'target': (0, 0, 0),
                       'scale': width+2, 'light_scale': 1.2})
    scene = bpy.context.scene
    scene.render.resolution_y = round(a.width*height/(width+2))
    scene['source_states_manifest_sha256'] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    scene['source_state_panel_offsets'] = json.dumps(placements)
    save_and_render(a, {'source_states': m, 'panel_offsets_nm': placements,
        'purpose': 'Source-state qualification, not animation; local complete-chain selections, no ribosomal subunit dissociation.',
        'no_molecular_cut': True, 'trajectory_modeled': False,
        'background_scope': 'Plain technical comparison; not a finished cellular environment.'})


if __name__ == '__main__': main()
