"""Structure-derived tRNA substrate cycle with intact molecules and authored motion.

9SPI provides the ribosome, tRNAs and seven modeled mRNA residues. The motion,
short peptide and layered presentation are teaching constructions, not an MD
trajectory or a set of experimentally observed intermediate structures.
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
from detail_scene import (Parts, collection, linear, load_chain, material,
                          path, save_and_render, studio)
from build_molecular_scene import nucleic_atoms


def smooth(t, start, end):
    x = min(1., max(0., (t-start)/(end-start)))
    return x*x*(3-2*x)


def entry_progress(t):
    x = min(1., max(0., (t-.10)/(.34-.10)))
    return x*x*x*(x*(6*x-15)+10)


def entry_pose(data, progress, right, up, front):
    # Approach from the open upper-left side; the former right-hand path
    # crossed the mRNA. Keep the entire tRNA and its amino acid rigid.
    controls = [right*(-4)+up*2.5+front*2,
                right*(-3.5)+up*4+front*2,
                right*(-1.5)+up+front*.25]
    p = progress
    offset = ((1-p)**3*controls[0] + 3*(1-p)**2*p*controls[1] +
              3*(1-p)*p*p*controls[2])
    turn = min(1., p/.85)
    turn = turn*turn*turn*(turn*(6*turn-15)+10)
    angle = math.radians(-16)*(1-turn)
    rotation = np.array(Quaternion(Vector(front), angle).to_matrix())
    center = data['atoms'].mean(0)
    return rotation, center+offset-rotation @ center


def atom(data, residue, name):
    ids = np.flatnonzero((data['residue_numbers'] == residue) &
                        (data['atom_names'] == name))
    if len(ids) != 1:
        raise ValueError(f'Expected one {residue}:{name}, got {len(ids)}')
    return np.asarray(data['atoms'][ids[0]], dtype=float)


def rigid_fit(source, target):
    """Fit shared C4' sites, retaining every atom of the moving tRNA."""
    def sites(d):
        return {int(r): p for r, n, p in
                zip(d['residue_numbers'], d['atom_names'], d['atoms']) if n == "C4'"}
    a, b = sites(source), sites(target)
    common = sorted(a.keys() & b.keys())
    x, y = np.array([a[k] for k in common]), np.array([b[k] for k in common])
    cx, cy = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-cx).T @ (y-cy))
    if np.linalg.det(vt.T @ u.T) < 0:
        vt[-1] *= -1
    rot = vt.T @ u.T
    translation = cy-rot @ cx
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rot.T+translation-y)**2, axis=1))))
    return rot, translation, rmsd


def move(obj, rot, translation, frame):
    obj.rotation_mode = 'QUATERNION'
    obj.rotation_quaternion = Matrix(rot).to_quaternion()
    obj.location = translation
    obj.keyframe_insert(data_path='rotation_quaternion', frame=frame)
    obj.keyframe_insert(data_path='location', frame=frame)


def trajectory(data, fit, progress, front, bow=.7):
    rotation, translation, _ = fit
    q = Quaternion((1, 0, 0, 0)).slerp(Matrix(rotation).to_quaternion(), progress)
    r = np.array(q.to_matrix())
    center = data['atoms'].mean(0)
    endpoint = rotation @ center + translation
    c = (1-progress)*center + progress*endpoint + front*bow*math.sin(math.pi*progress)
    return r, c-r @ center


def update_line(obj, points, frame):
    for p, xyz in zip(obj.data.splines[0].points, points):
        p.co = (*xyz, 1)
        p.keyframe_insert(data_path='co', frame=frame)


def emission(name, color):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    n = m.node_tree.nodes
    n.clear()
    shader = n.new('ShaderNodeEmission')
    shader.inputs['Color'].default_value = (*linear(color), 1)
    shader.inputs['Strength'].default_value = 1
    out = n.new('ShaderNodeOutputMaterial')
    m.node_tree.links.new(shader.outputs[0], out.inputs['Surface'])
    return m


def text_object(body, name, size, mat, group, font=None):
    data = bpy.data.curves.new(name, 'FONT')
    data.body, data.size = body, size
    if font:
        data.font = font
    data.materials.append(mat)
    obj = bpy.data.objects.new(name, data)
    obj.visible_shadow = False
    group.objects.link(obj)
    return obj


def mrna_frames(data):
    frames = {}
    for r in np.unique(data['residue_numbers']):
        r = int(r)
        c = atom(data, r, "C4'")
        x = atom(data, r, "C1'")-c
        x /= np.linalg.norm(x)
        y = atom(data, r, "C3'")-c
        y -= x*np.dot(x, y)
        y /= np.linalg.norm(y)
        frames[r] = (c, np.column_stack((x, y, np.cross(x, y))))
    return frames


def track_pose(frames, u):
    lo, hi = min(frames), max(frames)
    if u <= lo:
        return frames[lo][0]+(u-lo)*(frames[lo+1][0]-frames[lo][0]), frames[lo][1]
    if u >= hi:
        return frames[hi][0]+(u-hi)*(frames[hi][0]-frames[hi-1][0]), frames[hi][1]
    i = math.floor(u)
    x = u-i
    def point(k):
        if k < lo: return frames[lo][0]+(k-lo)*(frames[lo+1][0]-frames[lo][0])
        if k > hi: return frames[hi][0]+(k-hi)*(frames[hi][0]-frames[hi-1][0])
        return frames[k][0]
    p0, p1, p2, p3 = [point(k) for k in (i-1, i, i+1, i+2)]
    c = .5*((2*p1)+(-p0+p2)*x+(2*p0-5*p1+4*p2-p3)*x*x+
             (-p0+3*p1-3*p2+p3)*x*x*x)
    q = Matrix(frames[i][1]).to_quaternion().slerp(
        Matrix(frames[i+1][1]).to_quaternion(), x)
    return c, np.array(q.to_matrix())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--style-profile', type=Path,
                   default=Path(__file__).resolve().parent.parent/'assets/style-molecular-reference.json')
    p.add_argument('--font', type=Path, help='CJK font for Chinese labels; otherwise use English labels')
    p.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    p.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'])
    p.add_argument('--width', type=int, default=1920)
    p.add_argument('--samples', type=int, default=64)
    p.add_argument('--frames', type=int, default=192)
    p.add_argument('--fps', type=int, default=24)
    p.add_argument('--frame', type=int, default=80)
    p.add_argument('--render', action='store_true')
    p.add_argument('--animate', action='store_true')
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists(): p.error('Choose a fresh output directory.')
    if args.frames < 24 or not 1 <= args.frame <= args.frames: p.error('Check frames and preview frame.')
    source = json.loads((args.asset_dir/'manifest.json').read_text())
    if not source['source_file'].upper().startswith('9SPI'):
        p.error('This chain/sequence map is specific to 9SPI.')
    profile = json.loads(args.style_profile.read_text())
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    mats = {k: material(k, v, profile['roughness'], profile['bump'])
            for k, v in profile['palette'].items()}
    target = np.array((-2.4, -.3, -2.2))
    # View along the exposed tRNA face: the former -X view hid the CCA ends
    # behind the tRNA bodies even when the whole ribosomal subunit was offset.
    position = target+50*np.array((-.254887, -.951251, -.173648))
    front = position-target
    front /= np.linalg.norm(front)
    right = np.cross((0, 0, 1), front)
    right /= np.linalg.norm(right)
    up = np.cross(front, right)
    context_offset = np.zeros(3)
    context = collection('01 Whole ribosomal subunits — layered presentation')
    active = collection('02 Intact tRNAs — identity retained throughout')
    rna_group = collection('03 mRNA — native residues, schematic three-nt transport')
    peptides = collection('04 Short peptide residue symbols — schematic')
    annotation = collection('05 Explanation and binding sites')
    hidden = []
    maximum_context_depth = -float('inf')
    context_mats = {}
    for color in ['alpha', 'beta', 'heads']:
        mat = material('Depth context '+color, profile['palette'][color], .64, profile['bump'])
        shader = mat.node_tree.nodes.get('Principled BSDF')
        rgba = shader.inputs['Base Color'].default_value
        shader.inputs['Base Color'].default_value = tuple(v*.40 for v in rgba[:3])+(1,)
        context_mats[color] = mat
    for item in source['chains']:
        chain = item['chain']
        if chain in ['At', 'Pt', 'Et', 'Mr']: continue
        large = chain.startswith('L')
        color = ('alpha' if large else 'beta') if item['nucleic_acid'] else 'heads'
        obj, chain_data = load_chain(args.asset_dir, chain, context_mats[color], context,
                            f'{"60S" if large else "40S"} complete chain {chain}')
        if large:
            maximum_context_depth = max(maximum_context_depth,
                                        float((chain_data['vertices'] @ front).max()))
        else:
            obj.hide_render = obj.hide_viewport = True
            hidden.append(chain)
    data = {c: np.load(args.asset_dir/f'chain_{c}.npz', allow_pickle=False)
            for c in ['At', 'Pt', 'Et', 'Mr']}
    trnas = {}
    for ch, color in [('At', 'protein'), ('Pt', 'rna'), ('Et', 'accent')]:
        obj = nucleic_atoms(data[ch], 'tRNA identity '+ch, mats[color], active)
        obj['source_nucleic_chain'] = ch
        trnas[ch] = obj
    fits = {'At': rigid_fit(data['At'], data['Pt']),
            'Pt': rigid_fit(data['Pt'], data['Et'])}
    mr = data['Mr']
    residue_frames = mrna_frames(mr)
    mr_objects = {}
    for res in residue_frames:
        selected = np.flatnonzero(mr['residue_numbers'] == res)
        index = {int(old): new for new, old in enumerate(selected)}
        bonds = [(index[int(a)], index[int(b)]) for a, b in mr['bonds']
                 if int(a) in index and int(b) in index]
        subset = {'atoms': mr['atoms'][selected], 'elements': mr['elements'][selected], 'bonds': bonds}
        ob = nucleic_atoms(subset, f'mRNA residue {res} {mr["residue_names"][selected[0]]}',
                           mats['dna1' if res in (5, 6, 7) else 'dna2'], rna_group)
        ob['source_nucleic_chain'], ob['residue_number'] = 'Mr', res
        mr_objects[res] = ob
    bridges = {r: path(f'mRNA phosphodiester {r}-{r+1}', [(0,0,0)]*2, .047,
                       mats['dna1'], rna_group) for r in range(min(residue_frames), max(residue_frames))}
    peptide_mat = material('Peptide residue symbols', 'F8E8CD', .42, .006)
    aa_mat = material('Incoming Ile residue symbol', 'FE9E56', .42, .006)
    beads = []
    for i, label in enumerate(['Gly', 'Ala', 'Arg', 'Ile']):
        parts = Parts()
        parts.ball((0, 0, 0), .18, seg=20, rings=12)
        beads.append(parts.object('Peptide residue '+label, [aa_mat if i == 3 else peptide_mat], peptides))
    peptide_links = [path('Peptide connection '+str(i), [(0,0,0)]*2, .045,
                          peptide_mat, peptides) for i in range(3)]
    ester_p = path('P-tRNA ester connection', [(0,0,0)]*2, .038, peptide_mat, peptides)
    ester_a = path('A-tRNA ester connection', [(0,0,0)]*2, .038, aa_mat, peptides)
    tip_a, tip_p = atom(data['At'], 76, "O3'"), atom(data['Pt'], 76, "O3'")
    upstream = tip_p-tip_a
    upstream /= np.linalg.norm(upstream)
    initial_c = tip_p+front*.42
    ile = tip_a+front*.42
    initial_peptide = np.array([initial_c+upstream*.38*(2-i) for i in range(3)])
    transferred_peptide = np.array([ile+upstream*.38*(3-i) for i in range(3)])
    poses = []
    bridge_lengths = []
    minimum_actor_depth = float((mr['atoms'] @ front).min())
    for frame in range(1, args.frames+1):
        t = (frame-1)/(args.frames-1)
        dock = entry_progress(t)
        transfer = smooth(t, .37, .48)
        release = smooth(t, .49, .62)
        step = smooth(t, .62, .90)
        ra, ta = (entry_pose(data['At'], dock, right, up, front) if dock < 1 else
                  trajectory(data['At'], fits['At'], step, front))
        rp, tp = trajectory(data['Pt'], fits['Pt'], step, front)
        re, te = np.eye(3), release*(right*6-up*9+front*4)
        frame_poses = {}
        for ch, rotation, translation in [('At',ra,ta), ('Pt',rp,tp), ('Et',re,te)]:
            move(trnas[ch], rotation, translation, frame)
            minimum_actor_depth = min(minimum_actor_depth,
                float(((data[ch]['atoms'] @ rotation.T+translation) @ front).min()))
            frame_poses[ch] = {'rotation':rotation.tolist(), 'translation':translation.tolist()}
        mr_poses = {}
        for res, obj in mr_objects.items():
            center, orientation = track_pose(residue_frames, res-3*step)
            rotation = orientation @ residue_frames[res][1].T
            translation = center-rotation @ residue_frames[res][0]
            move(obj, rotation, translation, frame)
            mr_poses[res] = rotation, translation
        for res, obj in bridges.items():
            r1, d1 = mr_poses[res]
            r2, d2 = mr_poses[res+1]
            a = r1 @ atom(mr, res, "O3'")+d1
            b = r2 @ atom(mr, res+1, 'P')+d2
            update_line(obj, [a, b], frame)
            bridge_lengths.append(float(np.linalg.norm(b-a)))
        peptide = ((1-transfer)*initial_peptide+transfer*transferred_peptide) @ ra.T+ta
        # Peptide is attached to the stationary P-tRNA before transfer; at that
        # stage A is already docked, so its transform is identity.
        if transfer == 0:
            peptide = initial_peptide.copy()
        aa = ra @ ile+ta
        points = list(peptide)+[aa]
        for obj, xyz in zip(beads, points):
            obj.location = xyz
            obj.keyframe_insert(data_path='location', frame=frame)
        for i, obj in enumerate(peptide_links):
            update_line(obj, [points[i], points[i+1]], frame)
        peptide_links[2].data.bevel_depth = .045*transfer
        peptide_links[2].data.keyframe_insert(data_path='bevel_depth', frame=frame)
        update_line(ester_p, [rp @ tip_p+tp, points[2]], frame)
        ester_p.data.bevel_depth = .038*(1-transfer)
        ester_p.data.keyframe_insert(data_path='bevel_depth', frame=frame)
        update_line(ester_a, [ra @ tip_a+ta, aa], frame)
        poses.append({'frame':frame, 'dock':dock, 'peptide_transfer':transfer,
                      'old_E_release':release, 'codon_shift_nt':3*step,
                      'trnas':frame_poses, 'peptide_residue_positions':np.array(points).tolist()})
    # A labeled depth-separated view must actually expose the complete active
    # molecules. Derive the whole-60S offset from geometry instead of guessing.
    context_offset = -front*(maximum_context_depth-minimum_actor_depth+2.0)
    for obj in context.objects:
        if obj['source_chain'].startswith('L'):
            obj.location = context_offset
            obj['presentation_offset_nm'] = context_offset.tolist()
    profile['lighting'].append({'name':'Active site softbox','position':[-8,-25,10],
                                'power':13000,'color':'FFE9D4','size':16})
    camera = studio(profile, args, {'position':position.tolist(), 'target':target.tolist(),
                                   'scale':26., 'light_scale':1.})
    font = bpy.data.fonts.load(str(args.font)) if args.font else None
    light_background = np.dot(linear(profile['background']), (.2126,.7152,.0722)) > .45
    ink = emission('Readable annotation', '24323A' if light_background else 'E5EFF2')
    dim = emission('Secondary annotation', '465660' if light_background else 'AFC3CB')
    zh = bool(font)
    stages = [(0,.34, '氨酰 tRNA 进入 A 位点' if zh else 'Aminoacyl-tRNA enters the A site'),
              (.34,.49, '肽链从 P 位点转移至 A 位点 tRNA' if zh else 'Peptide transfers from P-tRNA to A-tRNA'),
              (.49,.62, '原 E 位点的空载 tRNA 离开' if zh else 'The previous E-site tRNA leaves'),
              (.62,.90, '转位：A → P，P → E；mRNA 前移一个密码子' if zh else 'Translocation: A to P, P to E; mRNA advances one codon'),
              (.90,1.01, '肽链增加一个残基，完成一次底物转位' if zh else 'One residue added; substrate translocation complete')]
    for i, (a,b,body) in enumerate(stages):
        ob = text_object(body, f'Stage {i+1}', .42, ink, annotation, font)
        ob.parent = camera
        ob.location = (-11.6, 7.35, -13)
        for frame in range(1,args.frames+1):
            t=(frame-1)/(args.frames-1)
            ob.hide_render = not a <= t < b
            ob.keyframe_insert(data_path='hide_render',frame=frame)
    footer = ('分层观察：40S 隐藏，60S 后移｜运动与短肽为教学示意' if zh else
              'Layered view: 40S hidden, 60S offset | Motion and short peptide are schematic')
    ob = text_object(footer, 'Presentation scope', .26, dim, annotation, font)
    ob.parent = camera
    ob.location = (-11.6, -7.5, -13)
    for ch, name in [('At','A'),('Pt','P'),('Et','E')]:
        center = data[ch]['atoms'].mean(0)
        location = target+right*5.5+up*np.dot(center-target,up)+front*8
        ob = text_object(name+(' 位点' if zh else ' site'), 'Binding site '+name, .43, ink, annotation, font)
        ob.location = location
        ob.rotation_euler = camera.rotation_euler
    for terminal in ['N', 'C']:
        ob = text_object(terminal, 'Peptide '+terminal+' terminus', .30, ink, annotation)
        ob.rotation_euler = camera.rotation_euler
        ob.data.align_x = ob.data.align_y = 'CENTER'
        for state in poses:
            points = np.array(state['peptide_residue_positions'])
            transfer = state['peptide_transfer']
            xyz = points[0] if terminal == 'N' else (1-transfer)*points[2]+transfer*points[3]
            direction = points[0]-points[1] if terminal=='N' else (
                (1-transfer)*(points[2]-points[1])+transfer*(points[3]-points[2]))
            screen_direction = right*np.dot(direction,right)+up*np.dot(direction,up)
            screen_direction /= np.linalg.norm(screen_direction)
            ob.location = xyz+front*20+screen_direction*.6
            ob.keyframe_insert(data_path='location',frame=state['frame'])
    for res, label, dx, dy in [(2,'5′',3.4,-.25),
                               (8,'mRNA 片段 3′' if zh else 'mRNA fragment 3′',-3.4,.25)]:
        ob = text_object(label, 'mRNA fragment polarity '+label, .29, ink, annotation, font)
        ob.rotation_euler = camera.rotation_euler
        ob.data.align_x = 'LEFT' if dx>0 else 'RIGHT'
        guide = path('Polarity guide '+str(res), [(0,0,0)]*2, .011, dim, annotation)
        for state in poses:
            center, _ = track_pose(residue_frames,res-state['codon_shift_nt'])
            ob.location = center+front*20+right*dx+up*dy
            ob.keyframe_insert(data_path='location',frame=state['frame'])
            update_line(guide,[center+front*20, np.array(ob.location)-right*math.copysign(.18,dx)],state['frame'])
    # Local glyph contrast matters: labels may cross a dark molecular surface
    # even in a light-background style. A small opposite-color halo stays legible.
    outline_mat = emission('Annotation halo', 'F6F8F8' if light_background else '10212B')
    for original in list(annotation.objects):
        if original.type != 'FONT': continue
        curve = original.data.copy()
        curve.offset = .018
        curve.materials.clear()
        curve.materials.append(outline_mat)
        outline = bpy.data.objects.new('Annotation outline '+original.name,curve)
        annotation.objects.link(outline)
        outline.parent = original
        outline.location = (0,0,-.025)
        outline.visible_shadow = False
        if original.name.startswith('Stage '):
            a,b,_ = stages[int(original.name.split()[-1])-1]
            for frame in range(1,args.frames+1):
                t=(frame-1)/(args.frames-1)
                outline.hide_render = not a<=t<b
                outline.keyframe_insert(data_path='hide_render',frame=frame)
    bpy.ops.file.pack_all()
    manifest = {'scene':'Human translation elongation substrate cycle', 'source_asset':source,
        'coordinate_units':'1 Blender unit = 1 nm',
        'presentation':{'hidden_complete_chains':hidden, 'large_subunit_offset_nm':context_offset.tolist(),
                        'molecules_sliced':False,'meaning':'Explicit layered view, not ribosome dissociation.'},
        'motion':'Authored rigid-body tRNA trajectories; substrate identities persist. No inferred atomistic pathway.',
        'entry_motion':{'path':'Cubic Bezier from the open upper-left foreground, with a 16-degree turn before docking.',
                        'easing':'Quintic smoothstep: zero endpoint velocity and acceleration.',
                        'control_points_camera_basis_nm':[[-4,2.5,2],[-3.5,4,2],[-1.5,1,.25],[0,0,0]],
                        'basis_order':['right','up','front'],
                        'rotation_axis':'camera front','initial_angle_degrees':-16,
                        'scope':'Authored presentation path; eEF1A-mediated delivery and accommodation are not modeled.'},
        'trna_endpoint_fit_rmsd_nm':{ch:fit[2] for ch,fit in fits.items()},
        'mrna':{'source_residues':list(residue_frames),'source_sequence':'CGUAUCC',
                'transport_nt':3,'A_codon_residues':[5,6,7], 'A_codon':'AUC',
                'P_codon_residues':[2,3,4],'P_codon':'CGU',
                'bridge_length_nm_range':[min(bridge_lengths),max(bridge_lengths)],
                'motion':'Residue poses move along an interpolated backbone track; bonds retained as display connectors.'},
        'peptide':{'residue_symbols':['Gly','Ala','Arg','Ile'],'initial_chain_length':3,'final_chain_length':4,
                   'evidence':'Illustrative short chain; source P tRNA decodes Arg and A tRNA decodes Ile. Peptide not deposited in 9SPI.'},
        'frame_poses':poses, 'stages':[{'start':a,'end':b,'text':body} for a,b,body in stages],
        'sources':['https://www.rcsb.org/structure/9SPI',
                   'https://doi.org/10.1016/j.cell.2015.03.052'],
        'limitations':['Not a complete translation cycle: initiation, termination and recycling are not shown.',
                       'Elongation factors, GTP chemistry, subunit rotations and hybrid tRNA states are omitted.',
                       'E-site release timing is an authored sequence, not a universal kinetic rule.',
                       'Source has only modeled residues; unresolved tRNA regions and mRNA flanks remain absent.',
                       'tRNA endpoint rigid fits are approximate; displayed peptide residues are symbols, not atomic chemistry.',
                       '60S is offset as an intact assembly and 40S is hidden for visibility; restore source view with build_complex_scene.py.']}
    save_and_render(args, manifest)


if __name__ == '__main__':
    main()
