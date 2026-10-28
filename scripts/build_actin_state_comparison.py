"""Render three source-derived actin nucleotide states, without cutting molecules.

This is a static state comparison, not a hydrolysis animation or full actin cycle.
Requires Blender 5.x and assets from prepare_actin_states.py.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import Parts, collection, material, mesh, save_and_render, studio


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--profile', type=Path, required=True)
    p.add_argument('--font', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    p.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'], default='clean')
    p.add_argument('--width', type=int, default=1800)
    p.add_argument('--samples', type=int, default=64)
    p.add_argument('--render', action='store_true')
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    args.fps, args.frames, args.frame, args.animate = 24, 1, 1, False
    if args.output.exists():
        p.error('Use a fresh output directory')
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile = json.loads(args.profile.read_text())
    assets = json.loads((args.asset_dir / 'manifest.json').read_text())
    group = collection('01 Whole F-actin subunits and retained ligands')
    insets = collection('02 Ligand-only insets — enlarged 3 times')
    labels = collection('03 Editable Chinese annotations')
    protein_mat = material('Source-derived F-actin envelope', profile['palette']['alpha'], .55, 0)
    # Element colours are a display convention, not a nucleotide-state measurement.
    colors = {'C': '87939D', 'N': '487DD0', 'O': 'D35151', 'P': 'D99B38',
              'Mg': '6DA65A', 'Be': 'B271B0', 'F': '5EA09A'}
    elements = list(colors)
    mats = [material('Element '+e, c, .48, 0) for e,c in colors.items()]
    font = bpy.data.fonts.load(str(args.font))
    font.pack()
    light_background = sum(int(profile['background'][i:i+2],16) for i in (0,2,4)) > 400
    ink = material('Annotation ink', '263F4A' if light_background else 'E0ECEF', .8, 0)
    shader = ink.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Emission Color'].default_value = ink.diffuse_color
    shader.inputs['Emission Strength'].default_value = .7

    def text(body, x, z, size=.36):
        data = bpy.data.curves.new(body, 'FONT')
        data.body, data.font, data.size, data.align_x = body, font, size, 'CENTER'
        obj = bpy.data.objects.new(body, data)
        labels.objects.link(obj)
        obj.location, obj.rotation_euler = (x, -6, z), (np.pi/2, 0, 0)
        data.materials.append(ink)
        return obj

    # The top row retains molecular scale. Only the explicitly marked ligand
    # insets below use a different scale; no protein is sectioned or exploded.
    for i, record in enumerate(assets['states']):
        code, x = record['id'], (i-1)*8.2
        array = np.load(args.asset_dir / (code+'.npz'), allow_pickle=False)
        offset = np.array([x, 0, 2.2])
        obj = mesh(code+' whole modeled chain C', (array['vertices']+offset).tolist(),
                   array['faces'].tolist(), [protein_mat], group)
        obj['source_pdb'], obj['author_chain'] = code, 'C'
        obj['state'], obj['evidence'] = record['state'], 'experimental-coordinate-derived envelope'
        obj['not_G_actin'] = True
        adp = next(r for r in record['ligands'] if r['ccd']=='ADP')
        center = np.mean(adp['positions_nm'], axis=0)
        for ligand in record['ligands']:
            points = np.array(ligand['positions_nm'])
            for inset in (False, True):
                parts = Parts()
                loc = (points-center)*3+np.array([x, -1.5, -5.0]) if inset else points+offset
                for pos, element in zip(loc, ligand['elements']):
                    radius = .13 if element=='Mg' else .095
                    parts.ball(pos, radius*(3 if inset else 1), elements.index(element), seg=16, rings=10)
                # Display bonds within one CCD residue only; no Mg coordination
                # or ADP-to-Pi covalent bond is invented.
                cov = {'C': .076, 'N': .071, 'O': .066, 'P': .107, 'F': .057, 'Be': .096}
                for a in range(len(points)):
                    for b in range(a):
                        ea, eb = ligand['elements'][a], ligand['elements'][b]
                        if ea not in cov or eb not in cov:
                            continue
                        if .04 < np.linalg.norm(points[a]-points[b]) < cov[ea]+cov[eb]+.035:
                            parts.tube([loc[a],loc[b]], .032*(3 if inset else 1), elements.index(ea), sides=8)
                target = parts.object(code+' '+ligand['ccd']+(' inset x3' if inset else ' native site'),
                                      mats, insets if inset else group)
                target['source_pdb'], target['ccd'] = code, ligand['ccd']
                target['source_residue'] = ligand['author_residue']
                target['representation_scale'] = 3 if inset else 1
                target['bond_display'] = 'distance-based within-residue covalent display only; not bond order'
        text(code, x, 7.8, .48)
        text(['Mg–ADP–BeF3', 'Mg–ADP·Pi', 'Mg–ADP'][i], x, 7.05, .49)
        text(['ATP 状态的类似物', '水解后，Pi 仍保留', 'Pi 已释放后的状态'][i], x, 6.4, .34)
        text(['ADP + BeF3 + Mg2+', 'ADP + Pi + Mg2+', 'ADP + Mg2+'][i], x, -7.8, .35)
    text('核苷酸状态：保留真实配体，区分水解与磷酸释放', 0, 9.25, .52)
    text('上：完整的 F-actin 结构亚基（兔骨骼肌 α-actin）   下：配体单独放大 ×3', 0, 8.65, .31)
    text('结构状态对照；不是反应轨迹。8A2R 含 ADP–BeF3，并非真实 ATP。', 0, -8.6, .31)
    text('来源：Oosterheert et al., Nature 2022 / PDB 8A2R · 8A2S · 8A2T', 0, -9.15, .28)
    studio(profile, args, {'position': (0,-55,0), 'target': (0,0,0), 'scale': 32, 'light_scale': .9})
    save_and_render(args, {'purpose': 'Static molecular-state asset validation for full-cycle planning',
                           'source_assets': assets, 'full_cycle_rendered': False,
                           'ligand_inset_magnification': 3,
                           'no_molecular_cut': True})


if __name__ == '__main__':
    main()
