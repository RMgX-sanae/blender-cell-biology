"""Detailed continuous-inner-membrane mitochondrion; selectable image-reference styles.

Schematic shapes and exaggerated proteins, not a calibrated tomography model.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
import organelle_geometry as geo
from detail_scene import collection, material, save_and_render, studio


def build_mitochondrion(profile):
    groups = {key: collection(name) for key, name in [
        ('membrane', '01 Outer and continuous inner membranes'),
        ('detail', '02 Bilayer cut edges and headgroups'),
        ('complexes', '03 Matrix-facing synthase dimer rows'),
        ('matrix', '04 Matrix molecules and circular mtDNA'),
        ('prototypes', '05 Hidden editable protein prototypes')]}
    colors = {'outer': '4E9194', 'outer_light': '91BCBA', 'inner': 'CF8D6C',
              'lipid': 'DBAA84', 'rim': 'E6D5B6', 'seam': '574A43',
              'protein': 'E3CE9A', 'gold': 'B99560', 'protein_dark': '4B797B',
              'rust': 'B46D55', 'dna': '8BC0C7', 'matrix_protein': 'A3B4A0',
              'carbon': '4B5457', 'oxygen': 'C65C4C', 'bond': 'A0ABA9'}
    colors.update(profile.get('organelle_palette', {}))
    mats = {name: material(name, color, profile['roughness'], .005)
            for name, color in colors.items()}
    geo.membranes(groups['membrane'], mats)
    geo.membrane_detail(groups['detail'], mats)
    dimers = geo.proteins(groups['complexes'], groups['prototypes'], mats)
    geo.matrix_contents(groups['matrix'], groups['prototypes'], mats)
    return groups, dimers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--style-profile', type=Path,
                        default=Path(__file__).resolve().parent.parent / 'assets/style-organelle.json')
    parser.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    parser.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'])
    parser.add_argument('--width', type=int, default=2400)
    parser.add_argument('--samples', type=int, default=64)
    parser.add_argument('--frames', type=int, default=96)
    parser.add_argument('--fps', type=int, default=24)
    parser.add_argument('--frame', type=int, default=1)
    parser.add_argument('--render', action='store_true')
    parser.add_argument('--animate', action='store_true', help='Slow camera inspection, not biological motion')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    if args.output.exists():
        raise SystemExit('Choose a fresh output directory.')
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile = json.loads(args.style_profile.read_text())
    groups, dimers = build_mitochondrion(profile)
    camera = studio(profile, args, {'position': (5, -12, 10), 'target': (0, 0, -.03),
                                   'scale': 12.6, 'light_scale': .34})
    if args.animate:
        camera.keyframe_insert(data_path='location', frame=1)
        camera.keyframe_insert(data_path='rotation_euler', frame=1)
        from detail_scene import aim
        camera.location = (3.5, -12.4, 9.2)
        aim(camera, (0, 0, -.03))
        camera.keyframe_insert(data_path='location', frame=args.frames)
        camera.keyframe_insert(data_path='rotation_euler', frame=args.frames)
    manifest = {'scene': 'Mitochondrion with a continuous inner membrane',
                'synthase_dimers': dimers,
                'evidence': 'Schematic hierarchical organelle geometry; not atomic or tomography data',
                'removed_quadrant': 'y < 0 and z > 0',
                'source_geometry': 'Adapted from user workspace mitochondrion_astra_20260905, original preserved',
                'limitations': ['Protein and lipid sizes/densities are exaggerated for the organelle view.',
                               'Crista junctions and local membrane thickness are not fitted to EM.',
                               'ATP synthases are procedural motifs oriented into matrix; no atomic conformation is claimed.',
                               'Optional animation is camera inspection only, not ATP production or mitochondrial dynamics.']}
    save_and_render(args, manifest)


if __name__ == '__main__':
    main()
