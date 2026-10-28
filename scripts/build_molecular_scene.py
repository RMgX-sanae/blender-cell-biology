"""Detailed molecular scenes using coordinate-derived assets and selectable style.

Example: blender -b --python build_molecular_scene.py -- --scene membrane
 --asset-dir /absolute/beta2 --output /absolute/new_scene --render
"""
import argparse
import json
import math
from pathlib import Path
import random
import sys

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from distributions import patch_points
from detail_scene import (Parts, collection, duplicate, load_chain, material,
                          path, save_and_render, studio)


def materials(profile):
    return {k: material(k, color, profile['roughness'], profile['bump'])
            for k, color in profile['palette'].items()}


def membrane_front(x, view):
    return -3.3 + (4.2*math.exp(-.5*(x/4.0)**2) if view=='arc' else 0)


def membrane(args, mats):
    c = collection('01 beta2AR and Gs — 3SN6 experimental coordinates')
    chain_names = {'R': ('receptor', 'beta2 adrenergic receptor'),
                   'A': ('alpha', 'Gs alpha'), 'B': ('beta', 'Gs beta'),
                   'G': ('gamma', 'Gs gamma')}
    points = []
    placements = [(0, 0, 0)]
    if args.environment:
        placements += [(-8.5, 11, -.45), (8.0, 18, .5)]
    for chain, (color, label) in chain_names.items():
        obj, data = load_chain(args.asset_dir, chain, mats[color], c, label, args.representation)
        for index,(px,py,angle) in enumerate(placements):
            rotation = np.array([[math.cos(angle),-math.sin(angle),0],
                                 [math.sin(angle),math.cos(angle),0],[0,0,1]])
            coords = data['atoms'] @ rotation.T + (px,py,0)
            points.extend(coords[np.abs(coords[:,2])<2.45].tolist())
            if index:
                copy=duplicate(obj,f'{label} background complex {index}',c)
                copy.location=(px,py,0)
                copy.rotation_euler.z=angle
    avoid = np.array(points)
    group = collection('02 Bilayer — schematic lipid headgroups and two tails')
    heads, tails = Parts(), Parts()
    rng = random.Random(42)
    lipid_placements = []
    width, height = (40, 30) if args.environment else (26.4, 13.6)
    for side, seed in [(-1,42),(1,43)]:
        candidates = patch_points(int(width*height*1.4), width, height,
                                  min_distance=.62, seed=seed)
        for x, y in candidates:
            yy = y + height/2 - 3.3
            if yy < membrane_front(x,args.membrane_view):continue
            if np.min(np.linalg.norm(avoid[:, :2] - (x, yy), axis=1)) < .62:continue
            h = .12 * math.sin(x * .32) * math.sin(yy * .35)
            z = side * 1.8 + h
            theta = rng.random() * math.tau
            u = np.array((math.cos(theta), math.sin(theta), 0))
            v = np.array((-math.sin(theta), math.cos(theta), 0))
            centre = np.array((x, yy, z))
            heads.ball(centre, (.31, .31, .27), seg=14, rings=9)
            lipid_placements.append([x, yy, z, side])
            for tail_side in (-1, 1):
                vertices = []
                bend=rng.uniform(-.21,.21)
                length=rng.uniform(1.36,1.52)
                for k in range(11):
                    t = k / 10
                    position = centre + u * (tail_side * .12 + .16*math.sin(math.pi*t))
                    position += v * (((k % 2) * 2 - 1) * .045 + bend*t*t)
                    position[2] -= side * (.25 + length*t)
                    vertices.append(tuple(position))
                tails.tube(vertices, .067, sides=6)
    head_object = heads.object('Two resolved hydrophilic leaflets', [mats['heads']], group)
    tail_object = tails.object('Paired hydrophobic chains in the bilayer core', [mats['tails']], group)
    for obj in (head_object, tail_object):
        obj['evidence'] = 'schematic lipid geometry; authored packing, no density measurement'
    environment_info = molecular_environment(args, mats) if args.environment else None
    camera = {'position': (22, -45, 16), 'target': (0, 2, -2.5), 'scale': 34}
    if args.environment:
        camera.update({'position':(17,-41,12),'target':(.5,3,-3.4),'type':'PERSP','lens':64,'fstop':.22})
    if args.membrane_view=='side':
        camera.update({'position':(8,-46,3),'target':(0,5,-3),'type':'PERSP','lens':62,'fstop':.3})
    return {'scene': 'beta2AR–Gs in membrane', 'lipids': len(lipid_placements),
            'lipid_positions_nm_and_leaflet': lipid_placements,
            'receptor_placements_xy_nm_and_rotation': placements,
            'environment': environment_info,
            'membrane_head_plane_separation_nm': 3.6,
            'membrane_view':args.membrane_view,
            'coordinate_units': '1 Blender unit = 1 nm',
            'limitations': ['Membrane composition, thickness, curvature and lipid packing are illustrative.',
                            'Side view has no receptor-centered notch; optional arc view is a local cutaway. Patch boundaries are display limits, not biological membrane edges.',
                            'The rigid receptor–Gs structure is a selected engineered active complex, not its dynamics.',
                            'Receptor orientation is inferred from TM annotation; lipids are excluded around projected protein coordinates.']}, \
           camera


def molecular_environment(args, mats):
    """Actual molecular shapes; authored counts/positions, with an empty hero region."""
    if args.crowder_dir is None or args.ribosome_dir is None:
        raise ValueError('Environment needs --crowder-dir (1UBQ) and --ribosome-dir (9SPI) prepared assets.')
    group=collection('03 Cytosolic depth layers — schematic placement of source-derived molecules')
    distant=material('Distant molecular environment','255667',.62,.006)
    nearer=material('Near peripheral molecular environment','437A84',.6,.006)
    ribo=json.loads((args.ribosome_dir/'manifest.json').read_text())
    if not ribo['source_file'].upper().startswith('9SPI'):
        raise ValueError('This depth layout is sized for 9SPI; prepare another layout for another structure.')
    for item in ribo['chains']:
        obj,_=load_chain(args.ribosome_dir,item['chain'],distant,group,f'Background ribosome {item["chain"]}')
        obj.location=(-22,18,-26)
        obj.rotation_euler.z=.7
        copy=duplicate(obj,f'Foreground ribosome {item["chain"]}',group)
        copy.location=(23,-15,-29)
        copy.rotation_euler=(.2,.5,-.7)
        copy.data=obj.data.copy()
        copy.data.materials.clear()
        copy.data.materials.append(nearer)
    ubq=json.loads((args.crowder_dir/'manifest.json').read_text())
    if not ubq['source_file'].upper().startswith('1UBQ'):
        raise ValueError('This crowding example is sized for ubiquitin 1UBQ.')
    molecule,_=load_chain(args.crowder_dir,'A',nearer,group,'Ubiquitin — source shape')
    rng=random.Random(611)
    positions=[]
    for i in range(24):
        for _ in range(1000):
            xyz=(rng.uniform(-18,20),rng.uniform(2,30),rng.uniform(-25,-8))
            # Reserve the hero silhouette and the two much larger ribosomes.
            if -7<xyz[0]<8 and xyz[1]<10:continue
            if any(np.linalg.norm(np.array(xyz)-(px,py,-4.5))<8.5
                   for px,py in [(0,0),(-8.5,11),(8,18)]):continue
            if np.linalg.norm(np.array(xyz)-(-22,18,-26))<17:continue
            if np.linalg.norm(np.array(xyz)-(23,-15,-29))<17:continue
            if any(np.linalg.norm(np.array(xyz)-p)<4.6 for p in positions):continue
            positions.append(xyz)
            break
        else:raise RuntimeError('Crowder placement exhausted; reduce the authored count.')
        obj=molecule if i==0 else duplicate(molecule,f'Ubiquitin depth instance {i+1}',group)
        obj.location=xyz
        obj.rotation_euler=tuple(rng.uniform(-math.pi,math.pi) for _ in range(3))
    return {'ribosome_source':ribo['entry_url'],'ubiquitin_source':ubq['entry_url'],
            'ribosome_centers_nm':[[-22,18,-26],[23,-15,-29]],
            'ubiquitin_centers_nm':positions,
            'interpretation':'Illustrative multiscale depth; source geometries retain nm scale, counts/positions are not measured.',
            'optics':'Art-directed shallow depth of field at scene units; f-number is not recovered from the reference image.'}


def nucleic_atoms(data, name, color, group):
    p = Parts()
    atoms, elements = data['atoms'], data['elements']
    radii = {'C': .105, 'N': .103, 'O': .10, 'P': .14, 'S': .13}
    for point, element in zip(atoms, elements):
        p.ball(point, radii.get(element, .105), seg=10, rings=7)
    for a, b in data['bonds']:
        p.tube([atoms[a], atoms[b]], .047, sides=6)
    obj = p.object(name, [color], group)
    obj['evidence'] = 'experimental atom coordinates; bond candidates inferred by distance, display radii are artistic'
    return obj


def transcription(args, mats):
    group = collection('01 Pol II elongation complex — 5FLM')
    for chain in 'ABCDEFGHIJKL':
        color = mats['protein'] if chain in 'AB' else mats['alpha']
        load_chain(args.asset_dir, chain, color, group, f'Pol II subunit {chain}', args.representation)
    nucleic = collection('02 Experimental DNA and RNA')
    for chain, color, name in [('N', 'dna1', 'DNA non-template chain N'),
                               ('T', 'dna2', 'DNA template chain T'),
                               ('P', 'rna', 'RNA chain P')]:
        data = np.load(args.asset_dir / f'chain_{chain}.npz', allow_pickle=False)
        nucleic_atoms(data, name, mats[color], nucleic)
    return {'scene': 'Mammalian transcribing RNA polymerase II',
            'coordinate_units': '1 Blender unit = 1 nm',
            'limitations': ['Static experimental elongation complex; no claim of simulated nucleotide addition.',
                            'Protein Gaussian envelopes and atom display radii are rendering representations.',
                            'Only modeled nucleic acid residues in 5FLM are displayed; missing regions are not invented.']}, \
           {'position': (23, -34, 21), 'target': (0, 0, 1.4), 'scale': 34, 'light_scale': .7}


def motor(args, mats):
    proto = collection('01 Microtubule — illustrative 13-protofilament lattice')
    originals = {}
    atom_data = {}
    for chain, color in [('A', 'alpha'), ('B', 'beta'), ('K', 'protein')]:
        obj, data = load_chain(args.asset_dir, chain, mats[color], proto, f'4HNA chain {chain}', args.representation)
        originals[chain], atom_data[chain] = obj, data
    # A local lattice example, not an assembly inferred from the crystal unit cell.
    for chain in 'AB':
        originals[chain].hide_render = True
        originals[chain].hide_viewport = True
        for row in range(13):
            angle = math.tau * row / 13
            for step in range(-3, 4):
                obj = duplicate(originals[chain], f'{chain} protofilament {row+1} dimer {step+4}', proto)
                obj.hide_render, obj.hide_viewport = False, False
                obj.rotation_euler.x = angle
                obj.location = (step * 8 + row * .25, -12 * math.sin(angle), 12 * math.cos(angle))
    head = originals['K']
    head.location = (0, 0, 12)
    head.name = 'Leading kinesin motor — K coordinate-derived shape'
    other = duplicate(head, 'Trailing motor — schematic duplicated head', proto)
    other.location.x = -8
    cargo_group = collection('02 Cargo and coiled-coil — schematic')
    # The moving assembly illustrates direction and discrete stepping, not ATP chemistry.
    cargo = bpy.data.objects.new('Transport assembly', None)
    cargo_group.objects.link(cargo)
    centroid = atom_data['K']['atoms'].mean(0)
    necks = [path('Flexible schematic neck linker', [(0, 0, 0)] * 20,
                  .13, mats['protein'], cargo_group) for _ in range(2)]
    trajectories = []
    for frame in range(1, args.frames + 1):
        phase = (frame - 1) / max(args.frames - 1, 1) * 2
        positions = []
        for index, obj in enumerate((other, head)):
            t = min(1, max(0, phase - index))
            ease = t * t * (3 - 2 * t)
            x = (-8 if index == 0 else 0) + 16 * ease
            lift = 2.1 * math.sin(math.pi * t) ** 2
            # The rigid coordinate-derived heads need an explicit side passage;
            # otherwise the mid-swing silhouettes and heavy atoms interpenetrate.
            lateral = 7.5 * math.sin(math.pi * t) ** 2
            obj.location = (x, lateral, 12 + lift)
            obj.keyframe_insert(data_path='location', frame=frame)
            positions.append((x, lateral, lift))
        body_x = sum(p[0] for p in positions) / 2
        cargo.location = (body_x, 0, 0)
        cargo.keyframe_insert(data_path='location', frame=frame)
        for index, ((x, lateral, lift), neck) in enumerate(zip(positions, necks)):
            a = Vector((centroid[0] + x, centroid[1] + lateral, 18.4 + lift))
            b = Vector((centroid[0] + body_x, centroid[1], 23))
            for j, point in enumerate(neck.data.splines[0].points):
                t = j / 19
                p = a.lerp(b, t)
                p.y += .5 * math.sin(math.pi * t) * (1 if index == 0 else -1)
                point.co = (*p, 1)
                point.keyframe_insert(data_path='co', frame=frame)
        trajectories.append({'frame': frame, 'head_offsets_xyz_nm': positions, 'cargo_x_nm': body_x})
    for sign in (-1, 1):
        pts = [(centroid[0] + .35 * math.cos(t * 9 + sign * math.pi / 2),
                centroid[1] + .35 * math.sin(t * 9 + sign * math.pi / 2),
                23 + t * 7) for t in np.linspace(0, 1, 120)]
        obj = path('Schematic coiled-coil strand', pts, .18, mats['protein'], cargo_group)
        obj.parent = cargo
    # Vesicle diameter reduced for the close view; explicitly not an intact measured cargo.
    bpy.ops.mesh.primitive_uv_sphere_add(segments=64, ring_count=32, radius=6.5,
                                      location=(centroid[0], centroid[1], 35))
    vesicle = bpy.context.object
    vesicle.name = 'Illustrative vesicle envelope — scale chosen for local view'
    vesicle.data.materials.append(mats['heads'])
    for p in vesicle.data.polygons:
        p.use_smooth = True
    vesicle.parent = cargo
    for old in list(vesicle.users_collection):
        old.objects.unlink(vesicle)
    cargo_group.objects.link(vesicle)
    lipids = Parts()
    for i in range(1800):
        z = 1 - 2 * (i + .5) / 1800
        angle = math.pi * (3 - math.sqrt(5)) * i
        r = math.sqrt(1 - z * z)
        p = (centroid[0] + 6.55 * r * math.cos(angle),
             centroid[1] + 6.55 * r * math.sin(angle), 35 + 6.55 * z)
        lipids.ball(p, .23, seg=8, rings=6)
    shell = lipids.object('Vesicle outer lipid headgroups — schematic packing', [mats['heads']], cargo_group)
    shell.parent = cargo
    return {'scene': 'Kinesin-associated cargo transport',
            'coordinate_units': '1 Blender unit = 1 nm for source assets',
            'trajectories': trajectories,
            'limitations': ['Microtubule lattice is a simplified 13-protofilament placement; seam, nucleotide heterogeneity and exact helical assembly are omitted.',
                            'Motor head duplication, coiled-coil, cargo size and motion are schematic.',
                            'Alternating stationary contacts and 16-nm head swings illustrate hand-over-hand transport; the lateral clearance path, linker deformation and timings are authored, not measured or mechanochemical simulation.']}, \
           {'position': (40, -88, 41), 'target': (4, 0, 13), 'scale': 106, 'light_scale': 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scene', choices=['membrane', 'transcription', 'motor'], required=True)
    parser.add_argument('--asset-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--style-profile', type=Path,
                        default=Path(__file__).resolve().parent.parent / 'assets/style-studio.json')
    parser.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    parser.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'])
    parser.add_argument('--width', type=int, default=1600)
    parser.add_argument('--samples', type=int, default=48)
    parser.add_argument('--fps', type=int, default=24)
    parser.add_argument('--frames', type=int, default=96)
    parser.add_argument('--frame', type=int, default=1)
    parser.add_argument('--render', action='store_true')
    parser.add_argument('--animate', action='store_true')
    parser.add_argument('--environment', action='store_true', help='Add authored molecular depth to the membrane scene')
    parser.add_argument('--crowder-dir', type=Path)
    parser.add_argument('--ribosome-dir', type=Path)
    parser.add_argument('--membrane-view',choices=['arc','side'],default='side')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    if args.environment and args.scene != 'membrane':
        parser.error('--environment is implemented for the membrane scene only.')
    if args.animate and args.scene != 'motor':
        parser.error('Only motor has a biological-process rig here. Static structures need a project-specific animation.')
    if args.output.exists():
        raise SystemExit('Choose a fresh output directory; existing work is not overwritten.')
    source = json.loads((args.asset_dir / 'manifest.json').read_text())
    required = {'membrane': 'beta2', 'transcription': 'transcription', 'motor': 'motor'}
    if source['preset'] != required[args.scene]:
        raise SystemExit('Asset preset does not match this scene.')
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile = json.loads(args.style_profile.read_text())
    args.representation = profile.get('representation', 'surface')
    mats = materials(profile)
    manifest, camera = globals()[args.scene](args, mats)
    manifest['source_asset'] = source
    studio(profile, args, camera)
    save_and_render(args, manifest)


if __name__ == '__main__':
    main()
