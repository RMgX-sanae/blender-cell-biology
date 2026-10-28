"""Animal-cell overview with connected ER/NE and detailed organelle instances.

The cell is an authored teaching model, not a reconstruction of a named cell type.
"""
import argparse
import json
import math
from pathlib import Path
import random
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import Parts, collection, material, mesh, path, save_and_render, studio
from build_organelle_scene import build_mitochondrion


def ellipsoid(name, center, radii, mat, group, cut=False, ripple=0):
    vertices, faces = [], []
    rows, columns = 80, 160
    for j in range(rows + 1):
        phi = math.pi * j / rows
        for i in range(columns):
            theta = math.tau * i / columns
            r = 1 + ripple * math.sin(theta * 4) * math.sin(phi * 3)
            vertices.append((center[0] + radii[0] * math.sin(phi) * math.cos(theta) * r,
                             center[1] + radii[1] * math.sin(phi) * math.sin(theta) * r,
                             center[2] + radii[2] * math.cos(phi) * r))
    for j in range(rows):
        for i in range(columns):
            indices = (j * columns + i, j * columns + (i + 1) % columns,
                       (j + 1) * columns + (i + 1) % columns, (j + 1) * columns + i)
            centroid = np.mean([vertices[k] for k in indices], axis=0)
            if cut and centroid[1] < 0 and centroid[2] > 0:
                continue
            faces.append(indices)
    obj = mesh(name, vertices, faces, [mat], group)
    return obj


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--style-profile', type=Path,
                   default=Path(__file__).resolve().parent.parent/'assets/style-organelle.json')
    p.add_argument('--engine', choices=['cycles', 'eevee'], default='cycles')
    p.add_argument('--compositor', choices=['clean', 'reference', 'monochrome'])
    p.add_argument('--width', type=int, default=2000)
    p.add_argument('--samples', type=int, default=64)
    p.add_argument('--fps', type=int, default=24)
    p.add_argument('--frames', type=int, default=96)
    p.add_argument('--frame', type=int, default=1)
    p.add_argument('--render', action='store_true')
    p.add_argument('--animate', action='store_true', help='Reserved for project-specific cell motion')
    args = p.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    if args.animate:
        raise SystemExit('This overview is static. Use a process-specific rig for animation.')
    if args.output.exists():
        raise SystemExit('Choose a fresh output directory.')
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile = json.loads(args.style_profile.read_text())
    rng = random.Random(47)
    colors = {'cell': '9BBAB5', 'cell_inner': 'BFD0BE', 'envelope': 'B7B1CB',
              'chromatin': '807A9E', 'nucleolus': 'A997BD', 'golgi': 'C6AD82',
              'ribosome': 'D9CBB0', 'lysosome': 'BD8D73', 'peroxisome': 'A4B58C',
              'microtubule': '95B1BB', 'centrosome': 'B6A488', 'actin': '9DB7AA'}
    colors.update(profile.get('cell_palette', {}))
    mats = {k: material(k, v, .51, .004) for k, v in colors.items()}
    cell = collection('01 Cell membrane — deliberate quarter cutaway')
    outer = ellipsoid('Plasma membrane outer face', (0, 0, 0), (7.25, 5.25, 3.5), mats['cell'], cell, cut=True)
    shell = outer.modifiers.new('Visible cut edge; membrane thickness exaggerated', 'SOLIDIFY')
    shell.thickness = .045
    shell.offset = -1
    er = collection('02 Connected nuclear-envelope and ER lumen boundary')
    data = np.load(args.asset_dir/'er_ne.npz', allow_pickle=False)
    membrane = mesh('Shared ER and nuclear-envelope membrane boundary', data['vertices'].tolist(),
                    data['faces'].tolist(), [mats['envelope']], er)
    smooth = membrane.modifiers.new('Sub-voxel contour smoothing', 'SMOOTH')
    smooth.iterations, smooth.factor = 2, .35
    manifest = json.loads((args.asset_dir/'manifest.json').read_text())
    nucleus = collection('03 Nucleoplasm — chromatin and unbounded nucleolus')
    center, radii = np.array(manifest['nuclear_center']), np.array(manifest['nuclear_radii'])
    ellipsoid('Nucleolus — no lipid boundary', center + (.0, -.25, .3), (.78, .62, .65),
              mats['nucleolus'], nucleus, ripple=.08)
    for strand in range(25):
        q = np.array([rng.uniform(-.6,.6) for _ in range(3)])
        direction = np.array([rng.uniform(-1,1) for _ in range(3)])
        pts = []
        for _ in range(260):
            direction = .91*direction + np.array([rng.uniform(-.35,.35) for _ in range(3)])
            direction /= np.linalg.norm(direction)
            candidate = q + .025*direction
            if np.linalg.norm(candidate) > .84:
                direction -= 2*np.dot(direction,q)*q/np.dot(q,q)
                candidate = q + .025*direction
            # Keep the coarse chromatin outside the unbounded nucleolus.
            if np.linalg.norm(candidate*np.array([2.45,2,1.65])-(0,-.25,.3)) < .8:
                direction = -direction
                candidate = q + .025*direction
            q = candidate
            pts.append(tuple(center + radii*q))
        path('Interphase chromatin path — coarse representation', pts, .024, mats['chromatin'], nucleus)
    ribosomes = Parts()
    # Exterior leaflet only: selected cisternal patches, away from the nucleus.
    for _ in range(260):
        x, y = rng.uniform(-1.8, 2.4), rng.uniform(1.7, 4.1)
        if ((x-.15)/2.5)**2 + ((y-2.8)/1.55)**2 > .86:
            continue
        z = 1.16 + .09 * (x-.15)**2 + .17
        pos = np.array((x, y, z))
        if np.linalg.norm((pos-center)/radii) < 1.1:
            continue
        ribosomes.ball(pos, (.055, .045, .04), seg=8, rings=6)
        ribosomes.ball(pos + (0, 0, .04), (.04, .035, .028), seg=8, rings=6)
    ribosomes.object('Ribosome motifs on ER cytosolic side; exaggerated', [mats['ribosome']], er)
    golgi = collection('04 Polarized stack of separate Golgi cisternae')
    for i in range(6):
        obj = ellipsoid(f'Golgi cisterna {i+1}', (2.65, -.55, -.6 + i*.32),
                         (.72 + i*.045, 1.6 - i*.08, .10), mats['golgi'], golgi)
        for vertex in obj.data.vertices:
            vertex.co.z += .20 * (vertex.co.y + .55)**2
            vertex.co.x += .24 * (vertex.co.y + .55)**2 + .035*math.sin(vertex.co.y*5+i)
    # Reuse the same detailed model at a scene-appropriate level of magnification.
    groups, _ = build_mitochondrion(profile)
    template = bpy.data.collections.new('Mitochondrion source geometry — collection instances')
    for group in groups.values():
        bpy.context.scene.collection.children.unlink(group)
        template.children.link(group)
    mito_group = collection('05 Mitochondria — separate illustrative local cutaways')
    for i, (pos, scale, angle) in enumerate([
        ((-3.6,-2.25,-.9), .27, -.12), ((3.9,2.0,-.5), .25, 1.1),
        ((2.6,-2.8,-1.1), .26, -.45)]):
        obj = bpy.data.objects.new(f'Mitochondrion {i+1}', None)
        obj.instance_type, obj.instance_collection = 'COLLECTION', template
        mito_group.objects.link(obj)
        obj.location, obj.scale = pos, (scale,)*3
        obj.rotation_euler.z = angle
    vesicles = collection('06 Endomembrane vesicles and single-membrane organelles')
    for i, (pos, radius, color) in enumerate([
        ((4.3,-1.7,-.4),.48,'lysosome'), ((-4.7,.1,-.6),.38,'peroxisome'),
        ((1.45,-3.5,-.7),.24,'golgi'), ((3.7,-.45,1.2),.22,'golgi'),
        ((4.1,.3,.7),.25,'golgi')]):
        ellipsoid(f'{color} {i+1}', pos, (radius,)*3, mats[color], vesicles, ripple=.025)
    cytoskeleton = collection('07 Selected microtubules, cortical actin and centrosome')
    for i in range(16):
        theta = -1.35 + i * 2.7 / 15
        end = np.array((6.6*math.cos(theta), 4.5*math.sin(theta), -.9 + .4*math.cos(theta)))
        start = np.array((1.8, 1.1, .65))
        pts = []
        for t in np.linspace(0, 1, 80):
            pnt = start*(1-t) + end*t
            pnt[2] -= .6 * math.sin(math.pi*t)
            pts.append(tuple(pnt))
        path('Microtubule — selected radial tracks', pts, .017, mats['microtubule'], cytoskeleton)
    for band in range(9):
        pts = []
        phi = 1.67 + .08*band
        for theta in np.linspace(math.pi, math.tau, 200):
            pts.append((7.03*math.sin(phi)*math.cos(theta), 5.04*math.sin(phi)*math.sin(theta),
                        3.28*math.cos(phi) + .035*math.sin(theta*14+band)))
        path('Cortical actin bundle — illustrative subset', pts, .018, mats['actin'], cytoskeleton)
    centrioles = Parts()
    for orientation in (0,1):
        for k in range(9):
            angle = k * math.tau/9
            for offset in (-.035,0,.035):
                pnt = np.array((1.8+.13*math.cos(angle)+offset,
                                1.1+.13*math.sin(angle), .5))
                direction = np.array((.42,0,0)) if orientation else np.array((0,0,.42))
                if orientation: pnt += (.2,.26,.1)
                centrioles.tube([pnt, pnt+direction], .018, sides=7)
    centrioles.object('Two centrioles — nine triplet motifs each', [mats['centrosome']], cytoskeleton)
    studio(profile,args,{'position':(10,-19,16),'target':(0,0,-.25),'scale':18.5,'light_scale':.52})
    save_and_render(args,{'scene':'Animal secretory-cell overview, schematic',
                          'coordinate_units':'1 Blender unit = 1 um; small structures exaggerated',
                          'er_ne':manifest,
                          'limitations':['Not a reconstruction of a named cell, cell line or species.',
                                         'Organelle sizes, counts, positions and cell outline are authored teaching choices.',
                                         'Lipid/protein motifs on organelles are enlarged; whole-cell and molecular scenes use separate scales.',
                                         'Nuclear pores, many organelles and most cytoskeletal filaments are omitted from this overview.',
                                         'Each mitochondrion has its own explanatory cutaway; these are not one physical section plane.']})


if __name__ == '__main__':
    main()
