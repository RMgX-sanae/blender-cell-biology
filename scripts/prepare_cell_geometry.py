"""Connected ER lumen/nuclear-envelope model for a schematic animal cell.

Runs outside Blender, with numpy/scipy/scikit-image. No microscopy reconstruction
is claimed. A shared implicit volume makes lumen continuity inspectable.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label
from skimage.measure import marching_cubes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--spacing', type=float, default=.045)
    args = p.parse_args()
    if args.output.exists():
        raise SystemExit('Choose a fresh output directory.')
    spacing = args.spacing
    origin = np.array([-4.9, -3.5, -2.2])
    stop = np.array([3.3, 4.9, 2.7])
    axes = [np.arange(a, b + spacing, spacing, dtype=np.float32) for a, b in zip(origin, stop)]
    x, y, z = np.meshgrid(*axes, indexing='ij', sparse=True)
    center = np.array([-1.5, .65, .2])
    radii = np.array([2.45, 2.0, 1.65])
    nuclear_distance = (np.sqrt(((x-center[0])/radii[0])**2 +
                                ((y-center[1])/radii[1])**2 +
                                ((z-center[2])/radii[2])**2) - 1) * 1.8
    field = np.abs(nuclear_distance) - .09
    # Flattened, slightly curved cisternae sit behind and beside the nucleus.
    cisternae = []
    for i in range(5):
        cx, cy, cz = .15, 2.8, -.8 + .49 * i
        bend = .09 * (x - cx) ** 2
        d = (np.sqrt(((x-cx)/2.5)**2 + ((y-cy)/1.55)**2 +
                     ((z-cz-bend)/.17)**2) - 1) * .17
        field = np.minimum(field, d)
        cisternae.append([cx, cy, cz])
    # Two real branch lumens connect every sheet to the nuclear-envelope space.
    for xx in (-.65, .8):
        # Curved tube axis sampled densely; union is resolved on the common grid.
        for t in np.linspace(0, 1, 70):
            yy = 2.0 + 1.0 * t
            zz = -.7 + 2.1 * t
            distance = np.sqrt((x-xx)**2 + (y-yy)**2 + (z-zz)**2) - .16
            field = np.minimum(field, distance)
    _, count = label(field < 0)
    if count != 1:
        raise RuntimeError(f'ER/NE lumen must be connected; found {count} components.')
    # Deliberate front upper cut exposes the nuclear interior.
    # A face-removed cut, not a cap across the nucleus or an invented membrane.
    vertices, faces, _, _ = marching_cubes(field, level=0, spacing=(spacing,)*3,
                                          gradient_direction='ascent')
    vertices += origin
    centroids = vertices[faces].mean(1)
    keep = ~((centroids[:,1] < center[1]) & (centroids[:,2] > center[2]+.4))
    faces = faces[keep]
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'er_ne.npz', vertices=vertices.astype(np.float32), faces=faces)
    manifest = {'scene': 'Schematic connected ER and nuclear-envelope lumen',
                'units': 'um', 'grid_spacing_um': spacing, 'lumen_components_before_cut': count,
                'nuclear_center': center.tolist(), 'nuclear_radii': radii.tolist(),
                'cisternae_centers': cisternae,
                'limitations': ['Authored teaching geometry, not segmented microscopy.',
                               'Nuclear-envelope gap and ER luminal thickness exaggerated for visibility.',
                               'Nuclear pores omitted in this overview; use a dedicated pore close-up.',
                               'Front upper cutaway intentionally removes membrane faces; it is not a biological opening.']}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'vertices':len(vertices),'triangles':len(faces), 'lumen_components_before_cut':count}))


if __name__ == '__main__':
    main()
