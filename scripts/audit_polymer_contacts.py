"""Flag severe inter-unit heavy-atom overlaps from evaluated Blender pose tracks.

Outside Blender: NumPy / SciPy. Threshold is a display-path diagnostic, not
an energy model. Source resolution and known endpoint contacts still matter.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--poses', type=Path, required=True)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threshold-nm', type=float, default=.12)
    args = p.parse_args()
    if not 0 < args.threshold_nm <= .2:
        p.error('Use a positive diagnostic threshold <= 0.2 nm')
    tracks = np.load(args.poses, allow_pickle=False)
    atoms = [np.load(args.asset_dir/k/f'chain_{c}.npz', allow_pickle=False)['atoms']
             for k, c in zip(tracks['kinds'], tracks['chains'])]
    identity = np.concatenate([np.full(len(a), i) for i, a in enumerate(atoms)])
    clashes = []
    for index, matrices in enumerate(tracks['matrices']):
        coords = np.concatenate([a@m[:3,:3].T+m[:3,3] for a,m in zip(atoms,matrices)])
        pairs = cKDTree(coords).query_pairs(args.threshold_nm, output_type='ndarray')
        pairs = pairs[identity[pairs[:,0]] != identity[pairs[:,1]]]
        if len(pairs):
            distances = np.linalg.norm(coords[pairs[:,0]]-coords[pairs[:,1]], axis=1)
            worst = pairs[np.argmin(distances)]
            clashes.append({'frame': index+1, 'pairs': len(pairs), 'minimum_nm': float(distances.min()),
                            'units': [str(tracks['ids'][identity[i]]) for i in worst]})
    report = {'frames_checked': len(tracks['matrices']), 'threshold_nm': args.threshold_nm,
              'clashing_frames': clashes, 'passed': not clashes,
              'limits': 'Detects severe cross-unit heavy-atom overlaps only, not vdW energies, hydrogens, ligands, omitted context or surface-envelope overlap.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'frames_checked': report['frames_checked'], 'clashing_frames': len(clashes),
                      'first_clash': clashes[:1]}))
    if clashes:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
