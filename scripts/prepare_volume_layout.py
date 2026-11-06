"""Place one exported complete Blender body in a qualified voxel-region bundle."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from volume_distribution import sample,validate


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--region',type=Path,required=True,help='Directory from prepare_volume_case.py')
    p.add_argument('--body',type=Path,required=True,help='Directory from blender_volume.py export')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--count',type=int,required=True);p.add_argument('--seed',type=int,default=23)
    p.add_argument('--boundary-gap-um',type=float,default=0.);p.add_argument('--pair-gap-um',type=float,default=0.)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    region=json.loads((a.region/'manifest.json').read_text());body=json.loads((a.body/'manifest.json').read_text())
    if region.get('kind')!='voxel-volume-placement-case' or region.get('units')!='micrometer':
        p.error('Expected a qualified micrometer source-mask bundle')
    if body.get('kind')!='evaluated-volume-body' or body.get('units')!='micrometer':
        p.error('Expected a complete evaluated micrometer body')
    if body.get('source_geometry_reopened_and_compared') is not True:
        p.error('Re-export the body with the saved-source geometry check')
    if digest(a.region/'source-masks.npz')!=region['data_sha256']['source-masks.npz']:
        p.error('Region data differs from its source manifest')
    if digest(a.body/'body.npz')!=body['data_sha256']:p.error('Body geometry differs from its export manifest')
    with np.load(a.region/'source-masks.npz',allow_pickle=False) as data:
        mask=data['available'];spacing=data['spacing_xyz_um'];origin=data['origin_xyz_um']
    with np.load(a.body/'body.npz',allow_pickle=False) as data:points=data['vertices']
    layout=sample(mask,spacing,origin,points,a.count,seed=a.seed,
                  boundary_gap=a.boundary_gap_um,pair_gap=a.pair_gap_um)
    check=validate(mask,spacing,origin,points,layout['transforms'],
                   boundary_gap=a.boundary_gap_um,pair_gap=a.pair_gap_um)
    a.output.mkdir(parents=True)
    shutil.copyfile(a.region/'source-masks.npz',a.output/'source-masks.npz')
    shutil.copyfile(a.body/'body.npz',a.output/'body.npz')
    np.savez_compressed(a.output/'layout.npz',**layout)
    record={'kind':'complete-body-volume-layout','version':1,'units':'micrometer',
            'region':region,'body':body,'seed':a.seed,'count':a.count,
            'boundary_gap_um':a.boundary_gap_um,'pair_gap_um':a.pair_gap_um,
            'files':{name:digest(a.output/name) for name in ['source-masks.npz','body.npz','layout.npz']},
            'audit':check,'limits':'Static copies of the exported geometry in the declared voxel union; not measured abundance or dynamics.'}
    (a.output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':True,'count':a.count,'body_radius_um':layout['body_radius'],'attempts':layout['attempts']}))


if __name__=='__main__':main()
