"""Rebuild a bounded Allen-cell numerical placement case; do not render a movie.

This adapter requires a qualified isotropic ZCYX segmentation and metadata row.
The generic volume_distribution API also accepts explicitly calibrated anisotropy.
"""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import Request,urlopen

import numpy as np
import tifffile

from volume_distribution import sample,validate


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_file(root, spec, fetch):
    if Path(spec['filename']).name!=spec['filename']:raise ValueError('Use a filename inside the input directory')
    path=root/spec['filename']
    if not path.exists():
        if not fetch:raise ValueError('Missing source; supply the cache or explicitly use --fetch: '+str(path))
        headers={'Range':spec['range']} if 'range' in spec else {}
        with urlopen(Request(spec['url'],headers=headers),timeout=30) as response:
            if 'range' in spec and response.status!=206:raise ValueError('Server did not honor the bounded metadata Range request')
            data=response.read(spec['maximum_bytes']+1)
        if len(data)>spec['maximum_bytes'] or hashlib.sha256(data).hexdigest()!=spec['sha256']:
            raise ValueError('Source size or hash changed; requalify rather than overriding the fingerprint')
        root.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    if digest(path)!=spec['sha256']:raise ValueError('Cached source hash differs: '+str(path))
    return path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'assets/volume-cell-230745.json')
    parser.add_argument('--fetch',action='store_true',help='Fetch only missing, size-limited, hash-pinned public inputs')
    args=parser.parse_args()
    if args.output.exists():parser.error('Use a fresh output directory')
    config=json.loads(args.config.read_text())
    if config['version']!=1:parser.error('Unsupported case configuration')
    image=source_file(args.inputs,config['image'],args.fetch);metadata=source_file(args.inputs,config['metadata'],args.fetch)
    with metadata.open() as f:rows=[r for r in csv.DictReader(f) if r['CellId']==config['cell_id']]
    if len(rows)!=1:parser.error('Expected one matching cell in the qualified metadata excerpt')
    row=rows[0]
    if (row['cell_stage']!=config['expected_cell_stage'] or row['structure_name']!=config['expected_structure'] or
        not config['image']['url'].endswith('/'+row['crop_seg'])):parser.error('Source path, cell state or structure does not match the qualified row')
    names=ast.literal_eval(row['name_dict'])['crop_seg'];spacing=np.asarray(ast.literal_eval(row['scale_micron']),float)
    if spacing.shape!=(3,) or not np.isfinite(spacing).all() or np.min(spacing)<=0 or not np.all(spacing==spacing[0]):
        parser.error('This adapter only qualifies isotropic metadata; establish axis order explicitly for other inputs')
    with tifffile.TiffFile(image) as f:
        if f.series[0].axes!='ZCYX':parser.error('Expected explicitly ordered ZCYX source axes')
        data=f.asarray()
    if data.dtype!=np.uint8 or not set(np.unique(data)).issubset({0,255}) or data.shape[1]!=len(names):
        parser.error('Expected declared binary 0/255 segmentation channels')
    masks={role:data[:,names.index(name)]!=0 for role,name in
           [('cell','membrane_segmentation'),('dna','dna_segmentation'),('structure','struct_segmentation')]}
    available=masks['cell']&~masks['dna']&~masks['structure'];origin=np.zeros(3)
    points=np.asarray(config['probe_vertices_xyz_um'],float)
    layout=sample(available,spacing,origin,points,config['count'],seed=config['seed'],
                  boundary_gap=config['boundary_gap_um'],pair_gap=config['pair_gap_um'])
    audit=validate(available,spacing,origin,points,layout['transforms'],
                   boundary_gap=config['boundary_gap_um'],pair_gap=config['pair_gap_um'])
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'source-masks.npz',**masks,available=available,spacing_xyz_um=spacing,origin_xyz_um=origin)
    np.savez_compressed(args.output/'layout.npz',body_vertices_xyz_um=points,**layout)
    record={'kind':'voxel-volume-placement-case','passed':audit['passed'],'config':config,
            'source_sha256':{'image':digest(image),'metadata':digest(metadata)},'selected_metadata':row,
            'data_sha256':{n:digest(args.output/n) for n in ['source-masks.npz','layout.npz']},
            'units':'micrometer','coordinate_frame':'ZYX arrays; XYZ positions; local crop voxel0 center at origin0. No global FOV registration.',
            'mask_policy':'cell AND NOT dna AND NOT structure; original source masks retained without repair',
            'mask_voxels':{k:int(v.sum()) for k,v in {**masks,'available':available}.items()},
            'source_conflicts':{'structure_in_dna_voxels':int(np.count_nonzero(masks['structure']&masks['dna'])),
                               'dna_outside_cell_voxels':int(np.count_nonzero(masks['dna']&~masks['cell'])),
                               'structure_outside_cell_voxels':int(np.count_nonzero(masks['structure']&~masks['cell']))},
            'audit':audit,'attempts':layout['attempts'],
            'limits':'Technical whole-body placement probe; not source segmentation validation, complete cytoplasmic occupancy, measured abundance, motion or a saved Blender scene.'}
    (args.output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['passed','mask_voxels','source_conflicts','audit','attempts']}))


if __name__=='__main__':main()
