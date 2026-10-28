"""Prepare coarse public-structure envelopes for a schematic cytosolic backdrop."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from prepare_molecular_asset import atom_table,xyz,surface
from prepare_filament_assets import rigid_fit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['ubiquitin','ribosome','actin']:p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh environment asset directory')
    sources={};prepared={}
    for name,code in [('ubiquitin','1UBQ'),('ribosome','9SPI'),('actin','8F8R')]:
        path=getattr(a,name);record=json.loads(path.with_suffix(path.suffix+'.source.json').read_text())
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=record['sha256'] or record['id'].upper()!=code:raise ValueError('Source identity/hash mismatch: '+name)
        table=atom_table(path)
        rows=table['D'] if name=='actin' else [row for chain in table.values() for row in chain]
        points=xyz(rows);origin=points.mean(0)
        spacing,sigma=(2.5,2.2) if name=='ribosome' else (1.2,1.5)
        vertices,faces=surface(points-origin,spacing,sigma,.6)
        data={'vertices':vertices*.1,'faces':faces}
        if name=='actin':
            step,rmsd=rigid_fit(table['D'],table['E'])
            nm=np.eye(4);nm[:3,:3]*=.1;nm[:3,3]=-origin*.1
            data['repeat_nm']=nm@step@np.linalg.inv(nm)
            data['repeat_fit_rmsd_angstrom']=np.array(rmsd)
        prepared[name]=data
        sources[name]={'id':code,'sha256':digest,'entry_url':record['entry_url'],
            'selected_chains':['D'] if name=='actin' else list(table),
            'spacing_angstrom':spacing,'sigma_angstrom':sigma,'atoms':len(rows)}
    a.output.mkdir(parents=True)
    for name,data in prepared.items():np.savez_compressed(a.output/(name+'.npz'),**data)
    manifest={'kind':'schematic-cytosolic-environment','units':'nm','sources':sources,
        'limits':['Public experimental structures guide the shapes, not the spatial concentration or arrangement.',
            'Human ubiquitin and a human 80S ribosome are background examples. The 8F8R rabbit actin shape is an explicitly selected mammalian filament proxy.',
            'Only resolved polymer atoms contribute to the coarse display envelopes; these are not EM maps.',
            'The filament repeats the intact internal D subunit using its fitted D-to-E transform; it is not a newly observed filament.']}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'output':str(a.output),'assets':list(prepared)}))


if __name__=='__main__':main()
