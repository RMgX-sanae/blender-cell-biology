"""Fetch and prepare the Allen Visual Guide's public background display meshes.

Use for noncommercial research/education under the provider's original terms.
These display surfaces retain source-resolution limits; they are not atomic models.
Downloaded OBJ and derived NPZ files remain outside the skill repository.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

import numpy as np


def parse_obj(path):
    vertices,faces=[],[]
    for line in path.read_text().splitlines():
        fields=line.split()
        if not fields:continue
        if fields[0]=='v':vertices.append([float(x) for x in fields[1:4]])
        elif fields[0]=='f':
            indices=[int(x.split('/')[0]) for x in fields[1:]]
            indices=[x-1 if x>0 else len(vertices)+x for x in indices]
            faces.extend((indices[0],indices[i],indices[i+1]) for i in range(1,len(indices)-1))
    v=np.asarray(vertices,dtype=np.float32);f=np.asarray(faces,dtype=np.int32)
    if not len(v) or not len(f) or not np.isfinite(v).all():raise ValueError('Empty or invalid OBJ geometry: '+str(path))
    if f.min()<0 or f.max()>=len(v):raise ValueError('Invalid OBJ face index: '+str(path))
    return v,f


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source-dir',type=Path,help='Reuse previously downloaded named OBJ files after hash checks.')
    a=p.parse_args()
    if a.output.exists():p.error('Choose a fresh output directory.')
    catalog=json.loads((Path(__file__).resolve().parents[1]/'assets/allen-environment-sources.json').read_text())
    a.output.mkdir(parents=True);raw=a.output/'sources';raw.mkdir()
    records=[]
    for entry in catalog['assets']:
        name=entry['name'];cached=a.source_dir/(name+'.obj') if a.source_dir else None
        if cached and cached.is_file():data=cached.read_bytes()
        else:
            request=urllib.request.Request(entry['url'],headers={'User-Agent':'Blender-cell-biology-source-fetcher'})
            with urllib.request.urlopen(request,timeout=60) as response:data=response.read()
        digest=hashlib.sha256(data).hexdigest()
        if digest!=entry['sha256']:raise ValueError('Source changed; inspect before replacing recorded hash: '+name)
        path=raw/(name+'.obj');path.write_bytes(data)
        vertices,faces=parse_obj(path)
        np.savez_compressed(a.output/(name+'.npz'),vertices=vertices,faces=faces)
        records.append({**entry,'vertices':len(vertices),'triangles':len(faces),
                        'bounds':[vertices.min(0).tolist(),vertices.max(0).tolist()]})
    catalog.update(kind='allen-visual-guide-environment',assets=records,
                   coordinates='Original mesh coordinates retained; the scene adapter uses recorded rigid placement and uniform scale.')
    (a.output/'manifest.json').write_text(json.dumps(catalog,indent=2))
    print(json.dumps({'assets':len(records),'output':str(a.output)},indent=2))


if __name__=='__main__':main()
