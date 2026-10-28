"""Source-derived whole-chain envelopes for a continuous V-ATPase teaching loop."""
import hashlib
import json
from pathlib import Path
import argparse
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_molecular_asset import atom_table, surface, xyz

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--structure',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
source=args.structure.resolve();out=args.output.resolve()
if source.stem.upper()!='6VQ6':parser.error('This selection is specific to 6VQ6.')
if out.exists():parser.error('Choose a fresh asset output directory.')
if out.is_relative_to(Path(__file__).resolve().parents[1]):parser.error('Prepared assets must be outside the skill directory.')
table = atom_table(source)
excluded = ['Q', 'R', 'S']  # SidK purification binder, not part of the native pump.
ring = ['b'] + list('ghijklmno')
rotor = ['H', 'L', 'd', 'c', 'p'] + ring
carriers = []
for chain in ring:
    residue = 98 if chain == 'b' else 139
    rows = [r for r in table[chain] if r[3] == residue and r[2] in ['OE1', 'OE2']]
    assert len(rows) == 2
    carriers.append(xyz(rows).mean(0))
carriers = np.array(carriers)
centroid = carriers.mean(0)
_, _, vh = np.linalg.svd(carriers-centroid)
z = vh[-1]
head = np.concatenate([xyz(table[c]) for c in 'ABCDEF']).mean(0)
if np.dot(z, head-centroid) < 0:
    z *= -1
arg = xyz([r for r in table['a'] if r[3] == 741 and r[2] == 'CZ'])[0]
x = arg-centroid
x -= z*np.dot(x,z)
x /= np.linalg.norm(x)
y = np.cross(z,x)
basis = np.array([x,y,z]).T
p = (carriers-centroid)@basis
circle = np.linalg.lstsq(np.c_[2*p[:,:2],np.ones(10)], (p[:,:2]**2).sum(1),rcond=None)[0]
origin = centroid+basis@np.r_[circle[:2],0.]
out.mkdir(parents=True,exist_ok=True)
info = {
    'source':'6VQ6', 'url':'https://www.rcsb.org/structure/6VQ6',
    'citation':'https://doi.org/10.1126/science.aaz2924',
    'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'organism':'Rattus norvegicus', 'source_state':'ADP- and SidK-bound composite rotational state 1',
    'excluded_chains':excluded, 'excluded_reason':'Three Legionella SidK purification/stabilization binders; no native H subunit modeled in source.',
    'nonpolymers':'Source ADP/Mg not displayed as active turnover substrate; incoming ATP and products are explicit graphical tokens.',
    'units':'nm', 'origin_angstrom':origin.tolist(), 'basis_columns':basis.tolist(),
    'membrane_orientation':'Inferred from c-ring acidic residue plane, positive Z toward V1/cytosol; membrane plane and lipids schematic.',
    'rotor_chains':rotor, 'c_ring_chains':ring,
    'carriers':[], 'chains':[],
    'motion':'Clockwise viewed from V1 toward Vo; authored rigid rotation, not a measured trajectory or molecular dynamics.',
    'surface':{'spacing_angstrom':.95,'sigma_angstrom':1.15,'isovalue':.55,'type':'Gaussian envelope, not cryo-EM density'},
}
for c, pos in zip(ring,carriers):
    info['carriers'].append({'chain':c,'residue':98 if c=='b' else 139,'site_nm':((pos-origin)@basis*.1).tolist()})
for chain, rows in table.items():
    if chain in excluded:
        continue
    coords = (xyz(rows)-origin)@basis
    v,f = surface(coords,.95,1.15,.55)
    np.savez_compressed(out/f'chain_{ord(chain):02x}.npz',vertices=v*.1,faces=f,atoms=coords*.1,
                        atom_names=np.array([r[2] for r in rows]),residue_numbers=np.array([r[3] for r in rows]))
    item={'chain':chain,'file':f'chain_{ord(chain):02x}.npz','atoms':len(rows),'vertices':len(v),'faces':len(f),
          'centroid_nm':(coords.mean(0)*.1).tolist(),'bounds_nm':(np.array([coords.min(0),coords.max(0)])*.1).tolist(),
          'role':'rotor' if chain in rotor else 'stator'}
    info['chains'].append(item)
    print(chain,len(v),flush=True)
(out/'manifest.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
