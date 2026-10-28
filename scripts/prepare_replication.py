"""Prepare intact nucleotide-pair prototypes from 1BNA for a schematic fork."""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--structure',type=Path,default=Path(__file__).resolve().parent.parent/'assets/1BNA.pdb')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():p.error('Choose a fresh output directory.')
    if args.structure.stem.upper()!='1BNA':p.error('This residue map is specific to 1BNA.')
    structure=gemmi.read_structure(str(args.structure))
    residues={(c.name,r.seqid.num):r for c in structure[0] for r in c
              if r.entity_type==gemmi.EntityType.Polymer}
    def atom(chain,res,name):
        return np.array(list(next(a for a in residues[chain,res] if a.name==name).pos))*.1
    def pair_center(i):return (atom('A',i,"C1'")+atom('B',25-i,"C1'"))/2
    args.output.mkdir(parents=True)
    records=[]
    for base,i in [('C',3),('G',4),('A',5),('T',7)]:
        origin=pair_center(i)
        x=pair_center(i+1)-pair_center(i-1);x/=np.linalg.norm(x)
        y=atom('A',i,"C1'")-atom('B',25-i,"C1'")
        y-=x*np.dot(x,y);y/=np.linalg.norm(y)
        basis=np.column_stack((x,y,np.cross(x,y)))
        for chain,res in [('A',i),('B',25-i)]:
            atoms=[a for a in residues[chain,res] if not a.element.is_hydrogen and a.occ>0]
            points=(np.array([list(a.pos) for a in atoms])*.1-origin)@basis
            bonds=cKDTree(points).query_pairs(.185,output_type='ndarray').astype(np.int32)
            names=[a.name for a in atoms]
            if not {'P',"O3'","C4'"}.issubset(names):raise ValueError('Incomplete prototype backbone')
            np.savez_compressed(args.output/f'{base}_{chain}.npz',atoms=points.astype(np.float32),
                elements=np.array([a.element.name for a in atoms]),atom_names=np.array(names),bonds=bonds)
        records.append({'base_A':base,'source_residues':{'A':i,'B':25-i},
                        'origin_nm':origin.tolist(),'basis_columns':basis.tolist()})
    (args.output/'manifest.json').write_text(json.dumps({'source_file':args.structure.name,
        'source_sha256':hashlib.sha256(args.structure.read_bytes()).hexdigest(),
        'url':'https://www.rcsb.org/structure/1BNA','units':'nm','pairs':records,
        'scope':'Complete internal nucleotide prototypes. Repetition, helix pitch, curved fork and new chains are authored, not a measured replication fork.',
        'limitations':['Atomic geometry is reused across a constructed sequence; it is not sequence-specific structural prediction.',
                       'Inter-residue connectors in an animated curved fork are display guides, not force-field bonds.']},indent=2)+'\n')
    print('Prepared four complementary 1BNA nucleotide pairs.')


if __name__=='__main__':main()
