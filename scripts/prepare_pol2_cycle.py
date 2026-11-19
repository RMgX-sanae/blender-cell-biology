"""Prepare the source envelopes and common viewing frames for a mammalian Pol II teaching route."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import numpy as np

from prepare_binding_candidates import checked_source, rigid_fit_points, sha, transform
from prepare_molecular_asset import atom_table, surface, xyz


CONFIG = Path(__file__).resolve().parents[1] / 'assets/pol2-transcription-cycle.json'


def source(directory, spec):
    path = checked_source(directory, spec)
    structure = gemmi.read_structure(str(path))
    assembly = next(a for a in structure.assemblies if a.name == '1')
    if len(assembly.generators) != 1:
        raise ValueError('This recipe requires the reviewed single-generator assembly')
    gen = assembly.generators[0]
    if len(gen.operators) != 1 or not gen.operators[0].transform.is_identity():
        raise ValueError('The deposited assembly no longer matches the identity copy')
    labels = set(gen.subchains)
    companions = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.entity_type == gemmi.EntityType.Polymer and residue.subchain not in labels:
                raise ValueError('A polymer lies outside the reviewed assembly')
            if residue.subchain in labels and residue.name in {'ZN','MG','SO4'}:
                for atom in residue:
                    if atom.occ > 0 and atom.altloc in ('\x00','A'):
                        companions.setdefault(chain.name, []).append((list(atom.pos), residue.name))
    return atom_table(path, include_insertion_codes=True), companions


def fit_core(moving, fixed):
    def sites(rows):
        return {(r[3],r[5],r[4]): r[0] for r in rows if r[2] == 'CA'}
    a, b = sites(moving), sites(fixed)
    ids = sorted(a.keys() & b.keys())
    if len(ids) < 900:
        raise ValueError('Insufficient corresponding RPB1 residues for a cross-specimen view')
    pa, pb = np.array([a[i] for i in ids]), np.array([b[i] for i in ids])
    first, _ = rigid_fit_points(pa,pb)
    kept = np.linalg.norm(transform(pa,first)-pb,axis=1) < 4.
    if kept.sum() < .75*len(ids):
        raise ValueError('The conserved RPB1 core does not support this registration')
    matrix, rmsd = rigid_fit_points(pa[kept],pb[kept])
    return matrix, {'corresponding':len(ids), 'retained':int(kept.sum()), 'rmsd_angstrom':rmsd}


def backbone(rows):
    by_residue = {}
    for r in rows:
        by_residue.setdefault((r[3],r[5]),{})[r[2]] = r[0]
    sites, bases, ids = [], [], []
    for (number,_), atoms in sorted(by_residue.items()):
        if "C4'" not in atoms:
            continue
        ring = [v for k,v in atoms.items() if k in {'N1','C2','N3','C4','C5','C6'}]
        if not ring:
            raise ValueError('Nucleic residue has no resolved base landmark')
        sites.append(atoms["C4'"]); bases.append(np.mean(ring,axis=0)); ids.append(number)
    return np.array(sites), np.array(bases), np.array(ids)


def view_basis(right, toward):
    right = right/np.linalg.norm(right)
    towards_viewer = toward - right*np.dot(right,toward)
    towards_viewer /= np.linalg.norm(towards_viewer)
    depth = -towards_viewer
    return np.column_stack((right,depth,np.cross(right,depth)))


def prepare(directory, output):
    if output.exists():
        raise ValueError('Choose a new preparation directory')
    recipe = json.loads(CONFIG.read_text())
    data = {role:source(directory,spec) for role,spec in recipe['coordinates'].items()}
    closed, initial, elong = (data[k][0] for k in ('closed','initial','elongation'))
    initial_to_closed, fit1 = fit_core(initial['A'],closed['A'])
    elong_to_initial, fit2 = fit_core(elong['A'],initial['A'])
    strand, bases, ids = backbone(elong['P'])
    selected = (ids>=15)&(ids<=20)
    short = strand[selected]
    if len(short)!=6: raise ValueError('Missing 5FLM RNA backbone guide')
    guide_raw = transform(transform(short,elong_to_initial),initial_to_closed)
    origin = guide_raw.mean(axis=0)
    dna, _, _ = backbone(closed['X'])
    basis = view_basis(dna[-1]-dna[0], origin-xyz(closed['A']).mean(axis=0))
    def display(points, matrix):
        return (transform(np.asarray(points),matrix)-origin) @ basis / 10
    guide = (guide_raw-origin) @ basis / 10
    arrays = {'short_rna_guide':guide,
              'short_rna_bases':display(transform(bases[selected],elong_to_initial),initial_to_closed)}
    meta = {'kind':'pol2-transcription-route-v1','config_sha256':sha(CONFIG),'fits':{'5IY8_to_5IY6':fit1,'5FLM_to_5IY8':fit2},'states':{}}
    protein_ids = 'ABCDEFGHIJKLMNOPQRSTUVW0123'
    for state,matrix in [('closed',np.eye(4)),('initial',initial_to_closed)]:
        rows,comp = data[state]
        chains,occluders = [],[]
        for chain in protein_ids:
            points = display(xyz(rows[chain]),matrix)
            key = state+'_'+chain
            arrays[key+'_v'], arrays[key+'_f'] = surface(points,.21,.255,.80)
            chains.append(chain)
            if any(np.any((np.linalg.norm(points[:,[0,2]]-v[[0,2]],axis=1)<.82)&(points[:,1]<v[1])) for v in guide):
                occluders.append(chain)
            if chain in comp:
                arrays[key+'_ions'] = display([a for a,_ in comp[chain]],matrix)
        for chain in 'XY':
            backbone_points,bases,ids = backbone(rows[chain])
            key = state+'_'+chain
            arrays[key+'_backbone'] = display(backbone_points,matrix)
            arrays[key+'_base'] = display(bases,matrix)
            arrays[key+'_ids'] = ids
        meta['states'][state]={'chains':chains, 'occluders_at_rna':occluders}
    # An independent, direction-preserving frame for the schematic downstream shot.
    template,_,ids = backbone(elong['T'])
    late_origin = short[-1]
    late_basis = view_basis(template[ids==3][0]-template[ids==37][0],short.mean(axis=0)-xyz(elong['A']).mean(axis=0))
    anchor = max((r for r in elong['A'] if r[2]=='CA'),key=lambda r:r[3])
    meta['elongation_frame']={'origin_angstrom':late_origin.tolist(),'basis_columns':late_basis.tolist(),
        'rpb1_last_ca':anchor[3],'ctd_attachment_nm':((np.array(anchor[0])-late_origin)@late_basis/10).tolist()}
    cleavage,comp = data['cleavage']
    sulfate = [p for p,label in comp['A'] if label=='SO4']
    zinc = [p for p,label in comp['A'] if label=='ZN']
    if len(sulfate)!=5 or len(zinc)!=4: raise ValueError('The CPSF73 companion inventory changed')
    c = np.mean(sulfate,axis=0)
    arrays['cpsf_v'],arrays['cpsf_f'] = surface((xyz(cleavage['A'])-c)/10,.15,.19,.7)
    arrays['cpsf_zn']=(np.array(zinc)-c)/10
    arrays['cpsf_so4']=(np.array(sulfate)-c)/10
    output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temp:
        dest=Path(temp)/'prepared'; dest.mkdir()
        np.savez_compressed(dest/'geometry.npz',**arrays)
        meta['geometry_sha256']=sha(dest/'geometry.npz')
        (dest/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
        dest.rename(output)
    print(json.dumps(meta['fits']))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    prepare(args.sources,args.output)
