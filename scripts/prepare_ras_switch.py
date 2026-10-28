"""Prepare H-Ras G-domain/SOS/p120GAP states with explicit missing-atom inference.

This is a source-guided teaching model, not a molecular trajectory. The native
membrane anchor and upstream receptor pathway are outside the selected constructs.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np

from prepare_molecular_asset import atom_table, xyz, surface
from prepare_filament_assets import rigid_fit
from prepare_phosphotransfer import residue_atoms, ccd_file, peptide_bonds, repair_inferred_geometry
from interaction_motion import plan


def fit_points(first, second):
    x, y = np.asarray(first), np.asarray(second)
    u, _, vt = np.linalg.svd((x-x.mean(0)).T @ (y-y.mean(0)))
    r = vt.T @ np.diag([1, 1, np.linalg.det(vt.T@u.T)]) @ u.T
    return r, y.mean(0)-x.mean(0)@r.T


def cif_string(value):
    value=str(value)
    return value[1:-1] if len(value)>1 and value[0]==value[-1] and value[0] in "\"'" else value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    tables={};models={};sources={}
    for code in ['4Q21','1BKD','1QRA','1WQ1']:
        path=a.sources/(code+'.cif');meta=json.loads(path.with_suffix('.cif.source.json').read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']:raise ValueError('Source hash mismatch: '+code)
        if meta['id'].upper()!=code:raise ValueError('Source identity mismatch')
        sources[code]=meta;tables[code]=atom_table(path);models[code]=gemmi.read_structure(str(path))[0]
    reference=[r for r in tables['4Q21']['A'] if r[3]<=166]
    ids=[(r[3],r[2]) for r in reference]
    core=lambda rows:[r for r in rows if r[3] not in range(25,42) and r[3] not in range(58,77) and r[3]<=166]
    origin=xyz(tables['1WQ1']['R']).mean(0)
    transforms={};observed={};ligands={};maps={}
    for code in sources:
        chain='A' if code in ['4Q21','1QRA'] else 'R'
        matrix,rmsd=rigid_fit(core(tables[code][chain]),core(tables['1WQ1']['R']))
        transforms[code]={'matrix_angstrom':matrix.tolist(),'core_ca_rmsd_angstrom':rmsd}
        def tr(points):return (np.asarray(points)@matrix[:3,:3].T+matrix[:3,3]-origin)*.1
        observed[code]={ (r[3],r[2]):tr(r[0]) for r in tables[code][chain] if r[3]<=166 }
        ligands[code]={}
        for c in models[code]:
            for residue in c:
                if residue.name not in ['GDP','GTP','MG','AF3']:continue
                ligands[code][residue.name]={n:(tr(v[0]),v[1]) for n,v in residue_atoms(residue).items()}
        if code=='1BKD':sos=tr(xyz(tables[code]['S']))
        if code=='1WQ1':gap=tr(xyz(tables[code]['G']))
    states={};canonical=np.array([observed['4Q21'][key] for key in ids])
    for code,mapping in observed.items():
        points=[];missing=[]
        for key in ids:
            if key in mapping:points.append(mapping[key]);continue
            anchors=[(key[0],n) for n in ['N','CA','C','CB'] if (key[0],n) in mapping and (key[0],n) in observed['4Q21']]
            if len(anchors)<3:raise ValueError('Insufficient residue anchors for '+str(key))
            r,t=fit_points([observed['4Q21'][v] for v in anchors],[mapping[v] for v in anchors])
            points.append(observed['4Q21'][key]@r.T+t);missing.append(list(key))
        state=np.array(points)
        state,repair=repair_inferred_geometry(canonical,state,reference,[i for i,key in enumerate(ids) if key in mapping])
        states[code]=state;maps[code]={'observed_atoms':len(ids)-len(missing),'inferred_atoms':missing,'geometry_repair':repair}
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'ras.npz',atoms=states['4Q21'],gef=states['1BKD'],gtp=states['1QRA'],gap=states['1WQ1'],
        bonds=peptide_bonds(reference,canonical),residue_numbers=np.array([k[0] for k in ids]),atom_names=np.array([k[1] for k in ids]))
    for name,points,rows in [('sos',sos,tables['1BKD']['S']),('gap',gap,tables['1WQ1']['G'])]:
        vertices,faces=surface(points*10,.75,1.1,.55)
        np.savez_compressed(a.output/(name+'.npz'),atoms=points,vertices=vertices*.1,faces=faces,
            atom_names=np.array([r[2] for r in rows]),residue_numbers=np.array([r[3] for r in rows]),bonds=peptide_bonds(rows,points))
    chemical={};ccd_records={}
    for name,code,ccd_name in [('old_gdp','4Q21','GDP'),('new_gtp','1QRA','GTP'),('gap_gdp','1WQ1','GDP')]:
        ccd,record=ccd_file(a.sources,ccd_name);ccd_records[ccd_name]=record
        ligand=ligands[code][ccd_name];names=list(ligand);lookup={n:i for i,n in enumerate(names)}
        bonds=np.array([[lookup[cif_string(r[0])],lookup[cif_string(r[1])]] for r in ccd.find('_chem_comp_bond.',['atom_id_1','atom_id_2'])
            if cif_string(r[0]) in lookup and cif_string(r[1]) in lookup],dtype=np.int32)
        expected={cif_string(r[0]) for r in ccd.find('_chem_comp_atom.',['atom_id','type_symbol']) if str(r[1]) not in ['H','D']}
        if set(names)!=expected:raise ValueError('Incomplete heavy-atom nucleotide: '+name)
        reached={0}
        for _ in names:
            reached.update(int(j) for i,j in bonds if int(i) in reached)
            reached.update(int(i) for i,j in bonds if int(j) in reached)
        if len(reached)!=len(names):raise ValueError('Disconnected nucleotide graph; preserve atom-name primes in CIF strings')
        chemical[name]={'atoms':np.array([ligand[n][0] for n in names]),'atom_names':np.array(names),
            'elements':np.array([ligand[n][1] for n in names]),'bonds':bonds,'magnesium':np.array([v[0] for v in ligands[code]['MG'].values()])}
        np.savez_compressed(a.output/(name+'.npz'),**chemical[name])
    # Reusable source docking paths are independent of the coupled reaction.
    # A failed path remains diagnostic; it does not authorize a scene build.
    paths={}
    for name,points,fixed in [('sos',sos,states['1BKD']),('gap',gap,states['1WQ1'])]:
        _,report=plan(points,fixed,frames=240,distance=10)
        paths[name]=report
        print(name,report['collision_screen_passed'],report['max_new_or_worsened_overlap_atoms'],flush=True)
    manifest={'kind':'ras-gef-gap-switch','units':'nm','sources':sources,'ccd':ccd_records,
        'selected_constructs':{'Ras':'H-Ras G-domain residues 1-166; 4Q21 residues 167-169 and 1WQ1 terminal OXT excluded to define a common construct.',
            'SOS':'1BKD author S, modeled residues within 568-1044; unresolved loops remain absent.',
            'GAP':'1WQ1 author G, modeled residues 718-1037 of p120GAP.'},
        'transforms':transforms,'ras_mapping':maps,'paths':paths,
        'structure_limits':['1BKD is nucleotide-free. A nucleotide must not be parked in its phosphate/Mg-displacing interface.',
            '1QRA contains true GTP. Missing switch-II side-chain positions are inferred, not observed.',
            '1WQ1 contains GDP and AlF3, a transition-state analogue. AlF3 is not converted into a phosphate atom.',
            'Atomic-center screening is a geometry diagnostic, not energetics or kinetics.'],
        'mechanism_sources':['https://doi.org/10.1038/28548','https://doi.org/10.1016/S0969-2126(00)80021-0',
            'https://doi.org/10.1126/science.277.5324.333','https://doi.org/10.1073/pnas.0712095105'],
        'timing_scope':'Qualitative teaching sequence only; 2008 FTIR work includes neurofibromin and 260 K experiments, not a measured p120GAP animation trajectory.'}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'output':str(a.output),'ras_atoms':len(ids),'mapping':maps},ensure_ascii=False))


if __name__=='__main__':main()
