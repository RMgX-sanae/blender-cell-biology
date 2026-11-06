"""Reconstruct the qualified clathrin cage from source operators and shared assets."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree
from prepare_molecular_asset import atom_table,xyz


def rigid(matrix):
    m=np.asarray(matrix,float)
    if m.shape!=(4,4) or not np.isfinite(m).all() or not np.allclose(m[3],[0,0,0,1],atol=1e-8,rtol=0) or not np.allclose(m[:3,:3].T@m[:3,:3],np.eye(3),atol=2e-7,rtol=0) or not np.isclose(np.linalg.det(m[:3,:3]),1,atol=2e-7,rtol=0):
        raise ValueError('Assembly operator must be a finite proper rigid transform')
    return m


def compare_reference(table,reference,units):
    """Match by position, then independently compare all site identities and coordinates."""
    names=list(reference);tree=cKDTree([xyz(reference[c]).mean(0)*.1 for c in names]);used=set();errors=[]
    for unit in units:
        matrix=rigid(unit['matrix_nm']);unit['reference_chains']={}
        for chain in unit['chains']:
            points=xyz(table[chain])*.1@matrix[:3,:3].T+matrix[:3,3]
            distance,index=tree.query(points.mean(0));name=names[index];rows=reference[name]
            if distance>1e-4 or name in used or [r[1:] for r in rows]!=[r[1:] for r in table[chain]]:
                raise ValueError('Assembly chain identity, occupancy or independent reference mapping differs')
            error=float(abs(points-xyz(rows)*.1).max())
            if error>1e-4:raise ValueError('Assembly coordinates differ from independently downloaded reference')
            used.add(name);errors.append(error);unit['reference_chains'][chain]=name
    if used!=set(reference):raise ValueError('Assembly omits reference chains')
    return {'matched_chains':len(used),'maximum_coordinate_error_nm':max(errors),'tolerance_nm':1e-4}


def check_groups(table,groups):
    # Source-chain names alone did not identify intact triskelia in 1XI4.
    # 3IYV reassigns the fitted segments to whole modeled heavy chains. Its
    # low-resolution CA trace still contains long steps; retain and report them.
    trace={}
    for group in groups:
        for c in group['heavy']+group['light']:
            chain=table[c];steps=np.linalg.norm(np.diff(xyz(chain),axis=0),axis=1)
            if any(b[3]!=a[3]+1 for a,b in zip(chain,chain[1:])) or any(r[2]!='CA' for r in chain):
                raise ValueError('Require the reviewed consecutive CA identity sequence')
            if not np.isfinite(steps).all() or np.any(steps<=0) or np.any(steps>10):
                raise ValueError('Source chain has a spatial discontinuity; do not treat fitted sections as one triskelion')
            trace[c]={'maximum_CA_step_angstrom':float(steps.max()),
                      'long_steps':[{'residues':[chain[i][3],chain[i+1][3]],'distance_angstrom':float(steps[i])}
                                    for i in np.flatnonzero(steps>4.5)]}
    heavy=[c for g in groups for c in g['heavy']]
    centers=np.array([xyz([r for r in table[c] if r[3]>=1590]).mean(0) for c in heavy])*.1
    distances=np.linalg.norm(centers[:,None]-centers[None,:],axis=-1)
    trees={c:cKDTree(xyz(table[c])*.1) for c in heavy};links={}
    for group in groups:
        for c in group['heavy']:
            neighbors={heavy[i] for i in np.argsort(distances[heavy.index(c)])[:3]}
            if neighbors!=set(group['heavy']):raise ValueError('Heavy chains do not share the source terminal hub')
        for c in group['light']:
            scores={h:float(tree.query(xyz(table[c])*.1)[0].mean()) for h,tree in trees.items()}
            closest=min(scores,key=scores.get)
            if closest not in group['heavy']:raise ValueError('Modeled light chain belongs to another source hub')
            links[c]={'nearest_heavy':closest,'mean_CA_distance_nm':scores[closest]}
    return {'passed':True,'light_heavy_links':links,'trace_check':{'passed':True,'chains':trace,
            'maximum_permitted_step_angstrom':10.,'limits':'A source-specific coarse-trace gate against disconnected sections, not a normal peptide-bond or atomic-geometry certificate.'},
            'limits':'Source-specific trace and grouping consistency; not an affinity or universal grouping predictor.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--config',type=Path,default=Path(__file__).resolve().parents[1]/'assets/clathrin-cage.json')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    config=json.loads(a.config.read_text());source=a.sources/(config['source_id']+'.cif')
    reference=a.sources/(config['source_id']+'-assembly'+config['assembly_id']+'.cif')
    for file,key in [(source,'source_sha256'),(reference,'reference_sha256')]:
        if hashlib.sha256(file.read_bytes()).hexdigest()!=config[key]:raise ValueError('Requalify changed source: '+file.name)
    table=atom_table(source);ref=atom_table(reference)
    selected=[c for g in config['groups'] for c in g['heavy']+g['light']]
    if len(selected)!=len(set(selected)) or set(selected)!=set(table):raise ValueError('Group selection must retain each source chain once')
    for group in config['groups']:
        if len(group['heavy'])!=3 or len(group['light'])!=3:raise ValueError('Require three heavy and three modeled light chains per triskelion')
        for role,span in [('heavy',range(1,1631)),('light',range(95,165))]:
            for c in group[role]:
                if [r[3] for r in table[c]]!=list(span) or any(r[2]!='CA' for r in table[c]):raise ValueError('Unexpected source CA chain identity')
    grouping=check_groups(table,config['groups'])
    assembly=next(x for x in gemmi.read_structure(str(source)).assemblies if x.name==config['assembly_id'])
    if len(assembly.generators)!=1:raise ValueError('This qualified cage uses one assembly generator')
    units=[]
    for operator in assembly.generators[0].operators:
        matrix=np.eye(4);matrix[:3,:3]=np.array(operator.transform.mat);matrix[:3,3]=np.array(list(operator.transform.vec))*.1
        rigid(matrix)
        for group in config['groups']:
            hub=np.concatenate([xyz([r for r in table[c] if r[3]>=1590]) for c in group['heavy']]).mean(0)*.1
            units.append({'id':operator.name+'-'+group['id'],'group':group['id'],'operator':operator.name,
                          'chains':group['heavy']+group['light'],'matrix_nm':matrix.tolist(),'hub_nm':(hub@matrix[:3,:3].T+matrix[:3,3]).tolist()})
    comparison=compare_reference(table,ref,units)
    if len(units)!=36 or comparison['matched_chains']!=216:raise ValueError('Qualified cage inventory differs')
    # Reuse the existing project-selection and surface preparation entry point.
    a.output.mkdir(parents=True)
    selection={'id':config['id'],'source_id':config['source_id'],'source_sha256':config['source_sha256'],
       'chains':[{'chain':c} for c in selected],'selection_note':config['selection_note'],
       'orientation':{'origin_angstrom':[0,0,0],'basis_columns':np.eye(3).tolist(),'note':'Retain deposited point frame for assembly operators.','evidence_urls':config['evidence_urls']},'limitations':[config['limits']]}
    path=a.output/'source-selection.json';path.write_text(json.dumps(selection,indent=2))
    surface=config['surface'];script=Path(__file__).resolve().parent
    subprocess.run([sys.executable,str(script/'prepare_molecular_asset.py'),str(source),'--config',str(path),'--output',str(a.output/'prototypes'),
       '--spacing',str(surface['spacing_angstrom']),'--sigma',str(surface['sigma_angstrom']),'--level',str(surface['isovalue'])],check=True)
    files=[p for p in (a.output/'prototypes').iterdir() if p.suffix in {'.json','.npz'}]+[path]
    manifest={'version':1,'units':'nm','source_config':config,'instances':units,'independent_reference_check':comparison,'grouping_check':grouping,
       'data_sha256':{str(p.relative_to(a.output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
       'state':'Static source assembly; no growth trajectory or membrane attachment inferred',
       'display_geometry':'CA-derived Gaussian envelopes; not all-atom surfaces or experimental density'}
    (a.output/'assembly.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'units':len(units),'comparison':comparison}))


if __name__=='__main__':main()
