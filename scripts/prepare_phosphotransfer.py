"""Prepare a PKA/SP20 teaching reaction, with explicit analog-to-ATP reconstruction.

4IAC is an AMP-PCP complex, not an observed ATP reactant. 4IAD supplies the
ADP/phosphopeptide endpoint. Motion and chemical interpolation are illustrative.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial import cKDTree

from fetch_structure import download
from interaction_motion import plan
from prepare_filament_assets import rigid_fit
from prepare_molecular_asset import atom_table, surface, xyz


def residue_atoms(residue):
    return {a.name:(np.array(list(a.pos)),a.element.name) for a in residue
            if not a.element.is_hydrogen and a.occ>0 and a.altloc in ('\x00','A')}


def ccd_file(directory,code):
    path=directory/(code+'.ccd.cif');record=path.with_suffix('.cif.source.json')
    url='https://files.rcsb.org/ligands/download/'+code+'.cif'
    if not path.exists():
        data=download(url);block=gemmi.cif.read_string(data.decode()).sole_block()
        if block.name!=code:raise ValueError('Wrong CCD component')
        path.write_bytes(data)
        record.write_text(json.dumps({'component':code,'url':url,'sha256':hashlib.sha256(data).hexdigest(),
            'license':'wwPDB structural data CC0; cite the Chemical Component Dictionary'},indent=2))
    meta=json.loads(record.read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']:raise ValueError('CCD hash mismatch')
    return gemmi.cif.read_file(str(path)).sole_block(),meta


def peptide_bonds(rows,positions):
    radii={'C':.076,'N':.071,'O':.066,'S':.105,'P':.107}
    result=[]
    for i,a in enumerate(rows):
        for j in range(i):
            b=rows[j]
            same=a[3]==b[3]
            adjacent=(a[3]==b[3]+1 and a[2]=='N' and b[2]=='C') or (b[3]==a[3]+1 and b[2]=='N' and a[2]=='C')
            if (same or adjacent) and .07<np.linalg.norm(positions[i]-positions[j])<radii[a[1]]+radii[b[1]]+.035:
                result.append([j,i])
    return np.array(result,dtype=np.int32)


def displacement_field(points,anchors,delta):
    distances,indices=cKDTree(anchors).query(points,k=4)
    weights=1/np.maximum(distances,.01)**4;weights/=weights.sum(1,keepdims=True)
    return (delta[indices]*weights[:,:,None]).sum(1)


def repair_inferred_geometry(reference,target,rows,matched):
    """Keep reconstructed gaps connected without warping the surface mesh.

    Observed coordinates outside affected residues stay fixed. Bond and
    adjacent-angle distances regularize the reconstructed atoms, not energies.
    """
    bonds=peptide_bonds(rows,reference)
    original=np.linalg.norm(reference[bonds[:,0]]-reference[bonds[:,1]],axis=1)
    proposed=np.linalg.norm(target[bonds[:,0]]-target[bonds[:,1]],axis=1)
    bad=bonds[(proposed<.08)|(proposed>.22)|(np.abs(proposed/original-1)>.12)]
    if not len(bad):return target,{'repaired_atoms':0,'max_bond_nm':float(proposed.max())}
    affected={rows[i][3] for i in bad.ravel()}
    affected|={r+offset for r in list(affected) for offset in [-1,1]}
    mobile=np.array([i for i,r in enumerate(rows) if r[3] in affected]);lookup={int(i):j for j,i in enumerate(mobile)}
    neighbors=[set() for _ in rows];pairs={tuple(edge) for edge in bonds}
    for i,j in bonds:neighbors[i].add(j);neighbors[j].add(i)
    for connected in neighbors:
        for i in connected:
            for j in connected:
                if i<j:pairs.add((i,j))
    edges=np.array([edge for edge in sorted(pairs) if edge[0] in lookup or edge[1] in lookup])
    lengths=np.linalg.norm(reference[edges[:,0]]-reference[edges[:,1]],axis=1)
    known=set(matched);anchor=np.array([.5 if int(i) in known else .05 for i in mobile])[:,None]
    count=mobile.size*3;sparsity=lil_matrix((count+len(edges),count),dtype=int);sparsity.setdiag(1)
    for k,(i,j) in enumerate(edges):
        for atom in [i,j]:
            if atom in lookup:
                col=3*lookup[atom];sparsity[count+k,col:col+3]=1
    repaired=reference.copy();evaluations=0
    for fraction in np.linspace(.125,1,8):
        stage=reference*(1-fraction)+target*fraction
        def residual(flat):
            pos=stage.copy();pos[mobile]=flat.reshape(-1,3)
            distance=np.linalg.norm(pos[edges[:,0]]-pos[edges[:,1]],axis=1)
            return np.r_[((pos[mobile]-stage[mobile])*anchor).ravel(),(distance-lengths)*20]
        initial=repaired[mobile]+(target[mobile]-reference[mobile])*.125
        fit=least_squares(residual,initial.ravel(),jac_sparsity=sparsity.tocsr(),max_nfev=300,ftol=1e-5)
        repaired=stage.copy();repaired[mobile]=fit.x.reshape(-1,3);evaluations+=fit.nfev
    after=np.linalg.norm(repaired[bonds[:,0]]-repaired[bonds[:,1]],axis=1)
    relative_error=float(np.max(np.abs(after/original-1)))
    if not fit.success or after.min()<.08 or after.max()>.22 or relative_error>.15:
        raise ValueError(f'Inferred kinase gap reconstruction failed: converged={fit.success}, evaluations={evaluations}, bond range={after.min():.6f}..{after.max():.6f} nm, maximum fractional bond error={relative_error:.4f}, solver={fit.message}')
    return repaired,{'repaired_atoms':len(mobile),'affected_residues':sorted(affected),
        'before_max_bond_nm':float(proposed.max()),'max_bond_nm':float(after.max()),'min_bond_nm':float(after.min()),
        'max_atom_adjustment_nm':float(np.linalg.norm(repaired-target,axis=1).max()),
        'maximum_fractional_bond_error':relative_error,
        'method':'Bond and adjacent-angle distance regularization near inferred gaps; observed outside atoms fixed; not energy minimization.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh asset directory')
    tables={};provenance={};structures={}
    for code in ['4IAC','4IAD','1J3H']:
        path=a.sources/(code+'.cif');meta=json.loads(path.with_suffix('.cif.source.json').read_text())
        if meta['id']!=code or hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']:
            raise ValueError('Source identity/hash mismatch: '+code)
        tables[code]=atom_table(path);structures[code]=gemmi.read_structure(str(path))[0];provenance[code]=meta
    transform,rmsd=rigid_fit(tables['4IAD']['A'],tables['4IAC']['A'])
    origin=xyz(tables['4IAC']['A']).mean(0)
    def coordinates(points,product=False):
        points=np.asarray(points)
        if product:points=points@transform[:3,:3].T+transform[:3,3]
        return (points-origin)*.1
    def ligand(code,name):
        found=[r for c in structures[code] for r in c if r.name==name]
        if len(found)!=1:raise ValueError(f'Expected one {name} in {code}')
        return residue_atoms(found[0])
    substrate=tables['4IAC']['S'];product={(r[3],r[2]):r for r in tables['4IAD']['S']}
    peptide0=coordinates(xyz(substrate))
    peptide1=coordinates(np.array([product[r[3],r[2]][0] for r in substrate]),True)
    extras=set(product)-{(r[3],r[2]) for r in substrate}
    if extras!={(621,n) for n in ['P','O1P','O2P','O3P']}:raise ValueError('Unexpected phosphopeptide atom map')
    analog=ligand('4IAC','ACP');adp=ligand('4IAD','ADP')
    ccd,ccd_meta=ccd_file(a.sources,'ATP');unquote=gemmi.cif.as_string
    rows=list(ccd.find('_chem_comp_atom.',['atom_id','type_symbol','pdbx_model_Cartn_x_ideal','pdbx_model_Cartn_y_ideal','pdbx_model_Cartn_z_ideal']))
    ideal={unquote(r[0]):(r[1],np.array([float(r[k]) for k in (2,3,4)])*.1) for r in rows if r[1]!='H'}
    names=list(ideal);idx={n:i for i,n in enumerate(names)}
    analog_map={n:('C3B' if n=='O3B' else n) for n in names}
    if set(analog_map.values())!=set(analog):raise ValueError('Unexpected ATP/AMP-PCP correspondence')
    start=coordinates(np.array([analog[analog_map[n]][0] for n in names]))
    bonds=[(unquote(r[0]),unquote(r[1])) for r in ccd.find('_chem_comp_bond.',['atom_id_1','atom_id_2'])
           if unquote(r[0]) in idx and unquote(r[1]) in idx]
    edges=np.array([[idx[x],idx[y]] for x,y in bonds])
    target_lengths=np.array([np.linalg.norm(ideal[x][1]-ideal[y][1]) for x,y in bonds])
    def residual(flat):
        pos=flat.reshape(-1,3)
        lengths=np.linalg.norm(pos[edges[:,0]]-pos[edges[:,1]],axis=1)
        return np.r_[(pos-start).ravel(),(lengths-target_lengths)*8]
    fit=least_squares(residual,start.ravel(),max_nfev=100)
    atp=fit.x.reshape(-1,3)
    if not fit.success:raise ValueError('ATP geometric reconstruction did not converge')
    gamma=['PG','O1G','O2G','O3G'];core=[n for n in names if n not in gamma]
    core0=np.array([atp[idx[n]] for n in core]);core1=coordinates(np.array([adp[n][0] for n in core]),True)
    phosphate0=np.array([atp[idx[n]] for n in gamma])
    phosphate1=coordinates(np.array([product[621,n][0] for n in ['P','O1P','O2P','O3P']]),True)
    kinase_rows=tables['4IAC']['A'];kinase=coordinates(xyz(kinase_rows))
    product_rows={(r[3],r[2]):r for r in tables['4IAD']['A']}
    kinase_product=coordinates(np.array([product_rows[r[3],r[2]][0] for r in kinase_rows]),True)
    # Preserve every atom of the closed reference. The opening is explicitly
    # source-guided, including inferred positions where apo atoms are absent.
    apo=tables['1J3H']['A']
    core_rows=lambda rows:[r for r in rows if 128<=r[3]<=300]
    apo_transform,apo_rmsd=rigid_fit(core_rows(apo),core_rows(kinase_rows))
    normal=lambda name:{'SEP':'SER','TPO':'THR','CME':'CYS'}.get(name,name)
    apo_map={(r[3],r[2]):r for r in apo}
    ca=[r for r in kinase_rows if r[2]=='CA' and (r[3],'CA') in apo_map and normal(r[4])==normal(apo_map[r[3],'CA'][4])]
    closed_ca=coordinates(xyz(ca))
    def apo_coordinates(points):
        return coordinates(np.asarray(points)@apo_transform[:3,:3].T+apo_transform[:3,3])
    open_ca=apo_coordinates(np.array([apo_map[r[3],'CA'][0] for r in ca]))
    kinase_open=kinase+displacement_field(kinase,closed_ca,open_ca-closed_ca)
    matched=[]
    for i,r in enumerate(kinase_rows):
        other=apo_map.get((r[3],r[2]))
        if other is not None and normal(r[4])==normal(other[4]):
            kinase_open[i]=apo_coordinates(other[0]);matched.append(i)
    kinase_open,gap_repair=repair_inferred_geometry(kinase,kinase_open,kinase_rows,matched)
    # Retain the observed nucleotide-bound Glu127 side-chain rotamer while the
    # nucleotide is present. Let it adopt the apo reference only during exit.
    # This staged interpolation is an authored accommodation, not an observed
    # induced-fit trajectory; no atoms are removed to open a passage.
    gated0=kinase_open.copy();gated1=kinase_open.copy()
    rotamer=[i for i,r in enumerate(kinase_rows) if r[3]==127 and r[2] in ['CB','CG','CD','OE1','OE2']]
    gated0[rotamer]=kinase[rotamer];gated1[rotamer]=kinase_product[rotamer]
    metals=[]
    for code in ['4IAC','4IAD']:
        rs=sorted([r for r in structures[code]['A'] if r.name=='MG'],key=lambda r:r.seqid.num)
        if [r.seqid.num for r in rs]!=[401,402]:raise ValueError('Unexpected metal site identifiers')
        metals.append(coordinates(np.array([list(r[0].pos) for r in rs]),code=='4IAD'))
    a.output.mkdir(parents=True)
    vertices,faces=surface(kinase*10,.75,1.1,.55)
    vertices=vertices*.1
    np.savez_compressed(a.output/'kinase.npz',atoms=kinase,vertices=vertices,faces=faces,
        open_atoms=kinase_open,product_atoms=kinase_product,
        bound_open_atoms=gated0,product_bound_open_atoms=gated1,
        bonds=peptide_bonds(kinase_rows,kinase),
        residue_numbers=np.array([r[3] for r in kinase_rows]),atom_names=np.array([r[2] for r in kinase_rows]))
    np.savez_compressed(a.output/'peptide.npz',atoms=peptide0,product_atoms=peptide1,
        atom_names=np.array([r[2] for r in substrate]),residue_numbers=np.array([r[3] for r in substrate]),
        elements=np.array([r[1] for r in substrate]),bonds=peptide_bonds(substrate,peptide0))
    np.savez_compressed(a.output/'nucleotide_core.npz',atoms=core0,product_atoms=core1,atom_names=np.array(core),
        elements=np.array([ideal[n][0] for n in core]),bonds=np.array([[core.index(x),core.index(y)] for x,y in bonds if x in core and y in core]))
    np.savez_compressed(a.output/'phosphoryl.npz',atoms=phosphate0,product_atoms=phosphate1,atom_names=np.array(gamma),
        elements=np.array(['P','O','O','O']),bonds=np.array([[0,1],[0,2],[0,3]]))
    np.savez_compressed(a.output/'metals.npz',atoms=metals[0],product_atoms=metals[1],site_names=np.array(['M1','M2']))
    plans={}
    def obstacles(closed,opened,*others):
        return np.vstack([closed,*others]),np.vstack([opened,*others])
    specifications={
        'atp_in':(np.vstack([core0,phosphate0,metals[0][1:]]),obstacles(kinase,kinase_open)),
        'peptide_in':(peptide0,obstacles(kinase,gated0,core0,phosphate0,metals[0])),
        'peptide_out':(np.vstack([peptide1,phosphate1]),obstacles(kinase_product,gated1,core1,metals[1])),
        'adp_out':(np.vstack([core1,metals[1][1:]]),obstacles(gated1,kinase_open)),
        'm1':(metals[0][:1],obstacles(kinase,kinase,core0,phosphate0,metals[0][1:])),
        'm1_out':(metals[1][:1],obstacles(gated1,gated1,core1,metals[1][1:]))}
    state_pairs={'atp_in':['closed','apo'],'peptide_in':['closed','bound_open'],
        'peptide_out':['product','product_bound_open'],'adp_out':['product_bound_open','apo'],
        'm1':['closed','closed'],'m1_out':['product_bound_open','product_bound_open']}
    for name,(mobile,(fixed,opened)) in specifications.items():
        report=None
        for closure in [.8,.95,.6]:
            _,candidate=plan(mobile,fixed,frames=240,distance=7,obstacle_open=opened,closure_start=closure)
            if report is None or (candidate['max_new_or_worsened_overlap_atoms'],len(candidate['flagged_frames'])) < (report['max_new_or_worsened_overlap_atoms'],len(report['flagged_frames'])):
                report=candidate
            if report['collision_screen_passed']:break
        report['enzyme_states_closed_open']=state_pairs[name]
        plans[name]=report;np.savez_compressed(a.output/(name+'-path.npz'),moving=mobile,obstacle=fixed,obstacle_open=opened)
        print(name,report['collision_screen_passed'],report['max_new_or_worsened_overlap_atoms'],flush=True)
    manifest={'kind':'pka-phosphotransfer','sources':provenance,'ccd_atp':ccd_meta,'units':'nm',
        'alignment_rmsd_angstrom':rmsd,'alignment_product_to_reference':transform.tolist(),
        'reaction':'ATP + SP20(Ser21-OH) -> ADP + pSP20(Ser21-O-phosphoryl)',
        'source_model':'Mouse PKAc with engineered SP20 derived from human PKI, not a native full-length substrate',
        'reactant_reconstruction':'ATP CCD atom identity/bonds fitted geometrically to 4IAC AMP-PCP. C3B scaffold correspondence is replaced by oxygen O3B before animation; ATP is not observed in 4IAC.',
        'reconstruction_max_atom_shift_nm':float(np.linalg.norm(atp-start,axis=1).max()),
        'enzyme_opening':{'source':'1J3H author chain A','core_alignment_rmsd_angstrom':apo_rmsd,
            'source_atom_matches':len(matched),'inferred_atoms':len(kinase)-len(matched),
            'gap_geometry_repair':gap_repair,
            'nucleotide_bound_accommodation':'Glu127 CB/CG/CD/OE1/OE2 retain their bound-state positions until nucleotide exit; this timing is authored to avoid a clash, not experimentally measured.',
            'method':'Matched source atoms plus CA-displacement interpolation for unmatched atoms; the displayed envelope is rebuilt from the complete atom set in each frame.',
            'limits':'Apo construct/modified residues differ. No missing source atoms are added; unmatched positions and all intermediate conformations are inferred, not observed.'},
        'chemical_identity':{'peptide_core_atoms':len(peptide0),'nucleotide_core_atoms':len(core),
            'transferred_group_atoms':4,'heavy_atoms_initial':len(peptide0)+len(atp)+2,
            'heavy_atoms_final':len(peptide1)+len(core)+4+2,'transferred_phosphorus_id':'donor_PG',
            'broken_bond':['nucleotide:O3B','donor_PG'],'formed_bond':['peptide:621:OG','donor_PG'],
            'hydrogens_protonation_and_solvent':'omitted; this is a heavy-atom teaching diagram, not a charge/proton-balanced simulation'},
        'metal_sites':{'M1':'author A:401; release illustration before ADP','M2':'author A:402; ADP-associated release illustration'},
        'paths':plans,'path_screen_passed':all(r['collision_screen_passed'] for r in plans.values()),
        'sources_for_mechanism':['https://doi.org/10.1021/bi400066a','https://pmc.ncbi.nlm.nih.gov/articles/PMC4505467/'],
        'limits':['The endpoint structures are snapshots of different crystals, not one measured trajectory.',
            'Reaction interpolation and release order are teaching choices; no transition-state, proton-transfer or kinetic simulation.',
            'Whole kinase is retained in the overview; an active-site detail may omit that whole chain, never slice it.',
            'The intact kinase envelope is regenerated from source-guided 1J3H atom positions. Passing path checks does not validate an induced-fit trajectory.',
            'Source metals and the engineered substrate identity must remain explicit in project notes.']}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
