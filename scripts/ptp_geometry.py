"""PTP1B net-state geometry on one persistent source construct, in Angstrom."""
from pathlib import Path
import hashlib, itertools, json
import gemmi
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from prepare_molecular_asset import atom_table, xyz
from prepare_filament_assets import rigid_fit
from prepare_ras_switch import fit_points
from prepare_phosphotransfer import peptide_bonds
from prepare_ras_hydrolysis_motion import chemical_samples
from ubiquitin_geometry import rotate, index
from reaction_ledger import validate_sequence

SOURCE_HASHES = {
    '3I7Z':'1fc549b2c683239d9944e801982e08f8ae43294246d737970dec72ef4691cbfc',
    '1A5Y':'eb0039339cbf9e1e5fcc938ced3be4161c2a144faf2bc470c3ca07a3bf64ef83',
    '3I80':'7985be81d53995c7170dbc7710df246f2fed5cf7b0fc3e31e569f658f610d20d',
}


def unit(vector): return vector/np.linalg.norm(vector)


def tetrahedral(center, apical, reference, length=1.53):
    axis=unit(apical-center)
    tangent=reference[0]-center; tangent=unit(tangent-axis*(tangent@axis))
    cross=np.cross(axis,tangent)
    def endpoint(phi,order):
        angles=phi+np.array(order)*2*np.pi/3
        directions=-axis/3+np.sqrt(8/9)*(np.cos(angles)[:,None]*tangent+np.sin(angles)[:,None]*cross)
        return np.vstack([center,center+length*directions])
    best=None
    for order in itertools.permutations(range(3)):
        fit=minimize_scalar(lambda phi:float(np.sum((endpoint(phi,order)[1:]-reference)**2)),bounds=(-np.pi,np.pi),method='bounded')
        if best is None or fit.fun<best[0]:best=(fit.fun,fit.x,order)
    return endpoint(best[1],best[2])


class PtpGeometry:
    def __init__(self, directory):
        directory=Path(directory); models={}; tables={}; records={}
        for code,expected in SOURCE_HASHES.items():
            path=directory/(code+'.cif')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
                raise ValueError('Source changed; requalify '+code)
            record=json.loads(path.with_suffix('.cif.source.json').read_text())
            if record['id'].upper()!=code or record['sha256']!=expected:
                raise ValueError('Source provenance mismatch: '+code)
            models[code]=gemmi.read_structure(str(path))[0];tables[code]=atom_table(path)
            records[code]={k:record[k] for k in ['id','url','entry_url','sha256','license','license_url']}
        self.rows={'enzyme':tables['3I7Z']['A'],'peptide':tables['3I7Z']['B']}
        self.source={k:xyz(rows) for k,rows in self.rows.items()}
        self.indices={k:{(r[3],r[2]):i for i,r in enumerate(rows)} for k,rows in self.rows.items()}
        for key,rows in self.rows.items():
            if len(rows)!=len(self.indices[key]):raise ValueError('Ambiguous author residue/atom selection: '+key)
        if len(self.rows['enzyme'])!=2426 or len(self.rows['peptide'])!=51:
            raise ValueError('Unexpected resolved construct atom inventory')
        for number,name in [(215,'CYS'),(262,'GLN')]:
            if {r[4] for r in self.rows['enzyme'] if r[3]==number}!={name}:
                raise ValueError('Require native catalytic identities in the persistent construct')
        def atom(code,chain,residue,name):
            residues=[r for r in models[code][chain] if r.seqid.num==residue]
            matches=[a for r in residues for a in r if a.name==name and not a.element.is_hydrogen and a.occ>0 and a.altloc in ('\x00','A')]
            if len(matches)!=1:raise ValueError('Ambiguous source reaction atom')
            return np.array(list(matches[0].pos))
        def points(code,chain,residue,names):return np.array([atom(code,chain,residue,n) for n in names])
        e=self.source['enzyme'];self.sg=e[self.indices['enzyme'][215,'SG']]
        self.oh=self.source['peptide'][self.indices['peptide'][992,'OH']]
        analog=atom('3I7Z','A',322,'V');oxygen=points('3I7Z','A',322,['O1','O2','O3'])
        self.first=tetrahedral(self.oh+unit(analog-self.oh)*1.60,self.oh,oxygen)
        anchor=['N','CA','CB','SG'];old=points('1A5Y','A',215,anchor);new=points('3I7Z','A',215,anchor)
        rotation,shift=fit_points(old,new)
        middle=(points('1A5Y','A',215,['P','O1P','O2P','O3P'])-old[-1])@rotation.T+self.sg
        permutation=min(itertools.permutations(range(1,4)),key=lambda order:float(np.sum((middle[list(order)]-self.first[1:])**2)))
        self.middle=middle[[0,*permutation]]

        excluded={n for begin,end in [(47,49),(119,121),(179,187),(215,222),(261,262)] for n in range(begin,end+1)}
        core=lambda rows:[r for r in rows if r[3] not in excluded]
        matrix,rmsd=rigid_fit(core(tables['3I80']['A']),core(self.rows['enzyme']))
        native=xyz(tables['3I80']['A'])@matrix[:3,:3].T+matrix[:3,3]
        self.q_atoms=np.array([index(self.rows['enzyme'],262,n) for n in ['CG','CD','OE1','NE2']])
        self.q_axis=[index(self.rows['enzyme'],262,n) for n in ['CA','CB','CG']]
        native_ids=[index(tables['3I80']['A'],262,n) for n in ['CG','CD','OE1','NE2']]
        self.q_angles=np.zeros(2)
        def residual(angles):return (self.q_pose(angles)[self.q_atoms]-native[native_ids]).ravel()
        fit=least_squares(residual,[0.,np.radians(-120.)])
        if not fit.success:raise ValueError('Source-guided Q262 torsion fit failed')
        self.q_angles=fit.x;self.final_enzyme=self.q_pose(self.q_angles)
        rotation_w,shift_w=fit_points(native[native_ids[1:]],self.final_enzyme[self.q_atoms[1:]])
        water_reference=atom('3I80','A',322,'O4')@matrix[:3,:3].T+matrix[:3,3]
        self.water=water_reference@rotation_w.T+shift_w
        self.product=np.vstack([tetrahedral(self.water-unit(self.water-self.middle[0])*1.53,self.water,self.middle[1:]),self.water])
        self.bonds={key:peptide_bonds(rows,self.source[key]*.1) for key,rows in self.rows.items()}
        self.neighbors={}
        for key,edges in self.bonds.items():
            neighbors=[set() for _ in self.rows[key]]
            for a,b in edges:neighbors[a].add(int(b));neighbors[b].add(int(a))
            self.neighbors[key]=neighbors
        self.record={'sources':records,'persistent_construct':'3I7Z resolved author A residues 2-298 and B peptide 988-993, including terminal OXT. Native termini, ER anchor and regulatory pathway not modeled.',
            'first':'Modeled phosphate, tetrahedral P-O=1.53 A and Tyr OH-P=1.60 A; 3I7Z analog directions guide orientation. No V-to-P atom substitution.',
            'intermediate':'Only the CSP215 phosphoryl group from Q262A 1A5Y is transplanted in the local Cys frame. The mutant enzyme is not imported.',
            'CSP_local_anchor_rmsd_A':float(np.sqrt(np.mean(np.sum((old@rotation.T+shift-new)**2,axis=1)))),
            'CSP_source_oxygen_order':[ ['P','O1P','O2P','O3P'][i] for i in [0,*permutation]],
            'second_source_core_CA_rmsd_A':rmsd,'Q262_angles_degrees':np.degrees(self.q_angles).tolist(),
            'Q262_model_to_source_rmsd_A':float(np.sqrt(np.mean(np.sum(residual(fit.x).reshape(-1,3)**2,axis=1)))),
            'water':'Modeled water O placed by the 3I80 apical analog O4 relative to the Q262 amide frame; it is not an observed solvent-water atom.',
            'water_site_displacement_from_aligned_analog_A':float(np.linalg.norm(self.water-water_reference)),
            'product':'Idealized Pi with persistent phosphate/water identities; not the tungstate product analog.',
            'oxygen_identity':'Equivalent nonbridging oxygen correspondence fixed once at preparation by geometric matching; not isotope-resolved experimental trajectories.'}

    def q_pose(self,angles):
        points=self.source['enzyme'].copy();ca,cb,cg=self.q_axis
        rotate(points,ca,cb,self.q_atoms,angles[0]);rotate(points,cb,cg,self.q_atoms[1:],angles[1])
        return points

    def chemistry(self,first_fraction,second_fraction):
        first=chemical_samples(self.first,self.water,list(range(4)),np.vstack([self.middle,self.water]),first_fraction)
        second=chemical_samples(self.middle,self.water,list(range(4)),self.product,second_fraction)
        return np.where(np.asarray(second_fraction)[:,None,None]>0,second,first)

    def ledger(self):
        atoms=[[f'{key}:{r[3]}:{r[2]}',r[1]] for key,rows in self.rows.items() for r in rows]
        atoms += [[f'phosphoryl:{name}',element] for name,element in [('P','P'),('O1','O'),('O2','O'),('O3','O')]]+[['water:O','O']]
        base=[['phosphoryl:P','phosphoryl:'+name] for name in ['O1','O2','O3']]
        start=['peptide:992:OH','phosphoryl:P']; middle=['enzyme:215:SG','phosphoryl:P']; end=['water:O','phosphoryl:P']
        def stage(before,after,label):
            return {'version':1,'scope':'All selected resolved enzyme/peptide heavy atoms, modeled phosphoryl group and one water oxygen; waiting and released objects persist. '+label,
                'bond_scope':'P-O nonbridging and transferred-group carrier bonds only; no hydrogens, charges, bond order, proton relay or full enzyme covalent graph.',
                'reactant_atoms':atoms,'product_atoms':atoms,'atom_map':[[identity,identity] for identity,_ in atoms],
                'reactant_bonds':base+[before],'product_bonds':base+[after],'break':[before],'form':[after]}
        book={'version':1,'stages':[stage(start,middle,'Covalent enzyme intermediate.'),stage(middle,end,'Intermediate hydrolysis.')],
              'net':stage(start,end,'Independent net substrate dephosphorylation; Cys carrier bond absent at both net endpoints.')}
        validate_sequence(book)
        return book
