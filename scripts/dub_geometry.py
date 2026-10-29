"""Source-specific AMSH-LP net-hydrolysis geometry; inferred display poses, not MD."""
from pathlib import Path

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

from prepare_molecular_asset import atom_table, xyz
from prepare_phosphotransfer import peptide_bonds
from ubiquitin_geometry import index, rotate, digest

SOURCE_HASHES = {
    '2ZNV': 'f04dee313268542779080a2cae58b6a9a0ba3fece0af79d226a254ddaf0d1774',
    '2ZNR': '610a84b621d933d43c0db1b882b8cd366b8b5c70dd4e775db0cc2aea019b83a7',
}


def fit(first, second):
    first, second = np.asarray(first), np.asarray(second)
    u, _, vt = np.linalg.svd((first-first.mean(0)).T@(second-second.mean(0)))
    r = vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
    return r, second.mean(0)-r@first.mean(0)

def angle(a,b,c):
    a,c=a-b,c-b
    return np.arccos(np.clip(a@c/np.linalg.norm(a)/np.linalg.norm(c),-1,1))

class DubGeometry:
    def __init__(self, directory):
        sources = Path(directory)
        for name, expected in SOURCE_HASHES.items():
            if digest(sources/(name+'.cif')) != expected:
                raise ValueError('Source changed; requalify the construct: '+name)
        bound=atom_table(sources/'2ZNV.cif'); apo=atom_table(sources/'2ZNR.cif')
        source={r[2]:r for r in apo['A'] if r[3]==292}
        target={r[2]:r for r in bound['A'] if r[3]==292}
        anchors=['N','CA','C','CB']
        r,t=fit([source[n][0] for n in anchors],[target[n][0] for n in anchors])
        enzyme=[list(row) for row in bound['A']]
        for row in enzyme:
            if row[3]==292:row[4]='GLU'
        enzyme += [[(r@np.asarray(source[n][0])+t).tolist(),source[n][1],n,292,'GLU'] for n in ['CG','CD','OE1','OE2']]
        enzyme.sort(key=lambda row:(row[3],row[2]))
        self.rows={'enzyme':enzyme,'distal':bound['B'],'proximal':bound['C']}
        self.source={k:xyz(v) for k,v in self.rows.items()}
        ligands=[(347,'NE2'),(349,'NE2'),(360,'OD2')]
        aa=np.array([apo['A'][index(apo['A'],i,n)][0] for i,n in ligands])
        bb=np.array([enzyme[index(enzyme,i,n)][0] for i,n in ligands])
        r,t=fit(aa,bb)
        a=gemmi.read_structure(str(sources/'2ZNR.cif'))[0]['A']
        z=np.array(list(next(res for res in a if res.name=='ZN' and res.seqid.num==2)[0].pos))
        w=np.array(list(next(res for res in a if res.name=='HOH' and res.seqid.num==517)[0].pos))
        lengths=np.linalg.norm(aa-z,axis=1)
        z0=r@z+t
        zinc_fit=least_squares(lambda point:np.linalg.norm(bb-point,axis=1)-lengths,z0,gtol=1e-12,ftol=1e-12,xtol=1e-12)
        self.zinc=zinc_fit.x;self.water=r@(w-z)+self.zinc
        self.record={'metal_target_lengths_A':lengths.tolist(),'metal_anchor_fit_rmsd_A':float(np.sqrt(np.mean(np.sum((aa@r.T+t-bb)**2,axis=1)))),
            'metal_adjustment_A':float(np.linalg.norm(self.zinc-z0)), 'metal_fit_error_A':float(np.max(abs(zinc_fit.fun))),
            'source_water':'2ZNR A:517, rigidly oriented by the catalytic ligand triangle and translated with the fitted Zn; inferred in the 2ZNV substrate frame.'}
        self.axes=[]
        for residue in [74,75,76]:
            rows=self.rows['distal']
            self.axes += [('distal',index(rows,residue,'N'),index(rows,residue,'CA'),[i for i,v in enumerate(rows) if v[3]>residue or (v[3]==residue and v[2]!='N')]),
                          ('distal',index(rows,residue,'CA'),index(rows,residue,'C'),[i for i,v in enumerate(rows) if v[3]>residue or (v[3]==residue and v[2] in ['O','OXT'])])]
        rows=self.rows['proximal']
        for a,b,mobile in [('CA','CB',['CG','CD','CE','NZ']),('CB','CG',['CD','CE','NZ']),('CG','CD',['CE','NZ']),('CD','CE',['NZ'])]:
            self.axes.append(('proximal',index(rows,63,a),index(rows,63,b),[index(rows,63,n) for n in mobile]))
        rows=self.rows['distal']
        self.arg_axes=[(index(rows,74,a),index(rows,74,b),[index(rows,74,n) for n in mobile])
            for a,b,mobile in [('CA','CB',['CG','CD','NE','CZ','NH1','NH2']),('CB','CG',['CD','NE','CZ','NH1','NH2'])]]
        self.ids={n:index(self.rows['distal'],76,n) for n in ['C','O','CA']}
        self.ids.update({n:index(self.rows['proximal'],63,n) for n in ['NZ','CE']})
        self.graphs={};self.self_pairs={}
        for key,rows in self.rows.items():
            edges=peptide_bonds(rows,self.source[key]*.1);self.graphs[key]=edges
            neighbors=[set() for _ in rows]
            for a,b in edges:neighbors[a].add(int(b));neighbors[b].add(int(a))
            mobile=set(i for k,a,b,selected in self.axes if k==key for i in selected)
            if key=='enzyme':mobile={i for i,row in enumerate(rows) if row[3]==292 and row[2] in ('CG','CD','OE1','OE2')}
            pairs=set()
            for i in mobile:
                exclude={i}|neighbors[i]
                for n in neighbors[i]:exclude|=neighbors[n]
                pairs.update(tuple(sorted([i,j])) for j in range(len(rows)) if j not in exclude)
            self.self_pairs[key]=np.array(sorted(pairs))
        self.excluded={(self.ids['C'],self.ids['NZ']),(self.ids['O'],self.ids['NZ']),
                       (self.ids['CA'],self.ids['NZ']),(self.ids['C'],self.ids['CE'])}
        structural = [list(atom.pos) for chain in gemmi.read_structure(str(sources/'2ZNV.cif'))[0]
                      if chain.name == 'A' for residue in chain if residue.name == 'ZN' for atom in residue]
        if len(structural) != 1: raise ValueError('Expected one selected structural zinc')
        self.metals = np.vstack([structural, self.zinc])
        self.metal_ligands = [[index(enzyme, res, atom) for res, atom in pairs] for pairs in
                             [[(362,'NE2'),(402,'SG'),(408,'NE2'),(410,'NE2')], ligands]]
        self.record.update(restoration='2ZNV author A E292A scaffold; restore only Glu292 side chain from 2ZNR A local N/CA/C/CB fit.',
            retained_constructs='Distal B K63R; proximal C D77 extension. These are experimental constructs, not wild-type ubiquitin.',
            omitted='Other crystallographic copies, bulk solvent, praseodymium, ethylene glycol, hydrogens and unresolved native enzyme residues.')


    def pose(self,x):
        p={k:v.copy() for k,v in self.source.items()}
        for (key,a,b,indices),q in zip(self.axes,x[:10]):rotate(p[key],a,b,indices,q)
        for args,q in zip(self.arg_axes,x[22:24]):rotate(p['distal'],*args,q)
        for key,offset in [('distal',10),('proximal',16)]:
            center=self.source[key].mean(0)
            p[key]=(p[key]-center)@Rotation.from_rotvec(x[offset:offset+3]).as_matrix().T+center+x[offset+3:offset+6]
        return p

    def values(self,p):
        c,o,ca=[p['distal'][self.ids[n]] for n in ['C','O','CA']]
        n,ce=[p['proximal'][self.ids[n]] for n in ['NZ','CE']]
        normal=np.cross(o-c,ca-c);normal/=np.linalg.norm(normal)
        return c,o,ca,n,ce,normal

    def report(self,x):
        p=self.pose(x);c,o,ca,n,ce,normal=self.values(p)
        return {'C_N_A':float(np.linalg.norm(c-n)),'C_water_A':float(np.linalg.norm(c-self.water)),
            'O_C_N_deg':float(np.degrees(angle(o,c,n))),'C_N_CE_deg':float(np.degrees(angle(c,n,ce))),
            'N_plane_A':float(abs((n-c)@normal)),'CE_plane_A':float(abs((ce-c)@normal)),
            'water_plane_A':float(abs((self.water-c)@normal)),'O_C_water_deg':float(np.degrees(angle(o,c,self.water))),
            'source_shift_A':{k:float(np.linalg.norm(v-self.source[k],axis=1).max()) for k,v in p.items()},'source':self.record}

def terminal_frame(points, ids):
    c, o, ca = [points[ids[name]] for name in ('C','O','CA')]
    u = (ca-c)/np.linalg.norm(ca-c); v = o-c
    v -= u*(u@v); v /= np.linalg.norm(v)
    return c, np.array([u, v, np.cross(u,v)])


def oxygen_transport(model, product):
    """Keep the same water-derived oxygen covalently attached during release."""
    c, frame = terminal_frame(product['distal'], model.ids)
    coefficients = frame@(model.water-c)
    def position(points):
        c, frame = terminal_frame(points, model.ids)
        return c+coefficients@frame
    return position
