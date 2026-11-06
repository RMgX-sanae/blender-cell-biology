"""Mutation checks for actual compartment leaks, inventory loss and mesh faults."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
import prepare_vesicle_topology as m


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True,help='Prepared, qualified topology keys')
    a=p.parse_args();grid=m.prepare_grid(.08);checks=[]
    def rejects(name,call,message):
        try:call()
        except ValueError as error:
            if message not in str(error):raise AssertionError((name,str(error))) from error
            checks.append({'case':name,'rejected':True,'message':str(error)})
        else:raise AssertionError('Accepted '+name)
    def check(name,s,phi,message):rejects(name,lambda:m.inspect(s,phi,.08),message)
    s=m.state(.36);modified=copy.deepcopy(s);modified['t']=.34
    check('residual-neck-after-scission',s,m.field(modified,*grid),'connectivity')
    s=m.state(.77);check('premature-fusion-pore',s,m.field(m.state(.8),*grid),'connectivity')
    s=m.state(.8);check('apposition-mislabeled-as-fusion',s,m.field(m.state(.799),*grid),'connectivity')
    s=m.state(.77);phi=m.field(s,*grid);s['cargo'][0]=m.TARGET
    check('premature-cargo-arrival-in-disconnected-endosome',s,phi,'connected compartment')
    s=m.state(.5);phi=m.field(s,*grid);s['cargo'][0]=s['center']+[.88,0,0]
    check('center-inside-but-cargo-edge-through-membrane',s,phi,'Whole cargo crosses')
    s=m.state(.5);phi=m.field(s,*grid);s['normal']*=-1;s['ligand']=s['receptor']-s['normal']*.28
    check('flipped-receptor-domains',s,phi,'connected compartment')
    s=m.state(.5);phi=m.field(s,*grid);s['cargo_ids'][0]='cargo-2'
    check('duplicated-cargo-identity',s,phi,'inventory')
    s=m.state(.5);phi=m.field(s,*grid);s['cargo']=s['cargo'][:-1]
    check('lost-cargo',s,phi,'inventory')
    s=m.state(.5);phi=m.field(s,*grid);s['cargo'][0]=[-30,0,1]
    check('negative-index-grid-wrap',s,phi,'outside the grid')
    s=m.state(.5);phi=m.field(s,*grid);x,y,z=grid
    phi=np.minimum(phi,np.sqrt((x+3)**2+y**2+(z+6)**2)-.5)
    check('unaccounted-fluid-pocket',s,phi,'water pocket')
    s=m.state(.35);modified=copy.deepcopy(s);modified['center'][2]=-1.2
    check('separate-water-but-touching-membranes-at-scission',s,m.field(modified,*grid),'Membrane components')
    s=m.state(.77);phi=np.minimum(m.field(s,*grid),np.sqrt((x-3.2)**2+y*y+(z+6)**2)-2.)
    check('separate-water-but-touching-membranes-at-docking',s,phi,'Membrane components')
    with np.load(a.input/'state_0.800.npz',allow_pickle=False) as data:
        vertices=data['lumen_surface_vertices'];faces=data['lumen_surface_faces']
    s=m.state(.8)
    m.inspect_mesh(s,vertices,faces,-m.HALF,.08)
    mutations={'all-faces-reversed':(faces[:,::-1],'Reversed component'),
               'duplicated-face':(np.concatenate([faces,faces[:1]]),'Non-manifold')}
    flipped=faces.copy();flipped[0]=flipped[0,::-1]
    mutations['one-face-reversed']=(flipped,'adjacent face winding')
    internal=int(np.argmin(np.linalg.norm(vertices[faces].mean(1)-s['center'],axis=1)))
    mutations['internal-membrane-hole']=(np.delete(faces,internal,axis=0),'Euler')
    for name,(changed,message) in mutations.items():
        rejects(name,lambda:m.inspect_mesh(s,vertices,changed,-m.HALF,.08),message)
    print(json.dumps({'passed':True,'negative_cases':len(checks),'checks':checks},indent=2))


if __name__=='__main__':main()
