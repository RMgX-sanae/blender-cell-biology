"""Prepare the registered ribosome actors and a schematic peptide exit path for Blender."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import gemmi
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from skimage.graph import MCP_Geometric

from register_translation_states import prepare as register


def exit_guide(sources,registered):
    config=json.loads((registered/'5LZT.json').read_text())['orientation']
    center=np.asarray(config['origin_angstrom']);axes=np.asarray(config['basis_columns'])
    model=gemmi.read_structure(str(sources/'5LZT.cif'))[0]
    def coords(chain,atom_name=None):
        return np.asarray([list(a.pos) for res in chain for a in res
                           if a.occ>0 and a.altloc in ('\x00','A') and not a.element.is_hydrogen
                           and (atom_name is None or a.name==atom_name)])
    def nm(x):return (x-center)@axes/10
    peptide=model['1'].get_polymer()
    if [(r.seqid.num,r.name) for r in peptide][-1]!=(68,'VAL') or len(peptide)!=15:
        raise ValueError('Expected the resolved 5LZT nascent peptide 54–68')
    trace=nm(coords(peptide,'CA'))[::-1]
    carbon=nm(coords([peptide[-1]],'C'))[0]
    trace=np.vstack([carbon,trace])
    rim=nm(np.vstack([coords(model[ch].get_polymer(),'CA').mean(axis=0) for ch in ('X','Y','h')])).mean(axis=0)
    start=trace[-1];direction=rim-start;direction/=np.linalg.norm(direction)
    target=rim+4.5*direction
    wall=cKDTree(nm(np.vstack([coords(ch) for ch in model if ch.name!='1'])))
    # Between the last deposited residue and the exit, infer only an open route
    # for a teaching trace. The PDB does not supply coordinates for this motion.
    step=.18
    lower=np.minimum(start,target)-3.1;upper=np.maximum(start,target)+3.1
    shape=np.ceil((upper-lower)/step).astype(int)+1
    grid=np.indices(shape).reshape(3,-1).T*step+lower
    clearance=wall.query(grid)[0]
    delta=grid-start;offset=delta-np.outer(delta@direction,direction)
    costs=1+.18*np.sum(offset**2,axis=1)+.35/np.maximum(clearance-.17,.04)
    costs[clearance<.29]=np.inf
    first=tuple(np.rint((start-lower)/step).astype(int));last=tuple(np.rint((target-lower)/step).astype(int))
    graph=MCP_Geometric(costs.reshape(shape));distance,_=graph.find_costs([first],[last])
    if not np.isfinite(distance[last]):raise ValueError('No open exit route in the reviewed tunnel neighborhood')
    path=lower+step*np.array(graph.traceback(last))
    path[0]=start;path[-1]=target
    smooth=gaussian_filter1d(path,1.,axis=0,mode='nearest');smooth[0]=start;smooth[-1]=target
    def gap(points):
        samples=np.concatenate([a+(b-a)*np.linspace(0,1,8)[:,None] for a,b in zip(points,points[1:])])
        return float(wall.query(samples)[0].min())
    if gap(smooth)>=.26:path=smooth
    minimum=gap(path)
    if minimum<.25:raise ValueError('Teaching trace intersects the tunnel wall')
    length=np.linalg.norm(np.diff(trace,axis=0),axis=1).sum()
    guide=np.vstack([trace,path[1:],target+direction*(length+5)])
    return {'source':'5LZT','units':'nm','C_to_N_trace_nm':trace.tolist(),
            'world_path_nm':guide.tolist(),'mouth_landmark_nm':rim.tolist(),
            'outside_nm':target.tolist(),'clearance_after_resolved_segment_nm':minimum,
            'rim_source_chains':{'X':'uL23','Y':'uL24','h':'uL29'},
            'scope':'Resolved C68 carbonyl and 15 C-alpha positions followed by a clearance-guided display path. Only the resolved C-terminal fragment is tracked. This is not a measured exit trajectory or a folded protein.'}


def prepare(sources,donor,output):
    if output.exists():raise ValueError('Choose a fresh preparation output')
    output.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temp:
        pending=Path(temp)/'route';pending.mkdir()
        registered=pending/'registered'
        register(sources,donor,registered)
        for pdb,directory,spacing,sigma in [('6YAL','48s-initiation',2.2,2.6),
                                            ('5LZT','5LZT',1.5,1.85),('5LZU','5LZU',1.5,1.85),
                                            ('5LZV','5LZV',1.5,1.85)]:
            subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_molecular_asset.py')),
                            str(sources/(pdb+'.cif')),'--config',str(registered/(pdb+'.json')),
                            '--spacing',str(spacing),'--sigma',str(sigma),'--level','.8',
                            '--output',str(pending/directory)],check=True)
        guide=exit_guide(sources,registered)
        (pending/'exit-guide.json').write_text(json.dumps(guide,indent=2)+'\n')
        pending.rename(output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--elongation-manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();prepare(args.sources,args.elongation_manifest,args.output)
