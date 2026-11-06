"""Reuse the binding workflow for a source-oriented AP2/cargo membrane example."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import numpy as np
from scipy.spatial import cKDTree
from fetch_structure import download
from prepare_molecular_asset import atom_table
from interaction_membrane import plane_clearance,tail_paths

SOURCE_SHA='4cdfad6dfe3bd22c28aa117c70e97142319b7a75b94c129d892ceb52a63ea6cb'
MAP_SHA='b381446e1ae0fe9598068350baae6a00742c564a69b6462c09aca1da7bb6ff9b'
MAP_URL='https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-10748/map/emd_10748.map.gz'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--frames',type=int,default=240)
    p.add_argument('--mode',choices=['loop','teaching'],default='teaching');a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    source=a.sources/'6YAF.cif';map_file=a.sources/'emd_10748.map.gz'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=SOURCE_SHA:raise ValueError('Requalify the changed 6YAF source')
    if not map_file.exists():map_file.write_bytes(download(MAP_URL))
    if hashlib.sha256(map_file.read_bytes()).hexdigest()!=MAP_SHA:raise ValueError('Requalify the changed EMD-10748 map')
    raw=gzip.decompress(map_file.read_bytes())
    if struct.unpack('<4i',raw[:16])!=(128,128,128,2) or struct.unpack('<3i',raw[64:76])!=(1,2,3):
        raise ValueError('Expected the qualified map dimensions, mode and axes')
    spacing=np.array(struct.unpack('<3f',raw[40:52]))/np.array(struct.unpack('<3i',raw[28:40]))
    data=np.frombuffer(raw,dtype='<f4',offset=1024+struct.unpack('<i',raw[92:96])[0]).reshape(128,128,128)
    profile=data.mean((1,2));z=np.arange(128)*spacing[2];peaks=[]
    for low,high in [(35,55),(70,90)]:
        candidates=np.flatnonzero((z>=low)&(z<=high));peaks.append(float(z[candidates[np.argmax(profile[candidates])]]))
    script=Path(__file__).resolve().parent
    subprocess.run([sys.executable,str(script/'prepare_interaction.py'),'--examples',str(script.parent/'assets/ap2-capture.json'),
                    '--example','ap2-cargo-membrane','--sources',str(a.sources),'--output',str(a.output),
                    '--frames',str(a.frames),'--mode',a.mode],check=True)
    path=a.output/'manifest.json';m=json.loads(path.read_text());poses=np.load(a.output/'motion.npz')['matrices']
    atoms={};surfaces={}
    for role in ['moving','anchor']:
        bundles=[dict(np.load(a.output/actor['file'])) for actor in m['actors'] if actor['role']==role]
        atoms[role]=np.concatenate([b['atoms'] for b in bundles]);surfaces[role]=np.concatenate([b['vertices'] for b in bundles])
    origin=np.asarray(m['source_origin_angstrom']);mid=(peaks[0]+peaks[1])/2
    membrane={'kind':'source-oriented-planar-membrane','version':1,'normal':[0,0,1],
       'origin_nm':((np.array([113.92,113.92,mid])-origin)*.1).tolist(),
       'head_offset_nm':(peaks[1]-peaks[0])/20,'hydrophobic_half_width_nm':1.2,
       'head_radius_nm':.28,'tail_radius_nm':.055,'atom_guard_nm':.5,
       'view':{'azimuth_degrees':120,'elevation_degrees':5,
               'evidence':'Authored view selected by inspecting intact 6YAF cargo visibility; not an experimental camera.'},
       'map_url':MAP_URL,'map_sha256':MAP_SHA,'map_density_peak_z_angstrom':peaks,
       'plane_evidence':'Z orientation and bilayer density peaks from EMD-10748; a local planar display reference, not a universal bilayer thickness.',
       'limits':'Hydrophobic guard width, lipid symbols, packing and unresolved N-terminal cargo tether are authored. No atomically resolved PIP2, whole receptor, closed-to-open AP2 transition or uptake completion.'}
    checks={}
    for role in atoms:
        transform=poses if role=='moving' else np.eye(4)[None]
        atom_gap=plane_clearance(atoms[role],transform,membrane,.5)
        mesh_gap=plane_clearance(surfaces[role],transform,membrane)
        checks[role]={'atom_envelope_min_clearance_nm':float(atom_gap.min()),'surface_min_clearance_nm':float(mesh_gap.min())}
        if min(atom_gap.min(),mesh_gap.min())<0:raise ValueError('Protein or cargo crosses the authored membrane core')
    # Stable random packing, excluding the full trajectory of each intact actor.
    rng=np.random.default_rng(20261106);lipids=[];center=np.concatenate(list(atoms.values())).mean(0)
    candidates=rng.uniform([-10,-7],[10,7],(1600,2))+center[:2]
    trees={r:cKDTree(v) for r,v in atoms.items()}
    for side in [-1,1]:
        accepted=[]
        for xy in candidates:
            if accepted and np.min(np.linalg.norm(np.array(accepted)-xy,axis=1))<.63:continue
            accepted.append(xy)
        heads=np.column_stack([accepted,np.full(len(accepted),membrane['origin_nm'][2]+side*membrane['head_offset_nm'])])
        angles=rng.uniform(0,2*np.pi,len(heads));safe=np.ones(len(heads),bool)
        # Head clearance includes an illustrative molecular-envelope allowance.
        probes=np.stack([np.concatenate([h[None],*tail_paths(h,side,angle,membrane['head_offset_nm'])]) for h,angle in zip(heads,angles)])
        radii=np.r_[.68,np.full(18,.4)]
        for role,tree in trees.items():
            for pose in (poses if role=='moving' else np.eye(4)[None]):
                local=(probes-pose[:3,3])@pose[:3,:3]
                distances=tree.query(local.reshape(-1,3),k=1)[0].reshape(len(heads),-1)
                safe&=np.all(distances>=radii,axis=1)
        for head,angle in zip(heads[safe],angles[safe]):lipids.append([*head.tolist(),side,float(angle)])
    membrane['lipids']=lipids
    first_n=next(row[0] for row in atom_table(source)['P'] if row[2]=='N' and row[3]==1)
    first_n=(np.asarray(first_n)-origin)*.1
    # Attach to one actual displayed lipid, retaining that identity. The link is
    # unresolved in the coordinates and is not an invented atomic backbone.
    candidates=sorted((i for i,l in enumerate(lipids) if l[3]==1),key=lambda i:np.linalg.norm(np.asarray(lipids[i][:3])-first_n))
    for index in candidates:
        tether=np.linspace(first_n,lipids[index][:3],12)
        tether_gap=min(float(trees['moving'].query((tether-pose[:3,3])@pose[:3,:3])[0].min()) for pose in poses)
        if tether_gap>=.20:break
    else:raise ValueError('No clear schematic link to a membrane lipid')
    membrane['cargo_lipid_index']=index
    membrane['cargo_tether_atom']={'actor':'anchor_P','residue_number':1,'atom_name':'N'}
    membrane['cargo_tether_nm']=tether.tolist();membrane['trajectory_checks']=checks
    membrane['tether_min_ap2_atom_distance_nm']=tether_gap
    membrane['lipid_exclusion']={'checked_poses':len(poses),'head_to_actor_atom_centers_nm':.68,
       'sampled_tail_to_actor_atom_centers_nm':.4,'tether_to_ap2_atom_centers_nm':.20,
       'limits':'Point/envelope proximity screening; not atomistic packing, van der Waals or membrane energetics.'}
    m['membrane']=membrane;m['membrane_screen_passed']=True
    path.write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'membrane_screen_passed':True,'checks':checks,'lipids':len(lipids),'tether_min_nm':tether_gap}))


if __name__=='__main__':main()
