"""Join validated exchange and hydrolysis without resetting nucleotide identities.

This is one Ras state cycle with open substrate/product pools, not a closed
chemical loop. No previous movie or Blender scene is used as an input.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from interaction_motion import source_contacts
from prepare_phosphotransfer_motion import smooth
from prepare_ras_exchange_motion import apply, pair_check
from prepare_ras_hydrolysis_motion import merge_contacts


def caption_rows(exchange,hydrolysis,ratio):
    (a,b),(c,d)=exchange['intervals']
    release=a+(b-a)*exchange['plans']['old_gdp']['release_start']
    loaded=c+(d-c)*(1-exchange['plans']['new_gtp']['release_start'])
    rows=[(0,a,'Ras 与旧 GDP 结合'),(a,release,'SOS 接近 Ras–GDP'),
          (release,b,'旧 GDP 与 Mg2+ 释放'),(b,c,'无核苷酸的 Ras–SOS 复合物'),
          (c,loaded,'新的 GTP 与 Mg2+ 进入'),(loaded,.78,'Ras 转为 GTP 活性态，SOS 离开'),
          (.78,1,'SOS 离开，水中氧进入反应位点')]
    result=[(x*ratio,y*ratio,label) for x,y,label in rows if y>x]
    binding=hydrolysis['intervals']['binding'];reaction=hydrolysis['intervals']['hydrolysis'];end=hydrolysis['intervals']['release']
    departure=end[0]+.22*(end[1]-end[0])
    rows=[(0,binding[0],'Ras 与新 GTP 结合'),(binding[0],binding[1],'GAP 接近并稳定反应构象'),
          (binding[1],reaction[0],'催化位点对齐'),(reaction[0],reaction[1],'GTP 水解：水中的氧进入 Pi'),
          (reaction[1],departure,'新 GDP 与 Pi 暂时留在复合物中'),
          (departure,end[1],'Pi 释放，Ras 回到 GDP 关闭态'),(end[1],1,'Ras 结合新 GDP，等待下一轮交换')]
    return result+[(ratio+(1-ratio)*x,ratio+(1-ratio)*y,label) for x,y,label in rows if y>x]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exchange-dir',type=Path,required=True);p.add_argument('--hydrolysis-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    em=json.loads((a.exchange_dir/'manifest.json').read_text());hm=json.loads((a.hydrolysis_dir/'manifest.json').read_text())
    if not em['passed'] or not hm['passed']:p.error('Both independent motion checks must pass first')
    if em['mode']!=hm['mode'] or em['fps']!=hm['fps']:p.error('Use matching modes and fps')
    assets=Path(hm['asset_dir'])
    if assets.resolve()!=Path(em['asset_dir']).resolve():p.error('Prepare both stages from the same source asset directory')
    e=dict(np.load(a.exchange_dir/'exchange.npz'));h=dict(np.load(a.hydrolysis_dir/'hydrolysis.npz'))
    d={k:dict(np.load(assets/(k+'.npz'))) for k in ['ras','sos','gap','new_gtp','old_gdp','gap_gdp']}
    before=em['frames'];after=hm['frames'];count=before+after;ratio=before/count
    names=list(d['new_gtp']['atom_names']);core=[names.index(n) for n in hm['core_atom_names']];gamma=[names.index(n) for n in hm['gamma_atom_names']]
    handoff={key:float(np.max(np.abs(first-last))) for key,first,last in [
        ('ras',e['ras'][-1],h['ras'][0]),('core',e['new_gtp'][-1,core],h['nucleotide_core'][0]),
        ('gamma',e['new_gtp'][-1,gamma],h['phosphate'][0]),('magnesium',e['new_mg'][-1],h['magnesium'][0])]}
    if max(handoff.values())>1e-6:p.error('Stage endpoints do not agree; do not bridge a molecular coordinate jump')
    def join(first,second):return np.concatenate([first[:-1],second])
    def hold(points,length):return np.repeat(points[None],length,axis=0)
    # One catalytic water oxygen is explicitly brought from the omitted solvent
    # pool after GTP loading. Fixed early parking in the reaction site would
    # collide with the GDP/GEF Ras conformations.
    first=round(before*.78);progress=smooth((np.arange(before+1)-first)/(before-first))
    direction=np.array([-.6974195511633267,-.22249585873706484,.68125])
    lateral=np.array([-.0205326086,-.8842978830,-.4664715062])
    water=h['water_oxygen'][0]+(1-progress[:,None,None])*14*direction+14*(1-progress[:,None,None])**2*lateral
    # The common camera fits both enzymes. Park the approaching GAP farther out
    # until its phase starts so its distant tips do not leak into the SOS panel.
    route=json.loads((Path(__file__).resolve().parents[1]/'assets/ras-gap-route.json').read_text())
    time=np.linspace(0,1,after+1);release=hm['intervals']['release']
    exit_phase=(time-release[0])/(release[1]-release[0])
    gap_offset=10*(1-smooth(time/(hm['intervals']['binding'][1]*.55)))+10*smooth((exit_phase-.55)/.45)
    gap=h['gap']+gap_offset[:,None,None]*np.array(route['direction'])
    samples={
        'ras':join(e['ras'],h['ras']),
        'gap':join(hold(gap[0],before+1),gap),
        'sos_pose':join(e['sos_pose'],hold(e['sos_pose'][-1],after+1)),
        'nucleotide_core':join(e['new_gtp'][:,core],h['nucleotide_core']),
        'phosphate':join(e['new_gtp'][:,gamma],h['phosphate']),
        'magnesium':join(e['new_mg'],h['magnesium']),
        'water_oxygen':join(water,h['water_oxygen']),
        'old_gdp':join(e['old_gdp'],hold(e['old_gdp'][-1],after+1)),
        'old_mg':join(e['old_mg'],hold(e['old_mg'][-1],after+1)),
        'transfer':join(np.zeros(before+1),h['transfer']),
        'pi_release':join(np.zeros(before+1),h['pi_release']),
        'core_bonds':h['core_bonds'],'old_gdp_bonds':d['old_gdp']['bonds']}
    nt=d['new_gtp'];gd=d['gap_gdp'];old=d['old_gdp'];ras=d['ras']
    gap_indices=[list(gd['atom_names']).index(n) for n in hm['core_atom_names']]
    old_indices=[list(old['atom_names']).index(n) for n in hm['core_atom_names']]
    contacts=lambda x,y:source_contacts(x,y,.20)
    allow={
        'old_gdp_old_mg':contacts(old['atoms'],old['magnesium']),
        'old_gdp_ras':contacts(old['atoms'],ras['atoms']),'old_mg_ras':contacts(old['magnesium'],ras['atoms']),
        'nucleotide_core_phosphate':contacts(nt['atoms'][core],nt['atoms'][gamma]),
        'nucleotide_core_magnesium':merge_contacts(contacts(nt['atoms'][core],nt['magnesium']),contacts(gd['atoms'][gap_indices],gd['magnesium']),contacts(old['atoms'][old_indices],old['magnesium'])),
        'nucleotide_core_ras':merge_contacts(contacts(nt['atoms'][core],ras['gtp']),contacts(gd['atoms'][gap_indices],ras['gap']),contacts(old['atoms'][old_indices],ras['atoms'])),
        'phosphate_water_oxygen':{(0,0):.155},'phosphate_magnesium':contacts(nt['atoms'][gamma],nt['magnesium']),
        'phosphate_ras':contacts(nt['atoms'][gamma],ras['gtp']),
        'water_oxygen_ras':contacts(h['water_oxygen'][0],ras['gap']),
        'magnesium_ras':merge_contacts(contacts(nt['magnesium'],ras['gtp']),contacts(gd['magnesium'],ras['gap']),contacts(old['magnesium'],ras['atoms'])),
        'sos_ras':contacts(d['sos']['atoms'],ras['gef']),'gap_ras':contacts(d['gap']['atoms'],ras['gap'])}
    actors={k:samples[k] for k in ['old_gdp','old_mg','nucleotide_core','phosphate','water_oxygen','magnesium']}
    actors.update(sos=apply(d['sos']['atoms'],samples['sos_pose']),gap=samples['gap'],ras=samples['ras'])
    checks={}
    for first,second in itertools.combinations(actors,2):
        label=first+'_'+second;counts,distances=pair_check(actors[first],actors[second],allow.get(label))
        checks[label]={'maximum_overlap_atoms':int(counts.max()),'bad_frames':(np.flatnonzero(counts)+1).tolist(),'minimum_nm':float(distances.min())}
        print(label,checks[label]['maximum_overlap_atoms'],flush=True)
    bonds={}
    for key,edges in [('ras',ras['bonds']),('gap',d['gap']['bonds']),('sos',d['sos']['bonds']),('old_gdp',old['bonds']),('nucleotide_core',h['core_bonds'])]:
        distances=np.linalg.norm(actors[key][:,edges[:,0]]-actors[key][:,edges[:,1]],axis=2)
        bonds[key]=[float(distances.min()),float(distances.max())]
    report={'kind':'ras-state-cycle','asset_dir':str(assets.resolve()),'frames':count,'fps':em['fps'],'mode':em['mode'],
        'exchange_frames':before,'hydrolysis_frames':after,'handoff_index':before,'handoff_coordinate_errors_nm':handoff,
        'inputs':{name:{'directory':str(folder.resolve()),'motion_sha256':hashlib.sha256((folder/file).read_bytes()).hexdigest()}
                  for name,folder,file in [('exchange',a.exchange_dir,'exchange.npz'),('hydrolysis',a.hydrolysis_dir,'hydrolysis.npz')]},
        'intervals':{key:[ratio+(1-ratio)*v for v in values] for key,values in hm['intervals'].items()},
        'captions':caption_rows(em,hm,ratio),'core_atom_names':hm['core_atom_names'],'gamma_atom_names':hm['gamma_atom_names'],
        'old_gdp_atom_names':old['atom_names'].tolist(),'water_entry':{'start_frame':round(before*.78)+1,'end_frame':before+1,'distance_nm':14,'direction':direction.tolist(),'lateral_distance_nm':14,'lateral_direction':lateral.tolist(),'lateral_power':2},
        'gap_parking_translation':{'distance_nm':10,'direction':route['direction'],'entry_offset_zero_by_fraction':hm['intervals']['binding'][1]*.55,'exit_offset_release_phase':[.55,1],'scope':'Smooth rigid offset for off-screen arrival/departure, zero at docking; all combined contacts screened.'},
        'full_timeline_checks':checks,'covalent_length_ranges_nm':bonds,'inherited_gap_internal_checks':{k:v for k,v in hm['full_timeline_checks'].items() if k.startswith('gap_self_')},
        'heavy_atoms':63,'chemical_identity_scope':'GDP_A exits unchanged. GTP_B enters and becomes GDP_B plus Pi_B, whose fourth oxygen is water_O. Mg_A and Mg_B remain distinct.',
        'chemistry_limits':hm['chemistry_limits']+['GEF/GAP catalytic constructs, not full-length membrane-associated proteins or a full signaling pathway.',
            'One finite Ras state cycle with open substrate/product pools; the ending GDP molecule differs from the starting GDP molecule.',
            'Water entry is an authored screened path, not measured solvent diffusion.'],
        'passed':all(v['maximum_overlap_atoms']==0 for v in checks.values()) and all(.10<v[0] and v[1]<.22 for v in bonds.values())}
    a.output.mkdir(parents=True);np.savez_compressed(a.output/'state-cycle.npz',**samples)
    (a.output/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    if not report['passed']:raise SystemExit('Combined timeline failed; inspect the complete partner screen')
    print(json.dumps({'passed':True,'frames':count,'heavy_atoms':63,'output':str(a.output)}),flush=True)


if __name__=='__main__':main()
