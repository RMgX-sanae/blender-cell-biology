"""Generate a bounded source-guided CA candidate from explicit project constraints.

No cached trajectory is read for a fresh solve. The result is an authored
geometric candidate, not measured dynamics or a general protein-motion model.
"""
import argparse,datetime,json,os,platform,signal,sys
from pathlib import Path
import numpy as np
import scipy
from ca_candidate_model import CandidateModel
from ca_candidate_constraints import Constraints
from ca_path_constraints import inspect
from source_residue_mapping import digest
from prepare_ca_conformation import prepare


def validate_stages(config):
    stages=config['stages']
    if not stages or stages[0].get('kind')!='guided' or any(s['kind']!='settle' for s in stages[1:]):
        raise ValueError('Use one finite guided stage followed by explicit endpoint-settling stages')
    for row in stages:
        if type(row['steps']) is not int or row['steps']<1 or not np.isfinite(row['weight']) or row['weight']<=0:
            raise ValueError('Stages need positive step counts and target weights')
        if row['kind']=='settle':
            if (any(type(row[k]) is not int or not 0<=row[k]<=row['steps'] for k in ('minimum_steps','stagnation_after'))
                    or any(not np.isfinite(row[k]) or row[k]<=0 for k in ('rms_stop_nm','step_stop_nm'))):
                raise ValueError('Invalid finite endpoint-settling criteria')
    if not np.isfinite(config['maximum_final_correction_nm']) or config['maximum_final_correction_nm']<=0:
        raise ValueError('Declare a positive numerical endpoint-correction bound')


def binding(config_path,source_path):
    root=Path(__file__).resolve().parent
    files=['generate_ca_candidate.py','ca_candidate_model.py','ca_candidate_constraints.py','ca_internal_coordinates.py',
           'prepare_ca_conformation.py','ca_path_constraints.py','source_residue_mapping.py','ca_conformation_assets.py']
    return {'candidate_config_sha256':digest(config_path),'source_config_sha256':digest(source_path),
            'code_sha256':{name:digest(root/name) for name in files},
            'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'machine':platform.machine()}}


def checkpoint(folder,points,keys,state,identity,status):
    # A single archive holds the path and the matching stage cursor together.
    temporary=folder/'checkpoint.tmp.npz'
    np.savez_compressed(temporary,points=np.array(points),keys=keys,metadata=json.dumps({'binding':identity,**state}))
    temporary.replace(folder/'checkpoint.npz')
    report={'status':status,'passed':False,'pid':os.getpid(),'updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage':state['stage'],'next_step':state['step'],'samples':len(points),'last_step':state['rows'][-1] if state['rows'] else None}
    (folder/'run-state.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


def solve(model,folder,identity,resume=False,stop_after=None):
    config=model.config;validate_stages(config);constraints=Constraints(model)
    observed=[(model.source,np.ones(len(model.source),bool)),(model.target,model.mapped)]
    for endpoint in (model.source,model.target):
        review=inspect(np.array([endpoint,endpoint]),model.edges,model.policy,observed)
        if not review['passed']:raise ValueError('Reviewed endpoint or inferred fragment fails CA geometry: '+str(review))
    state={'stage':0,'step':0,'rows':[]};points=[model.source.copy()]
    if resume:
        with np.load(folder/'checkpoint.npz') as data:
            restored=json.loads(str(data['metadata']));points=list(data['points']);keys=data['keys']
        if restored.pop('binding')!=identity:raise ValueError('Checkpoint source/config/code/runtime binding changed')
        if not np.array_equal(keys,model.keys) or not np.array_equal(points[0],model.source):
            raise ValueError('Checkpoint source identities or first coordinates changed')
        if len(restored['rows'])!=len(points)-1:raise ValueError('Checkpoint step records do not match the retained path')
        if not 0<=restored['stage']<=len(config['stages']):raise ValueError('Checkpoint stage cursor is invalid')
        if restored['stage']<len(config['stages']) and not 0<=restored['step']<=config['stages'][restored['stage']]['steps']:
            raise ValueError('Checkpoint step cursor is invalid')
        if len(points)>1 and not inspect(np.array(points),model.edges,model.policy,observed)['passed']:
            raise ValueError('Checkpoint prefix fails independent geometry')
        state=restored
    else:checkpoint(folder,points,model.keys,state,identity,'running')
    stopping=[False]
    def request_stop(*_):stopping[0]=True
    old_handlers={sig:signal.signal(sig,request_stop) for sig in (signal.SIGTERM,signal.SIGINT)}
    new_steps=0
    try:
        while state['stage']<len(config['stages']):
            stage=config['stages'][state['stage']];done=state['step']>=stage['steps']
            if stage['kind']=='settle' and state['rows']:
                last=state['rows'][-1]
                done|=(state['step']>=stage['minimum_steps'] and last['rms_to_goal_nm']<stage['rms_stop_nm'])
                done|=(state['step']>=stage['stagnation_after'] and last['step_nm']<stage['step_stop_nm'])
            if done:
                state['stage']+=1;state['step']=0;checkpoint(folder,points,model.keys,state,identity,'running');continue
            if stopping[0] or (stop_after is not None and new_steps>=stop_after):
                checkpoint(folder,points,model.keys,state,identity,'checkpointed');return None
            q=(state['step']+1)/stage['steps'] if stage['kind']=='guided' else 1.
            goal,phase=model.goal(q) if stage['kind']=='guided' else (model.target,np.ones(len(model.names)))
            result,row=constraints.step(points[-1],goal,phase,stage['weight'])
            row.update(stage=state['stage'],step=state['step'],q=q,weight=stage['weight'])
            points.append(result);state['rows'].append(row);state['step']+=1;new_steps+=1
            if new_steps%10==0:checkpoint(folder,points,model.keys,state,identity,'running')
        # Only a tiny, independently checked final correction is allowed.
        delta=float(np.linalg.norm(points[-1]-model.target,axis=1).max())
        if delta>config['maximum_final_correction_nm']:
            raise ValueError('Candidate has not reached the source-derived endpoint; refusing a reset: '+str(delta))
        final=inspect(np.array([points[-1],model.target]),model.edges,model.policy,observed)
        if not final['passed']:raise ValueError('Numerical endpoint correction fails geometry')
        candidate=np.array(points+[model.target.copy()]);review=inspect(candidate,model.edges,model.policy,observed)
        if not review['passed']:raise ValueError('Complete candidate fails independent linear-path checks')
        checkpoint(folder,points,model.keys,state,identity,'solved_pending_source_qualification')
        return candidate,{'path_checks':review,'final_correction_nm':delta,'steps':state['rows'],'binding':identity,'fits':model.fits}
    except Exception:
        checkpoint(folder,points,model.keys,state,identity,'failed');raise
    finally:
        for sig,handler in old_handlers.items():signal.signal(sig,handler)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true')
    p.add_argument('--stop-after',type=int,help='Checkpoint after N new steps; a partial result is not a qualified candidate')
    a=p.parse_args()
    if a.stop_after is not None and a.stop_after<1:p.error('Use a positive new-step limit')
    if a.output.exists() and not a.resume:p.error('Use a fresh output directory or explicitly resume its checkpoint')
    if a.resume and not (a.output/'checkpoint.npz').is_file():p.error('No checkpoint exists at the requested output')
    if (a.output/'candidate-review.json').exists():p.error('Candidate is already qualified; use a fresh directory for another solve')
    config=json.loads(a.config.read_text());source_path=a.config.parent/config['source_config']
    source_spec=json.loads(source_path.read_text());identity=binding(a.config,source_path)
    model=CandidateModel(a.sources,source_spec,config);a.output.mkdir(parents=True,exist_ok=a.resume)
    result=solve(model,a.output,identity,a.resume,a.stop_after)
    if result is None:return
    candidate,review=result
    np.savez_compressed(a.output/'motion.npz',points=candidate,keys=model.keys)
    # Bind the newly checked candidate rather than inheriting a previous path's hash.
    source_spec['path_sha256']=digest(a.output/'motion.npz')
    (a.output/'conformation.json').write_text(json.dumps(source_spec,ensure_ascii=False,indent=2)+'\n')
    manifest,_=prepare(a.sources,a.output/'motion.npz',source_spec)
    review.update(passed=True,endpoint_checks=manifest['endpoint_checks'],scope=source_spec['scope'],
                  raw_source_hashes=manifest['mapping']['sha256'],candidate_sha256=source_spec['path_sha256'],cached_trajectory_read=False)
    (a.output/'candidate-review.json').write_text(json.dumps(review,indent=2)+'\n')
    (a.output/'run-state.json').write_text(json.dumps({'status':'qualified','passed':True,'samples':len(candidate),
        'candidate_sha256':source_spec['path_sha256'],'scope':source_spec['scope']},indent=2)+'\n')
    print(json.dumps({'status':'qualified','samples':len(candidate),'points':len(model.keys),'review':str(a.output/'candidate-review.json')}),flush=True)


if __name__=='__main__':main()
