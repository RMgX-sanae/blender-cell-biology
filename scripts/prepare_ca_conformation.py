"""Qualify an authored linear CA path against two reviewed mmCIF/SIFTS sources.

The path is an explicit input, not predicted from endpoints. All observed CAs of
the selected first-state chains are retained. Unmatched endpoint regions and
omitted companions require an explicit project explanation.
"""
import argparse,json
from pathlib import Path
import gemmi
import numpy as np
from source_residue_mapping import correspondence,digest
from ca_path_constraints import inspect
from ca_conformation_assets import sample_frames


def source_edges(path, keys):
    """Use polymer label sequence positions, preserving author insertion codes."""
    index={tuple(k):i for i,k in enumerate(keys)};edges=[];seen=set()
    for chain in gemmi.read_structure(str(path))[0]:
        rows=[]
        for residue in chain:
            key=(chain.name,str(residue.seqid))
            if key not in index:continue
            if key in seen:raise ValueError('Repeated source connectivity identity')
            seen.add(key)
            if residue.label_seq is None:raise ValueError('Source polymer lacks label sequence connectivity')
            rows.append((residue.label_seq,index[key]))
        rows.sort()
        if len({n for n,_ in rows})!=len(rows):raise ValueError('Ambiguous source label sequence position')
        edges.extend((i,j) for (a,i),(b,j) in zip(rows,rows[1:]) if b==a+1)
    if seen!=set(index):raise ValueError('Connectivity source is missing retained identities')
    return np.array(edges,int).reshape(-1,2)


def qualified_endpoints(sources, spec):
    """Return independently mapped endpoints; unobserved target CAs stay NaN."""
    if spec.get('version')!=1 or spec.get('interpolation')!='linear_ca' or not spec.get('scope'):
        raise ValueError('Declare version 1, linear_ca interpolation and the biological scope')
    files={}
    for role in ('first','second'):
        for kind in ('cif','sifts'):
            row=spec['sources'][role][kind];file=sources/row['file']
            if digest(file)!=row['sha256']:raise ValueError('Source input changed: '+str(file))
            files[role+'_'+kind]=file
        if gemmi.cif.read(str(files[role+'_cif'])).sole_block().find_value('_entry.id')!=spec['sources'][role]['entry']:
            raise ValueError('Declared source entry differs from the mmCIF identity')
    observed,mapping=correspondence(**files,chains=spec['chains'],insertions=spec.get('insertions',[]))
    for role in ('first','second'):
        declared=spec['unmatched'][role]
        if declared['residues']!=mapping['unmatched'][role] or (declared['residues'] and not declared.get('reason')):
            raise ValueError('Unmatched '+role+' endpoint regions differ from the reviewed selection')
    omitted={}
    for role in ('first','second'):
        selected={r[role] for r in spec['chains']}
        omitted[role]=[c.name for c in gemmi.read_structure(str(files[role+'_cif']))[0] if c.name not in selected]
        declaration=spec['omitted_chains'][role]
        if declaration['chains']!=omitted[role] or (omitted[role] and not declaration.get('reason')):
            raise ValueError('Unselected source chains need an explicit scope and companion review')
    keys=list(observed[0]);lookup={key:i for i,key in enumerate(keys)}
    fit_rows=[r for r in mapping['pairs'] if r['first'][0] in spec['alignment_chains']]
    if set(spec['alignment_chains'])-{r['first'][0] for r in fit_rows} or len(fit_rows)<3:
        raise ValueError('Alignment chains must have at least three observed corresponding CAs')
    x=np.array([observed[0][tuple(r['first'])]['xyz_nm'] for r in fit_rows])
    y=np.array([observed[1][tuple(r['second'])]['xyz_nm'] for r in fit_rows])
    a,b=x.mean(0),y.mean(0)
    if np.linalg.matrix_rank(x-a)<2 or np.linalg.matrix_rank(y-b)<2:raise ValueError('Degenerate alignment')
    u,_,v=np.linalg.svd((x-a).T@(y-b));rotation=v.T@np.diag([1,1,np.linalg.det(v.T@u.T)])@u.T
    translation=b-rotation@a
    first=np.array([r['xyz_nm'] for r in observed[0].values()])@rotation.T+translation
    second=np.full_like(first,np.nan);covered=np.zeros(len(keys),bool)
    for row in mapping['pairs']:
        i=lookup[tuple(row['first'])];second[i]=observed[1][tuple(row['second'])]['xyz_nm'];covered[i]=True
    edges=source_edges(files['first_cif'],keys);policy=json.loads(json.dumps(spec['geometry']))
    for gap in policy.get('unmodeled_gaps',[]):
        if not gap.get('reason'):raise ValueError('Missing-region endpoint limit needs a reviewed rationale')
        gap['indices']=[lookup[tuple(key)] for key in gap.pop('authors')]
    ids=[f"{spec['sources']['first']['entry']}:{c}:{r}:{observed[0][c,r]['name']}:CA" for c,r in keys]
    manifest={'version':1,'units':'nm','interpolation':'linear_ca','scope':spec['scope'],'config':spec,
              'point_ids':ids,'chains':list(dict.fromkeys(k[0] for k in keys)),
              'endpoint_checks':None,'path_checks':None,'mapping':mapping,'omitted_chains':omitted,
              'first_to_second':{'rotation':rotation.tolist(),'translation_nm':translation.tolist()},
              'geometry':policy,'passed':False,
              'limits':'Authored fixed-inventory CA path, not a predicted all-atom transition or a full activation mechanism'}
    return manifest,dict(keys=np.array(keys),edges=edges,first=first,second=second,second_observed=covered)


def prepare(sources, path, spec):
    manifest,source=qualified_endpoints(sources,spec)
    if digest(path)!=spec['path_sha256']:raise ValueError('Authored path changed; review the new candidate explicitly')
    with np.load(path) as data:
        points=np.asarray(data['points'],float);keys=data['keys']
    if (not np.array_equal(keys,source['keys']) or points.shape!=(len(points),len(keys),3)
            or len(points)<2 or not np.isfinite(points).all()):
        raise ValueError('Candidate path must retain every selected first-source CA identity, in source order')
    first=source['first'];second=source['second'].copy();covered=source['second_observed'];edges=source['edges']
    second[~covered]=points[-1,~covered]
    errors={'first_max_nm':float(np.linalg.norm(points[0]-first,axis=1).max()),
            'second_observed_max_nm':float(np.linalg.norm(points[-1,covered]-second[covered],axis=1).max())}
    if errors['first_max_nm']>1e-6 or errors['second_observed_max_nm']>1e-5:
        raise ValueError('Path endpoints differ from source coordinates: '+str(errors))
    checks=inspect(points,edges,manifest['geometry'],[(first,np.ones(len(keys),bool)),(second,covered)])
    if not checks['passed']:raise ValueError('Candidate CA path fails geometry: '+json.dumps(checks))
    times=sample_frames(points,spec['presentation'])
    manifest.update(endpoint_checks=errors,path_checks=checks,passed=True)
    return manifest,dict(points=points,keys=np.array(keys),edges=edges,frames=times,first=first,second=second,second_observed=covered)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--path',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh prepared directory')
    manifest,data=prepare(a.sources,a.path,json.loads(a.config.read_text()))
    a.output.mkdir(parents=True);np.savez_compressed(a.output/'motion.npz',**data)
    manifest['data_sha256']={'motion.npz':digest(a.output/'motion.npz')}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':True,'samples':len(data['points']),'points':len(manifest['point_ids']),'path_checks':manifest['path_checks']}))


if __name__=='__main__':main()
