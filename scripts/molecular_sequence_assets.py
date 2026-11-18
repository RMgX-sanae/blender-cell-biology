"""Verify consumed files and fixed atom identities before a molecular display build."""
import hashlib,json
from pathlib import Path
import numpy as np
from reaction_ledger import validate,validate_sequence
from reaction_pose_geometry import check_pose_distances


def _indices(values,count,label,pairs=False):
    array=np.asarray(values)
    if pairs and array.size==0:return
    if (array.ndim!=(2 if pairs else 1) or not array.size
            or (pairs and array.shape[1]!=2) or not np.issubdtype(array.dtype,np.integer)
            or np.any(array<0) or np.any(array>=count)
            or (pairs and np.any(array[:,0]==array[:,1]))):
        raise ValueError('Require retained integer atom indices for '+label)


def display_edges(settings, folder):
    if 'bond_file' in settings:
        with np.load(folder/settings['bond_file']) as archive:
            return archive['bonds']
    return settings.get('bonds', [])


def coordinates_only(settings, edges):
    if settings['representation'] == 'coordinates':
        if {'bonds', 'bond_file', 'bond_radius_nm', 'radius_nm'} & settings.keys():
            raise ValueError('A coordinates-only display cannot specify drawn atoms or bonds')
        return True
    # Older display files used a zero atom radius and no bonds for the same view.
    return (settings['representation'] in {'sticks', 'spheres'}
            and settings.get('radius_nm') == 0 and len(edges) == 0)


def load_assets(folder,display):
    folder=Path(folder);source=json.loads((folder/'manifest.json').read_text())
    if source.get('passed') is not True or source.get('checks',{}).get('passed') is not True:
        raise ValueError('Qualify the mechanism before rendering')
    if any(type(source.get(key)) is not int or source[key]<=0 for key in ('frames','fps')):
        raise ValueError('Require positive integer frames and fps')
    if source.get('units') is None and not display.get('unit_note'):
        raise ValueError('Qualify the units of legacy prepared coordinates in the display specification')
    if source.get('units',display.get('units'))!='nm':raise ValueError('Prepared display coordinates must be in nm')
    views=display['views'];ids=source['atom_ids']
    if display['version']!=1 or len(views)!=2 or len({v['id'] for v in views})!=2:
        raise ValueError('Supply a version-1 display with two distinct views')
    if not ids or any(not values or any(not isinstance(v,str) or not v for v in values) for values in ids.values()):
        raise ValueError('Require nonempty physical atom identities for every component')
    physical_ids=[v for values in ids.values() for v in values]
    if len(physical_ids)!=len(set(physical_ids)):raise ValueError('Atom identities must be unique across components')
    if set(views[0]['actors'])!=set(ids):raise ValueError('The overview must retain every prepared component')
    if any(set(view['actors'])-set(ids) for view in views):raise ValueError('Display refers to an absent component')
    captions=source.get('captions',source.get('config',{}).get('captions'))
    if captions is not None and [r[:2] for r in display['captions']]!=[r[:2] for r in captions]:
        raise ValueError('Display caption intervals differ from prepared mechanism events')
    if len(display['background_components'])!=3:raise ValueError('This background layout requires three source shape entries')
    required={'motion.npz',*('background_'+name+'.npz' for name in display['background_components'])}
    required.update(settings['bond_file'] for view in views for settings in view['actors'].values() if 'bond_file' in settings)
    ledgers=[('reaction-sequence.json','reaction_sequence',validate_sequence),('reaction-ledger.json','reaction_ledger',validate)]
    for name,field,_ in ledgers:
        if field in source or (folder/name).exists():required.add(name)
    geometry_file = 'reaction-geometry.json'
    if 'reaction_geometry' in source or (folder/geometry_file).exists():
        required.add(geometry_file)
    hashes=source.get('data_sha256',{})
    if not required<=set(hashes):raise ValueError('Missing prepared file hashes: '+', '.join(sorted(required-set(hashes))))
    for name,expected in hashes.items():
        if Path(name).name!=name:raise ValueError('Prepared files must be directly inside the asset directory')
        path=folder/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Prepared data missing or changed: '+name)
    with np.load(folder/'motion.npz') as archive:data=dict(archive)
    for key,identities in ids.items():
        if key not in data or data[key].shape!=(source['frames']+1,len(identities),3) or not np.isfinite(data[key]).all():
            raise ValueError('Invalid prepared coordinate sequence: '+key)
    for view in views:
        for key,settings in view['actors'].items():
            edges=display_edges(settings,folder)
            _indices(edges,len(ids[key]),key+' display edges',pairs=True)
            coordinates_only(settings,edges)
    for marker in display.get('highlights',[]):
        _indices(marker['indices'],len(ids[marker['component']]),marker['id'])
    for marker in display.get('bond_annotations',[]):
        if len(marker['atoms'])!=2:raise ValueError('A bond annotation requires two atoms')
        for key,index in marker['atoms']:_indices([index],len(ids[key]),marker['id'])
    ledger_checks={}
    ledger_specs={}
    for name,_,check in ledgers:
        if name not in hashes:continue
        spec=json.loads((folder/name).read_text());ledger_checks[name]=check(spec)
        ledger_specs[name]=spec
        first=spec['stages'][0] if name=='reaction-sequence.json' else spec
        # Coordinate IDs stay in the initial reactant basis. The checked atom_map
        # carries product aliases (for example a water oxygen becoming OXT).
        if {identity for identity,_ in first['reactant_atoms']}!=set(physical_ids):
            raise ValueError('Ledger and display initial physical atom IDs differ')
    receipt = {'verified_files':sorted(hashes),'required_files':sorted(required),'ledger_checks':ledger_checks}
    if geometry_file in hashes:
        review = json.loads((folder/geometry_file).read_text())
        name = review['ledger']
        if name not in ledger_specs:
            raise ValueError('Reaction geometry must name a verified reaction ledger or sequence')
        if review['ledger_sha256'] != hashes[name] or review['motion_sha256'] != hashes['motion.npz']:
            raise ValueError('Reaction geometry names a different ledger or motion revision')
        poses = check_pose_distances(review, ledger_specs[name], ids, data, source['frames'])
        receipt['reaction_geometry'] = dict(ledger=name, ledger_sha256=hashes[name],
                                            motion_sha256=hashes['motion.npz'], poses=poses)
    return source,data,receipt
