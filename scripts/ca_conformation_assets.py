"""Load a qualified fixed-identity CA display and preserve all path samples."""
import hashlib,json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sample_frames(points, presentation):
    if presentation.get('mode')!='continuous':raise ValueError('This finite conformation view uses continuous presentation, not a loop')
    duration=presentation['frames'];fps=presentation['fps']
    if any(type(v) is not int or v<=0 for v in (duration,fps)):raise ValueError('Require positive integer frames and fps')
    travel=np.linalg.norm(np.diff(points,axis=0),axis=2).max(1)
    if not travel.sum()>0:raise ValueError('A conformation path must contain motion')
    knots=np.r_[0,np.cumsum(travel)];knots/=knots[-1]
    # Keep every solver knot distinct, including nearly stationary refinements.
    knots=.9*knots+.1*np.linspace(0,1,len(knots))
    u=np.linspace(0,1,10001);ease=u*u*u*(10+u*(-15+6*u))
    return 1+duration*np.interp(knots,ease,u)


def load_assets(folder):
    folder=Path(folder);manifest=json.loads((folder/'manifest.json').read_text())
    if (manifest.get('version')!=1 or manifest.get('units')!='nm' or manifest.get('interpolation')!='linear_ca'
            or manifest.get('passed') is not True or manifest['path_checks'].get('passed') is not True):
        raise ValueError('Qualify the source-bound linear CA path before rendering')
    if digest(folder/'motion.npz')!=manifest['data_sha256']['motion.npz']:raise ValueError('Prepared path changed')
    with np.load(folder/'motion.npz') as archive:data=dict(archive)
    points=data['points'];ids=manifest['point_ids'];keys=data['keys']
    if (points.shape!=(len(points),len(ids),3) or not np.isfinite(points).all() or len(ids)!=len(set(ids))
            or keys.shape!=(len(ids),2) or len({tuple(k) for k in keys})!=len(ids)):
        raise ValueError('Prepared CA coordinates must retain unique fixed identities')
    if not np.array_equal(data['frames'],sample_frames(points,manifest['config']['presentation'])):
        raise ValueError('Prepared sample timing differs from presentation settings')
    if manifest['chains']!=list(dict.fromkeys(keys[:,0])):raise ValueError('Prepared chain inventory changed')
    return manifest,data
