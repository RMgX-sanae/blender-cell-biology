"""Perturb real prepared assets to exercise biological/geometry failure checks."""
import argparse
import copy
import json
from pathlib import Path
import tempfile

import numpy as np
from audit_clathrin_growth import audit,digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sources','config','motion','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();record=json.loads((a.motion/'manifest.json').read_text());data=dict(np.load(a.motion/'motion.npz'));checks=[]
    def change(case,m,d):
        if case=='missing-modeled-light-chain':m['atom_ids'][m['config']['moving_unit']]=m['atom_ids'][m['config']['moving_unit']][:-70]
        elif case=='wrong-hinge-pivot':m['joints'][0]['pivot_nm'][0]+=.5
        elif case=='unclosed-motion':d['moving_ca'][-1,0,0]+=.05
        elif case=='changed-source-bound-endpoint':d['moving_ca'][m['frames']//2,0,0]+=.05
        elif case=='global-protein-stretch':d['rotations'][0]*=1.1
        elif case=='new-rigid-body-collision':
            shift=d['obstacle_ca'][0]-d['moving_ca'][0,0]
            for i in (0,-1):d['moving_ca'][i]+=shift;d['centers'][i]+=shift
        elif case=='wrong-local-deformation':d['moving_ca'][1,0,0]+=.1
    cases={
        'missing-modeled-light-chain':'Source identities/order changed',
        'wrong-hinge-pivot':'Hinge pivot does not follow raw CA',
        'unclosed-motion':'Motion does not close',
        'changed-source-bound-endpoint':'Bound endpoint lost source coordinates',
        'global-protein-stretch':'Global scale/reflection detected',
        'new-rigid-body-collision':'New/worsened contact outside exact source-pair allowance',
        'wrong-local-deformation':'Actual deformation differs from rigid-segment source motion'}
    for case,diagnostic in cases.items():
        m=copy.deepcopy(record);d={k:v.copy() for k,v in data.items()};change(case,m,d)
        with tempfile.TemporaryDirectory(prefix='clathrin-negative-') as directory:
            folder=Path(directory);(folder/'path.npz').symlink_to((a.motion/'path.npz').resolve())
            np.savez_compressed(folder/'motion.npz',**d);m['data_sha256']['motion.npz']=digest(folder/'motion.npz')
            (folder/'manifest.json').write_text(json.dumps(m))
            try:audit(a.sources,a.config,folder)
            except ValueError as error:
                if str(error)!=diagnostic:raise AssertionError(f'{case}: unexpected rejection {error}') from error
                checks.append({'case':case,'rejected':True,'diagnostic':str(error)})
            else:raise AssertionError('Accepted '+case)
    result={'passed':True,'negative_cases':len(checks),'checks':checks}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
