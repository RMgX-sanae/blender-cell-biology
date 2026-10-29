"""Scientific negative controls using freshly downloaded, qualified sources."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from prepare_ubiquitin_transfer import sample, audit
from ubiquitin_geometry import UbiquitinGeometry, index


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    if a.output.exists(): p.error('Use a new test receipt')
    m = UbiquitinGeometry(a.sources)
    config = json.loads((Path(__file__).resolve().parents[1]/'assets/ubiquitin-transfer.json').read_text())
    arrays, states, fractions = sample(m, config, np.linspace(0, 1, 121))
    result = audit(m, arrays, states, fractions); assert result['passed']
    r = m.rows; assert r['e2'][index(r['e2'], 87, 'SG')][1] == 'S'
    assert not any(v[3] == 87 and v[2] == 'NZ' for v in r['e2'])
    assert r['e2'][index(r['e2'], 92, 'NZ')][4] == 'LYS'
    product = m.pose(config['torsion_parameters'][-1]); valid = m.check(product, 'product')
    assert 1.30 < valid['C_N_A'] < 1.36
    assert min(valid['nonbonded_minima_A'].values()) >= 2.38
    # The short C--N and adjacent 1-3 distances are valid only for a formed bond.
    wrong_state = m.check(product, 'reactant')
    assert min(wrong_state['nonbonded_minima_A'].values()) < 1.36
    collision = {k: v.copy() for k, v in product.items()}
    collision['acceptor'][0] = collision['e2'][0]
    assert min(m.check(collision, 'product')['nonbonded_minima_A'].values()) < 1e-9
    wrong_site = {k: v.copy() for k, v in product.items()}
    k48 = index(r['acceptor'], 48, 'NZ')
    wrong_site['acceptor'][k48] = wrong_site['donor'][m.ids['c']]+[1.33, 0, 0]
    assert min(m.check(wrong_site, 'product')['nonbonded_minima_A'].values()) < 1.34
    # A stretched source bond must fail even if no partner is penetrated.
    stretched = {k: v.copy() for k, v in arrays.items()}
    stretched['donor'][80, index(r['donor'], 20, 'CA')] += [2, 0, 0]
    assert not audit(m, stretched, states, fractions)['passed']
    # Correct C--N length alone does not make a valid amide geometry.
    twisted = {k: v.copy() for k, v in arrays.items()}; frame = 100
    center = twisted['acceptor'][frame, m.ids['nz']]
    axis = center-twisted['donor'][frame, m.ids['c']]; axis /= np.linalg.norm(axis)
    twisted['acceptor'][frame] = (twisted['acceptor'][frame]-center)@Rotation.from_rotvec(axis*np.pi/2).as_matrix().T+center
    bad_amide = audit(m, twisted, states, fractions)
    assert any('amide_angles' in failure for failure in bad_amide['failures'])
    receipt = {'passed': True, 'samples': 121, 'negative_controls': ['wrong chemical state', 'new interprotein collision', 'wrong K48 reaction site', 'stretched peptide bond', 'nonplanar amide with correct bond length'],
        'native_C87_and_K92_identity': True, 'positive_geometry': result,
        'limits': 'Source-specific geometric controls; not evidence of physiological rates or an elementary reaction mechanism.'}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'passed': True, 'negative_controls': receipt['negative_controls']}))


if __name__ == '__main__': main()
