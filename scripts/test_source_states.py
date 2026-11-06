"""Targeted failures for a reviewed source-state comparison, using real sources."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np

from audit_source_states import check_sources, digest
from prepare_source_states import qualify
from prepare_binding_candidates import fit_template


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    config = json.loads(a.config.read_text())
    report = check_sources(a.asset_dir, a.sources, a.config)
    rejected = []
    def rejects(name, run):
        try: run()
        except ValueError as error: rejected.append({'case': name, 'error': str(error)})
        else: raise AssertionError('Invalid input accepted: '+name)
    # This fixture intentionally exercises two proteins sharing author chain jj,
    # an AAQ catalytic mutant and an explicitly retained nucleotide analogue.
    by_id = {s['pdb']: s for s in config['states']}
    for name, code, change in [
        ('wrong source hash', '5LZT', lambda s: s.update(sha256='0'*64)),
        ('wrong assembly', '5LZT', lambda s: s.update(assembly_id='2')),
        ('same chain wrong partner', '5LZT', lambda s: s['chains'][1].update(description='ABCE1')),
        ('wrong canonical protein', '5LZT', lambda s: s['chains'][1].update(uniprot='G1SG72')),
        ('expression host mistaken for protein species', '5LZT', lambda s: s['chains'][0].update(source_organism='Escherichia coli')),
        ('AAQ mislabeled native GGQ', '5LZV', lambda s: s['chains'][0]['sites'][0].update(name='GLY')),
        ('unresolved site', '5LZU', lambda s: s['chains'][0]['sites'][0].update(residue=999)),
        ('unresolved site atom', '5LZU', lambda s: s['chains'][0]['sites'][0].update(atoms=['FAKE'])),
        ('wrong stop-codon base', '5LZU', lambda s: s['chains'][-1]['sites'][0].update(name='G')),
        ('analogue mislabeled GTP', '5LZT', lambda s: s['ligand_instances'][1].update(component='GTP')),
        ('wrong Fe-S instance', '5LZV', lambda s: s['ligand_instances'][1].update(residue_number=602)),
        ('duplicate polymer', '5LZU', lambda s: s['chains'].append(copy.deepcopy(s['chains'][0]))),
    ]:
        state = copy.deepcopy(by_id[code]); change(state)
        rejects(name, lambda: qualify(a.sources, state))
    qualified = qualify(a.sources, by_id['5LZU'])[2]['ii']
    for name, pairs in [('missing fit residue', [[6,6],[7,7],[999,999]]),
                        ('duplicate fit residue', [[6,6],[6,7],[8,8]])]:
        rejects(name, lambda: fit_template(qualified, qualified, {'residue_pairs': pairs,
            'expected_substitutions': [], 'min_CA_coverage': 0, 'max_CA_rmsd_angstrom': 2}))
    with tempfile.TemporaryDirectory(prefix='source-state-negatives-') as tmp:
        root = Path(tmp)
        for name in ['changed atom identity', 'shifted polymer', 'dropped Fe-S cluster', 'shifted cofactor', 'stretched state']:
            target = root/name.replace(' ', '-')
            shutil.copytree(a.asset_dir, target)
            m = json.loads((target/'manifest.json').read_text())
            if name in ('changed atom identity', 'shifted polymer'):
                file = m['states'][0]['actors'][0]['file']
                with np.load(target/file, allow_pickle=False) as old:
                    d = {k: old[k].copy() for k in old.files}
                if name == 'changed atom identity': d['residue_numbers'][0] += 1
                else: d['atoms'][:,0] += .1
                np.savez_compressed(target/file, **d)
                m['data_sha256'][file] = digest(target/file)
            elif name == 'dropped Fe-S cluster': m['states'][2]['ligands'].pop()
            elif name == 'shifted cofactor': m['states'][0]['ligands'][1]['positions_nm'][0][0] += .1
            else: m['states'][0]['alignment_matrix_angstrom'][0][0] *= 1.1
            (target/'manifest.json').write_text(json.dumps(m))
            rejects(name, lambda: check_sources(target, a.sources, a.config))
    report['negative_cases'] = rejected
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': True, 'negative_cases': len(rejected)}))


if __name__ == '__main__': main()
