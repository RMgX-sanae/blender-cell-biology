"""Check reviewed polymer site relations in one identity-operation assembly.

This checks source identities, CCD atom types, distances and connection records.
It neither infers chemistry from distance nor certifies a biological mechanism.
"""
import argparse
import json
from pathlib import Path

import gemmi
import numpy as np

from prepare_binding_candidates import residue_id, sha
from prepare_phosphotransfer import ccd_file
from prepare_source_states import qualify, rows


def connection_key(record, partner):
    prefix = f'ptnr{partner}_'
    number = record[prefix+'auth_seq_id']
    insertion = record.get(f'pdbx_ptnr{partner}_PDB_ins_code') or ''
    return (record[prefix+'auth_asym_id'], str(number)+insertion,
            record[prefix+'label_comp_id'], record[prefix+'label_atom_id'])


def inspect(sources, spec):
    if not spec.get('interpretation') or not spec.get('evidence'):
        raise ValueError('Explain the intended relation and its evidence')
    path, provenance, table, chain_info, _, assembly = qualify(sources, spec['source'])
    block = gemmi.cif.read(str(path)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    if len(structure) != 1:
        raise ValueError('Select and qualify a single model before checking site links')
    components = {}
    for code, expected in spec['ccd_sha256'].items():
        ccd, meta = ccd_file(sources, code)
        if meta['sha256'] != expected:
            raise ValueError('Requalify the chemical component: '+code)
        atom_types = {r['atom_id']: r['type_symbol'] for r in rows(ccd, '_chem_comp_atom.')}
        bonds = {frozenset((r['atom_id_1'], r['atom_id_2'])): r['value_order']
                 for r in rows(ccd, '_chem_comp_bond.')}
        components[code] = (atom_types, bonds, meta)
    inventories = []
    seen_residues = set()
    for item in spec.get('residue_inventories', []):
        key = (item['chain'], item['residue'])
        absent = item['expected_absent']
        if (key in seen_residues or not item.get('note') or not isinstance(absent, dict)
                or any(not isinstance(reason, str) or not reason.strip() for reason in absent.values())):
            raise ValueError('Require a unique residue inventory and explanations for expected omissions')
        seen_residues.add(key)
        selected = [r for r in table.get(item['chain'], []) if residue_id(r) == item['residue']]
        if not selected or {r[4] for r in selected} != {item['component']}:
            raise ValueError('Residue inventory source identity differs: '+str(key))
        expected = {name: gemmi.Element(element).name for name, element in components[item['component']][0].items()
                    if element not in ('H', 'D')}
        if not set(absent) <= set(expected):
            raise ValueError('Expected omission is not a CCD heavy atom: '+str(key))
        actual = {r[2]: r[1] for r in selected}
        wanted = {name: element for name, element in expected.items() if name not in absent}
        errors = []
        if actual != wanted:
            errors.append('Actual heavy-atom inventory differs from CCD minus the reviewed omissions')
        inventories.append({**item, 'actual_heavy_atoms': actual,
                            'missing_from_ccd': sorted(set(expected)-set(actual)),
                            'unexpected_atoms': sorted(set(actual)-set(expected)),
                            'element_mismatches': sorted(name for name in actual.keys() & expected.keys()
                                                         if actual[name] != expected[name]),
                            'passed': not errors, 'errors': errors})
    sites = {}
    for name, site in spec['sites'].items():
        selected = [r for r in table.get(site['chain'], [])
                    if residue_id(r) == site['residue'] and r[2] == site['atom']]
        if len(selected) != 1:
            raise ValueError('Missing or ambiguous retained site: '+name)
        row = selected[0]
        if row[4] != site['component'] or row[1] != site['element']:
            raise ValueError('Reviewed site component/element differs: '+name)
        if components[site['component']][0].get(site['atom']) != site['element']:
            raise ValueError('CCD atom identity/element differs: '+name)
        actual = [(r, a) for c in structure[0] if c.name == site['chain'] for r in c
                  if str(r.seqid) == str(site['residue']) and r.name == site['component']
                  for a in r if a.name == site['atom'] and a.occ > 0 and a.altloc in ('\x00', 'A')]
        if len(actual) != 1: raise ValueError('Ambiguous source atom metadata: '+name)
        residue, atom = actual[0]
        sites[name] = {**site, 'position_angstrom': row[0],
                       'label_chain': residue.subchain, 'label_residue': residue.label_seq,
                       'altloc': atom.altloc.strip('\x00 '),
                       'key': (site['chain'], str(site['residue']), site['component'], site['atom'])}
    if not spec['relations'] or len({r['id'] for r in spec['relations']}) != len(spec['relations']):
        raise ValueError('Require uniquely named site relations')
    connections = rows(block, '_struct_conn.')
    relations = []
    for relation in spec['relations']:
        if relation['kind'] not in ('covalent', 'noncovalent'):
            raise ValueError('Choose a reviewed covalent or noncovalent relation')
        bounds = np.asarray(relation['distance_angstrom'], dtype=float)
        if bounds.shape != (2,) or not np.isfinite(bounds).all() or not 0 < bounds[0] < bounds[1]:
            raise ValueError('Require finite, positive, ordered distance bounds')
        if relation['annotation'] not in ('required', 'inspect'):
            raise ValueError('Choose required or inspect connection annotation')
        first, second = [sites[n] for n in relation['sites']]
        if first['key'] == second['key']:
            raise ValueError('A relation needs two different source atoms')
        distance = float(np.linalg.norm(np.asarray(first['position_angstrom'])-second['position_angstrom']))
        if not np.isfinite(distance):
            raise ValueError('Non-finite source coordinates')
        pair = {first['key'], second['key']}
        matched = [r for r in connections if {connection_key(r, 1), connection_key(r, 2)} == pair]
        # Different symmetry/alternate conformers do not describe these coordinates.
        direct, other = [], []
        actual_sites = {s['key']: s for s in (first, second)}
        for record in matched:
            retained = all(record.get(f'ptnr{i}_symmetry') == '1_555'
                           and (record.get(f'pdbx_ptnr{i}_label_alt_id') or '') == actual_sites[connection_key(record, i)]['altloc']
                           and record[f'ptnr{i}_label_asym_id'] == actual_sites[connection_key(record, i)]['label_chain']
                           and str(record[f'ptnr{i}_label_seq_id']) == str(actual_sites[connection_key(record, i)]['label_residue'])
                           for i in (1, 2))
            (direct if retained else other).append(record)
        ccd_order = None
        if first['key'][:3] == second['key'][:3]:
            ccd_order = components[first['component']][1].get(frozenset((first['atom'], second['atom'])))
        annotated_covalent = ccd_order is not None or any(r['conn_type_id'] in ('covale', 'disulf') for r in direct)
        errors = []
        if not bounds[0] <= distance <= bounds[1]:
            errors.append('Distance outside reviewed bounds')
        if relation['kind'] == 'noncovalent' and annotated_covalent:
            errors.append('Covalent source/CCD annotation conflicts with the reviewed noncovalent relation')
        if relation['kind'] == 'covalent' and relation['annotation'] == 'required' and not annotated_covalent:
            errors.append('Required covalent annotation is absent for these source atoms')
        if relation['kind'] == 'noncovalent' and relation['annotation'] == 'required' and not direct:
            errors.append('Required noncovalent connection annotation is absent')
        relations.append({**relation, 'measured_distance_angstrom': distance,
                          'ccd_bond_order': ccd_order, 'direct_connections': direct,
                          'other_symmetry_conformer_or_label_connections': other,
                          'passed': not errors, 'errors': errors})
    return {'passed': all(r['passed'] for r in relations+inventories), 'id': spec['id'],
            'interpretation': spec['interpretation'], 'evidence': spec['evidence'],
            'source': provenance, 'assembly': assembly, 'chains': chain_info,
            'ccd_sources': {code: data[2] for code, data in components.items()},
            'sites': sites, 'relations': relations, 'residue_inventories': inventories,
            'limits': 'Only selected heavy-atom sites in a single model, retained blank/A alternate conformers and an identity-operation assembly. No global contact audit, hydrogen-bond angle/protonation check, affinity, dynamics, reaction-path or experimental validation. A missing struct_conn record alone does not disprove a bond. Thresholds and relation kinds need independent review.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists(): p.error('Choose a fresh report path')
    try:
        cases = json.loads(a.config.read_text())['cases']
        matches = [s for s in cases if s['id'] == a.case]
        if len(matches) != 1: raise ValueError('Require one matching case')
        report = inspect(a.sources, matches[0])
    except (ValueError, KeyError, TypeError, OSError, RuntimeError) as error:
        report = {'passed': False, 'id': a.case, 'error': str(error)}
    report['config_sha256'] = sha(a.config)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'passed': report['passed'], 'case': a.case, 'output': str(a.output)}))
    if not report['passed']: raise SystemExit(1)


if __name__ == '__main__': main()
