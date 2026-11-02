"""Check endpoint atom correspondence and declared covalent graph changes.

The caller supplies scientifically qualified endpoints and the expected reaction.
This checks bookkeeping, not coordinates, bond order, protonation or mechanism.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def atoms(rows, label):
    result = {}
    for identity, element in rows:
        if not isinstance(identity, str) or not identity or not isinstance(element, str) or not element:
            raise ValueError(label+': atom IDs and elements must be nonempty strings')
        if identity in result: raise ValueError(label+': repeated atom '+identity)
        result[identity] = element
    if not result: raise ValueError(label+': no atoms')
    return result


def bonds(rows, identities, label):
    result = set()
    for first, second in rows:
        if first not in identities or second not in identities:
            raise ValueError(label+': bond refers to an absent atom')
        if first == second: raise ValueError(label+': self bond')
        edge = tuple(sorted((first, second)))
        if edge in result: raise ValueError(label+': repeated bond')
        result.add(edge)
    return result


def validate(spec):
    if spec['version'] != 1: raise ValueError('Unsupported reaction ledger version')
    before = atoms(spec['reactant_atoms'], 'reactant')
    after = atoms(spec['product_atoms'], 'product')
    mapping, reverse = {}, {}
    for first, second in spec['atom_map']:
        if first in mapping or second in reverse: raise ValueError('Atom map is not one-to-one')
        mapping[first] = second; reverse[second] = first
    if set(mapping) != set(before) or set(reverse) != set(after):
        raise ValueError('Atom map does not cover both complete endpoint atom sets')
    for first, second in mapping.items():
        if before[first] != after[second]:
            raise ValueError('Element changed across atom map: '+first+' -> '+second)
    initial = bonds(spec['reactant_bonds'], before, 'reactant')
    product = bonds(spec['product_bonds'], after, 'product')
    final = {tuple(sorted((reverse[a], reverse[b]))) for a, b in product}
    broken = bonds(spec['break'], before, 'expected break')
    formed = bonds(spec['form'], before, 'expected form')
    if broken & formed: raise ValueError('A net bond change cannot both break and form the same edge')
    if initial-final != broken:
        raise ValueError('Broken bonds do not match the declared reaction')
    if final-initial != formed:
        raise ValueError('Formed bonds do not match the declared reaction')
    return {'passed': True, 'scope': spec['scope'], 'bond_scope': spec['bond_scope'],
            'reactant_atoms': len(before), 'product_atoms': len(after), 'mapped_atoms': len(mapping),
            'element_counts': dict(sorted(Counter(before.values()).items())),
            'atom_map_sha256': hashlib.sha256(json.dumps(sorted(mapping.items())).encode()).hexdigest(),
            'broken_bonds_in_reactant_ids': sorted(broken), 'formed_bonds_in_reactant_ids': sorted(formed),
            'unchanged_declared_bonds': len(initial & final),
            'limits': 'Only the supplied atoms and bond subgraph are checked. A consistent declaration can still be biologically wrong; source suitability, omitted atoms, bond order/charge, geometry, intermediate motion and rendered visibility need separate checks.'}


def _endpoint(spec, side):
    identities = atoms(spec[side+'_atoms'], side)
    return identities, bonds(spec[side+'_bonds'], identities, side)


def validate_sequence(spec):
    """Check a fixed atom inventory through stages against an independent net ledger."""
    if spec['version'] != 1: raise ValueError('Unsupported reaction sequence version')
    stages = spec['stages']
    if len(stages) < 2: raise ValueError('A reaction sequence needs at least two stages')
    checks = [validate(stage) for stage in stages]
    net = spec['net']; net_check = validate(net)
    mapping = dict(stages[0]['atom_map'])
    for index, (first, second) in enumerate(zip(stages, stages[1:]), 1):
        before, before_bonds = _endpoint(first, 'product')
        after, after_bonds = _endpoint(second, 'reactant')
        if before != after:
            raise ValueError(f'Handoff {index}: atom identities/elements differ')
        if before_bonds != after_bonds:
            raise ValueError(f'Handoff {index}: intermediate bond graph differs')
        next_map = dict(second['atom_map'])
        mapping = {start: next_map[current] for start, current in mapping.items()}
    if _endpoint(stages[0], 'reactant') != _endpoint(net, 'reactant'):
        raise ValueError('Net reactant state differs from the first stage')
    if _endpoint(stages[-1], 'product') != _endpoint(net, 'product'):
        raise ValueError('Net product state differs from the last stage')
    if mapping != dict(net['atom_map']):
        raise ValueError('Composed stage atom map differs from the independent net map')
    return {'passed': True, 'stages': checks, 'handoffs_checked': len(stages)-1,
            'mapped_atoms': len(mapping), 'net': net_check,
            'composed_atom_map_sha256': net_check['atom_map_sha256'],
            'limits': 'Fixed declared atom inventory and bond subgraph only. Shared IDs must mean the same atoms. Include waiting/released molecules in every stage; this does not validate open-pool exchange, coordinates, protonation, reaction order evidence or a closed catalytic cycle.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('ledger', type=Path); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sequence', action='store_true', help='Check stages, exact handoffs and an independent net ledger')
    a = p.parse_args()
    if a.output.exists(): p.error('Use a fresh audit output')
    try:
        check = validate_sequence if a.sequence else validate
        result = check(json.loads(a.ledger.read_text()))
    except (KeyError, TypeError, ValueError) as error:
        result = {'passed': False, 'error': str(error)}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False))
    if not result['passed']: raise SystemExit(1)


if __name__ == '__main__': main()
