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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('ledger', type=Path); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists(): p.error('Use a fresh audit output')
    try:
        result = validate(json.loads(a.ledger.read_text()))
    except (KeyError, TypeError, ValueError) as error:
        result = {'passed': False, 'error': str(error)}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False))
    if not result['passed']: raise SystemExit(1)


if __name__ == '__main__': main()
