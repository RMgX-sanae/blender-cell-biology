"""Check endpoint atom correspondence and declared covalent graph changes.

The caller supplies the endpoint graphs and independently expected changes.
Version 1 checks connectivity. Version 2 also checks explicit wwPDB bond orders.
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


CCD_ORDERS = frozenset({'arom', 'delo', 'doub', 'pi', 'poly', 'quad', 'sing', 'trip'})


def _bond_table(rows, identities, label, order_fields):
    result = {}
    if isinstance(rows, (str, bytes, dict)):
        raise ValueError(f'{label}: supply a list of bond rows')
    for row in rows:
        if isinstance(row, (str, bytes, dict)):
            raise ValueError(f'{label}: supply a row of atom IDs and bond orders')
        first, second, *orders = row
        if len(orders) != order_fields:
            raise ValueError(f'{label}: expected two atom IDs and {order_fields} bond-order fields')
        if not isinstance(first, str) or not isinstance(second, str) or first not in identities or second not in identities:
            raise ValueError(label+': bond refers to an absent atom')
        if first == second: raise ValueError(label+': self bond')
        edge = tuple(sorted((first, second)))
        if edge in result: raise ValueError(label+': repeated bond')
        if any(not isinstance(order, str) or order.lower() not in CCD_ORDERS for order in orders):
            raise ValueError(f'{label}: unsupported or unknown bond order on {edge}: {orders}')
        result[edge] = tuple(order.lower() for order in orders)
    return result


def bonds(rows, identities, label):
    return set(_bond_table(rows, identities, label, 0))


def _rows(table):
    return sorted((*edge, *orders) for edge, orders in table.items())


def _match(actual, declared, message, ordered):
    if actual == declared:
        return
    if ordered:
        edge = min(edge for edge in actual.keys() | declared.keys()
                   if actual.get(edge) != declared.get(edge))
        message += f': {edge}, endpoint graph {actual.get(edge)}, declared {declared.get(edge)}'
    raise ValueError(message)


def validate(spec):
    version = spec['version']
    if type(version) is not int or version not in (1, 2):
        raise ValueError('Unsupported reaction ledger version')
    ordered = version == 2
    if not ordered and 'change_order' in spec:
        raise ValueError('change_order requires a version-2 reaction ledger')
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
    fields = 1 if ordered else 0
    initial = _bond_table(spec['reactant_bonds'], before, 'reactant', fields)
    product = _bond_table(spec['product_bonds'], after, 'product', fields)
    final = {tuple(sorted((reverse[a], reverse[b]))): order for (a, b), order in product.items()}
    broken = _bond_table(spec['break'], before, 'expected break', fields)
    formed = _bond_table(spec['form'], before, 'expected form', fields)
    if broken.keys() & formed.keys():
        raise ValueError('A net bond change cannot both break and form the same edge')
    _match({edge: initial[edge] for edge in initial.keys() - final.keys()}, broken,
           'Broken bonds do not match the declared reaction', ordered)
    _match({edge: final[edge] for edge in final.keys() - initial.keys()}, formed,
           'Formed bonds do not match the declared reaction', ordered)
    retained = initial.keys() & final.keys()
    changed = {edge: initial[edge] + final[edge] for edge in retained if initial[edge] != final[edge]}
    if ordered:
        if 'change_order' not in spec:
            raise ValueError('Version 2 requires change_order, including an empty list for no order changes')
        declared = _bond_table(spec['change_order'], before, 'expected order change', 2)
        _match(changed, declared, 'Bond-order changes do not match the declared reaction', ordered)
    result = {'passed': True, 'scope': spec['scope'], 'bond_scope': spec['bond_scope'],
            'reactant_atoms': len(before), 'product_atoms': len(after), 'mapped_atoms': len(mapping),
            'element_counts': dict(sorted(Counter(before.values()).items())),
            'atom_map_sha256': hashlib.sha256(json.dumps(sorted(mapping.items())).encode()).hexdigest(),
            'broken_bonds_in_reactant_ids': _rows(broken), 'formed_bonds_in_reactant_ids': _rows(formed),
            'unchanged_declared_bonds': len(retained) - len(changed),
            'limits': 'Only the supplied atoms and bond subgraph are checked. A consistent declaration can still be biologically wrong; source suitability, omitted atoms, bond order/charge, geometry, intermediate motion and rendered visibility need separate checks.'}
    if ordered:
        result.update(version=2, changed_bond_orders_in_reactant_ids=_rows(changed))
        for side, table in [('reactant', initial), ('product', product)]:
            payload = json.dumps(_rows(table), separators=(',', ':')).encode()
            result[side + '_bonds_sha256'] = hashlib.sha256(payload).hexdigest()
        result['limits'] = ('Checks the supplied atom map, covalent subgraph and categorical bond orders. '
                            'Source assignments, protonation, charge, full valence, coordinates and motion require separate evidence.')
    return result


def _endpoint(spec, side):
    identities = atoms(spec[side+'_atoms'], side)
    fields = 1 if spec['version'] == 2 else 0
    return identities, _bond_table(spec[side+'_bonds'], identities, side, fields)


def validate_sequence(spec):
    """Check a fixed atom inventory through stages against an independent net ledger."""
    version = spec['version']
    if type(version) is not int or version not in (1, 2):
        raise ValueError('Unsupported reaction sequence version')
    stages = spec['stages']
    if len(stages) < 2: raise ValueError('A reaction sequence needs at least two stages')
    net = spec['net']
    if any(stage['version'] != version for stage in [*stages, net]):
        raise ValueError('Stages and net ledger must use the same version as the reaction sequence')
    checks = [validate(stage) for stage in stages]
    net_check = validate(net)
    mapping = dict(stages[0]['atom_map'])
    for index, (first, second) in enumerate(zip(stages, stages[1:]), 1):
        before, before_bonds = _endpoint(first, 'product')
        after, after_bonds = _endpoint(second, 'reactant')
        if before != after:
            raise ValueError(f'Handoff {index}: atom identities/elements differ')
        _match(before_bonds, after_bonds, f'Handoff {index}: intermediate bond graph differs', version == 2)
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
