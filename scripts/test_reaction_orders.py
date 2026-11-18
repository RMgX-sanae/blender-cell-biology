"""Technical bond-graph regressions for single reactions and stage handoffs."""
import unittest

from reaction_ledger import bonds, validate, validate_sequence


ATOMS = [('c', 'C'), ('o', 'O'), ('in', 'O'), ('out', 'O'), ('anchor', 'C')]
INITIAL = [('c', 'o', 'doub'), ('c', 'out', 'sing'), ('anchor', 'c', 'sing')]
MIDDLE = [('o', 'c', 'sing'), ('c', 'in', 'sing'), ('c', 'out', 'sing'), ('anchor', 'c', 'sing')]
FINAL = [('o', 'c', 'doub'), ('in', 'c', 'sing'), ('anchor', 'c', 'sing')]


def stage(initial, final, *, broken=(), formed=(), changed=()):
    return {'version': 2, 'scope': 'five technical sites', 'bond_scope': 'declared covalent edges',
            'reactant_atoms': ATOMS, 'product_atoms': ATOMS, 'atom_map': [(i, i) for i, _ in ATOMS],
            'reactant_bonds': list(initial), 'product_bonds': list(final),
            'break': list(broken), 'form': list(formed), 'change_order': list(changed)}


def sequence():
    add = stage(INITIAL, MIDDLE, formed=[('in', 'c', 'SING')], changed=[('c', 'o', 'DOUB', 'SING')])
    leave = stage(MIDDLE, FINAL, broken=[('out', 'c', 'sing')], changed=[('o', 'c', 'sing', 'doub')])
    net = stage(INITIAL, FINAL, broken=[('out', 'c', 'sing')], formed=[('c', 'in', 'sing')])
    return {'version': 2, 'stages': [add, leave], 'net': net}


class BondOrderChecks(unittest.TestCase):
    def test_transient_change_and_net_cancellation(self):
        receipt = validate_sequence(sequence())
        first, second = receipt['stages']
        self.assertEqual(first['changed_bond_orders_in_reactant_ids'], [('c', 'o', 'doub', 'sing')])
        self.assertEqual(second['changed_bond_orders_in_reactant_ids'], [('c', 'o', 'sing', 'doub')])
        self.assertEqual(first['formed_bonds_in_reactant_ids'], [('c', 'in', 'sing')])
        self.assertEqual(second['broken_bonds_in_reactant_ids'], [('c', 'out', 'sing')])
        self.assertEqual(first['unchanged_declared_bonds'], 2)
        self.assertEqual(receipt['net']['changed_bond_orders_in_reactant_ids'], [])
        self.assertEqual(first['product_bonds_sha256'], second['reactant_bonds_sha256'])

    def test_product_names_can_change(self):
        spec = sequence()['stages'][0]
        spec['product_atoms'] = [('p/' + atom, element) for atom, element in ATOMS]
        spec['atom_map'] = [(atom, 'p/' + atom) for atom, _ in ATOMS]
        spec['product_bonds'] = [('p/' + b, 'p/' + a, order.upper()) for a, b, order in reversed(MIDDLE)]
        self.assertEqual(validate(spec)['changed_bond_orders_in_reactant_ids'], [('c', 'o', 'doub', 'sing')])

    def test_graph_fingerprint_uses_canonical_edges_and_orders(self):
        endpoint = stage(INITIAL, INITIAL)
        fingerprint = validate(endpoint)['product_bonds_sha256']
        endpoint['product_bonds'] = [(b, a, o.upper()) for a, b, o in reversed(INITIAL)]
        self.assertEqual(validate(endpoint)['product_bonds_sha256'], fingerprint)
        endpoint['product_bonds'][0] = ('c', 'anchor', 'delo')
        endpoint['change_order'] = [('c', 'anchor', 'SING', 'DELO')]
        self.assertNotEqual(validate(endpoint)['product_bonds_sha256'], fingerprint)

    def test_wrong_orders_and_unreported_changes_fail(self):
        for field, replacement in [
                ('product_bonds', [('c', 'o', 'doub'), *MIDDLE[1:]]),
                ('change_order', []), ('change_order', [('c', 'o', 'sing', 'doub')]),
                ('change_order', [('c', 'o', 'doub', 'doub')]),
                ('form', [('in', 'c', 'doub')])]:
            with self.subTest(field=field, replacement=replacement):
                spec = sequence()['stages'][0]
                spec[field] = replacement
                with self.assertRaisesRegex(ValueError, 'do not match'):
                    validate(spec)
        spec = sequence()['stages'][1]
        spec['break'] = [('c', 'out', 'trip')]
        with self.assertRaisesRegex(ValueError, 'Broken bonds'):
            validate(spec)

    def test_ambiguous_or_invalid_rows_fail(self):
        for row in [('c', 'o'), ('c', 'o', '?'), ('c', 'o', 2), ('c', 'o', None),
                    ('c', 'o', 'sing', 'doub'), ('o', 'o', 'sing'), ('o', 'absent', 'sing')]:
            with self.subTest(row=row):
                spec = sequence()['stages'][0]
                spec['reactant_bonds'][0] = row
                with self.assertRaises(ValueError):
                    validate(spec)
        for field, duplicate in [('reactant_bonds', ('o', 'c', 'sing')),
                                 ('change_order', ('o', 'c', 'sing', 'doub'))]:
            spec = sequence()['stages'][0]
            spec[field].append(duplicate)
            with self.assertRaisesRegex(ValueError, 'repeated bond'):
                validate(spec)

    def test_retained_order_change_cannot_be_disguised_as_break_and_form(self):
        spec = sequence()['stages'][0]
        spec['break'] = [('c', 'o', 'doub')]
        spec['form'].append(('c', 'o', 'sing'))
        spec['change_order'] = []
        with self.assertRaisesRegex(ValueError, 'break and form'):
            validate(spec)

    def test_equal_connectivity_at_handoff_is_insufficient(self):
        spec = sequence()
        second = spec['stages'][1]
        second['reactant_bonds'][0] = ('c', 'o', 'doub')
        second['change_order'] = []
        self.assertTrue(validate(second)['passed'])
        with self.assertRaisesRegex(ValueError, 'Handoff 1.*intermediate bond graph'):
            validate_sequence(spec)

    def test_independent_net_must_match_endpoint_orders(self):
        for side in ('reactant', 'product'):
            spec = sequence()
            net = spec['net']
            net[side + '_bonds'][0] = ('o', 'c', 'sing')
            orders = ('sing', 'doub') if side == 'reactant' else ('doub', 'sing')
            net['change_order'] = [('c', 'o', *orders)]
            self.assertTrue(validate(net)['passed'])
            with self.assertRaisesRegex(ValueError, 'Net ' + side + ' state differs'):
                validate_sequence(spec)

    def test_schema_versions_cannot_drop_orders(self):
        spec = sequence()
        for entry in [spec['stages'][0], spec['net']]:
            entry['version'] = 1
            with self.assertRaisesRegex(ValueError, 'same version'):
                validate_sequence(spec)
            entry['version'] = 2
        stage_one = dict(spec['stages'][0])
        del stage_one['change_order']
        with self.assertRaisesRegex(ValueError, 'requires change_order'):
            validate(stage_one)
        with self.assertRaisesRegex(ValueError, 'version-2'):
            validate(dict(spec['stages'][0], version=1))
        self.assertEqual(bonds([('b', 'a')], {'a', 'b'}, 'legacy'), {('a', 'b')})

    def test_bond_labels_follow_the_deposited_representation(self):
        for order in ('arom', 'delo', 'doub', 'pi', 'poly', 'quad', 'sing', 'trip'):
            with self.subTest(order=order):
                same = stage([('c', 'o', order.upper())], [('o', 'c', order)])
                self.assertEqual(validate(same)['unchanged_declared_bonds'], 1)
        spec = stage([('c', 'o', 'arom')], [('c', 'o', 'delo')])
        with self.assertRaisesRegex(ValueError, 'Bond-order changes'):
            validate(spec)


if __name__ == '__main__':
    unittest.main()
