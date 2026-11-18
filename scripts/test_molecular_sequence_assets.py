"""Small technical fixtures for new-mechanism display handoffs; no biological models."""
import copy,hashlib,json,tempfile,unittest
from pathlib import Path
import numpy as np
from molecular_sequence_assets import coordinates_only, load_assets


class PreparedSequence(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.source={'passed':True,'checks':{'passed':True},'frames':2,'fps':24,'units':'nm',
                     'atom_ids':{'body':['body:C','body:O']},'data_sha256':{}}
        np.savez(self.root/'motion.npz',body=np.zeros((3,2,3)))
        np.savez(self.root/'bonds.npz',bonds=np.array([[0,1]]))
        for name in 'abc':np.savez(self.root/f'background_{name}.npz',vertices=np.zeros((3,3)),faces=[[0,1,2]])
        for path in self.root.iterdir():self.hash_file(path.name)
        actors={'body':{'representation':'sticks','bond_file':'bonds.npz'}}
        self.display={'version':1,'views':[{'id':'whole','actors':actors},{'id':'detail','actors':copy.deepcopy(actors)}],
                      'background_components':list('abc'),'highlights':[{'id':'site','component':'body','indices':[0]}],
                      'bond_annotations':[{'id':'link','atoms':[['body',0],['body',1]]}]}

    def hash_file(self,name):
        self.source['data_sha256'][name]=hashlib.sha256((self.root/name).read_bytes()).hexdigest()

    def load(self):
        (self.root/'manifest.json').write_text(json.dumps(self.source));return load_assets(self.root,self.display)

    def test_complete_inputs(self):
        source,data,receipt=self.load()
        self.assertEqual(data['body'].shape,(3,2,3));self.assertEqual(source,self.source)
        self.assertEqual(len(receipt['required_files']),5)

    def test_undrawn_coordinates_keep_the_same_source_atoms(self):
        self.display['views'][0]['actors']['body']={'representation':'coordinates'}
        source, data, _=self.load()
        self.assertEqual(source['atom_ids']['body'],['body:C','body:O'])
        self.assertEqual(data['body'].shape,(3,2,3))
        self.assertTrue(coordinates_only({'representation':'spheres','radius_nm':0},[]))
        self.assertFalse(coordinates_only({'representation':'sticks','radius_nm':0},[[0,1]]))

    def test_coordinates_cannot_declare_visible_geometry(self):
        for field,value in [('radius_nm',.05),('bond_radius_nm',.02),('bonds',[[0,1]]),('bond_file','bonds.npz')]:
            self.display['views'][0]['actors']['body']={'representation':'coordinates',field:value}
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'coordinates-only'):
                self.load()

    def test_consumed_files_cannot_be_omitted_from_hash_list(self):
        original=self.source['data_sha256'].copy()
        for name in original:
            self.source['data_sha256']=original.copy();del self.source['data_sha256'][name]
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,'Missing prepared file hashes'):self.load()

    def test_changed_motion_cannot_pass(self):
        np.savez(self.root/'motion.npz',body=np.ones((3,2,3)))
        with self.assertRaisesRegex(ValueError,'missing or changed'):self.load()

    def test_wrong_point_inventory_cannot_pass_even_with_current_hash(self):
        np.savez(self.root/'motion.npz',body=np.zeros((3,3,3)));self.hash_file('motion.npz')
        with self.assertRaisesRegex(ValueError,'coordinate sequence'):self.load()

    def test_negative_or_fractional_display_indices_cannot_wrap_or_truncate(self):
        for indices in [[-1],[2],[.5],[True],[]]:
            self.display['highlights'][0]['indices']=indices
            with self.subTest(indices=indices),self.assertRaisesRegex(ValueError,'integer atom indices'):self.load()

    def test_fractional_bond_file_cannot_be_cast_to_integer(self):
        np.savez(self.root/'bonds.npz',bonds=[[0.,1.5]]);self.hash_file('bonds.npz')
        with self.assertRaisesRegex(ValueError,'integer atom indices'):self.load()

    def test_bond_annotation_must_reference_existing_atoms(self):
        self.display['bond_annotations'][0]['atoms'][1][1]=-1
        with self.assertRaisesRegex(ValueError,'integer atom indices'):self.load()

    def test_duplicate_physical_identity_cannot_pass(self):
        self.source['atom_ids']['body'][1]='body:C'
        with self.assertRaisesRegex(ValueError,'unique'):self.load()

    def test_unqualified_or_wrong_units_cannot_pass(self):
        self.source['passed']='true'
        with self.assertRaisesRegex(ValueError,'Qualify the mechanism'):self.load()
        self.source['passed']=True;self.source['units']='angstrom'
        with self.assertRaisesRegex(ValueError,'must be in nm'):self.load()

    def test_legacy_units_require_an_explicit_conversion_note(self):
        del self.source['units'];self.display['units']='nm'
        with self.assertRaisesRegex(ValueError,'legacy prepared coordinates'):self.load()
        self.display['unit_note']='Technical fixture coordinates are already nm';self.load()

    def test_declared_ledger_cannot_be_omitted(self):
        self.source['reaction_ledger']={'passed':True}
        with self.assertRaisesRegex(ValueError,'reaction-ledger.json'):self.load()

    def test_prepared_config_caption_boundaries_are_not_ignored(self):
        self.source['config']={'captions':[[0,.5,'source before'],[.5,1,'source after']]}
        self.display['captions']=[[0,.5,'translated before'],[.5,1,'translated after']]
        self.load()
        self.display['captions'][0][1]=.6;self.display['captions'][1][0]=.6
        with self.assertRaisesRegex(ValueError,'caption intervals differ'):self.load()

    def test_product_aliases_preserve_initial_coordinate_identity(self):
        ledger={'version':1,'scope':'technical fixture','bond_scope':'technical fixture',
                'reactant_atoms':[['body:C','C'],['body:O','O']],
                'product_atoms':[['renamed:C','C'],['body:O','O']],
                'atom_map':[['body:C','renamed:C'],['body:O','body:O']],
                'reactant_bonds':[],'product_bonds':[],'break':[],'form':[]}
        path=self.root/'reaction-ledger.json';path.write_text(json.dumps(ledger));self.hash_file(path.name)
        self.assertTrue(self.load()[2]['ledger_checks'][path.name]['passed'])
        ledger['reactant_atoms'][0][0]='unknown:C';ledger['atom_map'][0][0]='unknown:C'
        path.write_text(json.dumps(ledger));self.hash_file(path.name)
        with self.assertRaisesRegex(ValueError,'initial physical atom IDs differ'):self.load()

    def test_display_rechecks_order_changes_in_hashed_assets(self):
        atoms = [['body:C', 'C'], ['body:O', 'O']]
        ledger = {'version': 2, 'scope': 'technical sites', 'bond_scope': 'one edge',
                  'reactant_atoms': atoms, 'product_atoms': atoms, 'atom_map': [[i, i] for i, _ in atoms],
                  'reactant_bonds': [['body:C', 'body:O', 'doub']], 'product_bonds': [['body:C', 'body:O', 'sing']],
                  'break': [], 'form': [], 'change_order': [['body:O', 'body:C', 'doub', 'sing']]}
        path = self.root/'reaction-ledger.json'
        path.write_text(json.dumps(ledger))
        self.hash_file(path.name)
        self.assertTrue(self.load()[2]['ledger_checks'][path.name]['passed'])
        ledger['product_bonds'][0][2] = 'doub'
        path.write_text(json.dumps(ledger))
        self.hash_file(path.name)
        with self.assertRaisesRegex(ValueError, 'Bond-order changes'):
            self.load()

    def pose_inputs(self):
        initial = [['body:C', 'C'], ['body:O', 'O']]
        middle = [['mid:C', 'C'], ['mid:O', 'O']]
        final = [['released:C', 'C'], ['released:O', 'O']]

        def stage(before, after, old, new):
            a, b = [row[0] for row in before]
            x, y = [row[0] for row in after]
            return dict(version=2, scope='technical', bond_scope='one edge', reactant_atoms=before, product_atoms=after,
                        atom_map=[[a, x], [b, y]], reactant_bonds=[[a, b, old]], product_bonds=[[y, x, new]],
                        change_order=[[a, b, old, new]] if old != new else [], form=[], **{'break': []})

        ledger = {'version': 2, 'stages': [stage(initial, middle, 'doub', 'sing'), stage(middle, final, 'sing', 'doub')],
                  'net': stage(initial, final, 'doub', 'doub')}
        name = 'reaction-sequence.json'
        (self.root/name).write_text(json.dumps(ledger))
        self.hash_file(name)
        xyz = np.zeros((3, 2, 3))
        xyz[:, 1, 0] = [.122, .14, .124]
        np.savez(self.root/'motion.npz', body=xyz)
        self.hash_file('motion.npz')
        poses = []
        for i, (atoms, order, window) in enumerate([(initial, 'doub', [1.21, 1.23]), (middle, 'sing', [1.39, 1.41]),
                                                    (final, 'doub', [1.23, 1.25])]):
            rule = dict(atoms=[atom for atom, _ in atoms], order=order, length_A=window, basis='fixed technical fixture')
            poses.append(dict(state=i, sample=i, bonds=[rule]))
        return dict(version=1, ledger=name, ledger_sha256=self.source['data_sha256'][name],
                    motion_sha256=self.source['data_sha256']['motion.npz'], poses=poses)

    def save_pose_review(self, review):
        path = self.root/'reaction-geometry.json'
        path.write_text(json.dumps(review))
        self.hash_file(path.name)

    def test_pose_checks_follow_product_aliases_and_use_angstrom(self):
        review = self.pose_inputs()
        review['poses'].reverse()
        self.save_pose_review(review)
        rows = self.load()[2]['reaction_geometry']['poses']
        self.assertEqual([row['bonds'][0]['length_A'] for row in rows], [1.22, 1.40, 1.24])
        self.assertEqual(rows[-1]['bonds'][0]['atoms'], ['released:C', 'released:O'])
        self.assertEqual(rows[-1]['bonds'][0]['physical_atoms'], ['body:C', 'body:O'])

    def test_pose_checks_bind_the_ledger_and_motion_revisions(self):
        review = self.pose_inputs()
        for key in ('ledger_sha256', 'motion_sha256', 'ledger'):
            self.save_pose_review(dict(review, **{key: 'stale'}))
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'ledger|revision'):
                self.load()
        self.save_pose_review(review)
        del self.source['data_sha256']['reaction-geometry.json']
        with self.assertRaisesRegex(ValueError, 'Missing prepared file hashes'):
            self.load()

    def test_pose_checks_refuse_wrong_sites_states_windows_or_order(self):
        original = self.pose_inputs()
        for target, key, value, error in [
                ('pose', 'state', True, 'Select state'), ('pose', 'state', 3, 'Select state'),
                ('pose', 'sample', 3, 'in-range'), ('pose', 'sample', 0, 'unique'),
                ('pose', 'bonds', [], 'at least one bond'),
                ('bond', 'order', 'doub', 'absent'), ('bond', 'atoms', ['body:C', 'body:O'], 'absent'),
                ('bond', 'length_A', [1.2, 1.3], 'outside'), ('bond', 'length_A', [1.6, 1.2], 'interval'),
                ('bond', 'length_A', [True, 2], 'interval'), ('bond', 'basis', '', 'basis')]:
            review = copy.deepcopy(original)
            pose = review['poses'][1]
            (pose if target == 'pose' else pose['bonds'][0])[key] = value
            self.save_pose_review(review)
            with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, error):
                self.load()
        review = copy.deepcopy(original)
        review['poses'][0]['sample'], review['poses'][1]['sample'] = 1, 0
        self.save_pose_review(review)
        with self.assertRaisesRegex(ValueError, 'backwards'):
            self.load()
        rule = original['poses'][1]['bonds'][0]
        original['poses'][1]['bonds'].append(dict(rule, atoms=rule['atoms'][::-1]))
        self.save_pose_review(original)
        with self.assertRaisesRegex(ValueError, 'Repeated bond'):
            self.load()


if __name__=='__main__':unittest.main()
