"""Small technical fixtures for new-mechanism display handoffs; no biological models."""
import copy,hashlib,json,tempfile,unittest
from pathlib import Path
import numpy as np
from molecular_sequence_assets import load_assets


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


if __name__=='__main__':unittest.main()
