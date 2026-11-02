"""Project selections and rigid source frames; biological qualification is separate."""
import copy
import unittest
import numpy as np
from prepare_molecular_asset import project_selection, xyz


class ProjectAssetTests(unittest.TestCase):
    def setUp(self):
        self.table = {'A': [([1.,2.,3.], 'C', 'CA', 10, 'ALA'),
                            ([2.,4.,6.], 'N', 'N', 11, 'GLY')],
                      'A-2': [([3.,2.,1.], 'C', 'CA', 10, 'ALA')]}
        self.config = {'id':'project-a', 'source_id':'1J4N', 'source_sha256':'a'*64,
            'chains':[{'chain':'A'}, {'chain':'A-2'}], 'selection_note':'Test fixture only',
            'orientation':{'origin_angstrom':[1.,1.,1.],
                'basis_columns':[[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]],
                'note':'Test frame only', 'evidence_urls':['https://www.rcsb.org/structure/1J4N']}}

    def select(self, config=None):
        return project_selection(self.table, config or self.config, 'a'*64)

    def test_selection_preserves_atoms_order_and_rigid_distances(self):
        selected,origin,basis = self.select()
        self.assertEqual(list(selected), ['A','A-2'])
        self.assertEqual(selected, self.table)
        points = np.concatenate([xyz(rows) for rows in selected.values()])
        local = (points-origin)@basis
        self.assertTrue(np.allclose(np.linalg.norm(points[:,None]-points,axis=-1),
                                    np.linalg.norm(local[:,None]-local,axis=-1)))

    def test_inclusive_ranges_do_not_duplicate_atoms(self):
        self.config['chains'] = [{'chain':'A','residue_ranges':[[10,10],[10,10]]}]
        selected,_,_ = self.select()
        self.assertEqual(selected['A'], self.table['A'][:1])

    def test_wrong_source_cannot_reuse_frame(self):
        self.config['source_sha256'] = 'b'*64
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'): self.select()

    def test_reflection_stretch_shear_and_nonfinite_frames_fail(self):
        for matrix in (np.diag([1,1,-1]), np.diag([1,1,1.01]),
                       [[1,.01,0],[0,1,0],[0,0,1]], np.full((3,3),np.nan)):
            config=copy.deepcopy(self.config);config['orientation']['basis_columns']=matrix
            with self.subTest(matrix=matrix), self.assertRaisesRegex(ValueError, 'proper rigid'): self.select(config)
        self.config['orientation']['origin_angstrom'][0]=np.inf
        with self.assertRaisesRegex(ValueError, 'finite'): self.select()

    def test_missing_duplicate_and_unsafe_chains_fail(self):
        for chains in ([{'chain':'B'}], [{'chain':'A'},{'chain':'A'}], [{'chain':'../A'}], []):
            config=copy.deepcopy(self.config);config['chains']=chains
            with self.subTest(chains=chains), self.assertRaises(ValueError): self.select(config)

    def test_empty_or_malformed_residue_selection_fails(self):
        for ranges in ([[100,110]], [], [[11,10]], [[10.,11.]], [[10,11,12]]):
            config=copy.deepcopy(self.config);config['chains']=[{'chain':'A','residue_ranges':ranges}]
            with self.subTest(ranges=ranges), self.assertRaises(ValueError): self.select(config)

    def test_matrix_without_documented_evidence_is_insufficient(self):
        self.config['orientation']['evidence_urls']=[]
        with self.assertRaisesRegex(ValueError, 'evidence'): self.select()


if __name__ == '__main__': unittest.main()
