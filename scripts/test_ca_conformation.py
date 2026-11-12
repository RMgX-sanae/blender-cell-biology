"""Small behavioral regressions for linear CA paths and source connectivity."""
import tempfile,unittest
from pathlib import Path
import gemmi
import numpy as np
from ca_path_constraints import inspect,neighborhoods
from ca_conformation_assets import sample_frames
from prepare_ca_conformation import source_edges


class ConformationTests(unittest.TestCase):
    def check_path(self, points, edges=(), **changes):
        points=np.array(points,float)
        policy={'nonlocal_guard_nm':.32,'source_contact_allowance_nm':.01,
                'adjacent_distance_nm':[.35,.41],'maximum_step_nm':5.,**changes}
        return inspect(points,edges,policy,[(points[0],np.ones(points.shape[1],bool))])

    def test_interior_crossing_despite_clear_endpoints(self):
        report=self.check_path([[[-1,0,0],[1,0,0]],[[1,0,0],[-1,0,0]]])
        self.assertFalse(report['passed']);self.assertLess(report['nonlocal_CA_margin_nm'],-.3)

    def test_translation_preserves_connected_shape(self):
        source=np.array([[0,0,0],[.38,0,0],[.76,0,0]])
        self.assertTrue(self.check_path([source,source+[0,1,0]],[[0,1],[1,2]])['passed'])

    def test_stretched_backbone(self):
        self.assertFalse(self.check_path([[[0,0,0],[.38,0,0]],[[0,0,0],[.6,0,0]]],[[0,1]])['passed'])

    def test_connected_points_cannot_swap_through_each_other(self):
        self.assertFalse(self.check_path([[[0,0,0],[.38,0,0]],[[.38,0,0],[0,0,0]]],[[0,1]])['passed'])

    def test_source_contact_allowance_is_pair_specific(self):
        self.assertTrue(self.check_path([[[0,0,0],[.2,0,0]],[[0,0,0],[.195,0,0]]])['passed'])
        self.assertFalse(self.check_path([[[0,0,0],[.2,0,0],[2,0,0]],[[0,0,0],[.2,0,0],[0,.2,0]]])['passed'])

    def test_inferred_target_does_not_authorize_contact(self):
        points=np.array([[[0,0,0],[1,0,0]],[[0,0,0],[.2,0,0]]])
        policy={'nonlocal_guard_nm':.32,'source_contact_allowance_nm':.01,'adjacent_distance_nm':[.35,.41],'maximum_step_nm':5.}
        self.assertFalse(inspect(points,[],policy,[(points[0],np.ones(2,bool)),(points[-1],np.array([True,False]))])['passed'])

    def test_halfspace_and_unmodeled_span(self):
        points=[[[0,0,0],[2,0,0]],[[0,0,-.1],[2,0,0]]]
        self.assertFalse(self.check_path(points,halfspace={'normal':[0,0,1],'offset_nm':0})['passed'])
        self.assertFalse(self.check_path(points,unmodeled_gaps=[{'indices':[0,1],'maximum_nm':1.9}])['passed'])

    def test_source_connectivity_uses_label_positions_not_author_numbers(self):
        structure=gemmi.Structure();model=gemmi.Model('1');chain=gemmi.Chain('A')
        for author,seq in [('42',1),('42A',2),('42B',3),('100',5)]:
            residue=gemmi.Residue();residue.name='ALA'
            residue.seqid=gemmi.SeqId(int(author.rstrip('AB')),author[-1] if author[-1] in 'AB' else ' ')
            residue.label_seq=seq;residue.subchain='A'
            atom=gemmi.Atom();atom.name='CA';atom.element=gemmi.Element('C');atom.pos=gemmi.Position(seq*3.8,0,0);residue.add_atom(atom);chain.add_residue(residue)
        model.add_chain(chain);structure.add_model(model)
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'insertions.cif';structure.make_mmcif_document().write_file(str(file))
            edges=source_edges(file,[('A',n) for n in ('42','42A','42B','100')])
        self.assertEqual(edges.tolist(),[[0,1],[1,2]])
        self.assertNotIn((2,3),neighborhoods(edges,4))

    def test_nonfinite_path_rejected(self):
        with self.assertRaises(ValueError):self.check_path([[[0,0,0]],[[float('nan'),0,0]]])

    def test_finite_presentation_preserves_stationary_refinements(self):
        points=np.array([[[0,0,0]],[[1,0,0]],[[1,0,0]],[[1.1,0,0]]])
        frames=sample_frames(points,{'mode':'continuous','frames':240,'fps':24})
        self.assertTrue(np.all(np.diff(frames)>0));self.assertEqual((frames[0],frames[-1]),(1,241))
        with self.assertRaises(ValueError):sample_frames(points,{'mode':'loop','frames':240,'fps':24})


if __name__=='__main__':unittest.main()
