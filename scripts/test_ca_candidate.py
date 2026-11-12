"""Behavioral tests for a reviewed source-specific candidate configuration."""
import argparse,copy,json,sys,unittest
from pathlib import Path
import numpy as np
from ca_candidate_model import CandidateModel
from ca_candidate_constraints import Constraints,traces
from ca_internal_coordinates import parameters,reconstruct,basis
from ca_path_constraints import inspect
from generate_ca_candidate import validate_stages


class CandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config=json.loads(cls.config_path.read_text());cls.source=json.loads((cls.config_path.parent/cls.config['source_config']).read_text())
        cls.model=CandidateModel(cls.sources,cls.source,cls.config)

    def build(self,change):
        config=copy.deepcopy(self.config);change(config)
        return CandidateModel(self.sources,self.source,config)

    def test_goals_keep_both_qualified_endpoints(self):
        np.testing.assert_allclose(self.model.goal(0)[0],self.model.source,rtol=0,atol=1e-10)
        np.testing.assert_allclose(self.model.goal(1)[0],self.model.target,rtol=0,atol=1e-10)

    def test_virtual_backbone_reconstructs_nonplanar_fragment(self):
        points=np.array([[0,0,0],[.38,0,0],[.38,.38,0],[.38,.38,.38],[.76,.38,.38]])
        np.testing.assert_allclose(reconstruct(points[:3],parameters(points)),points,atol=1e-12)

    def test_collinear_anchors_are_not_silently_repaired(self):
        with self.assertRaisesRegex(ValueError,'Collinear'):basis(np.zeros(3),np.ones(3),np.ones(3)*2)

    def test_unassigned_partition_rejected(self):
        with self.assertRaisesRegex(ValueError,'Every retained'):self.build(lambda c:c['partitions'].pop(0))

    def test_overlapping_partitions_rejected(self):
        def change(c):c['partitions'][0]['select'][0]['label_ranges']=[[340,460]]
        with self.assertRaisesRegex(ValueError,'disjoint'):self.build(change)

    def test_inference_cannot_replace_observed_target(self):
        def change(c):c['inferred_targets'][0]['select'][0]['label_ranges']=[[135,141]]
        with self.assertRaisesRegex(ValueError,'Infer only'):self.build(change)

    def test_unexplained_missing_target_rejected(self):
        with self.assertRaisesRegex(ValueError,'explicit inference'):self.build(lambda c:c.update(inferred_targets=[]))

    def test_terminal_cannot_cross_source_gap(self):
        def change(c):c['terminals'][1]['label_range']=[160,130]
        with self.assertRaisesRegex(ValueError,'every source label'):self.build(change)

    def test_trace_order_is_connectivity_not_array_index(self):
        self.assertEqual(traces([[0,2],[2,1]],4),[[0,2,1],[3]])

    def test_branched_and_cyclic_inputs_need_another_solver(self):
        with self.assertRaisesRegex(ValueError,'unbranched'):traces([[0,1],[0,2],[0,3]],4)
        with self.assertRaisesRegex(ValueError,'Cyclic'):traces([[0,1],[1,2],[2,0]],3)

    def test_objective_gradient_against_finite_differences(self):
        model=self.model;constraints=Constraints(model);goal,phase=model.goal(.0025)
        loss=constraints.objective(model.source,goal,phase,1.)
        x=model.source.ravel()+np.sin(np.arange(model.source.size))*.0002
        _,gradient=loss(x);h=1e-6
        for index in np.linspace(0,len(x)-1,24).astype(int):
            first=x.copy();second=x.copy();first[index]+=h;second[index]-=h
            observed=(loss(first)[0]-loss(second)[0])/(2*h)
            self.assertAlmostEqual(float(gradient[index]),float(observed),delta=2e-4)

    def test_actual_bounded_step_passes_independent_geometry(self):
        model=self.model;goal,phase=model.goal(.0025);point,row=Constraints(model).step(model.source,goal,phase)
        check=inspect(np.array([model.source,point]),model.edges,model.policy,
                      [(model.source,np.ones(len(point),bool)),(model.target,model.mapped)])
        self.assertTrue(check['passed']);self.assertLessEqual(row['step_nm'],np.sqrt(3)*model.config['solver']['step_bound_nm']+1e-8)

    def test_negative_restraint_weight_rejected(self):
        model=self.build(lambda c:c['solver'].update(fold_weight=-1))
        with self.assertRaisesRegex(ValueError,'solver weights'):Constraints(model)

    def test_schedule_requires_finite_guided_then_settle(self):
        config=copy.deepcopy(self.config);config['stages'][1]['kind']='guided'
        with self.assertRaisesRegex(ValueError,'finite guided'):validate_stages(config)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--sources',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    args=parser.parse_args();CandidateTests.sources=args.sources;CandidateTests.config_path=args.config
    unittest.main(argv=[sys.argv[0],'-v'])
