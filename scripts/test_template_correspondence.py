"""Technical tests of shared CA correspondence, independent of a protein example."""
import copy,unittest
import numpy as np
from prepare_binding_candidates import fit_template,transform
from prepare_source_states import alignment_maps


class CorrespondenceChecks(unittest.TestCase):
    def setUp(self):
        points=np.array([[0.,0,0],[1,0,0],[0,2,0],[1,1,3]])
        self.matrix=np.eye(4);self.matrix[:3,:3]=[[0,-1,0],[1,0,0],[0,0,1]];self.matrix[:3,3]=[8,-3,2]
        names=['ALA','GLY','SER','ASP'];identities=[(1,''),(184,'A'),(184,''),(5,'')]
        self.first=[(point,'C','CA',number,name,code) for point,(number,code),name in zip(points,identities,names)]
        self.second=[(point,'C','CA',number+100,name,code) for point,(_,_,_,number,name,code) in zip(transform(points,self.matrix),self.first)]
        self.mapping={'residue_pairs':[[1,101],['184A','284A'],[184,284],[5,105]],'expected_substitutions':[],
                      'min_CA_coverage':1.,'max_CA_rmsd_angstrom':.001}

    def test_explicit_numbering_and_insertion_codes_preserve_rigid_fit(self):
        matrix,report=fit_template(self.first,self.second,self.mapping)
        np.testing.assert_allclose(matrix,self.matrix,atol=1e-12,rtol=0)
        self.assertEqual(report['CA_pairs'],self.mapping['residue_pairs'])
        self.assertEqual(report['source_CA_coverage'],1.)

    def test_float_bool_and_string_integer_cannot_alias_numbered_sites(self):
        for value in [1.,True,'1']:
            mapping=copy.deepcopy(self.mapping);mapping['residue_pairs'][0][0]=value
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'integer author numbers'):
                fit_template(self.first,self.second,mapping)

    def test_invalid_gates_cannot_disable_qualification(self):
        for field,value in [('min_CA_coverage',float('nan')),('min_CA_coverage',-1),('min_CA_coverage',1.1),
                            ('max_CA_rmsd_angstrom',float('nan')),('max_CA_rmsd_angstrom',float('inf')),
                            ('max_CA_rmsd_angstrom',-1),('max_CA_rmsd_angstrom',True)]:
            mapping={**self.mapping,field:value}
            with self.subTest(field=field,value=value),self.assertRaisesRegex(ValueError,'finite CA coverage'):
                fit_template(self.first,self.second,mapping)

    def test_duplicate_ca_is_not_silently_overwritten(self):
        with self.assertRaisesRegex(ValueError,'unique, finite'):
            fit_template(self.first+[self.first[0]],self.second,self.mapping)

    def test_missing_insertion_is_not_collapsed_to_neighbor(self):
        mapping=copy.deepcopy(self.mapping);mapping['residue_pairs'][1][0]=184
        with self.assertRaisesRegex(ValueError,'unique, present'):fit_template(self.first,self.second,mapping)

    def test_nonfinite_ca_is_rejected(self):
        rows=copy.deepcopy(self.first);rows[0][0][0]=float('nan')
        with self.assertRaisesRegex(ValueError,'unique, finite'):fit_template(rows,self.second,self.mapping)

    def test_substitutions_require_the_reviewed_identity(self):
        second=copy.deepcopy(self.second);row=list(second[0]);row[4]='VAL';second[0]=tuple(row)
        with self.assertRaisesRegex(ValueError,'substitutions'):fit_template(self.first,second,self.mapping)
        mapping=copy.deepcopy(self.mapping);mapping['expected_substitutions']=[{'source_residue':1,'template_residue':101,'source_name':'ALA','template_name':'VAL'}]
        matrix,_=fit_template(self.first,second,mapping);np.testing.assert_allclose(matrix,self.matrix,atol=1e-12)

    def test_sparse_map_still_obeys_coverage_gate(self):
        mapping={**self.mapping,'residue_pairs':self.mapping['residue_pairs'][:3]}
        with self.assertRaisesRegex(ValueError,'coverage/RMSD'):fit_template(self.first,self.second,mapping)

    def test_state_maps_can_use_different_author_chains(self):
        config={'reference':'second','states':[{'pdb':'first','chains':[{'chain':'S'}]},
                                              {'pdb':'second','chains':[{'chain':'R'}]}],
                'alignment':{'reference_chain':'R','maps':{
                    'first':{**self.mapping,'chain':'S','note':'Reviewed technical renumbering'},
                    'second':{**self.mapping,'chain':'R','note':'Fixed reference',
                              'residue_pairs':[[b,b] for a,b in self.mapping['residue_pairs']]}}}}
        chain,maps=alignment_maps(config);self.assertEqual(chain,'R')
        fit_template(self.first,self.second,maps['first'])
        for key,change in [('omitted state',lambda c:c['alignment']['maps'].pop('first')),
                           ('wrong chain',lambda c:c['alignment']['maps']['first'].update(chain='R')),
                           ('unreviewed map',lambda c:c['alignment']['maps']['first'].update(note='')),
                           ('conflicting formats',lambda c:c['alignment'].update(chain='R')),
                           ('moved reference identity',lambda c:c['alignment']['maps']['second']['residue_pairs'][0].__setitem__(0,1))]:
            changed=copy.deepcopy(config);change(changed)
            with self.subTest(case=key),self.assertRaises(ValueError):alignment_maps(changed)


if __name__=='__main__':unittest.main()
