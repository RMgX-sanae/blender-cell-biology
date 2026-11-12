"""Technical mmCIF fixtures: instance identity, composed operations, companions.

These invented coordinates test transfer behavior, not protein structures.
Operation composition follows the wwPDB dictionary's right-to-left convention.
"""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import gemmi
import numpy as np

from source_assembly import SourceAssembly,verify_reference


def fixture(path,expression='(1,2)',reference=False,mutate=None):
    doc=gemmi.cif.Document();block=doc.add_new_block('TEST')
    block.set_mmcif_category('_entity.',{'id':['1','2'],'type':['polymer','non-polymer'],
                                      'pdbx_description':['Technical polymer','Technical calcium']})
    block.set_mmcif_category('_entity_poly.',{'entity_id':['1'],'type':['polypeptide(L)'],'pdbx_strand_id':['A,B']})
    block.set_mmcif_category('_pdbx_struct_assembly.',{'id':['1'],'details':['author_defined_assembly'],
        'method_details':[None],'oligomeric_details':['dimeric'],'oligomeric_count':['2']})
    block.set_mmcif_category('_pdbx_struct_assembly_gen.',{'assembly_id':['1'],
        'oper_expression':[expression],'asym_id_list':['X,Y,Z']})
    matrices=[np.eye(4),np.eye(4)];matrices[0][0,3]=10
    matrices[1][:3,:3]=[[0,-1,0],[1,0,0],[0,0,1]]
    operators={'id':['1','2'],'type':['general operation','general operation']}
    for i in range(3):
        operators[f'vector[{i+1}]']=[str(m[i,3]) for m in matrices]
        for j in range(3):operators[f'matrix[{i+1}][{j+1}]']=[str(m[i,j]) for m in matrices]
    block.set_mmcif_category('_pdbx_struct_oper_list.',operators)
    sites=[]
    for label,auth,offset in [('X','A',0),('Y','B',6)]:
        for seq,(number,icode,point) in enumerate([(10,'?',(0,0,0)),(11,'A',(2,0,0)),(12,'?',(0,3,0))],1):
            sites.append(['ATOM','C','CA','ALA',label,'1',str(seq),icode,
                          np.asarray(point,float)+[offset,0,0],str(number),auth])
    sites.append(['HETATM','CA','CA','CA','Z','2','.','B',np.array([1.,1.,1.]),'99','A'])
    if reference:
        copies=[]
        for index in range(2):
            for site in sites:
                row=copy.deepcopy(site);x,y,z=row[8]
                row[8]=np.array([x+10,y,z]) if index==0 else np.array([-y,x,z])
                row[4]+=str(index+1);row[10]+=str(index+1);copies.append(row)
        sites=copies
    block.set_mmcif_category('_struct_asym.',{'id':list(dict.fromkeys(s[4] for s in sites)),
        'entity_id':[next(s[5] for s in sites if s[4]==label) for label in dict.fromkeys(s[4] for s in sites)]})
    columns=['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id',
        'label_entity_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','occupancy',
        'B_iso_or_equiv','auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id','pdbx_PDB_model_num']
    loop=block.init_loop('_atom_site.',columns)
    for index,s in enumerate(sites,1):
        loop.add_row([s[0],str(index),s[1],s[2],'.',s[3],s[4],s[5],s[6],s[7],
                      *map(str,s[8]),'1','10',s[9],s[3],s[10],s[2],'1'])
    if mutate:mutate(block)
    doc.write_file(str(path))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polymer(instance='A1',label='X',chain='A',operator='1'):
    return {'instance_id':instance,'label_chain':label,'chain':chain,'operator':operator}


def companion(instance='Ca1',operator='1',owner='A1'):
    return {'instance_id':instance,'label_chain':'Z','chain':'A','operator':operator,
            'residue_number':99,'insertion_code':'B','component':'CA','owner_instance':owner}


def qualified_fixture(directory,code):
    source=directory/(code+'.cif');checksum=fixture(source)
    reference=directory/(code+'-assembly1.cif');reference_sum=fixture(reference,reference=True)
    source.with_suffix('.cif.source.json').write_text(json.dumps({'id':code,'sha256':checksum}))
    reference.with_suffix('.cif.source.json').write_text(json.dumps({'id':code,'assembly':1,'sha256':reference_sum}))
    return {'pdb':code,'sha256':checksum,'assembly_id':'1',
            'selection_note':'Invented technical coordinates, no biological claim',
            'assembly_reference':{'sha256':reference_sum,'complete':False}}


class AssemblyTransfer(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory();self.addCleanup(self.temporary.cleanup)
        self.directory=Path(self.temporary.name);self.source=self.directory/'source.cif'
        fixture(self.source);self.source_assembly=SourceAssembly(self.source,'1')

    def test_copies_preserve_author_identity_and_insertion(self):
        a,b=self.source_assembly.select([polymer(),polymer('A2',operator='2')])
        self.assertEqual([r[1:] for r in a['rows']],[r[1:] for r in b['rows']])
        self.assertEqual(a['rows'][1][3:],(11,'ALA','A'))
        np.testing.assert_allclose(a['rows'][1][0],[12,0,0])
        np.testing.assert_allclose(b['rows'][1][0],[0,2,0])

    def test_compound_operation_order(self):
        fixture(self.source,expression='(1)(2)');assembly=SourceAssembly(self.source,'1')
        names={key[1] for key in assembly.operations};self.assertEqual(len(names),1)
        item=assembly.select([polymer(operator=next(iter(names)))])[0]
        np.testing.assert_allclose([r[0] for r in item['rows']],[[10,0,0],[10,2,0],[7,0,0]])

    def test_ranges_products_and_generator_scope(self):
        fixture(self.source,expression='(1-2)(1,2)');assembly=SourceAssembly(self.source,'1')
        self.assertEqual(len(assembly.operations),12)
        point=assembly.select([polymer(operator='(2)(1)')])[0]['rows'][0][0]
        np.testing.assert_allclose(point,[0,10,0])
        for expression in ['(2-1)','(1,1)','(1)(missing)','(1)2','(1,,2)','(1-3)','((1))']:
            fixture(self.source,expression=expression)
            with self.subTest(expression=expression),self.assertRaises(ValueError):SourceAssembly(self.source,'1')
        def mutate(block):
            block.set_mmcif_category('_pdbx_struct_assembly_gen.',{'assembly_id':['1','1'],
                'oper_expression':['1','2'],'asym_id_list':['X,Z','Y']})
        fixture(self.source,mutate=mutate);assembly=SourceAssembly(self.source,'1')
        assembly.select([polymer(),polymer('B2','Y','B','2')],[companion()])
        with self.assertRaisesRegex(ValueError,'outside'):assembly.select([polymer(operator='2')])

    def test_required_companion_follows_owner(self):
        members=self.source_assembly.select([polymer('A2',operator='2')],[companion(operator='2',owner='A2')])
        np.testing.assert_allclose(members[-1]['rows'][0][0],[-1,1,1])
        self.assertEqual(members[-1]['rows'][0][3:],(99,'CA','B'))
        for spec in [companion(owner='missing'),companion(operator='2'),{**companion(),'component':'MG'},
                     {**companion(),'insertion_code':''},{**companion(),'residue_number':99.0}]:
            with self.subTest(spec=spec),self.assertRaises(ValueError):self.source_assembly.select([polymer()],[spec])

    def test_wrong_label_chain_operation_and_duplicates(self):
        for specs in [[polymer(label='Y')],[polymer(operator='missing')],
                      [polymer(),polymer('alias')],[polymer(),polymer('a1',operator='2')],
                      [{**polymer(),'description':'Different protein'}]]:
            with self.subTest(specs=specs),self.assertRaises(ValueError):self.source_assembly.select(specs)

    def test_independent_reference_all_members(self):
        reference=self.directory/'reference.cif';checksum=fixture(reference,reference=True)
        polymers=[polymer('A'+op,operator=op) for op in ['1','2']]+[polymer('B'+op,'Y','B',op) for op in ['1','2']]
        ligands=[companion('Ca'+op,op,'A'+op) for op in ['1','2']]
        members=self.source_assembly.select(polymers,ligands)
        report=verify_reference(members,reference,checksum,complete=True)
        self.assertEqual(report['matched_instances'],6)
        self.assertEqual(report['maximum_coordinate_error_angstrom'],0)
        reordered=self.directory/'reordered.cif'
        def reorder(block):
            block.set_mmcif_category('_atom_site.',{k:list(reversed(v)) for k,v in block.get_mmcif_category('_atom_site.').items()})
        reordered_sum=fixture(reordered,reference=True,mutate=reorder)
        self.assertTrue(verify_reference(members,reordered,reordered_sum,complete=True)['passed'])
        with self.assertRaisesRegex(ValueError,'omits'):verify_reference(members[:-1],reference,checksum,complete=True)
        with self.assertRaisesRegex(ValueError,'hash changed'):verify_reference(members,reference,'0'*64)
        changed=copy.deepcopy(members);changed[0]['rows'][0][0][0]+=.01
        with self.assertRaisesRegex(ValueError,'identity/coordinates'):verify_reference(changed,reference,checksum)
        changed=copy.deepcopy(members);row=list(changed[0]['rows'][0]);row[2]='N';changed[0]['rows'][0]=tuple(row)
        with self.assertRaisesRegex(ValueError,'identity/coordinates'):verify_reference(changed,reference,checksum)

    def test_invalid_source_operation_rejected(self):
        for value in ['-1','1.01','nan']:
            def mutate(block):block.find_loop('_pdbx_struct_oper_list.matrix[1][1]')[0]=value
            fixture(self.source,mutate=mutate)
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'proper rigid'):SourceAssembly(self.source,'1')

    def test_state_handoff_and_independent_tamper_checks(self):
        from prepare_source_states import prepare
        from audit_source_states import check_sources, digest
        states=[]
        for code in ['F001','F002']:
            spec=qualified_fixture(self.directory,code);spec['assembly_reference']['complete']=True
            polymers=[polymer('A'+op,operator=op) for op in ['1','2']]+[polymer('B'+op,'Y','B',op) for op in ['1','2']]
            states.append({**spec,
                'chains':[{**p,'description':'Technical polymer','label':p['instance_id'],'color':'70C8E0',
                           'sites':[{'residue':'11A','name':'ALA','atoms':['CA']}]} for p in polymers],
                'ligand_instances':[companion('Ca'+op,op,'A'+op) for op in ['1','2']]})
        config={'states':states,'reference':'F001','origin_angstrom':[2,3,4],
            'basis_columns':[[0,-1,0],[1,0,0],[0,0,1]],'limitations':['Technical fixture, not a protein'],
            'alignment':{'reference_chain':'A','reference_instance_id':'A1','maps':{
                code:{'chain':'A','instance_id':instance,'residue_pairs':[[10,10],['11A','11A'],[12,12]],
                    'expected_substitutions':[],'min_CA_coverage':1,'max_CA_rmsd_angstrom':1e-8,'note':'Known technical correspondence'}
                for code,instance in [('F001','A1'),('F002','A2')]}}}
        config_path=self.directory/'config.json';config_path.write_text(json.dumps(config))
        output=self.directory/'prepared';prepare(self.directory,config_path,output)
        report=check_sources(output,self.directory,config_path)
        self.assertTrue(report['passed']);self.assertEqual(report['source_counts'][0]['polymer_atoms'],12)
        manifest=json.loads((output/'manifest.json').read_text())
        a=np.load(output/manifest['states'][0]['actors'][0]['file'])
        b=np.load(output/manifest['states'][0]['actors'][1]['file'])
        np.testing.assert_allclose(a['atoms'][0],[-.3,-.8,-.4],atol=1e-12)
        self.assertFalse(np.allclose(a['atoms'],b['atoms']))
        for case in ['swapped_copies','missing_copy','missing_companion','moved_companion','wrong_owner','changed_operation','wrong_fit_instance','missing_assembly_proof']:
            target=self.directory/case;shutil.copytree(output,target);changed=copy.deepcopy(manifest);state=changed['states'][0]
            if case=='swapped_copies':
                file=state['actors'][0]['file'];d={k:a[k].copy() for k in a.files};d['atoms']=b['atoms'].copy()
                np.savez_compressed(target/file,**d);changed['data_sha256'][file]=digest(target/file)
            elif case=='missing_copy':state['actors'].pop()
            elif case=='missing_companion':state['ligands'].pop()
            elif case=='moved_companion':state['ligands'][1]['positions_nm'][0][0]+=.1
            elif case=='wrong_owner':state['ligands'][0]['owner_instance']='A2'
            elif case=='changed_operation':state['actors'][0]['source_operation_angstrom'][0][3]+=1
            elif case=='missing_assembly_proof':state['assembly'].pop('instance_check')
            else:state['actors'][0]['instance_id']='A2'
            (target/'manifest.json').write_text(json.dumps(changed))
            with self.subTest(case=case),self.assertRaises(ValueError):check_sources(target,self.directory,config_path)
        # List order is presentation, while the identities and owner must agree.
        manifest['states'][0]['ligands'].reverse()
        (output/'manifest.json').write_text(json.dumps(manifest))
        self.assertTrue(check_sources(output,self.directory,config_path)['passed'])

    def test_multichain_candidate_handoff(self):
        from prepare_binding_candidates import main,load_partner,fit_partner_template,restraint_indices
        from audit_binding_candidates import check
        spec=qualified_fixture(self.directory,'F001');config={'id':'technical-assembly-test'}
        for role,op in [('anchor','1'),('moving','2')]:
            config[role]={**spec,'label':'Technical '+role,'chains':[polymer('X',operator=op),polymer('Y','Y','B',op)],
                          'ligand_instances':[companion('Ca',op,'X')]}
        config['template']={role:{**copy.deepcopy(config[role]),'maps':[
            {'instance_id':name,'template_instance_id':name,'residue_pairs':[[10,10],['11A','11A'],[12,12]],
             'expected_substitutions':[],'min_CA_coverage':1,'max_CA_rmsd_angstrom':1e-8} for name in ['X','Y']]}
            for role in ['anchor','moving']}
        config['restraints']=[{'anchor':res,'moving':res,'anchor_instance':instance,'moving_instance':instance,
            'anchor_atom':'CA','moving_atom':'CA','range_nm':bounds,'evidence':'Known technical coordinate separation, not biological contacts'}
            for res,instance,bounds in [(10,'X',[.97,1.03]),('11A','X',[1.18,1.25]),(12,'Y',[1.88,1.97])]]
        config['search']={'seed':17,'atom_guard_nm':.24,'max_template_displacement_nm':.6,'cluster_rms_nm':.08}
        config_path=self.directory/'binding.json';config_path.write_text(json.dumps(config));output=self.directory/'binding'
        with patch('sys.argv',['prepare_binding_candidates.py','--sources',str(self.directory),'--config',str(config_path),
                               '--output',str(output),'--method','template','--frames','12']):main()
        report=check(self.directory,output)
        self.assertEqual(report['checked_poses'],13)
        self.assertEqual(report['source_members']['anchor'],{'polymer_atoms':6,'companion_atoms':1,'members':3})
        members,_=load_partner(self.directory,config['anchor'],'anchor')
        changed=copy.deepcopy(members)
        for row in changed[1]['rows']:row[0][0]+=2
        with self.assertRaisesRegex(ValueError,'shared partner'):
            fit_partner_template(changed,members,config['template']['anchor']['maps'])
        with self.assertRaisesRegex(ValueError,'every retained'):
            fit_partner_template(members,members,config['template']['anchor']['maps'][:1])
        rows={role:[r for m in members if m['kind']=='polymer' for r in m['rows']] for role in ['anchor','moving']}
        ids={role:[m['instance_id'] for m in members if m['kind']=='polymer' for r in m['rows']] for role in rows}
        for instance in [None,'missing']:
            restraints=copy.deepcopy(config['restraints']);restraints[0]['anchor_instance']=instance
            with self.subTest(instance=instance),self.assertRaisesRegex(ValueError,'retained polymer instance'):
                restraint_indices(rows,restraints,ids)
        manifest=json.loads((output/'interaction/manifest.json').read_text())
        for case in ['missing_member','swapped_instance','wrong_owner','missing_assembly_proof','moved_member']:
            target=self.directory/case;shutil.copytree(output,target);changed=copy.deepcopy(manifest)
            if case=='missing_member':
                actor=changed['actors'].pop(1);changed['data_sha256'].pop(actor['file'])
            elif case=='swapped_instance':changed['actors'][0]['instance_id']='Y'
            elif case=='wrong_owner':changed['actors'][2]['owner_instance']='Y'
            elif case=='missing_assembly_proof':changed.pop('partner_assembly_checks')
            else:
                actor=changed['actors'][1];path=target/'interaction'/actor['file']
                with np.load(path) as old:data={k:old[k].copy() for k in old.files}
                data['atoms'][:,0]+=.1;np.savez_compressed(path,**data)
                changed['data_sha256'][actor['file']]=hashlib.sha256(path.read_bytes()).hexdigest()
            (target/'interaction/manifest.json').write_text(json.dumps(changed))
            with self.subTest(case=case),self.assertRaises(ValueError):check(self.directory,target)


if __name__=='__main__':unittest.main()
