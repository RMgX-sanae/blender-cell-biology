"""Real-source site-link regressions; no render or synthetic biology fixture."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import tempfile

import gemmi

from audit_site_links import inspect
from prepare_binding_candidates import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists(): p.error('Choose a fresh report path')
    cases = {s['id']: s for s in json.loads(a.config.read_text())['cases']}
    positive = {name: inspect(a.sources, spec) for name, spec in cases.items()}
    assert all(r['passed'] for r in positive.values())
    rejected = []
    def rejects(name, spec, sources=a.sources):
        try:
            result = inspect(sources, spec)
        except (ValueError, KeyError) as error:
            rejected.append({'case': name, 'error': str(error)})
        else:
            if result['passed']: raise AssertionError('Invalid input accepted: '+name)
            rejected.append({'case': name, 'relations': result['relations']})
    for name, change in [
        ('wrong structure hash', lambda s: s['source'].update(sha256='0'*64)),
        ('wrong CCD hash', lambda s: s['ccd_sha256'].update(MEQ='0'*64)),
        ('wrong source chain identity', lambda s: s['source']['chains'][0].update(description='eRF1')),
        ('native O3 absent in analogue', lambda s: s['sites']['tRNA_N3'].update(atom="O3'", element='O')),
        ('amide N relabeled ester O', lambda s: s['sites']['tRNA_N3'].update(element='O')),
        ('methylglutamine relabeled glutamine', lambda s: s['sites']['Q230_methyl'].update(component='GLN')),
        ('label chain used as author chain', lambda s: s['sites']['peptide_C'].update(chain='DB')),
        ('wrong insertion code', lambda s: s['sites']['peptide_C'].update(residue='7A')),
        ('covalent analogue called noncovalent', lambda s: s['relations'][0].update(kind='noncovalent')),
        ('nearby O2 called already bonded', lambda s: s['relations'][-1].update(kind='covalent', annotation='required', distance_angstrom=[1.2,1.6])),
        ('NaN distance threshold', lambda s: s['relations'][0].update(distance_angstrom=[1.2,float('nan')])),
        ('reversed distance threshold', lambda s: s['relations'][0].update(distance_angstrom=[1.5,1.2])),
        ('self relation', lambda s: s['relations'][0].update(sites=['peptide_C','peptide_C'])),
        ('duplicate relation ID', lambda s: s['relations'].append(copy.deepcopy(s['relations'][0]))),
    ]:
        spec = copy.deepcopy(cases['rf1-analogue']); change(spec); rejects(name, spec)
    # A genuine deposited annotation conflict must survive in the failed report.
    spec = copy.deepcopy(cases['rf1-analogue'])
    spec['source']['chains'].append({'chain':'1A','description':'23S Ribosomal RNA'})
    for code in ('ARG', '5MU'):
        meta = json.loads((a.sources/(code+'.ccd.cif.source.json')).read_text())
        spec['ccd_sha256'][code] = meta['sha256']
    spec['sites'].update(rRNA_O2={'chain':'1A','residue':1915,'component':'5MU','atom':'O2','element':'O'},
                         RF1_NH1={'chain':'1w','residue':286,'component':'ARG','atom':'NH1','element':'N'})
    spec['relations'] = [{'id':'unresolved_rRNA_RF1_contact', 'sites':['rRNA_O2','RF1_NH1'],
                         'kind':'noncovalent', 'distance_angstrom':[2.4,4.0], 'annotation':'inspect'}]
    rejects('deposited 1.30 A RNA-protein link is not auto-approved', spec)
    assert len(rejected[-1]['relations'][0]['errors']) == 2
    if 'rf1-deoxy-ester' in cases:
        for name, change in [
            ('nominal A silently treated as intact ribose', lambda s: s['residue_inventories'][0]['expected_absent'].pop("O2'")),
            ('present O3 silently treated as absent', lambda s: s['residue_inventories'][0]['expected_absent'].update({"O3'":'Incorrect omission'})),
            ('nonexistent CCD atom omitted', lambda s: s['residue_inventories'][0]['expected_absent'].update({'NOT_AN_ATOM':'Incorrect atom'})),
            ('missing construct explanation', lambda s: s['residue_inventories'][0]['expected_absent'].update({"O2'":''})),
            ('inventory assigned to wrong residue', lambda s: s['residue_inventories'][0].update(residue=75)),
            ('duplicate residue inventory', lambda s: s['residue_inventories'].append(copy.deepcopy(s['residue_inventories'][0]))),
        ]:
            spec = copy.deepcopy(cases['rf1-deoxy-ester']); change(spec); rejects(name, spec)
        assert positive['rf1-deoxy-ester']['residue_inventories'][0]['missing_from_ccd'] == ["O2'", 'OP3']
    # Saved/reloaded source mutations: update the receipt so geometry and
    # annotation checks must reject independently of the checksum gate.
    with tempfile.TemporaryDirectory(prefix='site-links-') as temporary:
        root = Path(temporary)
        for file in a.sources.glob('*.ccd.cif*'): shutil.copy2(file, root/file.name)
        original = a.sources/'2PTC.cif'
        meta = json.loads(original.with_suffix('.cif.source.json').read_text())
        for name in ('shifted inhibitor carbonyl', 'different symmetry disulfide',
                     'different conformer disulfide', 'inconsistent label chain',
                     'inconsistent label residue', 'multiple models'):
            block = gemmi.cif.read(str(original)).sole_block()
            if name == 'shifted inhibitor carbonyl':
                table = block.find('_atom_site.', ['auth_asym_id','auth_seq_id','label_atom_id','Cartn_x'])
                target = next(r for r in table if list(r)[:3] == ['I','15','C'])
                target[3] = str(float(target[3])+10)
            elif name != 'multiple models':
                field = {'different symmetry disulfide':'ptnr1_symmetry',
                         'different conformer disulfide':'pdbx_ptnr1_label_alt_id',
                         'inconsistent label chain':'ptnr1_label_asym_id',
                         'inconsistent label residue':'ptnr1_label_seq_id'}[name]
                table = block.find('_struct_conn.', ['conn_type_id','ptnr1_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id',field])
                target = next(r for r in table if list(r)[:4] == ['disulf','I','5','55'])
                target[4] = {'different symmetry disulfide':'2_555',
                             'different conformer disulfide':'A', 'inconsistent label chain':'NOT_I',
                             'inconsistent label residue':'999'}[name]
            else:
                loop = block.find_loop('_atom_site.pdbx_PDB_model_num').get_loop()
                model_column = list(loop.tags).index('_atom_site.pdbx_PDB_model_num')
                id_column = list(loop.tags).index('_atom_site.id')
                original_rows = [list(loop.values[i:i+loop.width()])
                                 for i in range(0, len(loop.values), loop.width())]
                for i, row in enumerate(original_rows, len(original_rows)+1):
                    row[model_column] = '2'; row[id_column] = str(i); loop.add_row(row)
            file = root/'2PTC.cif'; block.write_file(str(file))
            meta['sha256'] = sha(file)
            file.with_suffix('.cif.source.json').write_text(json.dumps(meta))
            spec = copy.deepcopy(cases['trypsin-bpti']); spec['source']['sha256'] = meta['sha256']
            rejects(name, spec, root)
            if name == 'multiple models':
                assert 'single model' in rejected[-1]['error']
    report = {'passed': True, 'positive': positive, 'negative_cases': rejected}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'passed': True, 'positive_cases': len(positive), 'negative_cases': len(rejected)}))


if __name__ == '__main__': main()
