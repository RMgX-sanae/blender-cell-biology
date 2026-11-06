"""Reject wrong grouping and source-assembly mappings using the downloaded case."""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from prepare_molecular_asset import atom_table
from prepare_clathrin_assembly import rigid,compare_reference,check_groups


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--assets',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--obsolete-source',type=Path,help='Optional cached 1XI4 for the real discontinuous-chain regression')
    a=p.parse_args()
    m=json.loads((a.assets/'assembly.json').read_text());config=m['source_config']
    table=atom_table(a.sources/(config['source_id']+'.cif'));ref=atom_table(a.sources/(config['source_id']+'-assembly'+config['assembly_id']+'.cif'))
    units=m['instances'];groups=m['source_config']['groups'];checks=[]
    compare_reference(table,ref,copy.deepcopy(units));check_groups(table,groups)
    def reject(name,call):
        try:call()
        except ValueError as error:checks.append({'case':name,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted '+name)
    for name,scale in [('reflection',-1.),('stretch',1.05),('nonfinite',float('nan'))]:
        matrix=np.eye(4);matrix[0,0]=scale;reject(name,lambda m=matrix:rigid(m))
    reject('missing-complete-unit',lambda:compare_reference(table,ref,copy.deepcopy(units[:-1])))
    duplicate=copy.deepcopy(units);duplicate[-1]=copy.deepcopy(duplicate[0])
    reject('duplicate-complete-unit',lambda:compare_reference(table,ref,duplicate))
    moved=copy.deepcopy(units);moved[0]['matrix_nm'][0][3]+=.1
    reject('wrong-source-placement',lambda:compare_reference(table,ref,moved))
    changed=dict(table);changed['A']=copy.deepcopy(table['A'])
    for index,offset in [(0,.02),(1,-.02)]:changed['A'][index][0][0]+=offset
    reject('same-centroid-wrong-sites',lambda:compare_reference(changed,ref,copy.deepcopy(units)))
    changed=dict(table);changed['A']=list(table['A']);changed['A'][0],changed['A'][1]=changed['A'][1],changed['A'][0]
    reject('same-centroid-wrong-identities',lambda:compare_reference(changed,ref,copy.deepcopy(units)))
    changed=copy.deepcopy(groups);changed[0]['heavy'][1],changed[1]['heavy'][1]=changed[1]['heavy'][1],changed[0]['heavy'][1]
    reject('mixed-terminal-hubs',lambda:check_groups(table,changed))
    changed=copy.deepcopy(groups);changed[0]['light'][0],changed[1]['light'][0]=changed[1]['light'][0],changed[0]['light'][0]
    reject('wrong-light-chain-owner',lambda:check_groups(table,changed))
    changed=copy.deepcopy(table);changed['A'][837][0][0]+=100
    reject('disconnected-source-sections',lambda:check_groups(changed,groups))
    changed=copy.deepcopy(table);changed['A'][837],changed['A'][838]=changed['A'][838],changed['A'][837]
    reject('nonconsecutive-source-identity',lambda:check_groups(changed,groups))
    if a.obsolete_source:
        reject('actual-1XI4-disconnected-author-chains',lambda:check_groups(atom_table(a.obsolete_source),groups))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps({'passed':True,'negative_cases':len(checks),'checks':checks},indent=2)+'\n')
    print(json.dumps({'passed':True,'negative_cases':len(checks)}))


if __name__=='__main__':main()
