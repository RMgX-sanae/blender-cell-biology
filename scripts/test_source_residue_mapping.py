"""Exercise real AP2 source correspondence and misleading mapping controls."""
import argparse
import copy
import gzip
import json
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

from source_residue_mapping import correspondence


CHAINS=[{'first':c,'second':c,'uniprot':u} for c,u in [('A','P18484'),('B','P63010'),('M','P84092'),('S','P62743')]]
INSERTIONS=[{'first':['A','272'],'second':['A','272'],'name':'GLU',
             'reason':'Reviewed shared insertion between canonical P18484 residues271 and272.'}]


def check(sources):
    paths=[sources/name for name in ('6YAE.cif','6YAF.cif','6yae.xml.gz','6yaf.xml.gz')]
    _,report=correspondence(*paths,CHAINS,INSERTIONS)
    pairs={tuple(p['second']):p for p in report['pairs']}
    assert report['observed_counts']==[1717,1762] and len(pairs)==1711
    assert pairs['M','350']['first']==['M','339'] and pairs['M','350']['name']=='LYS'
    assert pairs['M','367']['first']==['M','356'] and pairs['M','368']['first']==['M','357']
    assert {tuple(v['author']) for v in report['unmatched']['first']}=={('M',str(n)) for n in range(136,142)}
    assert all(('M',str(n)) not in pairs for n in range(237,248))
    _,automatic=correspondence(*paths,CHAINS)
    assert len(automatic['pairs'])==1710 and ('A','272') not in {tuple(p['second']) for p in automatic['pairs']}
    rejected=[]
    def reject(name,files=paths,chains=CHAINS,insertions=INSERTIONS):
        try:correspondence(*files,chains,insertions)
        except ValueError as error:rejected.append({'case':name,'reason':str(error)})
        else:raise AssertionError('Accepted invalid mapping: '+name)
    reject('swapped SIFTS entries',paths[:2]+[paths[3],paths[2]])
    wrong=copy.deepcopy(CHAINS);wrong[0]['uniprot']='P62743';reject('wrong protein accession',chains=wrong)
    wrong=copy.deepcopy(CHAINS);wrong[0]['second']='S';reject('duplicate target chain',chains=wrong)
    bad=copy.deepcopy(INSERTIONS);bad[0].update(first=['A','273'],second=['A','273'],name='ASP')
    reject('override canonical identity as an insertion',insertions=bad)
    bad=copy.deepcopy(INSERTIONS);bad[0]['reason']='';reject('unreviewed insertion',insertions=bad)
    with tempfile.TemporaryDirectory(prefix='ap2-sifts-controls-') as temp:
        for change in ('name','number'):
            tree=ET.fromstring(gzip.decompress(paths[3].read_bytes()))
            for res in tree.iter():
                if res.tag.split('}')[-1]!='residue':continue
                refs=[r for r in res if r.tag.split('}')[-1]=='crossRefDb']
                pdb=next((r for r in refs if r.attrib.get('dbSource')=='PDB' and r.attrib.get('dbChainId')=='M' and r.attrib.get('dbResNum')=='350'),None)
                if pdb is None:continue
                if change=='name':pdb.set('dbResName','ALA')
                else:
                    unp=next(r for r in refs if r.attrib.get('dbSource')=='UniProt');unp.set('dbResNum','340')
            modified=Path(temp)/(change+'.xml');modified.write_bytes(ET.tostring(tree))
            reject('changed SIFTS residue '+change,paths[:3]+[modified])
    return {'passed':True,'observed_counts':report['observed_counts'],'mapped_residues':len(pairs),
            'unmatched_counts':{k:len(v) for k,v in report['unmatched'].items()},'negative_controls':rejected}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path);a=p.parse_args();result=check(a.sources)
    if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
