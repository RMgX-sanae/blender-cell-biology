"""Register four independent ribosome specimens to the molecular frame of the 9SPI elongation scene."""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import numpy as np

from prepare_binding_candidates import checked_source, rigid_fit_points, transform


RECIPE=Path(__file__).resolve().parents[1]/'assets/translation-route-sources.json'


def read_model(directory,spec):
    structure=gemmi.read_structure(str(checked_source(directory,spec)))
    assembly=next(a for a in structure.assemblies if a.name=='1')
    if len(assembly.generators)!=1:
        raise ValueError(f"{spec['pdb']} requires a single assembly generator")
    gen=assembly.generators[0]
    if len(gen.operators)!=1 or not gen.operators[0].transform.is_identity():
        raise ValueError(f"{spec['pdb']} is outside the reviewed assembly frame")
    return structure,set(gen.subchains)


def landmarks(structure,chain):
    polymer=structure[0][chain].get_polymer()
    sequence=structure.get_entity_of(polymer).full_sequence
    bases=''.join(x if x in ('A','C','G','U') else 'N' for x in sequence)
    sites={res.label_seq:list(atom.pos) for res in polymer for atom in res
           if res.name in ('A','C','G','U') and atom.name=="C4'" and atom.occ>0 and atom.altloc in ('\x00','A')}
    return bases,sites


def register(source,reference):
    seq,sites=source
    ref,ref_sites=reference
    alignment=gemmi.align_string_sequences(list(seq),list(ref),[],gemmi.AlignmentScoring())
    i=j=0; pairs=[]
    for base,target in zip(alignment.add_gaps(seq,1),alignment.add_gaps(ref,2)):
        if base!='-':i+=1
        if target!='-':j+=1
        if base==target and base in 'ACGU' and i in sites and j in ref_sites:
            pairs.append((i,j))
    if len(pairs)<1000:
        raise ValueError('Insufficient corresponding 18S landmarks in these mammalian samples')
    moving=np.array([sites[i] for i,_ in pairs]); fixed=np.array([ref_sites[j] for _,j in pairs])
    initial,_=rigid_fit_points(moving,fixed)
    core=np.linalg.norm(transform(moving,initial)-fixed,axis=1)<5
    if core.sum()<.8*len(pairs):
        raise ValueError('The selected 18S cores disagree in this viewing frame')
    matrix,rmsd=rigid_fit_points(moving[core],fixed[core])
    return matrix,{'matched':len(pairs),'retained':int(core.sum()),'rmsd_angstrom':rmsd,
                   'sequence_cigar':alignment.cigar_str()}


def prepare(sources,manifest,output):
    recipe=json.loads(RECIPE.read_text())
    reference=recipe['reference']
    shot=json.loads(manifest.read_text())
    donor=shot['source_asset']
    if (shot.get('scene')!='Human translation elongation substrate cycle'
            or donor['source_sha256']!=reference['sha256'] or donor.get('chain_file_schema')!=1
            or len(set(row['file'].casefold() for row in donor['chains']))!=83):
        raise ValueError('Use the source-qualified 83-chain 9SPI elongation scene manifest')
    template,_=read_model(sources,reference)
    anchors=landmarks(template,reference['ssu_rna'])
    center=np.asarray(donor['source_origin_angstrom'])
    axes=np.asarray(donor['source_basis_columns'])
    if output.exists():raise ValueError('Choose a fresh frame-registration directory')
    output.parent.mkdir(parents=True,exist_ok=True)
    results={}
    with TemporaryDirectory(dir=output.parent) as temp:
        pending=Path(temp)/output.name;pending.mkdir()
        for spec in recipe['states']:
            entry=spec['pdb']
            structure,assembly=read_model(sources,spec)
            model=structure[0]
            move,fit=register(landmarks(structure,spec['ssu_rna']),anchors)
            rotation=move[:3,:3]
            origin=rotation.T@(center-move[:3,3]); basis=rotation.T@axes
            chains=([c.name for c in model if len(c.get_polymer())] if spec['chains']=='all' else spec['chains'])
            for ch in chains:
                polymer=model[ch].get_polymer()
                if not len(polymer) or any(r.subchain not in assembly for r in polymer):
                    raise ValueError(f'{entry} chain {ch} outside selected assembly')
            ch,first,codon=spec['codon']
            triplet=[r for r in model[ch].get_polymer() if first<=r.seqid.num<first+3]
            if ''.join(r.name for r in triplet)!=codon:
                raise ValueError(f'{entry} selected codon changed')
            for residue in triplet:
                required={"C4'",'N9' if residue.name in ('A','G') else 'N1'}
                resolved={a.name for a in residue if a.occ>0 and a.altloc in ('\x00','A')}
                if not required<=resolved:
                    raise ValueError(f'{entry} {ch}:{residue.seqid} lacks resolved codon landmarks')
            ions=[]
            for ch,num,component in spec['companions']:
                r=next(r for r in model[ch] if r.seqid.num==num and r.name==component and r.subchain in assembly)
                atoms=[{'atom':a.name,'element':a.element.name,'xyz_nm':((np.array(list(a.pos))-origin)@basis/10).tolist()}
                       for a in r if not a.element.is_hydrogen and a.occ>0 and a.altloc in ('\x00','A')]
                if not atoms:raise ValueError(f'{entry} empty companion {ch}:{num}')
                ions.append({'chain':ch,'number':num,'component':component,'atoms':atoms})
            config={'id':'translation_'+entry.lower(),'source_id':entry,'source_sha256':spec['sha256'],
                    'selection_note':spec['scope'],'chains':[{'chain':c} for c in chains],
                    'orientation':{'origin_angstrom':origin.tolist(),'basis_columns':basis.tolist(),
                                   'note':'Conserved 18S C4′ core in the same viewing frame as the separate 9SPI elongation scene.',
                                   'evidence_urls':['https://www.rcsb.org/structure/'+entry,'https://www.rcsb.org/structure/9SPI']}}
            (pending/(entry+'.json')).write_text(json.dumps(config,indent=2)+'\n')
            results[entry]={'fit':fit,'coordinate_unit':'nm','roles':spec['roles'],'codon':spec['codon'],
                            'polymer_chains':chains,'companions':ions,'config':entry+'.json'}
        (pending/'registration.json').write_text(json.dumps(results,indent=2)+'\n')
        pending.rename(output)
    return {k:v['fit'] for k,v in results.items()}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--elongation-manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(prepare(args.sources,args.elongation_manifest,args.output)))
