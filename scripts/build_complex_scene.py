"""Source-derived human ribosome/replisome structural views, not process animations."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from detail_scene import (collection, fit_orthographic, load_chain,
                          material, save_and_render, studio)
from build_molecular_scene import nucleic_atoms


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scene',choices=['ribosome','replisome'],required=True)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--style-profile',type=Path,
                   default=Path(__file__).resolve().parent.parent/'assets/style-molecular-reference.json')
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--compositor',choices=['clean','reference','monochrome'])
    p.add_argument('--width',type=int,default=1920)
    p.add_argument('--samples',type=int,default=64)
    p.add_argument('--render',action='store_true')
    p.add_argument('--reveal',action='store_true',help='Hide selected COMPLETE foreground subunits; never slice individual molecules')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    args.frames,args.fps,args.frame,args.animate=1,24,1,False
    if args.output.exists():p.error('Choose a fresh output directory.')
    source=json.loads((args.asset_dir/'manifest.json').read_text())
    required={'ribosome':'9SPI','replisome':'7PFO'}
    if not source['source_file'].upper().startswith(required[args.scene]):
        p.error('This chain map is specific to '+required[args.scene])
    profile=json.loads(args.style_profile.read_text())
    palette=profile['palette']
    mats={k:material(k,v,profile['roughness'],profile['bump']) for k,v in palette.items()}
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    proteins=collection('01 Structure-derived protein subunits')
    rna=collection('02 Structure-derived ribosomal RNA')
    active=collection('03 Modeled DNA, mRNA and tRNA residues')
    roles={}
    for item in source['chains']:
        chain=item['chain']
        if args.scene=='ribosome':
            role={'At':'A-site tRNA','Pt':'P-site tRNA','Et':'E-site tRNA','Mr':'mRNA'}.get(chain)
            if role:
                color={'At':'protein','Pt':'rna','Et':'accent','Mr':'dna1'}[chain]
                data=np.load(args.asset_dir/f'chain_{chain}.npz',allow_pickle=False)
                nucleic_atoms(data,role,mats[color],active)
            else:
                large=chain.startswith('L')
                role=('60S' if large else '40S')+(' rRNA' if item['nucleic_acid'] else ' protein')
                color=('alpha' if large else 'beta') if item['nucleic_acid'] else ('heads' if large else 'accent')
                load_chain(args.asset_dir,chain,mats[color],rna if item['nucleic_acid'] else proteins,
                           f'{role} chain {chain}',representation='surface')
        else:
            role={'M':'Leading-strand DNA','N':'Lagging-strand DNA',
                  'A':'Polymerase epsilon subunit 2','B':'Polymerase epsilon catalytic fragment',
                  'C':'CDC45','H':'AND-1','I':'AND-1','J':'AND-1',
                  'K':'TIMELESS','L':'TIPIN','Q':'Claspin'}.get(chain)
            if chain in ('M','N'):
                data=np.load(args.asset_dir/f'chain_{chain}.npz',allow_pickle=False)
                nucleic_atoms(data,role,mats['dna1' if chain=='M' else 'dna2'],active)
            else:
                if chain in '234567':role,color='MCM helicase subunit','alpha'
                elif chain in 'AB':color='protein'
                elif chain in 'DEFG':role,color='GINS subunit','beta'
                elif chain in 'HIJ':color='heads'
                else:color='accent'
                load_chain(args.asset_dir,chain,mats[color],proteins,
                           f'{role} chain {chain}',representation='surface')
        roles[chain]=role
    position=(-50,-12,-16) if args.scene=='ribosome' else (35,-45,35)
    if args.scene=='ribosome':
        profile['lighting'].append({'name':'Ribosome active-face softbox','position':[-28,-8,-12],
                                    'power':15000,'color':'FFE7CC','size':18})
    hidden=[]
    if args.reveal:
        # These are presentation omissions, not a dissociated biological state.
        # Keep all atomic coordinates and complete molecular meshes unchanged.
        for obj in list(proteins.objects)+list(rna.objects):
            chain=obj['source_chain']
            if (chain.startswith('S') if args.scene=='ribosome' else chain in ['K','L','Q']):
                obj.hide_render=obj.hide_viewport=True
                hidden.append(chain)
        if args.scene=='replisome':
            position=(-25,30,-45)
            profile['lighting'].append({'name':'DNA-facing softbox','position':[-14,18,-28],
                                        'power':12000,'color':'FFE7CC','size':20})
    camera=studio(profile,args,{'position':position,
                               'target':(0,0,0),'scale':48,'light_scale':1.2})
    subjects=list(proteins.objects)+list(rna.objects)+list(active.objects)
    fit_orthographic(camera,subjects)
    save_and_render(args,{'scene':args.scene,'source_asset':source,'chain_roles':roles,
                          'intact_subunit_reveal':{'hidden_chains':hidden,'coordinates_moved':False,
                                                   'molecular_meshes_sliced':False,
                                                   'note':'Selected whole subunits hidden for visibility; not biological dissociation.'},
                          'coordinate_units':'1 Blender unit = 1 nm',
                          'representation':'Gaussian display envelopes; active nucleic acids use atoms and distance-inferred bonds',
                          'limitations':['Static structure, not a modeled reaction cycle or a molecular dynamics trajectory.',
                                         'Only residues modeled in the deposited structure are shown; missing regions are not reconstructed.',
                                         'Subunits remain in deposited coordinates. Reveal view hides complete components and records them; it does not cut their meshes.',
                                         'Chain colors and Gaussian envelope parameters are artistic display choices.']})


if __name__=='__main__':main()
