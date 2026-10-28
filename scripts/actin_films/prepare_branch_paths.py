"""Prepare intact Arp2/3 recruitment and sequential daughter-actin teaching paths."""
import hashlib
import json

import numpy as np

from prepare_actin_paths import displacement,screen_translations,smooth
from prepare_contact_repairs import repair_path
from project_context import ROOT


def main():
    assets = ROOT/'assets';output = assets/'branch-path-repairs.json'
    record = {'passed':False,'motions':{},'intervals':{'arp':[.10,.30]},
        'daughter_reveal':.303,
        'source':{'entry':'7TPT','url':'https://www.rcsb.org/structure/7TPT',
            'sha256':hashlib.sha256((ROOT/'sources/7TPT.cif').read_bytes()).hexdigest()},
        'limits':'Mature stabilized branch geometry supplies the endpoint. Arp recruitment, daughter monomer order, WCA timing and paths are teaching choices, not observed activation or nucleation intermediates.'}
    output.write_text(json.dumps(record,indent=2)+'\n')
    atoms = {c:np.load(assets/'7TPT'/(c+'.npz'))['atoms'].astype(float) for c in 'ABCDEFGHIJKNOPQRS'}
    combine = lambda names:np.concatenate([atoms[c] for c in names])
    plans = json.loads((assets/'approach-plans.json').read_text())
    record['motions']['arp'] = repair_path(combine('ABCDEFG'),combine('NOPQRS'),plans['arp']['offset'])
    previous = 'ABCDEFGNOPQRS'
    for i,c in enumerate('HIJK'):
        key = 'branch_'+c
        record['motions'][key] = repair_path(atoms[c],combine(previous),np.array(plans[key]['offset'])*.5,[i*6,i*10,-i])
        record['intervals'][key] = [.31+i*.10,.43+i*.10]
        previous += c
        print('BRANCH_PATH_READY',key,record['motions'][key]['after_minimum_nm'],flush=True)

    def trajectory(fraction):
        first,last = record['intervals']['arp']
        shift = displacement(1-smooth((fraction-first)/(last-first)),record['motions']['arp'])
        shifts = {c:shift if c in 'ABCDEFG' else np.zeros(3) for c in 'ABCDEFGNOPQRS'}
        if fraction>=record['daughter_reveal']:
            for c in 'HIJK':
                first,last = record['intervals']['branch_'+c]
                shifts[c] = displacement(1-smooth((fraction-first)/(last-first)),record['motions']['branch_'+c])
        return shifts

    record['screen'] = screen_translations(atoms,trajectory)
    record['passed'] = record['screen']['passed']
    output.write_text(json.dumps(record,indent=2)+'\n')
    if not record['passed']:raise ValueError('Coupled branch recruitment still overlaps; inspect branch-path-repairs.json')
    print('BRANCH_RECRUITMENT_READY',flush=True)


if __name__=='__main__':main()
