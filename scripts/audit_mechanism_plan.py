"""Check coverage/dependencies and the tracked-actin chemical ledger in a plan.

This is deliberately a plan audit: it does not verify literature, coordinates,
Blender scenes, kinetics or rendered movies. Standard library only.
"""
import argparse
import json
from pathlib import Path


REQUIRED = {'input', 'activation', 'monomer_supply', 'nucleation', 'elongation',
            'hydrolysis', 'phosphate_release', 'functional_output', 'termination',
            'disassembly', 'recycling', 'reentry'}
TRANSITIONS = {
    'nucleotide_exchange': ('G_ADP', 'G_ATP'),
    'incorporation': ('G_ATP', 'F_ATP'),
    'ATP_hydrolysis': ('F_ATP', 'F_ADP_Pi'),
    'Pi_release': ('F_ADP_Pi', 'F_ADP'),
    'depolymerization': ('F_ADP', 'G_ADP'),
    'reentry': ('G_ATP', 'F_ATP'),
}
TAGS = {'nucleotide_exchange': {'monomer_supply','recycling'},
        'incorporation': {'elongation'}, 'ATP_hydrolysis': {'hydrolysis'},
        'Pi_release': {'phosphate_release'}, 'depolymerization': {'disassembly'},
        'reentry': {'reentry'}}


def audit(path):
    data = json.loads(path.read_text())
    errors = []
    def require(ok, message):
        if not ok:
            errors.append(message)
    require(data.get('kind') == 'actin_mechanism_plan', 'Not a full actin mechanism plan; a local motion template is insufficient.')
    require(data.get('route') in {'arp23','formin'}, 'Specify arp23 or formin route.')
    if errors:
        return {'passed': False, 'level': 'plan_only', 'errors': errors}
    catalog = json.loads((path.parent/data['source_catalog']).read_text())['sources']
    stages = data.get('stages', [])
    ids = [s['id'] for s in stages]
    require(len(ids) == len(set(ids)), 'Duplicate stage IDs.')
    by_id = {s['id']: s for s in stages}
    coverage = set()
    for s in stages:
        for field in ('name','actors','before','after','implementation','checks','sources','mechanism_limits'):
            require(bool(s.get(field)), s['id']+': missing '+field)
        coverage.update(s.get('coverage', []))
        for dep in s.get('depends_on', []):
            require(dep in by_id and dep != s['id'], s['id']+': invalid prerequisite '+dep)
        for ref in s.get('sources', []):
            require(ref in catalog, s['id']+': unknown source '+ref)
        require(s.get('motion_evidence') in {'schematic','inferred','measured'}, s['id']+': missing motion evidence label')
    require(REQUIRED <= coverage, 'Missing coverage: '+', '.join(sorted(REQUIRED-coverage)))
    visited, active = set(), set()
    def visit(sid):
        if sid in active:
            require(False, 'Causal dependency cycle at '+sid+'; use a new cycle instance instead.')
            return
        if sid in visited or sid not in by_id:
            return
        active.add(sid)
        for dep in by_id[sid].get('depends_on', []):
            visit(dep)
        active.remove(sid)
        visited.add(sid)
    for sid in ids:
        visit(sid)

    ledger = data.get('chemical_ledger', {})
    events = ledger.get('events', [])
    state, nucleotide, bound_pi, free_pi = ledger.get('initial_state'), None, set(), set()
    types = []
    protein = ledger.get('protein_id')
    require(bool(protein), 'Missing stable tracked protein ID.')
    for e in events:
        typ = e.get('type')
        types.append(typ)
        label = e.get('stage', '?')+': '
        require(e.get('stage') in by_id, label+'unknown chemistry stage')
        require(e.get('actor') == protein, label+'tracked protein identity changed')
        require(typ in TRANSITIONS, label+'unknown chemical event')
        require(e.get('before') == state, label+'state discontinuity')
        require((e.get('before'), e.get('after')) == TRANSITIONS.get(typ), label+'invalid state transition')
        if e.get('stage') in by_id:
            require(bool(set(by_id[e['stage']]['coverage']) & TAGS.get(typ,set())), label+'chemistry/coverage mismatch')
        old, new = e.get('nucleotide_before'), e.get('nucleotide_after')
        require(bool(old) and bool(new), label+'missing nucleotide identity')
        if nucleotide is not None:
            require(old == nucleotide, label+'nucleotide identity discontinuity')
        consumed = e.get('consumed_from_pool', [])
        released = e.get('released_to_pool', [])
        product = e.get('new_bound_product', [])
        require(not any('PPi' in x for x in consumed+released+product), label+'actin ATPase must not produce PPi')
        if typ == 'nucleotide_exchange':
            require(old != new, label+'exchange must replace the nucleotide, not phosphorylate bound ADP')
            require(len(consumed)==1 and consumed[0].startswith('ATP_'), label+'exchange needs a new pool ATP')
            require(len(released)==1 and released[0].startswith('ADP_'), label+'exchange must release old ADP')
            require(not product, label+'exchange does not create a bound Pi')
        else:
            require(old == new, label+'nucleotide moiety changed without exchange')
        if typ == 'ATP_hydrolysis':
            require(consumed == ['H2O'] and not released, label+'hydrolysis consumes water and retains Pi in this selected cycle')
            require(len(product)==1 and product[0].startswith('Pi_'), label+'hydrolysis must create retained Pi')
            require(not (set(product) & (bound_pi|free_pi)), label+'duplicate phosphate identity')
            bound_pi.update(product)
        elif typ == 'Pi_release':
            require(len(released)==1 and set(released) <= bound_pi, label+'Pi must be previously bound and released once')
            require(not consumed and not product, label+'release cannot create a second phosphate')
            bound_pi.difference_update(released)
            free_pi.update(released)
        elif typ != 'nucleotide_exchange':
            require(not consumed and not released and not product, label+'unexpected chemistry in a physical assembly event')
        state, nucleotide = e.get('after'), new
    require(set(TRANSITIONS) <= set(types), 'The tracked cycle omits a chemical/physical stage.')
    require(types.count('nucleotide_exchange') >= 2, 'Track initial preparation and post-disassembly reloading separately.')
    require(state == 'F_ATP' and types[-1:] == ['reentry'], 'Tracked original protein must re-enter an available filament end.')
    require(not bound_pi and len(free_pi)==1, 'Expected one retained-then-released phosphate in the tracked demonstration cycle.')
    return {'passed':not errors, 'level':'plan_only', 'route':data['route'],
            'stages':len(stages), 'coverage':sorted(coverage), 'chemical_events':len(events),
            'errors':errors, 'upstream_adapter_status':data.get('scope',{}).get('upstream_adapter_status'),
            'production_readiness':'Not assessed. No scene, movie, source-truth or contact validation performed.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    result = audit(args.plan)
    output = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    print(output, end='')
    raise SystemExit(0 if result['passed'] else 1)
