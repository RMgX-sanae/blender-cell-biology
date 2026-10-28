"""Find public experimental PDB assets before modeling; save candidates, not endorsements.

Standard-library only. Other asset repositories still require their own search
and item-level checks; this helper deliberately searches RCSB only.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


def get_json(url):
    for attempt in range(2):
        try:
            with urlopen(Request(url, headers={'User-Agent': 'blender-cell-biology/asset-search'}), timeout=45) as r:
                if r.status == 204:
                    return {}
                return json.load(r)
        except HTTPError:
            raise
        except URLError:
            if attempt:
                raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('query', help='English molecule / assembly / state terms')
    p.add_argument('--limit', type=int, default=8)
    p.add_argument('--title', action='store_true', help='Match words in the structure title, reducing full-text noise')
    p.add_argument('--taxon', type=int, help='Optional NCBI taxonomy lineage ID, including descendants; do not guess it')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if not args.query.strip():
        p.error('Supply a non-empty search query')
    if not 1 <= args.limit <= 25:
        p.error('--limit must be 1–25')
    if args.output.exists():
        p.error('Output exists; inspect/reuse it or choose another file')
    query = {'type': 'terminal', 'service': 'full_text', 'parameters': {'value': args.query}}
    if args.title:
        # A single contains_words value is OR-like; require every supplied word.
        query = {'type': 'group', 'logical_operator': 'and', 'nodes': [
            {'type': 'terminal', 'service': 'text', 'parameters': {
                'attribute': 'struct.title', 'operator': 'contains_words', 'value': word}}
            for word in args.query.split()]}
    if args.taxon:
        query = {'type': 'group', 'logical_operator': 'and', 'nodes': [query, {
            'type': 'terminal', 'service': 'text', 'parameters': {
                'attribute': 'rcsb_entity_source_organism.taxonomy_lineage.id',
                'operator': 'exact_match', 'value': str(args.taxon)}}]}
    request = {'query': query, 'return_type': 'entry', 'request_options': {
        'paginate': {'start': 0, 'rows': args.limit}, 'results_content_type': ['experimental']}}
    endpoint = 'https://search.rcsb.org/rcsbsearch/v2/query?'
    response = get_json(endpoint + urlencode({'json': json.dumps(request)}))
    candidates = []
    for hit in response.get('result_set', []):
        entry = hit['identifier']
        metadata = get_json(f'https://data.rcsb.org/rest/v1/core/entry/{entry}')
        info = metadata.get('rcsb_entry_info', {})
        candidate = {
            'id': entry, 'title': metadata.get('struct', {}).get('title'),
            'entry_url': f'https://www.rcsb.org/structure/{entry}',
            'coordinate_url': f'https://files.rcsb.org/download/{entry}.cif',
            'method': [e['method'] for e in metadata.get('exptl', [])],
            'resolution_angstrom': info.get('resolution_combined'),
            'assembly_ids': metadata.get('rcsb_entry_container_identifiers', {}).get('assembly_ids', []),
            'primary_doi': metadata.get('rcsb_primary_citation', {}).get('pdbx_database_id_DOI'),
            'citation': metadata.get('rcsb_primary_citation', {}),
            'polymers': [],
            'coordinate_license': 'CC0',
            'license_url': 'https://www.ebi.ac.uk/pdbe/about/public-data-access-statement',
            'decision': 'candidate — inspect species, construct, state, assembly and missing regions'}
        # Entity records carry actual species, unlike the entry title alone.
        for entity in metadata.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', []):
            data = get_json(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{entity}')
            candidate['polymers'].append({
                'entity': entity,
                'description': data.get('rcsb_polymer_entity', {}).get('pdbx_description'),
                'source_organisms': data.get('rcsb_entity_source_organism', []),
                'author_chains': data.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', []),
                'mutation': data.get('entity_poly', {}).get('pdbx_mutation')})
        candidates.append(candidate)
    result = {'query': args.query, 'request': request, 'total_hits': response.get('total_count', 0),
              'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'provider': 'RCSB PDB',
              'candidates': candidates,
              'limits': 'Search relevance is not biological suitability. Coordinates are not downloaded; website figures have separate rights.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'total_hits': result['total_hits'], 'ids': [c['id'] for c in candidates],
                      'output': str(args.output)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
