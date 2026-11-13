"""Appearance adapter for the existing authored replication fork."""
import json
from layered_molecular_appearance import shade
import layered_molecular_appearance as shared

ROLES = {'protein', 'nucleic_acid', 'primer_symbol'}
PALETTE = {'heads', 'protein', 'accent', 'dna1', 'dna2', 'ink'}


def read_profile(path):
    return shared.read_profile(path, roles=ROLES, palette=PALETTE, mechanism='replication')


def finish(scene, args, request, profile, context):
    record = shared.finish(scene, args, request, profile, context,
        record_key='replication_appearance', context_name='Whole CMG context',
        scope='Intact CMG is displaced behind an authored fork; nucleotide geometry, primer symbols and growth order are unchanged. Defocus is editorial, not measured optics.')
    record['language'] = 'zh' if args.font else 'en'
    scene['replication_appearance'] = json.dumps(record)
    return record
