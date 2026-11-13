"""Appearance adapter for the intact9SPI layered teaching view."""
from layered_molecular_appearance import digest, shade
import layered_molecular_appearance as shared

ROLES = {'protein', 'nucleic_acid', 'peptide_symbol'}
PALETTE = {'alpha', 'beta', 'heads', 'protein', 'rna', 'accent', 'dna1', 'dna2',
           'peptide', 'amino_acid', 'ink', 'secondary_ink'}



def read_profile(path):
    return shared.read_profile(path, roles=ROLES, palette=PALETTE, mechanism='translation')


def finish(scene, args, request, profile, context):
    return shared.finish(scene, args, request, profile, context,
        record_key='translation_appearance', context_name='Whole ribosomal context',
        scope='Appearance of the existing 9SPI layered teaching view. Whole 60S displacement, hidden 40S and all substrate motions are unchanged; defocus is editorial, not measured optics.')
