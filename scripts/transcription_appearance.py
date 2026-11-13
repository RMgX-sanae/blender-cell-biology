"""Appearance adapter for the existing5FLM C/U teaching reconstruction."""
import json
import bpy
import numpy as np
from compositor_stack import isolate_labels
from layered_molecular_appearance import digest, shade
import layered_molecular_appearance as shared

ROLES = {'protein', 'nucleic_acid', 'phosphate'}
PALETTE = {'heads', 'protein', 'rna', 'accent', 'dna1', 'dna2', 'ink'}


def read_profile(path):
    return shared.read_profile(path, roles=ROLES, palette=PALETTE, mechanism='transcription')


def finish(scene, args, request, profile, context):
    record = shared.finish(scene, args, request, profile, context,
        record_key='transcription_appearance', context_name='Whole Pol II context',
        scope='Intact Pol II is displaced for a layered teaching view; hidden non-template DNA and reaction geometry remain unchanged. Defocus is editorial, not measured optics.')
    record.update(language='zh' if args.font else 'en', cycles=1)
    scene['transcription_appearance'] = json.dumps(record)
    return record


def extend(scene, args, record, font_cjk, font_latin):
    """New captions join the inherited protected collection; keep its compositor."""
    scene.frame_set(1)
    bpy.context.view_layer.update()
    shared.apply_fonts(scene, font_cjk, font_latin)
    labels = isolate_labels(scene, record['studio']['background'], 26*args.width/1280)
    record['compositor']['protected_labels'] = labels
    scene.frame_set(args.frame)
    bpy.context.view_layer.update()
    record.update(profile_path=str(args.output.resolve()/'appearance.json'), cycles=2,
                  fonts=[{'object': o.name, 'body': o.data.body} for o in scene.objects if o.type == 'FONT'],
                  camera_matrix=np.array(scene.camera.matrix_world).tolist())
    record['scope'] += ' Two-cycle extension retains authored23-second timing and template/RNA motion; it is not a complete transcription cycle.'
    scene['transcription_appearance'] = json.dumps(record)
    return record
