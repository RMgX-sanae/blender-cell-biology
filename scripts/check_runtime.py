"""Check a selected molecular-video runtime before a long render.

Uses this Python, the supplied Blender/FFmpeg and local fonts. Installs nothing.
A tiny technical probe exercises shared geometry, typography and neutral
composition, then the existing video encoder/decoder. It is not a biology test.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


PROBE = r'''
import json, sys
from pathlib import Path
from types import SimpleNamespace
import bpy, numpy as np
from mathutils import Matrix
p=Path(sys.argv[sys.argv.index('--')+1]);spec=json.loads(p.read_text())
root=Path(spec['skill']);sys.path.insert(0,str(root/'scripts'))
import detail_scene as ds
import style_nodes as sn
from build_ubiquitin_transfer_scene import envelope
from compositor_stack import build_stack
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
a=SimpleNamespace(engine=spec['engine'],width=320,samples=4,fps=24,frames=3,compositor='clean')
profile=json.loads((root/'assets/style-molecular-reference.json').read_text())
profile.update(camera={},bump=0,background='#071f2b')
camera=ds.studio(profile,a,{'position':[0,0,10],'target':[0,0,0],'scale':5,'light_scale':1})
camera.matrix_world=Matrix.Translation((0,0,10));camera.data.dof.use_dof=False
scene=bpy.context.scene;scene.render.resolution_y=180
points=np.load(p.parent/'probe.npz')['atoms']
obj=sn.point_object('Technical geometry probe',points,ds.collection('Probe geometry'))
envelope(obj,ds.material('Probe material','#71bdc2',.5,0))
font=bpy.data.curves.new('Font probe','FONT');font.body='细胞 ATP 123'
font.font=bpy.data.fonts.load(spec['font_zh']);font.font_bold=bpy.data.fonts.load(spec['font_en'])
font.font.pack();font.font_bold.pack();font.align_x='CENTER';font.size=.35
for char,style in zip(font.body,font.body_format):style.use_bold=ord(char)<128
label=bpy.data.objects.new('Font probe',font);scene.collection.objects.link(label);label.location=(0,-.85,1)
recipe=next(r.copy() for r in json.loads((root/'assets/render-recipes.json').read_text())['recipes'] if r['id']=='neutral')
recipe.update(background=profile['background'],label_min_pixels=18)
compositor=build_stack(scene,recipe)
scene.render.filepath=str(p.parent/'probe.png');bpy.ops.render.render(write_still=True)
evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=evaluated.to_mesh()
surface_vertices=len(mesh.vertices);evaluated.to_mesh_clear()
if not surface_vertices:raise RuntimeError('Density probe produced an empty mesh')
if not font.font.packed_file or not font.font_bold.packed_file:raise RuntimeError('Probe fonts were not packed')
samples=scene.cycles.samples if a.engine=='cycles' else scene.eevee.taa_render_samples
if samples!=a.samples:raise RuntimeError('Requested render sample count was not applied')
result={'version':bpy.app.version_string,'python':sys.version,'numpy':np.__version__,'asset_atoms_loaded':len(points),
        'engine':scene.render.engine,'cycles_device':scene.cycles.device if a.engine=='cycles' else None,
        'requested_samples':a.samples,'actual_samples':samples,
        'geometry_nodes_executed':True,'surface_vertices':surface_vertices,'fonts_loaded_and_packed':True,
        'neutral_compositor':compositor,'resolution':[320,180]}
(p.parent/'blender-result.json').write_text(json.dumps(result,indent=2)+'\n')
'''


def executable(value):
    candidate = Path(value).expanduser()
    resolved = str(candidate.resolve()) if candidate.is_file() else shutil.which(value)
    if not resolved: raise ValueError('Executable not found: '+value)
    return resolved


def font_info(path, sample):
    from fontTools.ttLib import TTFont, TTLibError
    try:
        with TTFont(path, fontNumber=0) as font:
            family = sorted({record.toUnicode() for record in font['name'].names if record.nameID in (1,2,4)})
            cmap = font.getBestCmap() or {}
            missing = sorted({c for c in sample if not c.isspace() and ord(c) not in cmap})
    except (TTLibError, KeyError) as error:
        raise ValueError('Cannot read the selected font: '+str(error)) from error
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'family_and_style':family,'sample':sample,'missing_sample_characters':missing,
            'limits':'First font face and the listed sample only; review actual project labels separately.'}


def command(argv, seconds=120):
    result = subprocess.run(argv, capture_output=True, text=True, timeout=seconds)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}): '+str(argv[0])+'\n'+(result.stderr+'\n'+result.stdout)[-3000:])
    return result.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--blender', required=True); p.add_argument('--ffmpeg', required=True)
    p.add_argument('--font-zh', type=Path, required=True); p.add_argument('--font-en', type=Path, required=True)
    p.add_argument('--engine', choices=['eevee','cycles'], default='eevee')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); skill = Path(__file__).resolve().parents[1]; output = a.output.expanduser().resolve()
    if output.exists(): p.error('Use a new runtime report path')
    if output.is_relative_to(skill): p.error('Runtime output belongs in the project, outside the skill')
    report = {'passed':False,'python_executable':sys.executable,'python_version':sys.version,'python_prefix':sys.prefix,
              'packages':{},'fonts':{},'executables':{},'engine_requested':a.engine,'failures':[],
              'scope':'Prepare a bundled DNA-chain asset with the shared source/surface code, load it into the selected Blender density/typography/neutral compositor, and encode/decode three technical frames.',
              'limits':'No mechanism, source suitability, project-specific dependencies, long render, animation smoothness, alternate filter or GPU performance certification. Custom font choices are reported, not rejected by family name.'}
    for name in ('numpy','scipy','skimage','gemmi','fontTools'):
        try:
            module = importlib.import_module(name)
            report['packages'][name] = getattr(module,'__version__','unknown')
        except (ImportError, OSError) as error: report['failures'].append(name+': '+str(error))
    for key, value in [('blender',a.blender),('ffmpeg',a.ffmpeg)]:
        try: report['executables'][key] = executable(value)
        except ValueError as error: report['failures'].append(str(error))
    for key, path, sample in [('zh',a.font_zh,'细胞核膜囊泡'),('en',a.font_en,'ATP DNA RNA 0123456789+-()')]:
        try:
            info = font_info(path.expanduser().resolve(),sample); report['fonts'][key] = info
            if info['missing_sample_characters']: report['failures'].append('Font '+key+' lacks '+''.join(info['missing_sample_characters']))
        except (ImportError, OSError, ValueError) as error: report['failures'].append('Font '+key+': '+str(error))
    output.parent.mkdir(parents=True, exist_ok=True)
    if not report['failures']:
        try:
            # Import the actual preparation module, including compiled surface dependencies.
            import prepare_molecular_asset
            import numpy as np
            source = skill/'assets/1BNA.pdb'; table = prepare_molecular_asset.atom_table(source)
            if not table: raise ValueError('Bundled runtime source has no selected polymer atoms')
            chain, rows = next(iter(table.items())); points = prepare_molecular_asset.xyz(rows)
            vertices, faces = prepare_molecular_asset.surface(points,.8,1.15,.55)
            if not len(faces) or not np.isfinite(vertices).all(): raise ValueError('Invalid prepared probe surface')
            report['preparation'] = {'source':source.name,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                                     'author_chain':chain,'atoms':len(points),'surface_vertices':len(vertices),'surface_faces':len(faces)}
            with tempfile.TemporaryDirectory(prefix='blender-runtime-',dir=output.parent) as temp:
                temp = Path(temp); spec = {'skill':str(skill),'engine':a.engine,
                    'font_zh':report['fonts']['zh']['path'],'font_en':report['fonts']['en']['path']}
                (temp/'input.json').write_text(json.dumps(spec))
                np.savez_compressed(temp/'probe.npz',atoms=(points-points.mean(0))*.1)
                (temp/'probe.py').write_text(PROBE)
                command([report['executables']['blender'],'-b','--factory-startup','--python-exit-code','1',
                         '--python',str(temp/'probe.py'),'--',str(temp/'input.json')])
                report['blender'] = json.loads((temp/'blender-result.json').read_text())
                if report['blender']['asset_atoms_loaded'] != len(points): raise ValueError('Python/Blender asset handoff mismatch')
                from render_cache import valid_png
                if not valid_png(temp/'probe.png',320,180): raise ValueError('Invalid rendered probe PNG')
                frames = temp/'frames'; frames.mkdir()
                for index in range(1,4): shutil.copyfile(temp/'probe.png',frames/f'frame_{index:04d}.png')
                (temp/'manifest.json').write_text(json.dumps({'frames':3,'fps':24}))
                ffmpeg = report['executables']['ffmpeg']
                version = command([ffmpeg,'-version']).splitlines()
                if not version or not version[0].startswith('ffmpeg version'):
                    raise ValueError('Selected executable did not identify itself as FFmpeg')
                report['ffmpeg_version'] = version[0]
                command([sys.executable,str(skill/'scripts/encode_video.py'),str(temp),'--ffmpeg',ffmpeg])
                command([sys.executable,str(skill/'scripts/audit_video.py'),str(temp/'preview.mp4'),
                         '--output',str(temp/'video-check'),'--ffmpeg',ffmpeg,'--fps','24','--frames','3'])
                video = json.loads((temp/'video-check/video-audit.json').read_text())
                report['video_probe'] = {k:video[k] for k in ['decoded_frames','fps','resolution','all_pts_consecutive','decode_errors']}
                report['video_probe']['content'] = 'Three identical technical frames; no mechanism animation claimed.'
            report['temporary_media_removed'] = True
        except (ImportError, OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            report['failures'].append(str(error))
    report['passed'] = not report['failures']
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'report':str(output),'failures':report['failures']},ensure_ascii=False))
    if not report['passed']: raise SystemExit(1)


if __name__ == '__main__': main()
