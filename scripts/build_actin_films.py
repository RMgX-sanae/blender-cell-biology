"""Rebuild both complete actin teaching films without the original project.

Run with the scientific preparation Python (numpy/scipy/skimage/gemmi/fontTools/
Pillow). Blender is launched as a subprocess. Generated assets stay in --output.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

SKILL = Path(__file__).resolve().parents[1]
SCRIPTS = SKILL / 'scripts'
CHAPTERS = SCRIPTS / 'actin_films'
SOURCES = ['7USC', '7USE', '1K8K', '7TPT', '9B27', '2BTF', '6VAO',
           '6FM2-assembly1', '2BAP', '1FOE', '4JD2', '8F8R', '8F8Q',
           '7AQK', '3J6F', '8A2R', '8A2S', '8A2T', '1UBQ']
QUALITY = {'draft': {'width': 640, 'samples': 8, 'fps': 12},
           'review': {'width': 960, 'samples': 16, 'fps': 24},
           'final': {'width': 1600, 'samples': 32, 'fps': 24}}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(value, name):
    result = shutil.which(value or name)
    if not result:
        raise ValueError(f'{name} unavailable; pass its absolute executable path')
    return str(Path(result).resolve())


def font_file(value, candidates, destination, family):
    """Resolve installed fonts; extract one TTC face only into the project."""
    source = Path(value).expanduser() if value else next(
        (Path(p).expanduser() for p in candidates if Path(p).expanduser().is_file()), None)
    if source is None or not source.is_file():
        raise ValueError(f'{family} unavailable; specify an installed font file')
    if source.suffix.lower() != '.ttc':
        return str(source.resolve())
    from fontTools.ttLib import TTCollection
    collection = TTCollection(source)
    for font in collection.fonts:
        names = {n.toUnicode() for n in font['name'].names if n.nameID in (1, 2, 4, 6)}
        if any(family.lower() in n.lower() for n in names) and any('regular' in n.lower() for n in names):
            font.save(destination)
            collection.close()
            return str(destination)
    collection.close()
    raise ValueError(f'No {family} Regular face in {source}; supply a single-face TTF/OTF')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', default='prepare,build,render,encode,audit')
    parser.add_argument('--blender')
    parser.add_argument('--ffmpeg')
    parser.add_argument('--font-zh')
    parser.add_argument('--font-en')
    parser.add_argument('--quality', choices=QUALITY)
    parser.add_argument('--engine', choices=('eevee', 'cycles'))
    parser.add_argument('--compositor', choices=('clean', 'reference', 'monochrome'))
    parser.add_argument('--narration', choices=('none', 'macos', 'provided'))
    parser.add_argument('--voice')
    parser.add_argument('--shots', default='all')
    parser.add_argument('--preview', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    if root.is_relative_to(SKILL):
        parser.error('--output must be outside the skill directory')
    stages = args.steps.split(',')
    if set(stages) - {'prepare', 'build', 'render', 'encode', 'audit', 'preflight'}:
        parser.error('Unknown step')
    modules = ['numpy', 'scipy', 'skimage', 'gemmi', 'fontTools', 'PIL']
    missing = [m for m in modules if importlib.util.find_spec(m) is None]
    if missing:
        parser.error('Preparation Python lacks: ' + ', '.join(missing))
    root.mkdir(parents=True, exist_ok=True)
    for name in ['assets', 'sources', 'shots', 'movies', 'audio', 'validation', 'logs']:
        (root / name).mkdir(exist_ok=True)
    config = root / 'project.json'
    old = json.loads(config.read_text()) if config.exists() else {}
    options = dict(old)
    options.update(blender=executable(args.blender or old.get('blender'), 'blender'),
                   ffmpeg=executable(args.ffmpeg or old.get('ffmpeg'), 'ffmpeg'),
                   quality=args.quality or old.get('quality', 'review'),
                   engine=args.engine or old.get('engine', 'eevee'),
                   compositor=args.compositor or old.get('compositor', 'clean'),
                   narration=args.narration or old.get('narration', 'none'),
                   voice=args.voice or old.get('voice', 'Tingting'))
    options.update(QUALITY[options['quality']])
    options['font_zh'] = font_file(args.font_zh or old.get('font_zh'),
        ['/System/Library/Fonts/Supplemental/Songti.ttc'], root/'assets/Songti-Regular.ttf', 'Songti SC')
    options['font_en'] = font_file(args.font_en or old.get('font_en'),
        ['/System/Library/Fonts/Supplemental/Times New Roman.ttf', 'C:/Windows/Fonts/times.ttf'],
        root/'assets/Times-New-Roman-Regular.ttf', 'Times New Roman')
    if options['narration'] == 'macos' and not shutil.which('say'):
        parser.error('--narration macos requires say; use provided WAVs or none')
    build_keys = ('width', 'samples', 'fps', 'engine', 'compositor', 'font_zh', 'font_en', 'narration', 'voice')
    if old and any(old.get(k) != options[k] for k in build_keys):
        if 'build' not in stages and set(stages).intersection({'render', 'encode', 'audit'}):
            parser.error('Changed scene options require a build step before render/encode/audit')
    config.write_text(json.dumps(options, ensure_ascii=False, indent=2)+'\n')
    env = dict(os.environ, BIO_FILM_PROJECT=str(root))
    history = root/'validation/run-commands.json'
    record = json.loads(history.read_text()) if history.exists() else []
    previous_commands = len(record)

    def run(command, label):
        print(label, flush=True)
        log = root/'logs'/(f'{len(record):03d}-'+label+'.log')
        with log.open('w') as stream:
            result = subprocess.run([str(c) for c in command], stdout=stream, stderr=subprocess.STDOUT, env=env)
        record.append({'label': label, 'command': list(map(str, command)), 'exit_code': result.returncode, 'log': str(log)})
        (root/'validation/run-commands.json').write_text(json.dumps(record, ensure_ascii=False, indent=2))
        if result.returncode:
            raise RuntimeError(f'{label} failed; see {log}')

    def python(script, *extra):
        run([sys.executable, script, *extra], script.stem)

    def blender(script, *extra):
        run([options['blender'], '-b', '--factory-startup', '--python-exit-code', '1',
             '--python', script, '--', *extra], script.stem)

    def prepare_cached(script, folder, names, *extra):
        receipt = root/'validation'/(folder.name+'-preparation.json')
        signature = {str(p.relative_to(SKILL)): sha(p) for p in SCRIPTS.glob('prepare_*.py')}
        signature.update({name: sha(root/'sources'/(name+'.cif')) for name in names})
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if previous['inputs'] != signature:
                raise ValueError(f'Preparation changed for {folder}; use a new project directory')
            if all((folder/k).is_file() and sha(folder/k) == v for k,v in previous['outputs'].items()):
                print('CACHE_VERIFIED', folder.name, flush=True)
                return
        if folder.exists():
            folder.rename(folder.with_name(folder.name+'.incomplete-'+datetime.now().strftime('%Y%m%d%H%M%S')))
        python(script, *extra, '--output', folder)
        receipt.write_text(json.dumps({'inputs': signature, 'outputs': {str(p.relative_to(folder)): sha(p) for p in folder.rglob('*') if p.is_file()}}, indent=2))

    version = subprocess.check_output([options['blender'], '--version'], text=True).splitlines()[0]
    preflight = {'time_utc': datetime.now(timezone.utc).isoformat(), 'blender': version,
                 'python': sys.version, 'dependencies_available': modules, 'options': options,
                 'skill': str(SKILL), 'source_scripts_sha256': {str(p.relative_to(SKILL)): sha(p) for p in CHAPTERS.glob('*.py')},
                 'limits': 'Technical preflight; no mechanism or visual-quality claim.'}
    (root/'validation/preflight.json').write_text(json.dumps(preflight, ensure_ascii=False, indent=2))

    if 'prepare' in stages:
        for name in SOURCES:
            path = root/'sources'/(name+'.cif')
            provenance = path.with_suffix('.cif.source.json')
            if path.is_file() and provenance.is_file():
                if sha(path) != json.loads(provenance.read_text())['sha256']:
                    raise ValueError(f'Source hash mismatch: {path}')
                continue
            extra = ['--assembly', name.split('assembly')[1]] if 'assembly' in name else []
            run([sys.executable, SCRIPTS/'fetch_structure.py', name[:4], '--output', root/'sources', *extra], 'fetch-'+name)
        prepare_cached(SCRIPTS/'prepare_filament_assets.py', root/'assets/filaments', ['8F8R','8F8Q','7AQK','3J6F'], '--source-dir', root/'sources')
        prepare_cached(SCRIPTS/'prepare_actin_states.py', root/'assets/states', ['8A2R','8A2S','8A2T'], '--source-dir', root/'sources')
        prepare_cached(SCRIPTS/'prepare_molecular_asset.py', root/'assets/ubiquitin', ['1UBQ'], root/'sources/1UBQ.cif', '--preset', 'generic')
        python(CHAPTERS/'prepare_assets.py')
        python(CHAPTERS/'plan_clear_approaches.py')
        python(CHAPTERS/'prepare_contact_repairs.py')
        python(CHAPTERS/'prepare_fh2_seed.py')
        python(CHAPTERS/'prepare_fh2_growth.py')
    if 'build' in stages:
        python(CHAPTERS/'write_storyboard.py')
        blender(CHAPTERS/'build_shots.py', '--shots', args.shots, *(['--preview'] if args.preview else []))
    if 'render' in stages:
        blender(CHAPTERS/'render_chapters.py', '--shots', args.shots)
    if 'encode' in stages:
        python(CHAPTERS/'encode_movies.py', '--chapters', args.shots, *(['--routes'] if args.shots == 'all' else []))
    if 'audit' in stages:
        blender(CHAPTERS/'render_chapters.py', '--shots', args.shots, '--audit-only')
        blender(CHAPTERS/'audit_saved_assets.py')
        python(CHAPTERS/'check_contacts.py', *([] if args.shots=='all' else ['--shots',args.shots]))
        repaired = ['a05','a07','c05','f03','f04','f06','f07']
        strict = repaired if args.shots=='all' else [name for name in repaired if name in args.shots.split(',')]
        if strict:
            python(CHAPTERS/'check_contacts.py', '--shots', ','.join(strict),
                   '--threshold', '.20', '--fail-on-overlap', '--summary', 'strict-contacts-summary.json')
        if args.shots == 'all':
            python(CHAPTERS/'check_final_movies.py')
    print(json.dumps({'project': str(root), 'completed_commands': len(record)-previous_commands, 'steps': stages,
                      'visual_review_required': True}, ensure_ascii=False))


if __name__ == '__main__':
    main()
