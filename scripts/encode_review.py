"""Encode uniformly sampled review frames at their original source timing.

This makes a short or sparse review, not a complete-film verification. Use
encode_video.py for an explicitly requested complete contiguous deliverable.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from render_cache import valid_png


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,help='Directory made by render_saved.py')
    parser.add_argument('--ffmpeg',default=shutil.which('ffmpeg'))
    args=parser.parse_args();root=args.directory.resolve()
    report=json.loads((root/'render-report.json').read_text())
    signature=json.loads((root/'render-signature.json').read_text())
    frames=report['frames']
    if len(frames)<2 or any(type(f)!=int for f in frames):parser.error('Supply at least two sampled integer frames')
    step=frames[1]-frames[0]
    if step<=0 or any(b-a!=step for a,b in zip(frames,frames[1:])):
        parser.error('Use one uniformly sampled interval; keep irregular keyframes as stills')
    if report['linear_multilayer_exr']:parser.error('Review encoding requires display-referred PNG frames')
    if report['fps']!=signature['fps'] or report.get('fps_base',1)!=signature.get('fps_base',1):
        parser.error('Render report and signature timing differ')
    if hashlib.sha256(Path(report['source_blend']).read_bytes()).hexdigest()!=signature['scene_sha256']:
        parser.error('Source scene changed since these frames were rendered')
    fps=(Fraction(str(report['fps']))/Fraction(str(report.get('fps_base',1)))/step).limit_denominator(1000000)
    if fps<=0 or not args.ffmpeg:parser.error('Supply valid source timing and an FFmpeg executable')
    x,y,percent=signature['resolution'];width,height=x*percent//100,y*percent//100
    if width%2 or height%2:parser.error('H.264 review dimensions must be even')
    images=[root/f'frame_{f:04d}.png' for f in frames]
    if not all(valid_png(p,width,height) for p in images):parser.error('A selected PNG is absent, corrupt or has the wrong dimensions')
    target=root/'review.mp4';receipt=root/'review.json'
    if target.exists() or receipt.exists():parser.error('Review output already exists')
    with tempfile.TemporaryDirectory(prefix='review-input-',dir=root) as temporary:
        temporary=Path(temporary)
        for index,path in enumerate(images):
            (temporary/f'{index:06d}.png').symlink_to(path)
        subprocess.run([str(args.ffmpeg),'-v','error','-n','-framerate',str(fps),
                        '-start_number','0','-i',str(temporary/'%06d.png'),
                        '-c:v','libx264','-crf','20','-pix_fmt','yuv420p',
                        '-movflags','+faststart',str(target)],check=True)
    result={'kind':'sampled-animation-review','source_blend':report['source_blend'],
            'source_scene_sha256':signature['scene_sha256'],'source_frames':frames,
            'source_seconds':report.get('source_seconds'),'source_frame_step':step,
            'source_fps':str(fps*step),'review_fps':str(fps),'review_frames':len(frames),
            'duration_seconds':float(Fraction(len(frames))/fps),'resolution':[width,height],
            'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'limits':'Spacing between sampled source times is preserved; the last sample is held for one interval. No interpolation, audio, unseen-frame inspection or full-film validation.'}
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
