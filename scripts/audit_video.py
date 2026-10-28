"""Fully decode a movie, check timing, and extract actual encoded review frames."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('video',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ffmpeg',default=shutil.which('ffmpeg'));p.add_argument('--fps',type=int,required=True)
    p.add_argument('--frames',type=int,required=True);p.add_argument('--keyframes',default='')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh audit directory')
    a.output.mkdir(parents=True)
    digest=a.output/'frames.framemd5'
    process=subprocess.run([a.ffmpeg,'-v','error','-i',str(a.video),'-map','0:v:0','-f','framemd5',str(digest)],capture_output=True,text=True)
    (a.output/'decode.log').write_text(process.stderr)
    process.check_returncode()
    contents=digest.read_text();records=[r.split(',') for r in contents.splitlines() if r and not r.startswith('#')]
    tb=re.search(r'#tb 0: (\d+)/(\d+)',contents)
    dims=re.search(r'#dimensions 0: (\d+)x(\d+)',contents)
    assert len(records)==a.frames,(len(records),a.frames)
    assert [int(r[2]) for r in records]==list(range(a.frames))
    assert tb and int(tb[1])*a.fps==int(tb[2]),'Incorrect time base'
    indices=[int(x) for x in a.keyframes.split(',')] if a.keyframes else sorted({round(i*(a.frames-1)/8) for i in range(9)})
    assert all(0<=x<a.frames for x in indices)
    selection='+'.join(f'eq(n\\,{i})' for i in indices)
    subprocess.run([a.ffmpeg,'-v','error','-i',str(a.video),'-vf','select='+selection,'-fps_mode','vfr',str(a.output/'key-%03d.png')],check=True)
    result={'video':str(a.video.resolve()),'sha256':hashlib.sha256(a.video.read_bytes()).hexdigest(),
            'decoded_frames':len(records),'fps':a.fps,'seconds':a.frames/a.fps,
            'resolution':[int(dims[1]),int(dims[2])] if dims else None,
            'all_pts_consecutive':True,'decode_errors':process.stderr,'keyframe_indices_zero_based':indices,
            'visual_review_status':'required; successful decode does not establish biological or visual quality'}
    (a.output/'video-audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
