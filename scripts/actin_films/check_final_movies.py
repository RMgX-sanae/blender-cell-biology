"""Decode the delivered MP4 files and extract actual encoded-frame review sheets."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
from PIL import Image,ImageDraw,ImageFont

from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
from project_context import FFMPEG
data=json.loads((ROOT/'storyboard.json').read_text());shots={s['id']:s for s in data['shots']}
font=ImageFont.truetype(FONT,18)
reports=[]
for name,route in data['routes'].items():
    movie=ROOT/'movies'/('Arp23_完整循环.mp4' if name=='arp23' else 'Formin_完整循环.mp4')
    output=ROOT/'validation'/name;output.mkdir(exist_ok=True)
    md5=output/'decoded-frames.framemd5'
    command=[FFMPEG,'-hide_banner','-loglevel','error','-y','-i',str(movie),'-map','0:v:0','-f','framemd5',str(md5)]
    decoded=subprocess.run(command,capture_output=True,text=True,check=True)
    records=[line.split(',') for line in md5.read_text().splitlines() if line and not line.startswith('#')]
    expected=sum(shots[id]['frames'] for id in route['shots']);assert len(records)==expected
    assert [int(x[2]) for x in records]==list(range(expected))
    timebase=re.search(r'#tb 0: (\d+)/(\d+)',md5.read_text())
    assert timebase and int(timebase[1])*OPTIONS['fps']==int(timebase[2])
    maxvol=meanvol=None
    if OPTIONS['narration'] != 'none':
        sound=subprocess.run([FFMPEG,'-hide_banner','-i',str(movie),'-map','0:a:0','-af','volumedetect','-f','null','-'],capture_output=True,text=True,check=True)
        (output/'audio-check.log').write_text(sound.stderr)
        maxvol=float(re.search(r'max_volume: ([\-\d.]+) dB',sound.stderr).group(1))
        meanvol=float(re.search(r'mean_volume: ([\-\d.]+) dB',sound.stderr).group(1))
        assert maxvol<0 and meanvol>-60
    indices=[];labels=[];offset=0
    for id in route['shots']:
        s=shots[id]
        for fraction in [.05,.48,.85]:
            indices.append(offset+round((s['frames']-1)*fraction));labels.append(f"{id}  {round(fraction*100)}%  {s['title']}")
        offset+=s['frames']
    expression='+'.join(f'eq(n\\,{n})' for n in indices)
    subprocess.run([FFMPEG,'-hide_banner','-loglevel','error','-y','-i',str(movie),'-vf','select='+expression,
                    '-fps_mode','vfr',str(output/'frame-%03d.png')],check=True)
    assert len(list(output.glob('frame-*.png')))==len(indices)
    actual_resolution=list(Image.open(output/'frame-001.png').size)
    assert actual_resolution==[OPTIONS['width'],2*round(OPTIONS['width']*9/32)]
    for first in range(0,len(indices),12):
        group=list(range(first,min(first+12,len(indices))));rows=(len(group)+2)//3
        sheet=Image.new('RGB',(1440,rows*294),'#f2f2f2');draw=ImageDraw.Draw(sheet)
        for j,index in enumerate(group):
            row,col=divmod(j,3);img=Image.open(output/f'frame-{index+1:03d}.png').convert('RGB');img.thumbnail((480,270))
            sheet.paste(img,(480*col,294*row+24));draw.text((480*col+5,294*row+2),labels[index][:39],font=font,fill='#111111')
        sheet.save(output/f'review-{first//12+1}.jpg',quality=93)
    reports.append({'route':name,'file':str(movie),'sha256':hashlib.sha256(movie.read_bytes()).hexdigest(),
                    'bytes':movie.stat().st_size,'decoded_frames':len(records),'expected_frames':expected,
                    'duration_seconds':expected/OPTIONS['fps'],'resolution':actual_resolution,'fps':OPTIONS['fps'],
                    'video_decode_errors':decoded.stderr,'audio_mean_db':meanvol,'audio_peak_db':maxvol,
                    'review_images':len(indices),'all_pts_consecutive':True})
    focus='a03' if name=='arp23' else 'f04'
    start=sum(shots[id]['frames'] for id in route['shots'][:route['shots'].index(focus)])
    fractions=[.05,.12,.18,.24,.29,.33,.40,.47,.54,.61,.68,.75,.82,.90,.97]
    indices=[start+round((shots[focus]['frames']-1)*t) for t in fractions]
    expression='+'.join(f'eq(n\\,{n})' for n in indices)
    subprocess.run([FFMPEG,'-hide_banner','-loglevel','error','-y','-i',str(movie),'-vf','select='+expression,
                    '-fps_mode','vfr',str(output/'motion-%03d.png')],check=True)
    sheet=Image.new('RGB',(1440,5*294),'#f2f2f2');draw=ImageDraw.Draw(sheet)
    for j,t in enumerate(fractions):
        row,col=divmod(j,3);img=Image.open(output/f'motion-{j+1:03d}.png').convert('RGB');img.thumbnail((480,270))
        sheet.paste(img,(480*col,294*row+24));draw.text((480*col+5,294*row+2),f'{focus}  {t:.0%}',font=font,fill='#111111')
    sheet.save(output/'motion-review.jpg',quality=94)
(ROOT/'validation/final-movies.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(reports,ensure_ascii=False,indent=2))
