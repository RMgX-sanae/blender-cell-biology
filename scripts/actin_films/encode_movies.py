"""Encode actual rendered 24 fps chapters and assemble narrated route films."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
from project_context import FFMPEG

def run(args):subprocess.run([FFMPEG,'-hide_banner','-loglevel','error','-y']+args,check=True)
def quoted(path):return "'"+str(path).replace("'","'\\''")+"'"

def main():
    p=argparse.ArgumentParser();p.add_argument('--chapters',default='all');p.add_argument('--routes',action='store_true')
    args=p.parse_args();d=json.loads((ROOT/'storyboard.json').read_text());byid={s['id']:s for s in d['shots']}
    for s in d['shots']:
        if args.chapters!='all' and s['id'] not in args.chapters.split(','):continue
        folder=ROOT/'shots'/s['id'];idx=folder/'render-index.json'
        if not idx.exists():continue
        index=json.loads(idx.read_text());out=folder/'chapter.mp4'
        encoded=folder/'encoded-scene-sha256.txt'
        signature='image2-v3:'+str(index['fps'])+':'+OPTIONS['narration']+':'+index['scene_sha256']
        if out.exists() and encoded.exists() and encoded.read_text()==signature:continue
        seq=folder/'sequence';seq.mkdir(exist_ok=True)
        for i,source in enumerate(index['paths'],1):
            link=seq/f'{i:06d}.png'
            if link.exists() and not os.path.samefile(source,link):link.unlink()
            if not link.exists():os.link(source,link)
        run(['-framerate',str(index['fps']),'-i',str(seq/'%06d.png'),'-an','-frames:v',str(s['frames']),
             '-c:v','libx264','-preset','medium','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(out)])
        if OPTIONS['narration'] != 'none':
            run(['-i',str(ROOT/'audio'/(s['id']+'.wav')),'-af',f"adelay=400:all=1,apad=whole_dur={s['seconds']},atrim=duration={s['seconds']}",
                 '-ar','48000','-ac','1','-c:a','pcm_s16le',str(folder/'chapter.wav')])
        encoded.write_text(signature);print('ENCODED',s['id'],flush=True)
    if not args.routes:return
    for name,route in d['routes'].items():
        assert all((ROOT/'shots'/id/'chapter.mp4').exists() for id in route['shots']),name+' chapters missing'
        folder=ROOT/'movies';video_list=folder/(name+'-video.ffconcat');audio_list=folder/(name+'-audio.ffconcat')
        video_list.write_text('ffconcat version 1.0\n'+''.join('file '+quoted(ROOT/'shots'/id/'chapter.mp4')+'\n' for id in route['shots']))
        audio_list.write_text('ffconcat version 1.0\n'+''.join('file '+quoted(ROOT/'shots'/id/'chapter.wav')+'\n' for id in route['shots']))
        times=[];clock=0;meta=';FFMETADATA1\ntitle='+route['title']+'\ncomment=Source-derived structures; authored teaching motion. '+('Chinese synthetic narration.' if OPTIONS['narration']=='macos' else 'Optional supplied or no narration.')+'\n'
        subtitle=[];n=1
        def stamp(t):
            ms=round(t*1000);h,ms=divmod(ms,3600000);m,ms=divmod(ms,60000);sec,ms=divmod(ms,1000)
            return f'{h:02}:{m:02}:{sec:02},{ms:03}'
        for id in route['shots']:
            s=byid[id];end=clock+s['seconds'];times.append({'id':id,'start_seconds':clock,'end_seconds':end,'title':s['title']})
            meta+=f"[CHAPTER]\nTIMEBASE=1/1000\nSTART={clock*1000}\nEND={end*1000}\ntitle={s['title']}\n"
            for j,beat in enumerate(s['beats']):
                boundaries=s.get('beat_boundaries',[0,1/3,2/3,1])
                a=clock+boundaries[j]*s['seconds'];b=clock+boundaries[j+1]*s['seconds']
                subtitle.append(f'{n}\n{stamp(a)} --> {stamp(b)}\n{beat}\n');n+=1
            clock=end
        (folder/(name+'.chapters.json')).write_text(json.dumps(times,ensure_ascii=False,indent=2)+'\n')
        (folder/(name+'.srt')).write_text('\n'.join(subtitle))
        metadata=folder/(name+'.ffmetadata');metadata.write_text(meta)
        silent=folder/(name+'-picture.mp4')
        run(['-safe','0','-f','concat','-i',str(video_list),'-map','0:v','-c:v','copy','-an',str(silent)])
        out=folder/('Arp23_完整循环.mp4' if name=='arp23' else 'Formin_完整循环.mp4')
        if OPTIONS['narration'] == 'none':
            run(['-i',str(silent),'-i',str(metadata),'-map','0:v','-map_metadata','1','-map_chapters','1',
                 '-c:v','copy','-an','-movflags','+faststart',str(out)])
        else:
            voice=folder/(name+'-narration.wav')
            run(['-safe','0','-f','concat','-i',str(audio_list),'-c:a','pcm_s16le',str(voice)])
            run(['-i',str(silent),'-i',str(voice),'-i',str(metadata),'-map','0:v','-map','1:a','-map_metadata','2','-map_chapters','2',
                 '-c:v','copy','-c:a','aac','-b:a','160k','-shortest','-movflags','+faststart',str(out)])
        print('MOVIE_COMPLETE',name,clock,out,flush=True)

if __name__=='__main__':main()
