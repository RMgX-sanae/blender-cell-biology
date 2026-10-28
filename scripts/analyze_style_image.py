"""Measure a reference image to assist, not replace, visual style interpretation.

Uses Pillow and numpy. Crop away website UI before palette measurements. These
statistics cannot identify a renderer, material roughness or biological truth.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image


def measure(path, crop=None):
    image = Image.open(path).convert('RGB')
    if crop:
        image = image.crop(crop)
    dimensions = list(image.size)
    image.thumbnail((640, 640))
    a = np.asarray(image, dtype=float)/255
    luma = a @ np.array([.2126, .7152, .0722])
    saturation = (a.max(2)-a.min(2))/np.maximum(a.max(2), .001)
    gradient = np.hypot(np.gradient(luma, axis=0), np.gradient(luma, axis=1))
    quantized = image.quantize(colors=8)
    palette = quantized.getpalette()
    swatches = []
    for count, i in sorted(quantized.getcolors(), reverse=True):
        rgb = palette[3*i:3*i+3]
        swatches.append({'srgb_hex': '#'+''.join(f'{v:02x}' for v in rgb),
                         'fraction': round(count/(image.width*image.height), 4)})
    return {'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'crop_pixels': crop, 'cropped_dimensions': dimensions,
            'luma_percentiles': dict(zip(('p05','p50','p95'), np.quantile(luma, [.05,.5,.95]).tolist())),
            'mean_saturation': float(saturation.mean()), 'strong_edge_fraction': float((gradient>.08).mean()),
            'palette': swatches,
            'limitations': 'Image appearance statistics, not recovered shader parameters. UI, labels and highlights can bias these values.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path)
    p.add_argument('--crop',type=int,nargs=4,metavar=('LEFT','TOP','RIGHT','BOTTOM'))
    p.add_argument('--candidate',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():p.error('Choose a fresh report filename.')
    report={'reference':measure(args.reference,args.crop),
            'visual_review_required':['shape hierarchy and scale','flat versus shaded surfaces','edge treatment',
                                      'light direction and softness','depth layers and occlusion','camera and framing',
                                      'biology versus artistic assumptions'],
            'observations_to_be_filled_after_viewing':{},
            'implementation_hypotheses_to_be_tested':{}}
    if args.candidate:report['candidate']=measure(args.candidate)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(args.output)


if __name__=='__main__':main()
