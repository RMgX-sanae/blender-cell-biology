"""Verify that complex-scene surface subunits remain intact and unmoved."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    results=[]
    for obj in bpy.context.scene.objects:
        if obj.type!='MESH' or 'source_chain' not in obj:continue
        chain=obj['source_chain']
        data=np.load(args.asset_dir/f'chain_{chain}.npz',allow_pickle=False)
        coords=np.empty(len(obj.data.vertices)*3,dtype=np.float32)
        obj.data.vertices.foreach_get('co',coords)
        loops=np.empty(len(obj.data.loops),dtype=np.int32)
        obj.data.loops.foreach_get('vertex_index',loops)
        intact=(np.array_equal(coords.reshape(-1,3),data['vertices']) and
                np.array_equal(loops,data['faces'].ravel()))
        unmoved=bool(np.allclose(np.array(obj.matrix_world),np.eye(4)))
        results.append({'chain':chain,'intact_surface_topology':intact,
                        'unmoved':unmoved,'hidden_for_view':obj.hide_render})
    report={'blend':bpy.data.filepath,'subunits':results,
            'limits':'Checks raw source-derived surface meshes and transforms. Rendering may smooth envelopes; this does not validate a biological state.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2))
    if not results or not all(r['intact_surface_topology'] and r['unmoved'] for r in results):
        raise RuntimeError('Intact-molecule audit failed; inspect the report.')
    print(json.dumps({'intact_unmoved_surface_subunits':len(results),
                      'hidden_complete_subunits':[r['chain'] for r in results if r['hidden_for_view']]}))


if __name__=='__main__':main()
