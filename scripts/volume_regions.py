"""Calibrate one structure's volume fraction against physically spaced nuclear-distance shells."""
import numpy as np


def placement_mask(data,policy):
    if policy=='remaining-cytoplasm':return data['available']
    if policy=='replace-measured-structure':
        cell=data['cell'];dna=data['dna'];target=data['structure']
        if (cell.ndim!=3 or cell.shape!=dna.shape or cell.shape!=target.shape or
                any(v.dtype!=bool for v in (cell,dna,target))):
            raise ValueError('Use three co-registered boolean ZYX masks')
        return cell & ~dna
    raise ValueError('Unknown placement compartment '+str(policy))


def nuclear_shells(data,edges):
    from scipy.ndimage import distance_transform_edt
    edges=np.array(edges,dtype=float)
    if edges.ndim!=1 or not len(edges) or not np.isfinite(edges).all() or np.any(edges<=0) or np.any(np.diff(edges)<=0):
        raise ValueError('Supply strictly increasing positive shell edges in micrometers')
    mask=placement_mask(data,'replace-measured-structure')
    dna=data['dna'];observed=data['structure']
    if not dna.any() or not np.any(mask & observed):
        raise ValueError('The grid requires a nucleus and a measured target in cytoplasm')
    distance=distance_transform_edt(~dna,sampling=data['spacing_xyz_um'][::-1])
    shell=np.searchsorted(edges,distance,side='right').astype(np.int32)
    num=len(edges)+1
    volume=np.bincount(shell[mask],minlength=num)
    target=np.bincount(shell[mask & observed],minlength=num)
    fraction=np.divide(target,volume,out=np.zeros(num),where=volume>0)
    enrichment=fraction/(target.sum()/volume.sum())
    stats={'kind':'nuclear-distance volume fraction','edges_um':edges.tolist(),
           'cytoplasm_voxels':volume.tolist(),'target_voxels':target.tolist(),
           'target_voxels_outside_cytoplasm':int(np.sum(observed & ~mask)),
           'target_volume_fraction':[float(v) if n else None for v,n in zip(fraction,volume)],
           'fold_over_global_fraction':[float(v) if n else None for v,n in zip(enrichment,volume)],
           'distance_definition':'Voxel-center distance to nearest segmented nuclear voxel center. Left-closed shells, including an unbounded last shell.',
           'observable':'Single-cell segmented target volume per cytoplasmic volume. Rigid-body count, orientation and packing are modeled.'}
    return mask,shell,stats


def regional_counts(count,target_voxels):
    if type(count)!=int or count<1:raise ValueError('Use a positive model count')
    target=np.asarray(target_voxels,dtype=float)
    if target.ndim!=1 or not len(target) or not np.isfinite(target).all() or np.any(target<0) or target.sum()<=0:
        raise ValueError('Use positive total measured occupancy in disjoint regions')
    shares=count*target/target.sum();counts=np.floor(shares).astype(int)
    counts[np.argsort(-(shares-counts),kind='stable')[:count-counts.sum()]]+=1
    return counts



def save_component_surfaces(masks,spacing,origin,rank,output,region_hash):
    """Save registered context and one intact labeled component, all in source micrometers."""
    import hashlib
    import json
    from scipy.ndimage import gaussian_filter,label,map_coordinates
    from skimage.measure import marching_cubes
    pool=masks['structure'];components,total=label(pool)
    counts=np.bincount(components.ravel())
    border=np.concatenate([components[0].ravel(),components[-1].ravel(),components[:,0].ravel(),
                           components[:,-1].ravel(),components[:,:,0].ravel(),components[:,:,-1].ravel()])
    excluded=set(border)|set(np.unique(components[~masks['cell'] | masks['dna']]))
    eligible=sorted((i for i in range(1,total+1) if i not in excluded),key=lambda i:(-counts[i],i))
    if not 1<=rank<=len(eligible):raise ValueError(f'Choose a component rank from 1 to {len(eligible)}')
    chosen=eligible[rank-1];arrays={}
    for role,mask in [('cell',masks['cell']),('nucleus',masks['dna']),('measured',pool),
                      ('body',components==chosen)]:
        zyx,faces,_,_=marching_cubes(np.pad(mask,1).astype(np.uint8),level=.5,allow_degenerate=False)
        # Smooth the lighting normals only. Preserve every source isosurface vertex.
        blurred=gaussian_filter(np.pad(mask,1).astype(np.float32),sigma=1.4)
        grad=np.gradient(blurred,*spacing[::-1])
        normals=-np.column_stack([map_coordinates(axis,zyx.T,order=1) for axis in grad[::-1]])
        lengths=np.linalg.norm(normals,axis=1)
        normals=np.divide(normals,lengths[:,None],out=np.zeros_like(normals),where=lengths[:,None]>1e-8)
        arrays[role+'_normals']=normals.astype(np.float32)
        vertices=origin+(zyx[:,::-1]-1)*spacing
        # Converting ZYX to XYZ reverses handedness; use outward physical normals.
        faces=faces[:,::-1].copy()
        volume=np.einsum('ij,ij->i',vertices[faces[:,0]],np.cross(vertices[faces[:,1]],vertices[faces[:,2]])).sum()
        if volume<0:faces=faces[:,::-1].copy()
        arrays[role+'_vertices']=vertices.astype(np.float32);arrays[role+'_faces']=faces.astype(np.int32)
    path=output/'surfaces.npz';np.savez_compressed(path,**arrays)
    info={'shading_normals':'Outward gradient of a 1.4-voxel Gaussian for lighting only; source surface coordinates preserved.',
          'units':'micrometer','source_masks_sha256':region_hash,
          'surfaces_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'component_rank':rank,
          'eligible_complete_components':len(eligible),'selected_source_label':int(chosen),
          'selected_source_voxels':int(counts[chosen]),
          'geometry':'Registered 0.5-level voxel surfaces. Component ranked by volume among face-connected source components fully inside cytoplasm, away from crop faces. A component need not equal an entire organelle.',
          'context_crop':'Exterior padding closes the display surface at any crop boundary. Source masks and regional measurements are unchanged.'}
    (output/'surfaces.json').write_text(json.dumps(info,indent=2)+'\n')
