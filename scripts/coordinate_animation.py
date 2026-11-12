"""Preserve checked coordinate samples and their timing in Blender shape keys.

This is a display handoff, not a motion planner or a chemistry/contact check.
Caller-supplied frames are presentation times and are never silently retimed.
"""
import bpy
import numpy as np


def _inputs(coordinates, frames):
    coordinates=np.asarray(coordinates, dtype=float)
    if (coordinates.ndim!=3 or coordinates.shape[2]!=3 or min(coordinates.shape[:2])<1
            or not np.isfinite(coordinates).all()):
        raise ValueError('Require finite coordinates shaped (samples, points, 3) with fixed point identity')
    frames=np.arange(1,len(coordinates)+1,dtype=float) if frames is None else np.asarray(frames,dtype=float)
    if frames.shape!=(len(coordinates),) or not np.isfinite(frames).all() or np.any(np.diff(frames)<=0):
        raise ValueError('Supply one strictly increasing finite frame for each coordinate sample')
    with np.errstate(over='ignore'):
        points32=coordinates.astype(np.float32);frames32=frames.astype(np.float32)
    if not np.isfinite(points32).all() or not np.isfinite(frames32).all() or np.any(np.diff(frames32)<=0):
        raise ValueError('Coordinates or distinct sample frames cannot be represented by Blender float32 storage')
    return coordinates,frames,points32


def _curves(keys):
    animation=keys.animation_data
    if (animation is None or animation.action is None or animation.action_slot is None
            or animation.drivers or animation.nla_tracks):
        raise ValueError('Require an explicit coordinate action without drivers or NLA remapping')
    # Blender can group the object's visibility and its shape keys in one action.
    # Only this datablock's slot controls its coordinate evaluation time.
    return [curve for layer in animation.action.layers for strip in layer.strips
            for bag in strip.channelbags if bag.slot_handle==animation.action_slot.handle for curve in bag.fcurves]


def verify_point_animation(obj, coordinates, frames=None):
    """Check all stored keys and linear timing; evaluated scene geometry is separate."""
    coordinates,frames,points32=_inputs(coordinates,frames)
    keys=obj.data.shape_keys
    if keys is None or keys.use_relative or len(keys.key_blocks)!=len(coordinates):
        raise ValueError('Absolute shape keys must retain every coordinate sample: '+obj.name)
    observed=np.empty(points32.shape[1]*3,dtype=np.float32)
    if len(obj.data.vertices)!=len(points32[0]):raise ValueError('Point count changed: '+obj.name)
    obj.data.vertices.foreach_get('co',observed)
    if not np.array_equal(observed.reshape(-1,3),points32[0]):raise ValueError('Basis coordinates changed: '+obj.name)
    for key,expected in zip(keys.key_blocks,points32):
        if key.interpolation!='KEY_LINEAR':raise ValueError('Shape-key interpolation must remain linear')
        key.data.foreach_get('co',observed)
        if not np.array_equal(observed.reshape(-1,3),expected):
            raise ValueError('Coordinate sample changed or points were reordered: '+obj.name+'/'+key.name)
    curves=_curves(keys)
    if len(curves)!=1 or curves[0].data_path!='eval_time':raise ValueError('Unexpected coordinate animation channels')
    curve=curves[0]
    if len(curve.keyframe_points)!=len(frames):
        raise ValueError('Blender merged coordinate time keys; supply more widely spaced presentation frames. Do not discard path samples without rechecking the new path.')
    actual=np.array([tuple(key.co) for key in curve.keyframe_points],dtype=np.float32)
    expected=np.column_stack([frames,[key.frame for key in keys.key_blocks]]).astype(np.float32)
    if not np.array_equal(actual,expected):raise ValueError('Coordinate timing differs from the supplied sample frames')
    if (curve.mute or curve.modifiers or curve.extrapolation!='CONSTANT'
            or any(key.interpolation!='LINEAR' for key in curve.keyframe_points)):
        raise ValueError('Coordinate timing must use unmodified linear interpolation and endpoint holds')
    return {'samples':len(frames),'points':len(points32[0]),'first_frame':float(frames[0]),
            'last_frame':float(frames[-1]),'maximum_float32_error':float(np.max(abs(coordinates-points32))),
            'scope':'Stored coordinate keys and timing; not identity provenance, evaluated surfaces or biological validity'}


def keyed_points(name, coordinates, group, edges=(), *, frames=None):
    coordinates,frames,points32=_inputs(coordinates,frames)
    edges=np.asarray(edges)
    if edges.size:
        if (edges.ndim!=2 or edges.shape[1]!=2 or not np.issubdtype(edges.dtype,np.integer)
                or np.any(edges<0) or np.any(edges>=len(points32[0])) or np.any(edges[:,0]==edges[:,1])):
            raise ValueError('Display edges must index two distinct retained points')
    edges=edges.reshape(-1,2).astype(int)
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(points32[0].tolist(),edges.tolist(),[]);mesh.update()
    obj=bpy.data.objects.new(name,mesh);group.objects.link(obj)
    for i,points in enumerate(points32):
        key=obj.shape_key_add(name=f'frame_{i+1:04d}');key.data.foreach_set('co',points.ravel());key.interpolation='KEY_LINEAR'
    keys=mesh.shape_keys;keys.use_relative=False
    # Use each actual absolute-key time, including Blender's float32 offsets.
    for key,frame in zip(keys.key_blocks,frames):
        keys.eval_time=key.frame;keys.keyframe_insert('eval_time',frame=float(frame))
    for curve in _curves(keys):
        curve.extrapolation='CONSTANT'
        for key in curve.keyframe_points:key.interpolation='LINEAR'
        curve.update()
    verify_point_animation(obj,coordinates,frames)
    return obj
