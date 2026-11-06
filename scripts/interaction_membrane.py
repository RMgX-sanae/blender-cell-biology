"""Planar membrane constraints and an explicitly schematic lipid display."""
import numpy as np


def plane_clearance(points,poses,membrane,radius=0):
    normal=np.asarray(membrane['normal'],float);origin=np.asarray(membrane['origin_nm'],float)
    if normal.shape!=(3,) or not np.isfinite(normal).all() or not np.isclose(np.linalg.norm(normal),1):
        raise ValueError('Require a finite unit cytosolic membrane normal')
    if origin.shape!=(3,) or not np.isfinite(origin).all() or not np.isfinite(membrane['hydrophobic_half_width_nm']) or membrane['hydrophobic_half_width_nm']<=0:
        raise ValueError('Require a finite plane origin and positive core half-width')
    directions=poses[:,:3,:3].transpose(0,2,1)@normal
    offset=(poses[:,:3,3]-origin)@normal-membrane['hydrophobic_half_width_nm']-radius
    return np.min(points@directions.T,axis=0)+offset


def tail_paths(head,side,angle,half_width):
    u=np.array([np.cos(angle),np.sin(angle),0]);v=np.array([-u[1],u[0],0])
    paths=[]
    for branch in [-1,1]:
        t=np.linspace(0,1,9)
        points=np.asarray(head)+u[None]*(branch*.105+.09*np.sin(np.pi*t))[:,None]
        points+=v[None]*(.04*np.sin(4*np.pi*t))[:,None]
        points[:,2]-=side*(.22+(half_width-.4)*t)
        paths.append(points)
    return paths


def display_parts(membrane):
    # Imported here so preparation and plane checks also run outside Blender.
    import detail_scene as ds
    if not np.allclose(membrane['normal'],[0,0,1],atol=1e-9,rtol=0):
        raise ValueError('This lipid display uses the source Z-axis membrane frame')
    tether_points=np.asarray(membrane['cargo_tether_nm'])
    lipid=membrane['lipids'][membrane['cargo_lipid_index']]
    if lipid[3]!=1 or not np.allclose(tether_points[-1],lipid[:3],atol=1e-9,rtol=0):
        raise ValueError('Cargo tether lost its cytosolic lipid identity')
    heads,tails,tether=ds.Parts(),ds.Parts(),ds.Parts()
    for x,y,z,side,angle in membrane['lipids']:
        heads.ball((x,y,z),membrane['head_radius_nm'],seg=12,rings=8)
        for points in tail_paths((x,y,z),side,angle,membrane['head_offset_nm']):
            tails.tube(points,membrane['tail_radius_nm'],sides=6)
    tether.tube(membrane['cargo_tether_nm'],.035,sides=8)
    return {'Membrane heads':heads,'Membrane tails':tails,'Cargo tether':tether}


def build_display(membrane,group,profile):
    import style_nodes as sn
    colors={'Membrane heads':'#638fa0','Membrane tails':'#365f72','Cargo tether':profile['palette'][0]}
    objects=[]
    for name,parts in display_parts(membrane).items():
        mat=sn.material(name,{**profile,'roughness':.8},colors[name],'pbr')
        obj=parts.object(name,[mat],group);obj['membrane_display_role']=name
        objects.append(obj)
    return objects


def audit_display(scene,membrane):
    for name,parts in display_parts(membrane).items():
        obj=scene.objects.get(name)
        if obj is None or obj.hide_render or not obj.visible_get() or obj.modifiers:
            raise ValueError('Membrane display is missing, hidden or modified: '+name)
        actual=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',actual)
        if actual.size!=len(parts.v)*3 or not np.allclose(actual.reshape(-1,3),parts.v,atol=2e-6,rtol=0) or [tuple(f.vertices) for f in obj.data.polygons]!=parts.f or not np.allclose(np.asarray(obj.matrix_world),np.eye(4),atol=1e-8,rtol=0):
            raise ValueError('Saved membrane geometry differs from prepared display: '+name)
    return {'passed':True,'objects':3,'lipid_display_count':len(membrane['lipids']),
            'limits':'Authored lipid symbols and unresolved tether, not measured molecular packing.'}
