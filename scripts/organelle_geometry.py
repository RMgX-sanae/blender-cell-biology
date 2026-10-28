"""Editable, procedural mitochondrion cutaway for Blender 5.2.

The inner boundary and curved cristae are ONE connected parameterized mesh.
Molecular complexes and lipids are illustrative geometry, not atomic models.
Geometry adapted from the user workspace mitochondrion_astra_20260905.
Original source is preserved. Used by build_organelle_scene.py.
"""
import math
import random

import bpy
import bmesh
from mathutils import Vector

TAU = math.tau
RNG = random.Random(20260905)
FOLDS = [(-3.78, 1.61), (-2.83, 2.34), (-1.55, 2.77), (-.15, 2.48), (1.04, 2.71), (2.07, 2.19), (3.43, 1.88)]
PHI0, PHI1 = 3.30, 6.53


def coll(name):
    c = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(c)
    return c


def mesh_object(name, verts, faces, materials, collection, indices=None):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    for mat in materials:
        mesh.materials.append(mat)
    for p in mesh.polygons:
        p.use_smooth = True
    if indices:
        for p, i in zip(mesh.polygons, indices):
            p.material_index = i
    return obj


def material(name, color, roughness=.4, metallic=.0, noise=0):
    m = bpy.data.materials.new(name)
    m.diffuse_color = (*color, 1)
    m.use_nodes = True
    n = m.node_tree.nodes
    links = m.node_tree.links
    p = n.get('Principled BSDF')
    p.inputs['Base Color'].default_value = (*color, 1)
    p.inputs['Roughness'].default_value = roughness
    p.inputs['Metallic'].default_value = metallic
    p.inputs['Subsurface Weight'].default_value = .035
    if noise:
        tex = n.new('ShaderNodeTexNoise')
        tex.inputs['Scale'].default_value = noise
        tex.inputs['Detail'].default_value = 3
        tex.inputs['Roughness'].default_value = .68
        coords = n.new('ShaderNodeTexCoord')
        links.new(coords.outputs['Object'], tex.inputs['Vector'])
        bump = n.new('ShaderNodeBump')
        bump.inputs['Strength'].default_value = .24
        bump.inputs['Distance'].default_value = .023
        links.new(tex.outputs['Fac'], bump.inputs['Height'])
        links.new(bump.outputs['Normal'], p.inputs['Normal'])
        ramp = n.new('ShaderNodeValToRGB')
        ramp.color_ramp.elements[0].position = .18
        ramp.color_ramp.elements[0].color = (*(v*.68 for v in color), 1)
        ramp.color_ramp.elements[1].position = .83
        ramp.color_ramp.elements[1].color = (*(min(1, v*1.15+.025) for v in color), 1)
        links.new(tex.outputs['Fac'], ramp.inputs[0])
        links.new(ramp.outputs['Color'], p.inputs['Base Color'])
    return m


class Parts:
    """Accumulate small geometric components into one editable material mesh."""
    def __init__(self):
        self.v, self.f, self.mi = [], [], []
        self.ball_count = 0

    def ball(self, pos, scale, mat=0, seg=10, rings=6):
        self.ball_count += 1
        if isinstance(scale, (float, int)):
            scale = (scale,)*3
        start = len(self.v)
        self.v.append((pos[0], pos[1], pos[2]+scale[2]))
        for j in range(1, rings):
            a = math.pi*j/rings
            for i in range(seg):
                t = TAU*i/seg
                self.v.append((pos[0]+scale[0]*math.sin(a)*math.cos(t), pos[1]+scale[1]*math.sin(a)*math.sin(t), pos[2]+scale[2]*math.cos(a)))
        bottom = len(self.v)
        self.v.append((pos[0], pos[1], pos[2]-scale[2]))
        for i in range(seg):
            self.f.append((start, start+1+i, start+1+(i+1)%seg))
            self.mi.append(mat)
        for j in range(rings-2):
            for i in range(seg):
                k = start+1+j*seg+i
                k2 = start+1+j*seg+(i+1)%seg
                self.f.append((k, k+seg, k2+seg, k2))
                self.mi.append(mat)
        for i in range(seg):
            self.f.append((bottom, bottom-seg+(i+1)%seg, bottom-seg+i))
            self.mi.append(mat)

    def rod(self, a, b, radius, mat=0, seg=8):
        a, b = Vector(a), Vector(b)
        axis = (b-a).normalized()
        u = axis.cross(Vector((0, 0, 1)))
        if u.length < .01:
            u = axis.cross(Vector((0, 1, 0)))
        u.normalize()
        v = axis.cross(u)
        start = len(self.v)
        for p in (a, b):
            for i in range(seg):
                q = p+radius*(u*math.cos(TAU*i/seg)+v*math.sin(TAU*i/seg))
                self.v.append(tuple(q))
        for i in range(seg):
            self.f.append((start+i, start+(i+1)%seg, start+seg+(i+1)%seg, start+seg+i))
            self.mi.append(mat)
        self.f.extend([tuple(start+i for i in reversed(range(seg))), tuple(start+seg+i for i in range(seg))])
        self.mi.extend([mat, mat])

    def object(self, name, materials, collection):
        obj = mesh_object(name, self.v, self.f, materials, collection, self.mi)
        obj['sphere_component_count'] = self.ball_count
        return obj


def curve(name, points, radius, mat, collection, cyclic=False):
    d = bpy.data.curves.new(name, 'CURVE')
    d.dimensions = '3D'
    d.bevel_depth = radius
    d.bevel_resolution = 2
    s = d.splines.new('POLY')
    s.points.add(len(points)-1)
    for p, co in zip(s.points, points):
        p.co = (*co, 1)
    s.use_cyclic_u = cyclic
    d.materials.append(mat)
    o = bpy.data.objects.new(name, d)
    collection.objects.link(o)
    return o



# The top cut edge is one CLOSED, deeply indented bean-shaped contour.
# Every U-shaped infolding returns to the same peripheral boundary.
RING_FOLDS = [(-3.25, 2.00, 1), (-.70, 2.10, 1), (1.85, 2.00, 1),
              (-2.10, 1.60, -1), (.50, 1.90, -1), (3.10, 1.30, -1)]
FOLDS = [(x,d) for x,d,side in RING_FOLDS]
RIM_HEIGHT = 0.


def ring_point(t, q, outer=False):
    s = math.sin(q)
    length, width = (5.35,1.96) if outer else (5.08,1.68)
    # Elliptical ends, with a gentle asymmetric bean-like axial bend.
    x = length*math.cos(t)
    y = width*math.sin(t)
    x_base = x
    if not outer:
        side = 1 if math.sin(t)>=0 else -1
        for k,(center,depth,root_side) in enumerate(RING_FOLDS):
            if root_side!=side:
                continue
            sigma = .113+.010*math.sin(k*1.7)
            delta=(x_base-center)/sigma
            profile=math.exp(-.5*delta*delta)
            y -= side*depth*profile*s**4
            x += sigma*1.8*delta*profile*s**4
            x += .07*math.sin(k*1.8)*profile*s**4
    y += .10*math.sin(x_base*.72)
    return Vector((x*s,y*s,RIM_HEIGHT-(1.98 if outer else 1.72)*math.cos(q)))


def ring_normal(t,q,outer=False):
    dt=ring_point(t+.0005,q,outer)-ring_point(t-.0005,q,outer)
    dq=ring_point(t,q+.0005,outer)-ring_point(t,q-.0005,outer)
    n=dq.cross(dt).normalized()
    return -n if outer else n


def fold_t(index):
    x,depth,side=RING_FOLDS[index]
    angle=math.acos(x/5.08)
    return angle if side>0 else TAU-angle


def ring_parameters(count,outer=False):
    # Uniform arc-length at the cut edge resolves both long walls and tight
    # hairpin tips without adding excessive samples to broad empty periphery.
    import bisect
    dense=16000
    lengths=[0.]
    previous=ring_point(0,math.pi/2,outer)
    for i in range(1,dense+1):
        point=ring_point(TAU*i/dense,math.pi/2,outer)
        lengths.append(lengths[-1]+(point-previous).length)
        previous=point
    parameters=[]
    for i in range(count):
        target=lengths[-1]*i/count
        j=max(1,bisect.bisect_left(lengths,target))
        fraction=(target-lengths[j-1])/(lengths[j]-lengths[j-1])
        parameters.append(TAU*(j-1+fraction)/dense)
    return parameters



def retained(point):
    return not (point[1]<-1e-6 and point[2]>1e-6)


def membranes(c,mats):
    for outer in (True,False):
        nt,nq=(1000,104) if outer else (1900,136)
        parameters=ring_parameters(nt,outer)
        verts=[tuple(ring_point(0,0,outer))]
        for j in range(1,nq):
            q=math.pi*j/nq
            for t in parameters:
                verts.append(tuple(ring_point(t,q,outer)))
        top=len(verts)
        verts.append(tuple(ring_point(0,math.pi,outer)))
        faces=[]
        for i in range(nt):
            faces.append((0,1+(i+1)%nt,1+i))
        for j in range(nq-2):
            for i in range(nt):
                a=1+j*nt+i
                b=1+j*nt+(i+1)%nt
                faces.append((a,b,b+nt,a+nt))
        for i in range(nt):
            faces.append((top,1+(nq-2)*nt+i,1+(nq-2)*nt+(i+1)%nt))
        name='01 OUTER MEMBRANE — intact three-quarter ellipsoid shell' if outer else '02 INNER MEMBRANE — continuous serpentine invaginations, upper-left wedge removed'
        ob=mesh_object(name,verts,faces,[mats['outer'] if outer else mats['inner'],mats['rim']],c)
        bm=bmesh.new()
        bm.from_mesh(ob.data)
        # Split both planes first; remove only their upper-left intersection.
        # No cap face is generated: the opening remains an unobstructed wedge.
        for normal in ((0,0,1),(0,1,0)):
            bmesh.ops.bisect_plane(bm,geom=list(bm.verts)+list(bm.edges)+list(bm.faces),dist=1e-7,plane_co=(0,0,0),plane_no=normal,clear_inner=False,clear_outer=False)
        doomed=[face for face in bm.faces if not retained(face.calc_center_median())]
        bmesh.ops.delete(bm,geom=doomed,context='FACES')
        loose=[v for v in bm.verts if not v.link_faces]
        if loose:
            bmesh.ops.delete(bm,geom=loose,context='VERTS')
        bmesh.ops.recalc_face_normals(bm,faces=list(bm.faces))
        bm.normal_update()
        # A boundary walk follows the real cut mesh through both cutting planes.
        boundary=[e for e in bm.edges if e.is_boundary]
        adjacency={}
        for e in boundary:
            a,b=e.verts
            adjacency.setdefault(a,[]).append(b)
            adjacency.setdefault(b,[]).append(a)
        first=next(iter(adjacency))
        loop=[first]
        previous=None
        current=first
        for _ in range(len(boundary)+1):
            candidates=[v for v in adjacency[current] if v is not previous]
            nxt=candidates[0]
            if nxt is first:
                break
            loop.append(nxt)
            previous,current=current,nxt
        points=[v.co.copy() for v in loop]
        normals=[v.normal.copy() for v in loop]
        bm.to_mesh(ob.data)
        bm.free()
        ob.data.update()
        for face in ob.data.polygons:
            face.use_smooth=True
        solid=ob.modifiers.new('Membrane thickness only — open wedge without caps','SOLIDIFY')
        solid.thickness=.064 if outer else .025
        solid.offset=0
        solid.material_offset_rim=1
        ob['Structure']='Full closed ellipsoid deformed into continuous U-shaped inner-membrane folds before removing only y<0 AND z>0. Remaining three-quarter shell has no cut-face cap.'
        ob['Cutaway definition']='Long axis X. View from +X toward origin (view direction -X): screen-left is -Y, up is +Z. Removed quadrant: y<0 and z>0.'
        import bisect
        lengths=[0.]
        for i in range(len(points)):
            lengths.append(lengths[-1]+(points[(i+1)%len(points)]-points[i]).length)
        edges=Parts()
        seam=[]
        for i in range(1600):
            target=lengths[-1]*i/1600
            j=max(0,bisect.bisect_right(lengths,target)-1)
            f=(target-lengths[j])/(lengths[j+1]-lengths[j])
            p=points[j].lerp(points[(j+1)%len(points)],f)
            n=normals[j].lerp(normals[(j+1)%len(points)],f).normalized()
            seam.append(tuple(p))
            for sign in (-1,1):
                edges.ball(p+n*sign*(.036 if outer else .017),.021 if outer else .017,i%3,seg=6,rings=4)
        curve(('OM' if outer else 'IM')+' wedge-cut bilayer seam — no cap',seam,.014,mats['seam'],c,True)
        edges.object(('OM' if outer else 'IM')+' wedge cut-edge two-leaflet heads',[mats['outer_light'] if outer else mats['lipid'],mats['rim'],mats['lipid']],c)


def membrane_detail(c,mats):
    p=Parts()
    for i in range(3600):
        t=RNG.random()*TAU
        q=RNG.uniform(.35,math.pi-.35)
        pos=ring_point(t,q,True)+ring_normal(t,q,True)*.038
        if not retained(pos):
            continue
        p.ball(pos,RNG.uniform(.012,.022),RNG.randrange(3),seg=6,rings=4)
    p.object('OM retained surface lipid heads',[mats['outer_light'],mats['outer'],mats['rim']],c)
    p=Parts()
    for k in range(6):
        tc=fold_t(k)
        for i in range(1000):
            t=tc+RNG.uniform(-.055,.055)
            q=RNG.uniform(.65,math.pi-.65)
            pos=ring_point(t,q)+ring_normal(t,q)*.022
            if not retained(pos):
                continue
            p.ball(pos,RNG.uniform(.009,.016),RNG.randrange(3),seg=6,rings=4)
    p.object('IM retained infolding lipid heads',[mats['lipid'],mats['inner'],mats['rim']],c)


def synthase_template(c, mats):
    p = Parts()
    # Local +Z points to matrix. Two catalytic heads sit on tilted stalks,
    # above a membrane rotor and a peripheral stator on each monomer.
    for side in (-1, 1):
        bx = side*.094
        hx = side*.123
        for j in range(10):
            a = TAU*j/10
            p.rod((bx+.039*math.cos(a), .039*math.sin(a), -.041), (bx+.039*math.cos(a), .039*math.sin(a), .042), .012, 2, 7)
        p.rod((bx, 0, .02), (hx, 0, .148), .016, 1)
        for j in range(6):
            a = TAU*j/6
            pos = (hx+.049*math.cos(a), .049*math.sin(a), .183)
            p.ball(pos, (.037, .037, .054), j%2, seg=10, rings=6)
            for z in (-.025, .016):
                p.ball((pos[0]+.018*math.cos(a), pos[1]+.018*math.sin(a), pos[2]+z), .018, j%2, seg=7, rings=4)
        p.ball((hx, 0, .231), (.025,.025,.018), 1, seg=9, rings=5)
        p.rod((bx+side*.054, 0, .003), (hx+side*.075, 0, .158), .010, 1)
        p.rod((hx+side*.075, 0, .158), (hx+side*.029, 0, .208), .009, 1)
    p.rod((-.075,0,-.01), (.075,0,-.01), .017, 2)
    ob = p.object('ATP synthase DIMER — linked geometric prototype', [mats['protein'], mats['gold'], mats['protein_dark']], c)
    ob.hide_render = True
    ob.hide_viewport = True
    return ob


def linked(template, name, pos, normal, collection, scale=1., twist=0.):
    if not retained(pos):
        return None
    o = bpy.data.objects.new(name, template.data)
    collection.objects.link(o)
    o.location = pos
    o.rotation_mode = 'QUATERNION'
    o.rotation_quaternion = Vector(normal).to_track_quat('Z', 'Y')
    if twist:
        from mathutils import Quaternion
        o.rotation_quaternion = o.rotation_quaternion @ Quaternion((0,0,1), twist)
    o.scale = (scale,)*3
    return o



def proteins(c,prototypes,mats):
    atp=synthase_template(prototypes,mats)
    total=0
    for k in range(6):
        t=fold_t(k)
        for j in range(26):
            q=.64+((math.pi-1.28)*j/25)
            n=ring_normal(t,q)
            if not retained(ring_point(t,q)+n*.026):
                continue
            linked(atp,'ATP synthase dimer | U-fold %02d bend %02d'%(k+1,j+1),ring_point(t,q)+n*.026,n,c,.73+RNG.random()*.12)
            total+=1
    p=Parts()
    for j in range(17):
        a=TAU*j/17
        p.ball((.091*math.cos(a),.045*math.sin(a),.012+.015*math.sin(3*a)),(.026,.030,.065),j%3,seg=8,rings=5)
    for pos,r in [((-.07,0,.097),(.057,.051,.065)),((-.04,0,.157),(.037,.033,.055)),((.043,0,.087),(.052,.047,.050))]:
        p.ball(pos,r,1,seg=10,rings=6)
    resp=p.object('Respiratory supercomplex — illustrative prototype',[mats['protein_dark'],mats['rust'],mats['gold']],prototypes)
    resp.hide_render=resp.hide_viewport=True
    for k in range(6):
        for j in range(7):
            t=fold_t(k)+RNG.choice((-1,1))*RNG.uniform(.024,.055)
            q=RNG.uniform(.82,math.pi-.82)
            n=ring_normal(t,q)
            linked(resp,'Respiratory complex %02d.%02d'%(k,j),ring_point(t,q)+n*.024,n,c,.7+RNG.random()*.2,RNG.random()*TAU)
    p=Parts()
    for j in range(12):
        a=TAU*j/12
        p.rod((.076*math.cos(a),.076*math.sin(a),-.05),(.076*math.cos(a+.15),.076*math.sin(a+.15),.066),.019,j%2,7)
        p.ball((.076*math.cos(a+.15),.076*math.sin(a+.15),.066),.021,j%2,seg=7,rings=4)
    pore=p.object('Outer-membrane barrel channel — prototype',[mats['protein'],mats['protein_dark']],prototypes)
    pore.hide_render=pore.hide_viewport=True
    for j in range(31):
        t,q=RNG.random()*TAU,RNG.uniform(.8,math.pi-.8)
        n=ring_normal(t,q,True)
        linked(pore,'Outer membrane channel %02d'%j,ring_point(t,q,True)+n*.04,n,c,RNG.uniform(.67,1.))
    return total


def matrix_contents(c, prototypes, mats):
    # Three compact, intertwined circular double-stranded DNA loops.
    centers = [(-2.67,-.10,-.15),(.00,.08,-.18),(2.65,-.10,-.12)]
    for k, center in enumerate(centers):
        for strand in range(2):
            pts=[]
            for j in range(480):
                t=TAU*j/480
                helix=.026*math.cos(23*t+strand*math.pi)
                pts.append((center[0]+(.34+helix)*math.cos(t),center[1]+(.23+helix)*math.sin(t),center[2]+.13*math.sin(3*t)+.026*math.sin(23*t+strand*math.pi)))
            curve('mtDNA nucleoid %d — circular strand %d'%(k+1,strand+1),pts,.012,mats['dna'] if strand==0 else mats['gold'],c,True)
    p=Parts()
    for j in range(12):
        a=TAU*j/12
        p.ball((.057*math.cos(a),.057*math.sin(a),.025*math.sin(3*a)),(.048,.042,.043),j%3,seg=8,rings=5)
    enzyme=p.object('Matrix enzyme oligomer — geometric prototype',[mats['protein'],mats['matrix_protein'],mats['gold']],prototypes)
    enzyme.hide_render=enzyme.hide_viewport=True
    for j in range(78):
        x=RNG.uniform(-4.48,4.48)
        y=RNG.uniform(-.92,.49)
        z=RNG.uniform(-.52,.33)
        # Keep bodies in matrix troughs, away from thin crista lumina.
        if min(abs(x-base) for base,_ in FOLDS)<.24:
            x+=.29
        linked(enzyme,'Matrix enzyme assembly %03d'%j,(x,y,z),(RNG.random(),RNG.random(),1),c,RNG.uniform(.53,1.08),RNG.random()*TAU)
    # Illustrative pyruvate anion connectivity, C3H3O3−; not optimized geometry.
    p=Parts()
    atomdata=[((-.15,0,0),.036,0),((0,0,0),.036,0),((.15,0,0),.036,0),((0,.112,0),.034,1),((.215,.093,0),.034,1),((.215,-.093,0),.034,1),((-.202,.071,.047),.022,2),((-.202,-.071,.047),.022,2),((-.202,0,-.082),.022,2)]
    for pos,rad,idx in atomdata:
        p.ball(pos,rad,idx,seg=10,rings=6)
    for a,b in [(0,1),(1,2),(2,4),(2,5),(0,6),(0,7),(0,8)]:
        p.rod(atomdata[a][0],atomdata[b][0],.010,3,8)
    for dx in (-.012,.012):
        p.rod((dx,.022,0),(dx,.09,0),.006,3,7)
    pyruvate=p.object('Pyruvate anion C3H3O3− — geometric prototype',[mats['carbon'],mats['oxygen'],mats['protein'],mats['bond']],prototypes)
    pyruvate.hide_render=pyruvate.hide_viewport=True
    for j,pos in enumerate([(-2.68,-.78,-.28),(.03,-.86,-.24),(1.15,-.15,-.18),(3.94,-.49,-.25)]):
        linked(pyruvate,'Pyruvate anion %d'%(j+1),pos,(.3,.2,1),c,1.02 if j<2 else .80, j*1.1)


