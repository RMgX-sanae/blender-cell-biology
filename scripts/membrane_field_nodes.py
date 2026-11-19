"""Native Geometry Nodes for a membrane with a changing connection between two aqueous regions."""
import bpy


def membrane_field(obj, record, bounds, materials, spacing=.0625):
    tree = bpy.data.node_groups.new('Membrane from signed distance', 'GeometryNodeTree')
    tree.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
    controls = {}
    for name, kind in [('Center','Vector'), ('Radius','Float'), ('Neck radius','Float'),
                       ('Neck open','Bool'), ('Pore radius','Float'), ('Pore open','Bool')]:
        socket = tree.interface.new_socket(name=name, in_out='INPUT', socket_type='NodeSocket'+kind)
        socket.force_non_field = True
        controls[name] = socket.identifier

    def node(kind, label, column, row):
        item = tree.nodes.new(kind)
        item.label = label
        item.location = (column*215, -row*215)
        return item

    def plug(target, value):
        if isinstance(value, bpy.types.NodeSocket):
            tree.links.new(value, target)
        else:
            target.default_value = value

    def scalar(op, a, b, label, column, row):
        n = node('ShaderNodeMath', label, column, row)
        n.operation = op
        plug(n.inputs[0], a)
        if b is not None:
            plug(n.inputs[1], b)
        return n.outputs[0]

    group = node('NodeGroupInput','Tracked field controls',0,2)
    position = node('GeometryNodeInputPosition','Sample point',0,0).outputs[0]
    center = group.outputs['Center']

    def sphere(location, radius, label, row):
        distance = node('ShaderNodeVectorMath', label+' center',1,row)
        distance.operation = 'DISTANCE'
        plug(distance.inputs[0], position)
        plug(distance.inputs[1], location)
        return scalar('SUBTRACT',distance.outputs['Value'],radius,label+' surface',2,row)

    def bridge(base, end, radius, opened, label, column, row):
        # Distance to a finite edge equals distance to the axis of the capped tube.
        line = node('GeometryNodeMeshLine', label+' axis',1,row)
        line.mode = 'END_POINTS'
        line.inputs['Count'].default_value = 2
        plug(line.inputs['Start Location'], center)
        plug(line.inputs['Offset'], end)
        nearest = node('GeometryNodeProximity', label+' distance',2,row)
        nearest.target_element = 'EDGES'
        plug(nearest.inputs['Geometry'], line.outputs['Mesh'])
        plug(nearest.inputs['Sample Position'], position)
        tube = scalar('SUBTRACT',nearest.outputs['Distance'],radius,label+' radius',3,row)
        joined = scalar('MINIMUM',base,tube,label+' union',column,row)
        switch = node('GeometryNodeSwitch',label+' present',column+1,row)
        switch.input_type = 'FLOAT'
        plug(switch.inputs['Switch'], opened)
        plug(switch.inputs['False'], base)
        plug(switch.inputs['True'], joined)
        return switch.outputs[0]

    axes = node('ShaderNodeSeparateXYZ','Donor plane',1,0)
    plug(axes.inputs[0],position)
    plane = scalar('MULTIPLY',axes.outputs['Z'],-1.,'Extracellular side',2,0)
    carrier = sphere(center,group.outputs['Radius'],'Carrier',1)
    target = sphere(record['target_center'],record['target_radius'],'Target compartment',2)
    distance = scalar('MINIMUM',plane,carrier,'Donor and carrier',3,0)
    distance = scalar('MINIMUM',distance,target,'Aqueous compartments',4,0)
    distance = bridge(distance,record['neck_end'],group.outputs['Neck radius'],group.outputs['Neck open'],'Donor neck',5,3)
    distance = bridge(distance,record['target_center'],group.outputs['Pore radius'],group.outputs['Pore open'],'Fusion pore',7,5)

    # A density band makes both physical boundaries of the same shell. The crop
    # closes only the thin rim of the sheet, rather than the entire water volume.
    half = record['half_thickness']
    absolute = scalar('ABSOLUTE',distance,None,'Distance from membrane middle',9,2)
    band = scalar('SUBTRACT',2*half,absolute,'Membrane band',10,2)
    density = scalar('MAXIMUM',band,0.,'Empty background',11,2)
    volume = node('GeometryNodeVolumeCube','Local membrane field',12,2)
    plug(volume.inputs['Density'],density)
    lo, hi = bounds
    volume.inputs['Min'].default_value = lo
    volume.inputs['Max'].default_value = hi
    for a,b,axis in zip(lo,hi,'XYZ'):
        volume.inputs['Resolution '+axis].default_value = round((b-a)/spacing)+1
    mesh = node('GeometryNodeVolumeToMesh','Both sides of membrane',13,2)
    mesh.inputs['Threshold'].default_value = half
    plug(mesh.inputs['Volume'],volume.outputs['Volume'])
    smooth = node('GeometryNodeSetShadeSmooth','Smooth membrane',14,2)
    plug(smooth.inputs['Geometry'],mesh.outputs['Mesh'])
    exterior = node('GeometryNodeSetMaterial','Cytosolic boundary',15,2)
    plug(exterior.inputs['Geometry'],smooth.outputs['Geometry'])
    plug(exterior.inputs['Material'],materials[0])
    interior = node('GeometryNodeSetMaterial','Lumen boundary',16,1)
    plug(interior.inputs['Geometry'],exterior.outputs['Geometry'])
    plug(interior.inputs['Selection'],scalar('LESS_THAN',distance,0.,'Lumen side',14,0))
    plug(interior.inputs['Material'],materials[1])
    output = node('NodeGroupOutput','Membrane',17,1)
    plug(output.inputs['Geometry'],interior.outputs['Geometry'])
    modifier = obj.modifiers.new('Membrane field and connections','NODES')
    modifier.node_group = tree
    return modifier, controls
