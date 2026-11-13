"""Check the moving3-prime guide against independently read reaction atom sites."""
import bpy
import numpy as np
from mathutils import Vector


def check_growing_end(scene, primer_o3, c20_o3, u21_o3=None):
    if 'transcription_appearance' not in scene:
        return None
    original = bpy.context.window.view_layer; frame = scene.frame_current
    beauty = scene.view_layers[scene['compositor_beauty_layer']]
    text = scene.view_layers[scene['compositor_label_layer']]
    guide = scene.objects['Leader RNA growing 3′']; label = scene.objects['RNA growing 3′ end']
    core = scene.objects['Incoming CTP core → RNA C20']
    front = np.array(scene.camera.rotation_euler.to_matrix())[:,2]
    right = np.array(scene.camera.rotation_euler.to_matrix())[:,0]
    maximum_anchor = maximum_label = 0.
    try:
        for f in range(scene.frame_start, scene.frame_end+1):
            bpy.context.window.view_layer = beauty
            scene.frame_set(f); beauty.update(); graph = bpy.context.evaluated_depsgraph_get()
            time = (f-1)/scene.render.fps
            if u21_o3 is None or time < 11:
                join = core.data.shape_keys.evaluated_get(graph).key_blocks['Observed 5FLM C20 product'].value
                expected = (1-join)*primer_o3+join*c20_o3
            else:
                parent = scene.objects['Authored one-position scaffold translocation'].evaluated_get(graph)
                expected = np.array(parent.matrix_world@Vector(c20_o3))
                if time >= 17.5:
                    bond = scene.objects['Second RNA O3′–P bond'].evaluated_get(graph)
                    join = bond.data.bevel_depth/.035
                    expected = (1-join)*expected+join*u21_o3
            evaluated = guide.evaluated_get(graph)
            points = [np.array(evaluated.matrix_world@Vector(p.co[:3])) for p in evaluated.data.splines[0].points]
            anchor_error = float(np.max(np.abs(points[0]-(expected+front*20))))
            bpy.context.window.view_layer = text
            scene.frame_set(f); text.update(); graph = bpy.context.evaluated_depsgraph_get()
            label_position = np.array(label.evaluated_get(graph).matrix_world.translation)
            label_error = float(np.max(np.abs(points[1]-(label_position+right*.1))))
            maximum_anchor = max(maximum_anchor, anchor_error); maximum_label = max(maximum_label, label_error)
            if max(anchor_error, label_error) > 4e-6:
                raise ValueError('Growing3-prime guide differs from RNA endpoint or label at frame '+str(f))
    finally:
        bpy.context.window.view_layer = original; scene.frame_set(frame); bpy.context.view_layer.update()
    return {'passed':True, 'frames':scene.frame_end, 'maximum_anchor_error_nm':maximum_anchor,
            'maximum_label_error_nm':maximum_label,
            'scope':'Displayed growing-end guide follows source O3-prime sites and the existing authored joining/translocation; not an additional catalytic trajectory.'}
