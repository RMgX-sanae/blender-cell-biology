"""Run in background Blender to check scene-local compositor branches after reopening."""
from pathlib import Path
import sys
import tempfile

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import compositors
from render_saved import saved_compositor


def main():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    first = bpy.context.scene
    first.name = 'Overview'
    second = bpy.data.scenes.new('Detail')
    bpy.context.window.scene = second
    profile = {'background': '102030', 'compositor': {'preset': 'reference'}}
    compositors(first, profile)
    compositors(second, profile)
    second['detail_compositor_branches']['clean'].name = 'Renamed detail grade'
    legacy = bpy.data.node_groups.new('COMP_STYLE_ink', 'CompositorNodeTree')
    legacy.use_fake_user = True
    with tempfile.TemporaryDirectory() as tmp:
        path = str(Path(tmp) / 'camera-studies.blend')
        bpy.ops.wm.save_as_mainfile(filepath=path)
        bpy.ops.wm.open_mainfile(filepath=path)
        for scene in bpy.data.scenes:
            for mode in ('clean', 'reference', 'monochrome'):
                group = saved_compositor(scene, mode)
                assert group == scene['detail_compositor_branches'][mode]
                layers = [node for node in group.nodes if node.type == 'R_LAYERS']
                assert len(layers) == 1 and layers[0].scene == scene, (scene.name, mode)
            assert saved_compositor(scene, 'ink') == bpy.data.node_groups['COMP_STYLE_ink']
            assert saved_compositor(scene, 'absent') is None
        overview = bpy.data.scenes['Overview']
        del overview['detail_compositor_branches']
        assert saved_compositor(overview, 'clean') == bpy.data.node_groups['COMP_clean']
    print('Two scenes, six reopened render inputs, renamed grade and legacy lookup verified.')


if __name__ == '__main__':
    main()
