"""Blender integration checks for saved role-shader handoff and rejection."""
import copy
import sys
from pathlib import Path
import unittest

import bpy

sys.path.insert(0,str(Path(__file__).resolve().parent))
from detail_scene import material
from audit_surface_materials import check_materials
from apply_scale_environment import layer_compositor
from audit_interaction_scene import displayed_materials


class SurfaceHandoff(unittest.TestCase):
    def setUp(self):
        for mat in list(bpy.data.materials):bpy.data.materials.remove(mat)
        self.profile={'surface_roles':{'protein':{'roughness':.31,'bump':.002,'texture_scale':7,'subsurface_weight':.07}}}
        self.mat=material('Checked protein','DEA835',settings=self.profile['surface_roles']['protein'])
        self.mat['surface_role']='protein'
        self.shader=self.mat.node_tree.nodes.get('Principled BSDF')

    def check(self):return check_materials(self.profile,['protein'])

    def test_applied_and_transparent_surface(self):
        self.assertTrue(self.check()['passed'])
        tree=self.mat.node_tree
        transparent=tree.nodes.new('ShaderNodeBsdfTransparent')
        mix=tree.nodes.new('ShaderNodeMixShader');mix.inputs[0].default_value=.35
        tree.links.new(transparent.outputs[0],mix.inputs[1]);tree.links.new(self.shader.outputs[0],mix.inputs[2])
        tree.links.new(mix.outputs[0],tree.nodes.get('Material Output').inputs['Surface'])
        self.assertTrue(self.check()['passed'])

    def test_changed_node_rejected(self):
        self.shader.inputs['Roughness'].default_value=.8
        with self.assertRaisesRegex(ValueError,'Saved shader differs'):self.check()

    def test_unused_shader_rejected(self):
        tree=self.mat.node_tree;other=tree.nodes.new('ShaderNodeEmission')
        tree.links.new(other.outputs[0],tree.nodes.get('Material Output').inputs['Surface'])
        with self.assertRaisesRegex(ValueError,'active surface output'):self.check()

    def test_missing_normal_connection_rejected(self):
        self.mat.node_tree.links.remove(self.shader.inputs['Normal'].links[0])
        with self.assertRaisesRegex(ValueError,'Bump does not reach'):self.check()

    def test_displacement_rejected(self):
        tree=self.mat.node_tree;noise=next(n for n in tree.nodes if n.type=='TEX_NOISE')
        tree.links.new(noise.outputs['Fac'],tree.nodes.get('Material Output').inputs['Displacement'])
        with self.assertRaisesRegex(ValueError,'geometry displacement'):self.check()

    def test_profile_drift_and_missing_role_rejected(self):
        profile=copy.deepcopy(self.profile);profile['surface_roles']['protein']['roughness']=.8
        with self.assertRaisesRegex(ValueError,'Profile setting was not applied'):check_materials(profile,['protein'])
        del self.mat['surface_role']
        with self.assertRaisesRegex(ValueError,'Required role materials missing'):self.check()

    def test_texture_space_and_driven_frequency_rejected(self):
        tree=self.mat.node_tree;noise=next(n for n in tree.nodes if n.type=='TEX_NOISE')
        coords=next(n for n in tree.nodes if n.type=='TEX_COORD')
        tree.links.new(coords.outputs['Generated'],noise.inputs['Vector'])
        with self.assertRaisesRegex(ValueError,'Object texture coordinates'):self.check()
        tree.links.new(coords.outputs['Object'],noise.inputs['Vector'])
        value=tree.nodes.new('ShaderNodeValue');tree.links.new(value.outputs[0],noise.inputs['Scale'])
        with self.assertRaisesRegex(ValueError,'texture setting override'):self.check()

    def test_invalid_environment_compositor_rejected(self):
        for settings in [{'opacity':[.2,float('nan'),.2]},{'opacity':[.2,1.1,.2]},
                         {'opacity':[.2,.3]},{'blur_pixels_1280':[-1,10,20]},
                         {'blur_pixels_1280':[1,float('inf'),20]},{'blur_pixel_1280':[1,2,3]}]:
            with self.subTest(settings=settings),self.assertRaises(ValueError):
                layer_compositor(bpy.context.scene,[],'dark',1280,'cell',{'environment_compositor':{'dark':settings}})

    def test_actual_surface_and_atom_material_assignment(self):
        from detail_scene import mesh
        from style_nodes import point_object,instance_spheres
        other=material('Different protein','70C8E0')
        group=bpy.context.scene.collection
        surface=mesh('Assignment surface',[(0,0,0),(1,0,0),(0,1,0)],[(0,1,2)],[self.mat],group)
        atoms=point_object('Assignment atoms',[(0,0,0)],group)
        try:
            self.assertEqual(displayed_materials(surface),{self.mat})
            surface.material_slots[0].link='OBJECT';surface.material_slots[0].material=other
            self.assertEqual(displayed_materials(surface),{other})
            tree=instance_spheres(atoms,self.mat)
            setter=next(n for n in tree.nodes if n.type=='SET_MATERIAL')
            unused=tree.nodes.new('GeometryNodeSetMaterial');unused.inputs['Material'].default_value=other
            self.assertEqual(displayed_materials(atoms),{self.mat})
            setter.inputs['Material'].default_value=other
            self.assertEqual(displayed_materials(atoms),{other})
            for link in list(setter.outputs['Geometry'].links):tree.links.remove(link)
            self.assertNotEqual(displayed_materials(atoms),{self.mat})
        finally:
            bpy.data.objects.remove(surface,do_unlink=True);bpy.data.objects.remove(atoms,do_unlink=True)

    def test_molecular_depth_controls_and_rejections(self):
        from cellular_environment import composite
        layers=[bpy.context.scene.view_layers.new('Depth test '+name) for name in ('far','middle','near')]
        tree=bpy.data.node_groups.new('Depth control test','CompositorNodeTree')
        try:
            composite(tree,layers,'081821',640,settings={'blur_pixels_1280':[12,6,20],'opacity':[.3,.2,.1]})
            for layer,pixels,opacity in zip(layers,[6,3,10],[.3,.2,.1]):
                blur=next(n for n in tree.nodes if n.label==layer.name+' defocus')
                over=next(n for n in tree.nodes if n.label==layer.name+' background strength')
                self.assertEqual(tuple(blur.inputs['Size'].default_value),(pixels,pixels))
                self.assertAlmostEqual(over.inputs['Factor'].default_value,opacity)
                self.assertEqual(over.inputs['Foreground'].links[0].from_node,blur)
            for settings in [{'opacity':[.2,float('nan'),.2]},{'opacity':[.2,1.1,.2]},
                             {'blur_pixels_1280':[1,2]},{'blur_pixels_1280':[-1,2,3]},
                             {'blur_pixel_1280':[1,2,3]}]:
                with self.subTest(settings=settings),self.assertRaises(ValueError):
                    composite(tree,layers,'081821',640,settings=settings)
            with self.assertRaises(ValueError):composite(tree,layers[:2],'081821',640,settings={'opacity':[.2,.2,.2]})
        finally:
            for layer in layers:bpy.context.scene.view_layers.remove(layer)
            bpy.data.node_groups.remove(tree)

    def test_studio_look_follows_selected_view_transform(self):
        from types import SimpleNamespace
        from detail_scene import studio
        args=SimpleNamespace(engine='eevee',samples=1,width=64,fps=24,frames=1,compositor='clean')
        profile={'world':'839BAA','background':'081821','lighting':[],'compositor':{}}
        camera={'position':[0,0,10],'target':[0,0,0],'scale':10}
        for settings,transform,look in [({},'AgX','AgX - Medium High Contrast'),
                ({'view_transform':'Standard'},'Standard','None'),
                ({'view_transform':'Standard','look':'Medium High Contrast'},'Standard','Medium High Contrast')]:
            with self.subTest(settings=settings):
                studio({**profile,**settings},args,camera)
                self.assertEqual(bpy.context.scene.view_settings.view_transform,transform)
                self.assertEqual(bpy.context.scene.view_settings.look,look)
        with self.assertRaises(TypeError):
            studio({**profile,'view_transform':'Standard','look':'AgX - Medium High Contrast'},args,camera)


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SurfaceHandoff))
    if not result.wasSuccessful():raise SystemExit(1)
