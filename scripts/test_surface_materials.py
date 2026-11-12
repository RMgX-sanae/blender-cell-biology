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


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SurfaceHandoff))
    if not result.wasSuccessful():raise SystemExit(1)
