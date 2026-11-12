"""Check profile handoff, stale-scene rejection and immutable inputs."""
import json,tempfile,unittest
from pathlib import Path
from appearance_config import configure,load_saved,read_profile,chapter_profile


class AppearanceInput(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.path=self.root/'style.json'
        self.path.write_text('{}')

    def request(self,profile):
        self.path.write_text(json.dumps(profile));return read_profile(self.path)

    def test_supported_camera_and_chapter_override(self):
        profile=self.request({'camera':{'scale_multiplier':1.2},'chapters':{'a03':{'camera':{'center_offset_nm':[1,-2]}}}})
        self.assertEqual(chapter_profile(profile,'a03')['camera'],{'scale_multiplier':1.2,'center_offset_nm':[1,-2]})
        self.assertEqual(chapter_profile(profile,'f04')['camera'],{'scale_multiplier':1.2})

    def test_invalid_controls(self):
        for profile in [{'camera':{'type':'PERSP'}},{'camera':{'scale_multiplier':0}},
                        {'camera':{'center_offset_nm':[0,float('nan')]}},{'palette':{'unknown':'FFFFFF'}},
                        {'palette':{'teal':'not-a-color'}},{'surface_roles':{'unknown':{}}},
                        {'chapters':{'a99':{}}},{'chapters':{'a03':{'lighting':[]}}},
                        {'lighting':[]},{'exposure':float('inf')},{'environment':{'compartment':'nucleus'}}]:
            with self.subTest(profile=profile),self.assertRaises(ValueError):self.request(profile)

    def test_snapshot_and_changed_profile_need_rebuild(self):
        options={};configure(options,{},self.path,None,self.root,['build'])
        self.assertEqual(load_saved(options),{})
        self.request({'exposure':.5})
        with self.assertRaisesRegex(ValueError,'requires rebuilding'):configure({},options,self.path,None,self.root,['render'])
        updated={};configure(updated,options,self.path,None,self.root,['build'])
        self.assertEqual(load_saved(updated),{'exposure':.5})
        self.assertEqual(load_saved(options),{})
        Path(updated['style_profile']).write_text('{}')
        with self.assertRaisesRegex(ValueError,'snapshot changed'):configure({},updated,None,None,self.root,['audit'])
        with self.assertRaisesRegex(ValueError,'snapshot is corrupted'):configure({},updated,self.path,None,self.root,['build'])

    def test_environment_signature_and_removal(self):
        env=self.root/'environment';env.mkdir();manifest=env/'manifest.json'
        manifest.write_text(json.dumps({'kind':'schematic-cytosolic-environment','units':'nm'}));(env/'actin.npz').write_bytes(b'input fingerprint fixture')
        self.request({'environment':{'compartment':'cytosol','context_note':'reviewed illustrative layout','context_sources':['source']}})
        with self.assertRaisesRegex(ValueError,'environment-dir'):configure({}, {}, self.path,None,self.root,['build'])
        options={};configure(options,{},self.path,env,self.root,['build'])
        (env/'actin.npz').write_bytes(b'changed input')
        with self.assertRaisesRegex(ValueError,'requires rebuilding'):configure({},options,None,None,self.root,['audit'])
        self.request({});updated={};configure(updated,options,self.path,None,self.root,['build'])
        self.assertIsNone(updated['environment_dir'])
        manifest.write_text(json.dumps({'kind':'schematic-cytosolic-environment','units':'angstrom'}))
        self.request({'environment':{'compartment':'cytosol','context_note':'reviewed','context_sources':['source']}})
        with self.assertRaisesRegex(ValueError,'native-nm'):configure({}, {}, self.path,env,self.root,['build'])


if __name__=='__main__':unittest.main()
