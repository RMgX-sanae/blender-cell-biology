"""Behavioral checks for geometry sampling; run with ordinary Python."""
import math
import unittest
from distributions import patch_points, sphere_points, triangle_points


class SamplingTests(unittest.TestCase):
    def test_surface_area_uniformity(self):
        points=sphere_points(20000,radius=2,seed=7)
        self.assertLess(max(abs(sum(x*x for x in p)-4) for p in points),1e-10)
        # Equal-area hemisphere and second moment checks catch latitude sampling bias.
        self.assertLess(abs(sum(p[2]>0 for p in points)/len(points)-.5),.02)
        self.assertLess(abs(sum(p[2]**2 for p in points)/len(points)-4/3),.04)

    def test_reproducibility_and_exclusion(self):
        for mode in ("uniform","clustered","polarized"):
            a=patch_points(32,mode=mode,min_distance=.6,seed=9)
            self.assertEqual(a,patch_points(32,mode=mode,min_distance=.6,seed=9))
            self.assertTrue(all(-5<=x<=5 and -3<=y<=3 for x,y in a))
            self.assertGreaterEqual(min(math.dist(p,q) for i,p in enumerate(a) for q in a[i+1:]),.6)

    def test_impossible_packing_fails(self):
        with self.assertRaises(ValueError):
            patch_points(2,width=1,height=1,min_distance=2)

    def test_mesh_area_weighting(self):
        triangles=[((0,0,0),(1,0,0),(0,1,0)),((10,0,0),(13,0,0),(10,3,0))]
        sample=triangle_points(triangles,10000,seed=1)
        self.assertLess(abs(sum(p[0]>5 for p,n in sample)/len(sample)-.9),.02)
        self.assertTrue(all(abs(sum(v*v for v in n)-1)<1e-12 for p,n in sample))

    def test_polarization_changes_expected_location(self):
        u=patch_points(2000,mode="uniform",seed=4)
        p=patch_points(2000,mode="polarized",seed=4)
        self.assertGreater(sum(x for x,y in p)/len(p)-sum(x for x,y in u)/len(u),1)


if __name__=="__main__":
    unittest.main()
