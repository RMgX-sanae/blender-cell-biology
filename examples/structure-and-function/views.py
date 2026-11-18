"""Scene-local studio and editable, closed ring meshes for the three studies."""

from math import cos, sin, tau
import bpy
from mathutils import Vector
import studio as st
from detail_scene import compositors


def view(scene, target, eye, span, palette, *, roll=0, fstop=None, light=0.60):
    target = Vector(target)
    direction = Vector(eye).normalized()
    distance = span * 2.65
    c = st.camera(scene, target + direction * distance, target, span, perspective=True, fstop=fstop)
    c.rotation_euler.rotate_axis("Z", roll)
    st.studio(scene, c, target, span, palette, intensity=light)
    compositors(scene, {"background": palette[-1], "compositor": {"preset": "reference", "contrast": 1.8}})
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "VIEW_3D":
                area.spaces.active.shading.color_type = "MATERIAL"
    return c


def loft(name, sections, col, mat, *, sides=32, phase=0, ripple=0):
    """Horizontal elliptical profiles. Each row gives center xyz and radii xy."""
    verts, faces = [], []
    for j, (x, y, z, rx, ry) in enumerate(sections):
        for k in range(sides):
            a = tau * k / sides
            r = 1 + ripple * (cos(3 * a + 0.24 * j + phase) + 0.38 * sin(5 * a - 0.08 * j + phase))
            verts.append((x + rx * cos(a) * r, y + ry * sin(a) * r, z))
        if j:
            a0, a1 = (j - 1) * sides, j * sides
            faces.extend((a0 + k, a0 + (k + 1) % sides, a1 + (k + 1) % sides, a1 + k) for k in range(sides))
    faces.extend([tuple(reversed(range(sides))), tuple(range(len(verts) - sides, len(verts)))])
    return st.mesh(name, verts, faces, col, mat)


def bio_material(name, color, *, rough=0.43, sss=0.055, bump=0.019, scale=3.5):
    return st.material(
        name,
        color,
        settings=dict(
            roughness=rough,
            subsurface_weight=sss,
            subsurface_radius=[0.9, 0.46, 0.26],
            subsurface_scale=0.12,
            specular=0.34,
            bump=bump,
            bump_strength=0.19,
            texture_scale=scale,
            texture_detail=2.5,
        ),
    )
