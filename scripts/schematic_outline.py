"""Shape authored envelopes and their attached structures in local coordinates."""


def curved_point(point, radii, *, bend, taper):
    """Bend the long axis and vary cross-section width smoothly from end to end."""
    x, y, z = point
    along = x / radii[0]
    width = 1 - taper * along
    offset = bend * radii[1] * along * along
    return x, offset + width * y, width * z


def curve_surface(obj, radii, *, bend, taper):
    for vertex in obj.data.vertices:
        vertex.co = curved_point(vertex.co, radii, bend=bend, taper=taper)
    obj.data.update()
