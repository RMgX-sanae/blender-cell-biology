"""A geometrical interpretation of the pore hierarchy of Coscinodiscus granii."""

from pathlib import Path
import os
from collections import defaultdict
import json
import numpy as np
from scipy.spatial import Voronoi, cKDTree

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SHOWCASE_PROJECT", SOURCE_DIR))
OUT = ROOT / "prepared"
OUT.mkdir(exist_ok=True)
RINGS = 23
PITCH = 0.95
rng = np.random.default_rng(3401)
seeds = []
keys = []
for k in range(RINGS + 2):
    n = max(1, 6 * k)
    for j in range(n):
        t = 2 * np.pi * j / n
        r = PITCH * k
        # Six gently curved sectors, with local variation in the row spacing.
        rr = r * (1 + 0.023 * np.cos(6 * t + 0.29 * k / RINGS) + 0.011 * np.cos(3 * t + 0.6)) if k else 0
        ang = (
            t
            + 0.105 * (k / RINGS) ** 1.08
            + 0.022 * np.sin(3 * t + 0.085 * k)
            + 0.009 * np.sin(7 * t - 0.045 * k)
        )
        point = np.array([rr * np.cos(ang), rr * np.sin(ang)]) + rng.normal(0, 0.045, 2) * min(k, 1)
        seeds.append(point)
        keys.append((k, j))
seeds = np.array(seeds)
vor = Voronoi(seeds)


def ccw(p):
    return p if np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)) > 0 else p[::-1]


def dome(x, y):
    rr = (x * x + y * y) / (RINGS * PITCH) ** 2
    return (
        0.78 * (1 - rr)
        + 0.24 * (x * x - y * y) / (RINGS * PITCH) ** 2
        + 0.065 * np.sin(0.27 * x + 0.17 * y)
    )


def densify(p, n):
    return np.vstack([a + (b - a) * t for a, b in zip(p, np.roll(p, -1, axis=0)) for t in np.arange(n) / n])


def subdiv(poly, c, sites):
    # Clip each subcell against bisectors; every pore occupies a disjoint area.
    out = []
    for i, s in enumerate(sites):
        p = poly.copy()
        for j, s2 in enumerate(sites):
            if i == j:
                continue
            normal = s2 - s
            d = (p - (s + s2) * 0.5) @ normal
            q = []
            for v, w, dv, dw in zip(p, np.roll(p, -1, axis=0), d, np.roll(d, -1)):
                if dv <= 1e-10:
                    q.append(v)
                if (dv < -1e-10 and dw > 1e-10) or (dv > 1e-10 and dw < -1e-10):
                    q.append(v + (w - v) * dv / (dv - dw))
            p = np.array(q)
            if len(p) < 3:
                break
        if len(p) >= 3:
            out.append((ccw(p), s))
    return out


class Builder:
    def __init__(self):
        self.v = []
        self.f = []

    def strip(self, rings, reverse=False):
        n = len(rings[0][0])
        start = len(self.v)
        for pts, z in rings:
            self.v.extend((x, y, dome(x, y) + z) for x, y in pts)
        for j in range(len(rings) - 1):
            a = start + j * n
            b = a + n
            for i in range(n):
                h = (i + 1) % n
                fs = [(a + i, a + h, b + h), (a + i, b + h, b + i)]
                self.f.extend([tuple(reversed(f)) for f in fs] if reverse else fs)

    def arrays(self):
        return np.array(self.v, dtype="f"), np.array(self.f, dtype="i")


def aperture(b, poly, c, z, thick, fraction, segments):
    outer = densify(poly, segments)
    directions = outer - c
    # Minimum distance to a wall avoids crossing a neighbouring opening.
    e = np.roll(poly, -1, axis=0) - poly
    dist = np.abs(e[:, 0] * (poly[:, 1] - c[1]) - e[:, 1] * (poly[:, 0] - c[0])) / np.linalg.norm(e, axis=1)
    radius = dist.min() * fraction
    unit = directions / np.linalg.norm(directions, axis=1)[:, None]
    # Mild ellipticity follows the local radial row direction seen in the proximal reference.
    if thick >= 0.09:
        th = np.arctan2(unit[:, 1], unit[:, 0])
        phi = np.arctan2(c[1], c[0])
        radius = radius * (1 + 0.10 * np.cos(2 * (th - phi)) + 0.028 * np.cos(3 * th + 2 * phi))[:, None]
    rim = c + unit * radius
    b.strip(
        [
            (outer, z),
            (c + directions * 0.82 + unit * radius * 0.09, z + 0.002),
            (c + unit * radius * 1.08, z),
            (rim, z - thick * 0.26),
            (rim, z - thick * 0.74),
            (c + unit * radius * 1.07, z - thick),
            (outer, z - thick),
        ]
    )
    return radius


def sites7(c, r, a):
    theta = np.arange(6) * np.pi / 3 + a
    return np.vstack([c, c + np.stack([np.cos(theta), np.sin(theta)], axis=1) * r])


def make_cell(parts, poly, c, idx, micro):
    directions = poly - c
    pitch = np.linalg.norm(directions, axis=1).mean()
    parts["wall"].strip(
        [
            (c + directions * 0.994, 0.57),
            (c + directions * 0.90, 0.57),
            (c + directions * 0.88, 0.48),
            (c + directions * 0.89, -0.05),
            (c + directions * 0.93, -0.13),
            (c + directions * 0.994, -0.13),
            (c + directions * 0.994, 0.57),
        ]
    )
    aperture(parts["foramen"], poly, c, 0.67, 0.10, 0.61 + 0.044 * np.sin(idx * 1.788), 5)
    ang = np.arctan2(directions[0, 1], directions[0, 0])
    floor = -0.099
    for sp, ctr in subdiv(poly, c, sites7(c, pitch * 0.56, ang)):
        rad = aperture(parts["cribrum"], sp, ctr, floor, 0.071, 0.67, 3)
        if micro:
            circ = ctr + rad * 1.16 * np.stack(
                [np.cos(np.arange(30) * np.pi / 15), np.sin(np.arange(30) * np.pi / 15)], axis=1
            )
            for npoly, nctr in subdiv(circ, ctr, sites7(ctr, rad * 0.59, ang + 0.26)):
                aperture(parts["cribellum"], npoly, nctr, floor - 0.075, 0.018, 0.43, 3)


ids = [i for i, (k, j) in enumerate(keys) if k <= RINGS]
pick = keys.index((8, 14))
closest = set(cKDTree(seeds).query(seeds[pick], k=7)[1])
whole = defaultdict(Builder)
micro = defaultdict(Builder)
for i in ids:
    vi = vor.regions[vor.point_region[i]]
    assert -1 not in vi
    poly = ccw(vor.vertices[vi])
    c = seeds[i]
    make_cell(whole, poly, c, i, False)
    if i in closest:
        make_cell(micro, poly, c, i, True)
edges = []
for (i, j), ij in zip(vor.ridge_points, vor.ridge_vertices):
    if (keys[i][0] <= RINGS) != (keys[j][0] <= RINGS) and len(ij) == 2 and -1 not in ij:
        edges.extend(ij)
idx = np.unique(edges)
rim = vor.vertices[idx]
rim = rim[np.argsort(np.arctan2(rim[:, 1], rim[:, 0]))]
np.save(OUT / "rim.npy", rim)
unit = rim / np.linalg.norm(rim, axis=1)[:, None]
whole["mantle"].strip(
    [
        (rim, 0.665),
        (rim + unit * 0.17, 0.69),
        (rim + unit * 0.29, 0.62),
        (rim + unit * 0.30, 0.48),
        (rim + unit * 0.29, -0.47),
        (rim + unit * 0.23, -0.55),
        (rim + unit * 0.10, -0.54),
        (rim + unit * 0.11, -0.13),
        (rim, -0.13),
    ],
    reverse=True,
)
for j in range(17):
    t = 2 * np.pi * (j + 0.12 * np.sin(j * 2.1)) / 17
    r = RINGS * PITCH * (1 + 0.023 * np.cos(6 * t) + 0.011 * np.cos(3 * t + 0.6)) - 0.43
    c = np.array([r * np.cos(t), r * np.sin(t)])
    tang = np.array([-np.sin(t), np.cos(t)])
    radial = np.array([np.cos(t), np.sin(t)])
    a = np.arange(30) * 2 * np.pi / 30
    ellipse = 0.12 * np.cos(a)[:, None] * tang + 0.075 * np.sin(a)[:, None] * radial
    whole["rimoportula"].strip(
        [
            (c + ellipse * 1.55, 0.65),
            (c + ellipse * 1.15, 0.72),
            (c + ellipse + 0.10 * radial, 1.0),
            (c + ellipse * 0.72 + 0.10 * radial, 1.03),
            (c + ellipse * 0.42 + 0.10 * radial, 1.01),
            (c + ellipse * 0.42, 0.66),
        ]
    )
for name, data in [("valve", whole), ("window", micro)]:
    parts = {}
    for role, b in data.items():
        v, f = b.arrays()
        if name == "window":
            v -= np.array([*seeds[pick], dome(*seeds[pick])], dtype="f")
        parts[role + "_v"] = v
        parts[role + "_f"] = f
    np.savez_compressed(OUT / (name + ".npz"), **parts)
info = {
    "species_reference": "Coscinodiscus granii",
    "reference": "10.1371/journal.pone.0014300, figures 1 and 7",
    "model": "procedural illustration; initial dimensions in micrometres, representative, not morphometry",
    "seed": 3401,
    "areolae": len(ids),
    "radial_pitch_um": PITCH,
    "window_original_cell_ids": sorted(map(int, closest)),
    "levels": "proximal foramen -> chamber walls -> distal cribrum -> finer distal cribellum; fine sieve resolved in detail view",
    "cribrum_apertures_in_whole": len(ids) * 7,
    "window_areolae": len(closest),
    "fine_openings_in_window": len(closest) * 49,
}
(OUT / "design.json").write_text(json.dumps(info, indent=2))
print(info)
