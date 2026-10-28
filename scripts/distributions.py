"""Deterministic schematic positions. Units match the caller's scene, never assumed biological."""
import math
import random


def sphere_points(count, radius=1.0, seed=23):
    """Uniform surface area on a sphere, not uniform latitude."""
    if count < 0 or radius <= 0:
        raise ValueError("count must be nonnegative and radius positive")
    rng = random.Random(seed)
    points = []
    for _ in range(count):
        z = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2 * math.pi)
        r = math.sqrt(1 - z * z)
        points.append((radius * r * math.cos(phi), radius * r * math.sin(phi), radius * z))
    return points


def patch_points(count, width=10.0, height=6.0, mode="uniform", min_distance=0.0, seed=23):
    """Sample a bounded XY patch; clusters and polarization are explicit hypotheses."""
    if count < 0 or min_distance < 0 or width <= 0 or height <= 0:
        raise ValueError("Invalid count, distance or patch dimensions")
    if mode not in {"uniform", "clustered", "polarized"}:
        raise ValueError("mode must be uniform, clustered or polarized")
    rng = random.Random(seed)
    points = []
    centers = [(-0.28 * width, -0.15 * height), (0.24 * width, 0.19 * height)]
    for _ in range(max(1000, count * 1000)):
        if len(points) == count:
            return points
        if mode == "clustered":
            cx, cy = rng.choice(centers)
            x, y = rng.gauss(cx, width * 0.10), rng.gauss(cy, height * 0.12)
        else:
            x, y = rng.uniform(-width / 2, width / 2), rng.uniform(-height / 2, height / 2)
            if mode == "polarized" and rng.random() > 0.08 + 0.92 * (x / width + 0.5):
                continue
        if not (-width / 2 <= x <= width / 2 and -height / 2 <= y <= height / 2):
            continue
        if all((x - a) ** 2 + (y - b) ** 2 >= min_distance ** 2 for a, b in points):
            points.append((x, y))
    raise ValueError(f"Packing not reached: {len(points)}/{count}. Reduce density/footprint explicitly.")


def triangle_points(triangles, count, seed=23):
    """Area-weighted sampling of WORLD-SPACE triangles. Degenerate faces are ignored."""
    rng = random.Random(seed)
    valid, cumulative, total = [], [], 0.0
    for a, b, c in triangles:
        u, v = [b[i] - a[i] for i in range(3)], [c[i] - a[i] for i in range(3)]
        cross = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
        length = math.sqrt(sum(x*x for x in cross))
        if length == 0:
            continue
        total += length / 2
        valid.append((a, b, c, tuple(x/length for x in cross)))
        cumulative.append(total)
    if not valid or count < 0:
        raise ValueError("Need nondegenerate triangles and a nonnegative count")
    import bisect
    result = []
    for _ in range(count):
        a, b, c, normal = valid[bisect.bisect_left(cumulative, rng.random() * total)]
        u, v = math.sqrt(rng.random()), rng.random()
        p = tuple((1-u)*a[i] + u*(1-v)*b[i] + u*v*c[i] for i in range(3))
        result.append((p, normal))
    return result
