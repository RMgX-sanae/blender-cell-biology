"""Rebuild the separate native models against an installed or extracted Blender biology skill."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


PROJECTS = ("01_hemocyanin", "02_corti", "03_confocal")
FILES = (
    "views.py",
    "check_saved.py",
    "check_sources.py",
    "01_hemocyanin/prepare.py",
    "01_hemocyanin/build.py",
    "02_corti/prepare_field.py",
    "02_corti/build.py",
    "03_confocal/build.py",
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--blender", required=True)
    p.add_argument("--skill", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--source-cache", type=Path)
    p.add_argument("--width", type=int, default=1920)
    p.add_argument("--samples", type=int)
    a = p.parse_args()
    root = Path(__file__).resolve().parent
    out = a.out.resolve()
    skill = a.skill.resolve()
    if out.is_relative_to(skill) or out == root or (out.exists() and any(out.iterdir())):
        p.error("Choose a new or empty project directory outside the source and skill")
    if (
        a.width < 128
        or (a.samples is not None and a.samples < 1)
        or not (skill / "scripts/detail_scene.py").is_file()
    ):
        p.error("Invalid render settings or skill location")
    studio = root / "studio.py"
    if not studio.exists():
        studio = root.parent / "three-scale-studies/src/studio.py"
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(studio, out / "studio.py")
    for rel in FILES:
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / rel, out / rel)
    if a.source_cache:
        (out / "sources").mkdir()
        for n in ("6R83.cif.gz", "p21-labels.bin", "p21-calibration.xml"):
            source = a.source_cache / n
            if source.is_file():
                shutil.copyfile(source, out / "sources" / n)
    env = dict(
        os.environ,
        CELL_BIOLOGY_SKILL=str(skill),
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        PYTHONDONTWRITEBYTECODE="1",
    )

    def run(*args):
        subprocess.run(list(map(str, args)), check=True, cwd=out, env=env)

    run(sys.executable, out / "01_hemocyanin/prepare.py")
    run(sys.executable, out / "02_corti/prepare_field.py")
    run(sys.executable, out / "check_sources.py", out)
    for name in PROJECTS:
        proj = out / name
        samples = a.samples if a.samples is not None else (48 if name == "02_corti" else 64)
        run(
            a.blender,
            "-b",
            "-t",
            "5",
            "--factory-startup",
            "--python-exit-code",
            "1",
            "--python",
            proj / "build.py",
            "--",
            "--out",
            proj,
            "--width",
            a.width,
            "--samples",
            samples,
            "--save-only",
        )
        run(
            a.blender,
            "-b",
            "-t",
            "5",
            "--disable-autoexec",
            proj / "model.blend",
            "--python-exit-code",
            "1",
            "--python",
            out / "check_saved.py",
            "--",
            "--project",
            proj,
            "--out",
            proj,
            "--render",
        )


if __name__ == "__main__":
    main()
