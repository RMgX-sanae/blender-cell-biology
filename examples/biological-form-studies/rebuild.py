"""Rebuild the three independent teaching models from their raw inputs in an empty directory."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
PROJECTS = ("01_spliceosome", "02_foraminifer", "03_ultramicrotome")
FILES = (
    "check_saved.py",
    "check_sources.py",
    "01_spliceosome/prepare.py",
    "01_spliceosome/build.py",
    "02_foraminifer/prepare.py",
    "02_foraminifer/build.py",
    "03_ultramicrotome/build.py",
)
RAW = "https://download.pangaea.de/dataset/949585/files/GE_GAD_2_738x765x384-1.6827_um-ushort_BE.raw"
SHA = "c8874ba38ca0740a7c9b78008a17224f19ecd637446de9cb75bb0c70284c0906"


def get_scan(path):
    tmp = path.with_suffix(".part")
    # PANGAEA may respond 503 while the public specimen is being restored from tape.
    deadline = time.monotonic() + 900
    while True:
        try:
            with urlopen(RAW, timeout=120) as src, tmp.open("wb") as dst:
                shutil.copyfileobj(src, dst, 1024 * 1024)
            break
        except HTTPError as exc:
            if exc.code != 503 or time.monotonic() >= deadline:
                raise
            time.sleep(min(30, max(5, int(exc.headers.get("Retry-After", "15")))))
    with tmp.open("rb") as f:
        if hashlib.file_digest(f, "sha256").hexdigest() != SHA:
            raise ValueError("Downloaded CT specimen differs from the reviewed GE_GAD_2 volume")
    tmp.replace(path)


def main():
    q = argparse.ArgumentParser(description=__doc__)
    q.add_argument("--out", type=Path, required=True)
    q.add_argument("--blender", type=Path, required=True)
    q.add_argument("--skill", type=Path, required=True)
    q.add_argument("--font-en", type=Path, required=True)
    q.add_argument("--source-cache", type=Path)
    q.add_argument("--width", type=int, default=1920)
    q.add_argument("--samples", type=int, default=54)
    a = q.parse_args()
    out = a.out.resolve()
    skill = a.skill.resolve()
    if (out.exists() and any(out.iterdir())) or out.is_relative_to(skill) or out == ROOT:
        q.error("Choose an empty project outside the skill")
    if (
        not a.font_en.is_file()
        or not a.blender.is_file()
        or not (skill / "scripts/detail_scene.py").is_file()
    ):
        q.error("Provide Blender, the Times New Roman Regular font and the skill root")
    if a.width < 128 or a.samples < 1:
        q.error("Invalid resolution or sample count")
    out.mkdir(parents=True, exist_ok=True)
    for rel in FILES:
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / rel, out / rel)
    for name, sibling in [
        ("studio.py", "three-scale-studies/src/studio.py"),
        ("views.py", "structure-and-function/views.py"),
    ]:
        src = ROOT / name if (ROOT / name).exists() else ROOT.parent / sibling
        shutil.copyfile(src, out / name)
    data = out / "sources"
    data.mkdir()
    if a.source_cache:
        for name in ["7B9V.cif.gz", "GE_GAD_2.raw"]:
            src = a.source_cache / name
            if src.exists():
                try:
                    os.link(src, data / name)
                except OSError:
                    shutil.copyfile(src, data / name)
    if not (data / "GE_GAD_2.raw").exists():
        get_scan(data / "GE_GAD_2.raw")
    env = dict(
        os.environ,
        CELL_BIOLOGY_SKILL=str(skill),
        CELL_FONT_EN=str(a.font_en.resolve()),
        PYTHONDONTWRITEBYTECODE="1",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
    )
    for name in PROJECTS:
        d = out / name
        if (d / "prepare.py").is_file():
            subprocess.run([sys.executable, d / "prepare.py"], env=env, check=True)
        subprocess.run(
            [
                a.blender,
                "-b",
                "--factory-startup",
                "--python-exit-code",
                "1",
                "--python",
                d / "build.py",
                "--",
                "--save-only",
                "--width",
                str(a.width),
                "--samples",
                str(a.samples),
                "--out",
                d,
            ],
            env=env,
            check=True,
        )
        subprocess.run(
            [
                a.blender,
                "-b",
                "--disable-autoexec",
                d / "model.blend",
                "--python-exit-code",
                "1",
                "--python",
                out / "check_saved.py",
                "--",
                "--out",
                d,
                "--render",
            ],
            env=env,
            check=True,
        )
    subprocess.run([sys.executable, out / "check_sources.py", str(out)], env=env, check=True)
    print(json.dumps({"built": list(PROJECTS), "width": a.width, "directory": str(out)}))


if __name__ == "__main__":
    main()
