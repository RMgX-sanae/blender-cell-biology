"""Prepare, save and reopen each independent project. Run with the scientific Python environment."""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parent
p = argparse.ArgumentParser(description=__doc__)
p.add_argument("projects", nargs="+", choices=["all", "antenna", "diatom", "pipette"])
p.add_argument("--blender", required=True, help="Blender executable")
p.add_argument("--skill", type=Path, default=Path.home() / ".codex/skills/blender-cell-biology")
p.add_argument(
    "--fonts", type=Path, help="Directory containing Songti-SC-Regular.ttf and Times-New-Roman.ttf"
)
p.add_argument("--width", type=int, default=1920)
p.add_argument("--samples", type=int, default=48)
p.add_argument("--source-cache", type=Path, help="Directory with the reviewed 6KGX.cif.gz")
p.add_argument("--out", type=Path, required=True, help="Project directory outside the installed skill")
args = p.parse_args()
if args.out.resolve().is_relative_to(args.skill.resolve()):
    p.error("--out must be outside the installed skill")
if any(k in args.projects for k in ("pipette", "all")) and not (
    args.fonts or os.environ.get("SHOWCASE_FONTS")
):
    p.error("Provide --fonts for the pipette readout and labels")
projects = {"antenna": "01_phycobilisome", "diatom": "02_diatom", "pipette": "03_pipette"}
env = dict(
    os.environ,
    CELL_BIOLOGY_SKILL=str(args.skill.resolve()),
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1",
)
if args.fonts:
    env["SHOWCASE_FONTS"] = str(args.fonts.resolve())
for key in projects if "all" in args.projects else args.projects:
    work = root / projects[key]
    output = (args.out / projects[key]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if key == "antenna" and args.source_cache:
        source = output / "sources/6KGX.cif.gz"
        source.parent.mkdir(exist_ok=True)
        cached = args.source_cache.resolve() / "6KGX.cif.gz"
        if not source.exists():
            shutil.copyfile(cached, source)
    env["SHOWCASE_PROJECT"] = str(output)
    for name in ("prepare.py", "prepare_views.py"):
        if (work / name).exists():
            subprocess.run([sys.executable, str(work / name)], env=env, check=True)
    common = [args.blender, "-b", "-t", "5", "--factory-startup", "--python-exit-code", "1"]
    subprocess.run(
        [
            *common,
            "-P",
            str(work / "build.py"),
            "--",
            "--out",
            str(output),
            "--width",
            str(args.width),
            "--samples",
            str(args.samples),
            "--save-only",
        ],
        env=env,
        check=True,
    )
    subprocess.run(
        [
            args.blender,
            "--disable-autoexec",
            "-b",
            str(output / "model.blend"),
            "-t",
            "5",
            "--python-exit-code",
            "1",
            "-P",
            str(root / "src/render_saved.py"),
            "--",
            "--project",
            key,
            "--output",
            str(output),
        ],
        env=env,
        check=True,
    )
