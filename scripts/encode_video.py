"""Encode the scene's contiguous PNG frames without changing their playback timing."""
import argparse
import json
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--fps", type=int)
    parser.add_argument("--ffmpeg", default=shutil.which("ffmpeg"))
    args = parser.parse_args()
    manifest = args.directory / "manifest.json"
    metadata = json.loads(manifest.read_text()) if manifest.exists() else {}
    recorded_fps = metadata.get("fps", metadata.get("parameters", {}).get("fps"))
    recorded_frames = metadata.get("frames", metadata.get("parameters", {}).get("frames"))
    if args.fps is not None and recorded_fps is not None and args.fps != recorded_fps:
        parser.error("FPS differs from scene manifest; perform intentional retiming separately")
    args.fps = args.fps if args.fps is not None else recorded_fps
    if not args.ffmpeg or args.fps is None or args.fps <= 0:
        parser.error("Supply an existing --ffmpeg executable and positive --fps")
    frames = sorted((args.directory / "frames").glob("frame_*.png"))
    if not frames or [p.name for p in frames] != [f"frame_{i:04d}.png" for i in range(1, len(frames)+1)]:
        parser.error("Expected a complete contiguous sequence starting with frame_0001.png")
    if isinstance(recorded_frames, int) and len(frames) != recorded_frames:
        parser.error(f"Expected {recorded_frames} frames from manifest; found {len(frames)}")
    destination = args.directory / "preview.mp4"
    subprocess.run([str(args.ffmpeg), "-v", "error", "-n", "-framerate", str(args.fps),
                    "-start_number", "1", "-i", str(args.directory / "frames/frame_%04d.png"),
                    "-c:v", "libx264", "-crf", "19", "-pix_fmt", "yuv420p",
                    "-movflags", "+faststart", str(destination)], check=True)
    print(destination)


if __name__ == "__main__":
    main()
