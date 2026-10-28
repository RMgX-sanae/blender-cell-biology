"""Prepare a molecular style benchmark from coordinates, including annotated ribbons.

Run outside Blender with the same dependencies as prepare_molecular_asset.py.
Secondary structure comes from the deposited structure's annotations; unassigned
residues stay loops. This does not run DSSP or infer secondary structure by color.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import gemmi


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('structure', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.resolve().is_relative_to(Path(__file__).resolve().parents[1]):
        p.error('Prepared assets belong outside the skill directory.')
    structure = gemmi.read_structure(str(args.structure))
    names = [c.name for c in structure[0]]
    if len(set(n.casefold() for n in names)) != len(set(names)):
        p.error('Chain IDs differ only by case. Use a project-specific hex-encoded asset importer on case-insensitive filesystems.')
    subprocess.run([sys.executable, str(Path(__file__).with_name('prepare_molecular_asset.py')),
                    str(args.structure), '--output', str(args.output)], check=True)
    segments = []
    for helix in structure.helices:
        segments.append({'chain': helix.start.chain_name, 'start': helix.start.res_id.seqid.num,
                         'end': helix.end.res_id.seqid.num, 'kind': 'helix'})
    for sheet in structure.sheets:
        for strand in sheet.strands:
            segments.append({'chain': strand.start.chain_name, 'start': strand.start.res_id.seqid.num,
                             'end': strand.end.res_id.seqid.num, 'kind': 'sheet'})
    (args.output / 'secondary-structure.json').write_text(json.dumps({
        'source': args.structure.name, 'method': 'deposited HELIX/SHEET or mmCIF equivalent',
        'segments': segments,
        'limitations': 'Author numbering; insertion-code boundaries need project-specific handling. '
                       'No annotation means loops, not an inferred all-coil biological state.'}, indent=2))


if __name__ == '__main__':
    main()
