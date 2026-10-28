"""Download public RCSB coordinates and provenance, without executing remote content."""
import argparse
import hashlib
import json
import re
import time
from http.client import IncompleteRead
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


def download(url, limit=50_000_000):
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"User-Agent": "blender-cell-biology/1.0"}), timeout=45) as response:
                data = response.read(limit + 1)
            break
        except HTTPError as error:
            if error.code not in {429,500,502,503,504} or attempt==2:raise
        except (URLError,TimeoutError,ConnectionResetError,IncompleteRead):
            if attempt==2:raise
        time.sleep(2**(attempt+1))
    if len(data) > limit:
        raise ValueError("File exceeds 50 MB; select/download a bounded structure explicitly")
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdb_id")
    parser.add_argument("--format", choices=("cif", "pdb"), default="cif")
    parser.add_argument("--assembly", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    entry = args.pdb_id.upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", entry):
        parser.error("This small helper accepts legacy four-character PDB IDs; use official tools for extended IDs")
    if args.assembly is not None and (args.format != "cif" or args.assembly < 1):
        parser.error("Assembly download requires --format cif and an assembly number >= 1")
    name = f"{entry}-assembly{args.assembly}.cif" if args.assembly else f"{entry}.{args.format}"
    url = f"https://files.rcsb.org/download/{name}"
    path = args.output / name
    meta_path = args.output / (name + ".source.json")
    if path.exists() or meta_path.exists():
        parser.error("Destination exists; reuse it or choose a new output directory")
    coordinates = download(url)
    if (args.format == "cif" and not coordinates.lstrip().startswith(b"data_")) or (
        args.format == "pdb" and b"\nATOM  " not in coordinates):
        raise ValueError("Response is not the requested coordinate format")
    metadata_url = f"https://data.rcsb.org/rest/v1/core/entry/{entry}"
    metadata = json.loads(download(metadata_url))
    args.output.mkdir(parents=True, exist_ok=True)
    path.write_bytes(coordinates)
    provenance = {
        "id": entry, "url": url, "entry_url": f"https://www.rcsb.org/structure/{entry}",
        "metadata_url": metadata_url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": hashlib.sha256(coordinates).hexdigest(), "bytes": len(coordinates),
        "assembly": args.assembly, "coordinate_unit": "angstrom", "evidence": "experimental_structure_model",
        "license": "wwPDB structural data: CC0; cite entry and primary publication",
        "license_url": "https://www.ebi.ac.uk/pdbe/about/public-data-access-statement",
        "metadata": metadata,
    }
    meta_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"coordinates": str(path), "provenance": str(meta_path), "bytes": len(coordinates)}))


if __name__ == "__main__":
    main()
