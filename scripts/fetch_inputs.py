#!/usr/bin/env python3
"""Fetch public, commit-pinned files; fail closed on checksum mismatch."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]

def fetch(destination=None):
    destination = Path(destination or ROOT / '.cache' / 'inputs')
    manifest = json.loads((ROOT / 'data/input_manifest.json').read_text())
    for name, digest in manifest['files'].items():
        target = destination / name
        if not target.exists():
            url = f"https://raw.githubusercontent.com/EngreitzLab/CRISPR_comparison/{manifest['commit']}/{name}"
            with urlopen(url, timeout=90) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != digest:
                raise RuntimeError(f'Checksum mismatch: {name}')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
        if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            raise RuntimeError(f'Local input checksum mismatch: {name}')
    return destination

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path)
    args = parser.parse_args()
    fetch(args.destination)
    print('Verified all nine pinned input/config/license files.')
