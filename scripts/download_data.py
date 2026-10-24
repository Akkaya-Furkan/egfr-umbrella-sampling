#!/usr/bin/env python3
"""Download public 1M17; verify response and record immutable content hash."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen


def main():
    url = 'https://files.rcsb.org/download/1M17.pdb'
    with urlopen(url, timeout=60) as response:
        content = response.read()
    text = content.decode('ascii')
    if not text.startswith('HEADER') or 'AQ4' not in text or '\nATOM  ' not in text:
        raise ValueError('Unexpected PDB response')
    dest = Path('data/raw/1M17.pdb')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(content)
    manifest = dict(pdb_id='1M17', ligand='AQ4', source_url=url,
                    sha256=hashlib.sha256(dest.read_bytes()).hexdigest(), bytes=dest.stat().st_size)
    Path('results/data_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest))


if __name__ == '__main__':
    main()
