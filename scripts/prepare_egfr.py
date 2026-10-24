#!/usr/bin/env python3
"""Create geometric translations only: no parameterization, relaxation or MD."""
import hashlib
import json
from pathlib import Path
import numpy as np


def xyz(lines):
    return np.array([[float(line[i:i+8]) for i in (30, 38, 46)] for line in lines]) / 10


def main():
    raw = Path('data/raw/1M17.pdb').read_bytes()
    manifest = json.loads(Path('results/data_manifest.json').read_text())
    if hashlib.sha256(raw).hexdigest() != manifest['sha256']:
        raise ValueError('PDB checksum differs from download manifest')
    lines = raw.decode().splitlines()
    protein = [l for l in lines if l.startswith('ATOM  ') and l[21] == 'A' and l[16] in (' ', 'A')]
    ligand = [l for l in lines if l.startswith('HETATM') and l[17:20] == 'AQ4' and l[21] == 'A']
    if not protein or not ligand:
        raise ValueError('Expected protein chain A and AQ4 ligand')
    p, lig = xyz(protein), xyz(ligand)
    direction = lig.mean(axis=0) - p.mean(axis=0)
    direction /= np.linalg.norm(direction)
    out = Path('results/egfr')
    out.mkdir(parents=True, exist_ok=True)
    (out/'protein.pdb').write_text('\n'.join(protein)+'\nEND\n')
    (out/'ligand_AQ4.pdb').write_text('\n'.join(ligand)+'\nEND\n')
    records = []
    for w, displacement in enumerate(np.linspace(0, 1.6, 17)):
        moved = lig + displacement * direction
        translated = [l[:30]+''.join(f'{v:8.3f}' for v in point*10)+l[54:] for l, point in zip(ligand, moved)]
        dest = out/f'window_{w:02d}.pdb'
        dest.write_text('\n'.join(protein+translated)+'\nEND\n')
        reread = [l for l in dest.read_text().splitlines() if l.startswith('HETATM')]
        actual = xyz(reread)
        projection = float((actual.mean(axis=0)-lig.mean(axis=0)) @ direction)
        min_dist = float(np.linalg.norm(p[:, None, :] - actual[None, :, :], axis=-1).min())
        records.append(dict(window=w, target_displacement_nm=float(displacement),
                            actual_displacement_nm=projection, error_nm=abs(projection-displacement),
                            minimum_heavy_atom_distance_nm=min_dist))
    report = dict(status='GEOMETRIC PREPARATION ONLY; EGFR PMF not done', protein_atoms=len(p),
                  ligand_atoms=len(lig), direction=direction.tolist(), windows=records,
                  expected_windows=17, actual_windows=len(list(out.glob('window_*.pdb'))),
                  max_projection_error_nm=max(r['error_nm'] for r in records))
    (out/'preparation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'windows'}, indent=2))


if __name__ == '__main__':
    main()
