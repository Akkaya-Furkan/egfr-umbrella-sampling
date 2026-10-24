#!/usr/bin/env python3
"""Resolve terminal naming and audit AM1-BCC rounding before LEaP."""
import json
from pathlib import Path
import numpy as np

base = Path('data/atomistic')
lines = (base/'protein_amber.pdb').read_text().splitlines()
first = next(l for l in lines if l.startswith('ATOM'))[21:27]
renamed = 0
for i,l in enumerate(lines):
    if l.startswith('ATOM') and l[21:27] == first and l[12:16].strip() == 'H':
        lines[i] = l[:12]+' H1 '+l[16:]; renamed += 1
(base/'protein_amber.pdb').write_text('\n'.join(lines)+'\n')
lines = (base/'ligand.mol2').read_text().splitlines()
a = lines.index('@<TRIPOS>ATOM')+1
b = lines.index('@<TRIPOS>BOND')
rows = [l.split() for l in lines[a:b]]
charges = np.array([float(row[8]) for row in rows])
old = float(charges.sum())
if abs(old)>.005: raise ValueError(f'Unexpected ligand net charge {old}')
# Distribute only AM1-BCC text-rounding residual; preserve relative charges.
corrected = charges-charges.mean()
for i,(row,q) in enumerate(zip(rows,corrected)):
    row[8] = f'{q:.8f}'
    lines[a+i] = ' '.join(row)
(base/'ligand.mol2').write_text('\n'.join(lines)+'\n')
sdf = (base/'AQ4.sdf').read_text().splitlines()
nheavy = int(sdf[3][:3])
crystal = np.array([[float(l[i:i+10]) for i in (0,10,20)] for l in sdf[4:4+nheavy]])
ligand_heavy = np.array([[float(v) for v in row[2:5]] for row in rows if not row[5].startswith('h')])
max_change = float(np.linalg.norm(crystal-ligand_heavy,axis=1).max())
if max_change>1e-5: raise ValueError(f'Ligand crystal coordinates changed: {max_change} A')
report = dict(terminal_h_renamed=renamed, ligand_atoms=len(rows),
              original_charge_e=old, corrected_charge_e=float(sum(float(lines[a+i].split()[8]) for i in range(len(rows)))),
              uniform_charge_correction_e=float(-charges.mean()),
              max_heavy_coordinate_change_angstrom=max_change,
              charge_method='GAFF2 / sqm AM1-BCC, neutral AQ4; uniform text-rounding correction')
Path('results/atomistic/parameterization.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
