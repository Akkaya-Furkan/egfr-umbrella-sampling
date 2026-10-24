#!/usr/bin/env python3
"""Repair internal missing residues/sidechains; leave unresolved termini truncated."""
import json
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

out = Path('data/atomistic')
f = PDBFixer(filename='data/raw/1M17.pdb')
f.findMissingResidues()
missing = {str(k): v for k, v in f.missingResidues.items()}
chains = list(f.topology.chains())
f.missingResidues = {k: v for k, v in f.missingResidues.items()
                     if k[1] not in (0, len(list(chains[k[0]].residues())))}
modeled = {str(k): v for k, v in f.missingResidues.items()}
f.removeHeterogens(False)
f.findNonstandardResidues()
f.replaceNonstandardResidues()
f.findMissingAtoms()
atoms = {str(r): [a.name for a in aa] for r, aa in f.missingAtoms.items()}
f.addMissingAtoms(seed=23024)
f.addMissingHydrogens(7.4)
# Encode selected histidine tautomers explicitly for LEaP.
for res in f.topology.residues():
    if res.name == 'HIS':
        names = {a.name for a in res.atoms()}
        res.name = 'HIP' if {'HD1', 'HE2'} <= names else ('HID' if 'HD1' in names else 'HIE')
with (out/'protein_repaired.pdb').open('w') as handle:
    PDBFile.writeFile(f.topology, f.positions, handle, keepIds=True)
report = dict(pH=7.4, seed=23024, missing_residue_blocks=missing,
              modeled_internal_blocks=modeled, missing_heavy_atoms=atoms,
              residues=sum(1 for _ in f.topology.residues()),
              atoms=sum(1 for _ in f.topology.atoms()),
              limitations=['Internal loop is modeled, not experimentally resolved.',
                           'Unresolved N/C terminal residues omitted.',
                           'Default pH-based protonation; no pKa calculation.'])
Path('results/atomistic/repair.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report), flush=True)
