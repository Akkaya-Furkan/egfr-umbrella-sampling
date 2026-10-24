#!/usr/bin/env bash
set -euo pipefail
export AMBERHOME="${AMBERHOME:-$HOME/.conda/envs/md_env}"
export PATH="$AMBERHOME/bin:$PATH"
export OMP_NUM_THREADS=4
python_md="${OPENMM_PYTHON:-$HOME/.conda/envs/openmm_env/bin/python}"
"$python_md" scripts/repair_atomistic.py
cd data/atomistic
obabel AQ4.sdf -O ligand_h.mol2 -h
antechamber -i ligand_h.mol2 -fi mol2 -o ligand.mol2 -fo mol2 -at gaff2 -c bcc -nc 0 -rn AQ4 -s 2
pdb4amber -i protein_repaired.pdb -o protein_amber.pdb
(cd ../.. && "$python_md" scripts/finalize_atomistic_inputs.py)
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2
cat > leap.in <<'LEAP'
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
lig = loadmol2 ligand.mol2
prot = loadpdb protein_amber.pdb
complex = combine {prot lig}
check complex
saveamberparm complex dry.prmtop dry.inpcrd
solvatebox complex TIP3PBOX 20.0
addionsrand complex Na+ 0
addionsrand complex Cl- 0
saveamberparm complex solvated.prmtop solvated.inpcrd
savepdb complex solvated.pdb
quit
LEAP
tleap -f leap.in
