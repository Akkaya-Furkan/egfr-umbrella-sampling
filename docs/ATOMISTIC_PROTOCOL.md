# Atomistic pilot protocol

This continuation uses public 1M17 chain A and AQ4 (erlotinib). It is distinct
from the double-well validation and the earlier 17 clashing rigid translations.
See [SCOPE](SCOPE.md) for the gates specified before atomistic production.

Preparation repairs internal missing residues with PDBFixer at pH 7.4. The 12
residues 965–976 are modeled; six N-terminal and three C-terminal unresolved
residues are omitted. The retained protein has 324 residues. Default histidine
hydrogen choices are encoded as HID/HIE/HIP for Amber. The terminal amine's H
name is changed to H1 for LEaP. This is a naming fix, not a protonation change.
No missing sidechain heavy atoms were found in the resolved residues.

The RCSB ModelServer SDF provides crystal-instance heavy-atom coordinates and
bond orders. Open Babel adds hydrogens without moving the heavy atoms.
AmberTools 24.8 antechamber uses `-at gaff2 -c bcc -nc 0`; parmchk2 supplies
missing GAFF2 terms. The text charges summed to -0.001 e; a uniform correction
of +0.000019230769 e per atom makes the 52-atom ligand neutral to text precision.
This small, recorded charge correction is not a different charge method.
Read-back audits compare ligand coordinates against the instance SDF and check
net charges, group identities, topology sizes, and SHA256 hashes.

LEaP combines ff14SB protein and GAFF2 ligand, solvates with a rectangular TIP3P
box using 20 Å padding, and adds eight Na+ counterions. There is no added salt.
The system has 132590 atoms and 42439 waters. The initial LEaP build failed on
terminal H naming; after correction, regenerated dry and solvated topologies
were read back successfully. Those local inputs and full logs are ignored by git.

OpenMM uses PME with 1 nm real-space cutoff, H-bond constraints, rigid water,
LangevinMiddle at 300 K, friction 1/ps, and 2 fs steps. All 324 CA atoms have a
100 kJ/mol/nm² positional restraint, including during production. The pilot is
therefore conditional on a restrained protein; it is not an unrestrained binding
free energy. Energy minimization is capped at 2000 iterations. NVT lasts 20 ps;
NPT lasts 80 ps with a 1 bar Monte Carlo barostat every 25 steps. Production uses
the fixed equilibrated box. Temperature/density traces are retained as diagnostics;
100 ps of equilibration cannot establish equilibration of all protein modes.

The reaction coordinate is the minimum-image, mass-weighted distance between
all 2594 protein heavy atoms and the 29 ligand heavy atoms. A 1000 kJ/mol/nm²
harmonic bias is centered at 1.1649579038512015 to 2.7649579038512018 nm in 17
0.1 nm increments. It is a radial coordinate, not a fixed vector. After 20 ps of
bound-state umbrella settling, the center moves by 1.6 nm in 400 ps of SMD.
Every 25 ps supplies a candidate starting state; each receives a new seeded
Maxwell velocity draw, 20 ps window burn-in and a target 400 ps production.
Actual durations and completion, not these targets, determine acceptance.
Scalar COM distances are stored every 0.1 ps. No molecular trajectory is written.
MD seeds are recorded; PDB hydrogen placement and LEaP random ion placement may
differ on a fresh preparation, whose actual files must be audited again.
Each relaxed start must have a minimum protein–ligand heavy distance >=0.15 nm.
Radial distances must remain 0.1 nm below half the shortest box length.

The existing `umbrella.analysis` histogram, overlap, inefficiency, WHAM, RMSE and
circular-block-bootstrap functions are reused by `umbrella.atomistic`.
Histograms have 0.01 nm bins; WHAM uses effective counts and the actual harmonic
bias. Unsupported bins stay NaN. The assessment interval remains fixed between
the first and last centers. Full, half, quarter and bootstrap solves must all
satisfy the SCOPE gates; a numerical solution alone is insufficient. The radial
PMF includes radial entropy; no radial Jacobian or standard-state correction is
claimed. A failed reconstruction is named `UNACCEPTED_diagnostic.csv` and is
never shown as the EGFR result.

To reproduce on a compatible ROCm/HIP host, create the environment in
`environment-atomistic.yml`, set `OPENMM_PYTHON` and `AMBERHOME` to that environment,
and set `ANALYSIS_PYTHON` to a Python with this project installed. Use a clean run
directory (the runner refuses to overwrite a simulation record):

```bash
nohup systemd-inhibit --what=sleep:idle --why='P23 atomistic pilot' \
  bash scripts/run_atomistic.sh > atomistic-run.log 2>&1 < /dev/null &
```

The runner limits preparation plus MD to 9900 seconds, reserving time for analysis
within the overall 10800-second task budget. `simulation.json` records the actual
platform and device. This recorded attempt used existing separate `md_env` and
`openmm_env` environments; `preparation_audit.json` stores installation records,
while `simulation.json` supplies the live OpenMM version (8.6). A stale conda
OpenMM installation record (8.5.2) is not the live imported version. HIP support
was verified by creation and execution of the full solvated simulation context.

References: [RCSB 1M17](https://www.rcsb.org/structure/1M17),
[PDBFixer manual](https://github.com/openmm/pdbfixer/blob/master/Manual.html),
[OpenMM simulation guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html),
[OpenMM HIP installation](https://docs.openmm.org/latest/userguide/application/01_getting_started.html).
