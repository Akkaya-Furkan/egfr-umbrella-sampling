# Methods and interpretation

## Sampling
The executed model is a single 10 Da particle in a separable external potential.
Cartesian x is the reaction coordinate. Harmonic y/z confinement integrates to
an x-independent factor; no radial Jacobian or standard-state correction is
appropriate. The exact x marginal PMF equals the prescribed double-well
potential plus a constant. This makes a falsifiable method benchmark, not a
molecular ligand model. OpenMM Reference runs on the CPU; the available AMD GPU
is not used. Random seeds, integration settings and measured runtime are in
`results/model/simulation.json`. Production excludes the fixed burn-in.

## WHAM
Use equal-width histograms and harmonic bias energies evaluated at bin centers.
Solve the self-consistent binned WHAM equations in log space. The unknown common
base energy cancels from the reweighting equations. The analytical potential is
used for the accuracy check only, never as a prior or as substituted samples.
An empty bin is NaN; missing assessment bins cause analysis failure.

Estimate statistical inefficiency from both x and base potential energy using
FFT autocorrelations truncated at the first nonpositive consecutive-lag pair;
use the larger estimate for each replica/window. Divide its histogram and
sample count by this estimate before WHAM. These are estimated effective
counts, not a guarantee that every slow mode is detected.

Bootstrap circular contiguous blocks separately in each replica/window,
with block length rounded up from five times the largest inefficiency.
Reuse the original inefficiency estimates for bootstrap weights. Use a fixed
bootstrap seed, save all bootstrap PMFs, and report pointwise one-standard-
deviation uncertainty after alignment at the same reference bin. This is not
a simultaneous confidence band; finite sampling and bin-width bias remain.

Adjacent overlap is the sum over bins of the minimum of two normalized
histograms, pooling replicas of each center. Compare independently estimated
first/second-half and replica PMFs on the prespecified assessment interval.
These consistency checks and exact-model accuracy validate this benchmark
only. They do not establish convergence for protein–ligand simulations.

## EGFR coordinate preparation
Download the deposited PDB, check content, and record SHA256. Select chain A
protein ATOM records (blank or A alternate location) and AQ4 HETATM records.
Use the vector from the protein heavy-atom geometric centroid to the ligand
heavy-atom geometric centroid. Translate the intact ligand along this fixed
vector; retain the protein coordinates. Report projected translations after
rereading PDB files, and minimum protein–ligand heavy-atom distances.

These coordinates can contain clashes and are **not production-ready windows**.
The direction has not been shown to pass through a physical exit channel.
There is no force field, solvent, pulling trajectory or EGFR energy evaluation.
Before an EGFR PMF, repair/protonate the protein, parameterize erlotinib, solvate
and equilibrate, establish an exit route, run and extend atomistic windows,
and assess overlap and convergence for that molecular system.

## Sources
- [RCSB 1M17: EGFR with erlotinib](https://www.rcsb.org/structure/1M17)
- [OpenMM custom external forces](https://docs.openmm.org/7.4.0/userguide/application.html#using-a-customexternalforce)
- [Kumar et al., WHAM (1992)](https://doi.org/10.1002/jcc.540130812)
- [PyMBAR timeseries documentation: correlation and effective samples](https://pymbar.readthedocs.io/en/latest/timeseries.html)
  (methodological reference; this repository implements binned WHAM, not MBAR).
