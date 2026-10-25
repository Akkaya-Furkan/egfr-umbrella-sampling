# Same-day scope — 2026-10-24

## Original benchmark question and decision
What does the free-energy profile along an EGFR ligand unbinding path look like?
An atomistic, solvated EGFR unbinding PMF cannot be established within this
same-day scope. First validate the umbrella/WHAM method on a tractable model;
apply only the coordinate/window preparation to public EGFR–erlotinib structure
[1M17](https://www.rcsb.org/structure/1M17), ligand AQ4. No thesis ligand is used.
The model PMF must never be described as an EGFR PMF.

## Executable experiment (specified before production)
OpenMM CPU/Reference Langevin dynamics for one particle with Cartesian coordinate
x in nm and separable harmonic y/z confinement. Base potential:
U(x) = 8 * ((x / 0.5)^2 - 1)^2 kJ/mol.
Temperature 300 K, mass 10 Da, friction 10/ps, timestep 0.002 ps.
Use 17 harmonic umbrellas centered from -0.8 to +0.8 nm, spring constant
100 kJ/mol/nm². Two independent seeded replicas per window, each with 100 ps
burn-in followed by 1000 ps production; save x every 0.1 ps. Run under
`systemd-inhibit --what=sleep:idle`. No large molecular trajectories are needed.
Estimate a binned WHAM PMF, adjacent histogram overlap, autocorrelation,
block-bootstrap uncertainty, replica and split-half consistency. Compare with
the exact model potential only after estimation, never use it to construct PMF.
Report WHAM on [-1,1] nm, assess accuracy on [-0.65,0.65] nm.

## Numerical acceptance gates
All gates must pass to claim the model workflow is validated:
- 34/34 complete production series, 10000 finite samples each.
- All WHAM solves converge to free-offset change < 1e-9 kBT.
- Every adjacent pooled-window histogram intersection >= 0.10.
- Model PMF RMSE <= 1.0 kJ/mol versus the exact potential, with both curves
  zeroed at the bin nearest -0.5 nm, on the assessment interval.
- First-half versus second-half PMF RMSE <= 1.0 kJ/mol on that interval.
- Independent-replica PMF RMSE <= 1.0 kJ/mol on that interval.
- At least 100 estimated effective samples per replica/window; bootstrap block
  length at least 5 times the largest estimated statistical inefficiency.
- At least 100 successful block-bootstrap reconstructions; finite pointwise
  standard deviations on the assessment interval. These are uncertainty
  estimates, not proof of convergence of a molecular system.
- Public PDB download has a SHA256 manifest; extract protein and AQ4; generate
  17 translated ligand starting structures with measured projected displacement
  error < 1e-4 nm (PDB rounding). This is structural preparation, not MD or a validated exit path.
- Offline analysis tests pass in pytest and GitHub Actions (no MD in CI).

## Originally deferred (continuation below supersedes stage status)
Atomistic ligand parameterization, protein repair/protonation, explicit solvent,
equilibration, pulling, EGFR umbrella production, EGFR histogram overlap and
convergence, EGFR PMF and standard-state binding free energy: **not done**.
A straight geometric translation is not evidence of a viable ligand exit route.
No EGFR thermodynamic conclusion is permitted by this scope.

## Atomistic continuation — gates fixed before production (2026-10-24)

Attempt a bounded pilot with public 1M17/AQ4, not a thesis ligand. The previous
model validation remains separate. Repair internal missing residues and sidechains
with PDBFixer, omit unresolved terminal extensions, select default protonation at
pH 7.4. Parameterize neutral erlotinib with GAFF2/AM1-BCC (AmberTools), protein
ff14SB, TIP3P explicit water and neutralizing ions. Modeled loops and protonation
are assumptions, not validated structural observations.

Minimize, equilibrate at 300 K/1 bar, then use a moving harmonic restraint on the
mass-weighted protein-heavy/ligand-heavy COM distance to generate 17 starting
states at 0.1 nm increments over a 1.6 nm radial range. Physical steered MD replaces
the earlier rigid translations. A weak protein CA positional restraint defines
a conditional, restrained-protein pilot. Record its strength and all run settings.
Use PME/H-bond constraints, 2 fs Langevin dynamics. Prefer HIP GPU, test the actual
platform. Neutralizing ions only; no added salt. Save scalar coordinate time series
and small diagnostics in git; topology, restart coordinates and logs stay ignored.

Before any atomistic PMF may be called a result, require **all** of:
- 17/17 finite production series, each at least 100 ps, following at least 20 ps
  window burn-in. No severe protein–ligand heavy-atom clash (<0.15 nm) in any
  relaxed starting state. Check radial periodic-image safety during sampling.
- Every adjacent histogram intersection >=0.10 with fixed 0.01 nm bins.
- At least 100 estimated effective coordinate samples in every window.
- All full, half, quarter and bootstrap WHAM solves converge below 1e-9 kBT.
- All bins between the first and last umbrella centers supported in every
  full/half/quarter reconstruction (no interpolation or shrinking the interval).
- First/second-half and every pair of quarter PMFs have RMSE <=2 kJ/mol after
  alignment at the same fixed bin in that interval.
- 100/100 circular block-bootstrap reconstructions; block length >=5 times the
  largest estimated coordinate inefficiency, and no longer than any series.
- At least 10 complete blocks per production series; finite bootstrap uncertainty.

These are necessary pilot gates, not proof that hidden protein/ligand modes or
multiple routes equilibrated. A radial PMF includes radial entropy; do not label
it a standard-state binding free energy. Failure of any gate must be reported as
**PMF yakınsamadı/kanıtlanmadı**, with actual completed windows and failure reasons.
An unconverged numerical reconstruction may only be retained as explicitly
unaccepted diagnostic data, never as the EGFR result. Bound the entire attempt
to <=10800 wall seconds; use nohup, timeout and sleep/idle inhibition for MD.

## Round 2 (2026-10-25) — no gate changed
Thresholds above were not altered. Round 2 added one window (18 total, centers
unchanged otherwise, k=1000 kJ/mol/nm²) and 500 ps/window re-sampling, because
round 1 showed one overlap failure (pair 10–11) and autocorrelation too long for
400 ps. `assess` now takes `expected_windows` (=number of planned windows)
instead of the literal 17; the checks themselves are the same. Outcome: overlap
gate passes, effective-sample, ten-block, solver, support, half/quarter and
bootstrap gates fail → **PMF yakınsamadı/kanıtlanmadı**.
