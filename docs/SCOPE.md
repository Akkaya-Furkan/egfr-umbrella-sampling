# Same-day scope — 2026-10-24

## Question and decision
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

## Explicitly deferred
Atomistic ligand parameterization, protein repair/protonation, explicit solvent,
equilibration, pulling, EGFR umbrella production, EGFR histogram overlap and
convergence, EGFR PMF and standard-state binding free energy: **not done**.
A straight geometric translation is not evidence of a viable ligand exit route.
No EGFR thermodynamic conclusion is permitted by this scope.
