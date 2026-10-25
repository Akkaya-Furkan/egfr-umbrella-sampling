"""Fail-closed atomistic pilot diagnostics, reusing the model's WHAM primitives."""
from itertools import combinations
import numpy as np
from .analysis import R, align, histograms, inefficiency, overlap, resample_blocks, rmse, wham


def assess(series, centers, *, spring=1000., temperature=300., interval_ps=.1,
           burn_ps=20., starts=None, complete=False, bootstrap=100, seed=23024,
           expected_windows=17):
    """Return JSON-safe evidence plus an explicitly conditional diagnostic curve.

    A solver solution alone never implies physical convergence. Missing/invalid
    series fail closed; all statistical gates apply to the fixed center interval.
    """
    centers = np.asarray(centers, float)
    report = dict(status='PMF yakınsamadı/kanıtlanmadı', accepted=False,
                  expected_windows=expected_windows, actual_windows=len(series), failures=[], gates={})
    gates = report['gates']
    gates['complete_production'] = bool(complete and len(series)==len(centers)==expected_windows)
    gates['burn_in'] = bool(burn_ps >= 20)
    gates['relaxed_nonclashing_starts'] = bool(starts is not None and len(starts)==expected_windows and
                                              all(np.isfinite(s['min_heavy_nm']) and s['min_heavy_nm']>=.15 for s in starts))
    if len(series) != len(centers) or not series:
        report['failures'] = [k for k,v in gates.items() if not v]
        report['analysis_note'] = f'No complete {expected_windows}-window ensemble; WHAM and convergence are unassessed.'
        return report, None
    xs = [np.asarray(x,float) for x in series]
    gates['finite_samples'] = bool(all(x.ndim==1 and len(x)>=4 and np.isfinite(x).all() for x in xs))
    if not gates['finite_samples']:
        report['failures'] = [k for k,v in gates.items() if not v]
        return report, None
    report['samples_per_window'] = [len(x) for x in xs]
    report['production_ps_per_window'] = [len(x)*interval_ps for x in xs]
    gates['production_duration'] = bool(min(report['production_ps_per_window']) >=100)
    edges = np.arange(np.floor(min(x.min() for x in xs)*100)/100-.01,
                      np.ceil(max(x.max() for x in xs)*100)/100+.0201,.01)
    bins = (edges[1:]+edges[:-1])/2
    counts = histograms(xs,edges)
    ov = overlap(counts)
    report['adjacent_overlap'] = ov.tolist()
    report['minimum_overlap'] = float(ov.min())
    gates['adjacent_overlap'] = bool(ov.min()>=.10)
    gs = []
    for x in xs:
        try: gs.append(inefficiency(x))
        except ValueError: gs.append(float(len(x)))
    effective = np.array([len(x) for x in xs])/gs
    report['inefficiencies'] = gs
    report['effective_samples'] = effective.tolist()
    gates['effective_samples'] = bool(effective.min()>=100)
    block = int(np.ceil(5*max(gs)))
    report['bootstrap_block_samples'] = block
    report['complete_blocks_per_window'] = [len(x)//block for x in xs]
    gates['ten_blocks'] = bool(min(report['complete_blocks_per_window'])>=10)
    bias = .5*spring*(bins[None,:]-centers[:,None])**2/(R*temperature)
    mask = (bins>=centers[0]) & (bins<=centers[-1])
    report['assessment_interval_nm'] = [float(centers[0]),float(centers[-1])]
    report['assessment_bins'] = int(mask.sum())
    ref = int(np.argmin(abs(bins-centers[0])))
    solves = []
    report['reconstruction_support'] = []

    def solve(parts):
        hs = histograms(parts,edges)/np.asarray(gs)[:,None]
        f, diag = wham(hs,bias)
        solves.append(diag)
        report['reconstruction_support'].append(dict(solve=len(solves), unsupported_assessment_bins=int((~np.isfinite(f[mask])).sum()), reference_bin_sampled=bool(np.isfinite(f[ref]))))
        return align(f*R*temperature,ref)

    curve = None
    try:
        full = solve(xs)
        curve = np.column_stack((bins,full))
        halves = [solve([np.array_split(x,2)[h] for x in xs]) for h in range(2)]
        quarters = [solve([np.array_split(x,4)[q] for x in xs]) for q in range(4)]
        gates['fixed_interval_support'] = bool(mask.any() and all(np.isfinite(f[mask]).all() for f in [full]+halves+quarters))
        if gates['fixed_interval_support']:
            half_error = rmse(*halves,mask)
            quarter_errors = [rmse(a,b,mask) for a,b in combinations(quarters,2)]
            report['half_rmse_kj_mol'] = half_error
            report['quarter_pair_rmse_kj_mol'] = quarter_errors
            gates['half_agreement'] = bool(half_error<=2)
            gates['quarter_agreement'] = bool(max(quarter_errors)<=2)
        else:
            gates['half_agreement'] = gates['quarter_agreement'] = False
        boot = []
        if block <= min(len(x) for x in xs):
            rng = np.random.default_rng(seed)
            for _ in range(bootstrap):
                boot.append(solve([resample_blocks(x,block,rng) for x in xs]))
        report['bootstrap_successes'] = len(boot)
        gates['bootstrap'] = bool(len(boot)>=100 and np.isfinite(np.std(np.asarray(boot)[:,mask],axis=0,ddof=1)).all()) if boot else False
        gates['solver'] = bool(all(s['residual_kbt']<1e-9 for s in solves))
    except (ValueError,RuntimeError) as exc:
        report['solver_error'] = str(exc)
        gates['solver'] = False
    report['solves_completed'] = len(solves)
    if solves: report['maximum_solver_residual_kbt'] = max(s['residual_kbt'] for s in solves)
    # Required gates cannot be skipped by an early exception.
    for required in ('fixed_interval_support','half_agreement','quarter_agreement','bootstrap','solver'):
        gates.setdefault(required,False)
    report['failures'] = [k for k,v in gates.items() if not v]
    report['accepted'] = bool(all(gates.values()))
    if report['accepted']: report['status'] = 'Conditional radial pilot PMF accepted by prespecified numerical gates'
    return report,curve
