#!/usr/bin/env python3
"""Reconstruct PMF only from simulated samples, with independent diagnostics."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from umbrella.analysis import R, potential, histograms, overlap, wham, inefficiency, align, rmse, resample_blocks


def main():
    out = Path('results/model')
    meta = json.loads((out/'simulation.json').read_text())
    x = np.load(out/'samples.npz')['x_nm']
    c = meta['config']
    kbt = R*c['temperature_K']
    centers = np.array(c['centers_nm'])
    edges = np.linspace(-1, 1, 101)
    bins = (edges[1:]+edges[:-1])/2
    ref = int(np.argmin(abs(bins+0.5)))
    mask = abs(bins) <= 0.65
    bias = 0.5*c['spring_kj_mol_nm2']*(bins[None, :]-centers[:, None])**2/kbt
    # Guard against silent removal of observations beyond the histogram domain.
    if x.shape != (2, 17, 10000) or not np.all(np.isfinite(x)) or x.min() < edges[0] or x.max() > edges[-1]:
        raise ValueError('Incomplete, nonfinite or out-of-range production series')
    g = np.array([[max(inefficiency(s), inefficiency(potential(s))) for s in replica] for replica in x])
    counts = np.array([histograms(replica, edges) for replica in x])
    solves = []

    def reconstruct(series, factors):
        h = np.array([histograms(replica, edges) for replica in series]) / factors[:, :, None]
        f, info = wham(h.reshape(-1, len(bins)), np.tile(bias, (len(series), 1)))
        solves.append(info)
        return align(f*kbt, ref)

    pmf = reconstruct(x, g)
    first = reconstruct(x[:, :, :5000], g)
    second = reconstruct(x[:, :, 5000:], g)
    replicas = [reconstruct(x[r:r+1], g[r:r+1]) for r in range(2)]
    block = int(np.ceil(5*g.max()))
    rng = np.random.default_rng(230923)
    boot = []
    for _ in range(100):
        resampled = np.array([[resample_blocks(s, block, rng) for s in rep] for rep in x])
        boot.append(reconstruct(resampled, g))
    boot = np.array(boot)
    # Unsupported bins stay NaN; standard deviations computed only where all
    # bootstrap reconstructions have support (no nan-skipping optimism).
    std = np.std(boot, axis=0, ddof=1)
    exact = align(potential(bins), ref)
    adjacent = overlap(counts.sum(axis=0))
    neff = x.shape[-1]/g
    metrics = dict(expected_series=34, actual_series=int(np.prod(x.shape[:2])),
                   expected_samples_per_series=10000, actual_samples=int(x.size),
                   min_samples_per_series=int(x.shape[-1]),
                   min_adjacent_overlap=float(adjacent.min()),
                   exact_pmf_rmse_kj_mol=rmse(pmf, exact, mask),
                   split_half_rmse_kj_mol=rmse(first, second, mask),
                   replica_rmse_kj_mol=rmse(*replicas, mask),
                   min_effective_samples=float(neff.min()), max_inefficiency=float(g.max()),
                   bootstrap_block_samples=block, bootstrap_successful=len(boot),
                   max_wham_residual_kbt=max(s['residual_kbt'] for s in solves),
                   wham_solves=len(solves), assessment_bins=int(mask.sum()),
                   bootstrap_std_finite=bool(np.all(np.isfinite(std[mask]))),
                   simulation_elapsed_seconds=meta['elapsed_seconds'])
    gates = dict(complete_series=metrics['actual_series']==34 and x.shape[-1]==10000,
                 wham_converged=metrics['max_wham_residual_kbt'] < 1e-9,
                 overlap=metrics['min_adjacent_overlap']>=0.10,
                 exact_accuracy=metrics['exact_pmf_rmse_kj_mol']<=1,
                 split_consistency=metrics['split_half_rmse_kj_mol']<=1,
                 replica_consistency=metrics['replica_rmse_kj_mol']<=1,
                 effective_samples=metrics['min_effective_samples']>=100,
                 block_length=block>=5*g.max(), bootstrap=len(boot)>=100 and metrics['bootstrap_std_finite'])
    report = dict(system='one-particle double-well model, NOT EGFR', metrics=metrics,
                  gates={k: bool(v) for k, v in gates.items()}, model_gates_passed=all(gates.values()),
                  adjacent_overlap=adjacent.tolist(), inefficiencies=g.tolist(),
                  wham_solvers=solves)
    (out/'analysis.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(out/'histograms.npz', counts=counts, edges_nm=edges, inefficiencies=g,
                        bootstrap_pmf_kj_mol=boot)
    table = np.column_stack([bins, pmf, std, exact, first, second, *replicas])
    np.savetxt(out/'pmf.csv', table, delimiter=',', header='x_nm,pmf_kj_mol,bootstrap_sd_kj_mol,exact_kj_mol,first_half,second_half,replica_0,replica_1', comments='')
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bins, exact, 'k--', label='Exact model potential')
    ax.plot(bins, pmf, label='WHAM from OpenMM windows')
    ax.fill_between(bins, pmf-std, pmf+std, alpha=0.25, label='Pointwise bootstrap ±1 SD')
    ax.axvspan(-0.65, 0.65, color='gray', alpha=0.08, label='Acceptance interval')
    ax.set(xlabel='Model coordinate x (nm)', ylabel='Relative PMF (kJ/mol)', ylim=(-2, 18), title='Double-well method benchmark — NOT an EGFR PMF')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out/'pmf.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7))
    for w, h in enumerate(counts.sum(axis=0)):
        axes[0].plot(bins, h/h.sum(), label=f'{centers[w]:.1f}')
    axes[0].set(xlabel='Model x (nm)', ylabel='Probability per bin', title='Actual biased-window histograms (replicas pooled)')
    axes[0].legend(ncol=6, fontsize=7)
    axes[1].plot(np.arange(len(adjacent)), adjacent, 'o-')
    axes[1].axhline(0.10, color='r', linestyle='--', label='Acceptance threshold')
    axes[1].set(xlabel='Adjacent window pair index', ylabel='Histogram intersection', ylim=(0,1))
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(out/'window_histograms.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5))
    for f, label in [(first, 'First half'), (second, 'Second half'), (replicas[0], 'Replica 0'), (replicas[1], 'Replica 1')]:
        ax.plot(bins[mask], f[mask], label=label)
    ax.set(xlabel='Model x (nm)', ylabel='Relative PMF (kJ/mol)', title='Model consistency diagnostics (not proof for EGFR)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/'convergence.png', dpi=160)
    plt.close(fig)
    print(json.dumps(dict(metrics=metrics, gates={k: bool(v) for k, v in gates.items()}), indent=2))


if __name__ == '__main__':
    main()
