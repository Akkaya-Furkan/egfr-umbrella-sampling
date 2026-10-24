"""Binned WHAM in reduced units; empty bins remain NaN (never interpolated)."""
import numpy as np
from scipy.special import logsumexp

R = 0.00831446261815324  # kJ mol^-1 K^-1


def potential(x):
    return 8.0 * ((np.asarray(x) / 0.5) ** 2 - 1.0) ** 2


def histograms(series, edges):
    return np.array([np.histogram(x, edges)[0] for x in series], dtype=float)


def overlap(counts):
    h = np.asarray(counts, dtype=float)
    if h.ndim != 2 or np.any(h < 0) or np.any(h.sum(axis=1) == 0):
        raise ValueError("Nonempty nonnegative window histograms required")
    h = h / h.sum(axis=1, keepdims=True)
    return np.minimum(h[:-1], h[1:]).sum(axis=1)


def wham(counts, bias_kbt, tolerance=1e-9, max_iterations=50000):
    """Return PMF/kBT and solver diagnostics. Counts may be effective counts.

    p_b ∝ sum_k H_kb / sum_k N_k exp(f_k - u_kb),
    f_k = -log sum_b p_b exp(-u_kb). Equal-width bins only.
    """
    h, u = np.asarray(counts, float), np.asarray(bias_kbt, float)
    if h.ndim != 2 or h.shape != u.shape or not np.all(np.isfinite(u)):
        raise ValueError("Matching finite window-by-bin arrays required")
    if not np.all(np.isfinite(h)) or np.any(h < 0) or np.any(h.sum(axis=1) <= 0):
        raise ValueError("Each window must contain finite nonnegative counts")
    occupied = h.sum(axis=0) > 0
    log_h = np.log(h.sum(axis=0)[occupied])
    log_n = np.log(h.sum(axis=1))
    f = np.zeros(len(h))
    for iteration in range(1, max_iterations + 1):
        log_p = log_h - logsumexp(log_n[:, None] + f[:, None] - u[:, occupied], axis=0)
        log_p -= logsumexp(log_p)
        new = -logsumexp(log_p[None, :] - u[:, occupied], axis=1)
        new -= new[0]
        residual = float(np.max(np.abs(new - f)))
        f = new
        if residual < tolerance:
            # Recompute probability using converged offsets.
            log_p = log_h - logsumexp(log_n[:, None] + f[:, None] - u[:, occupied], axis=0)
            log_p -= logsumexp(log_p)
            pmf = np.full(h.shape[1], np.nan)
            pmf[occupied] = -log_p
            return pmf, {"iterations": iteration, "residual_kbt": residual}
    raise RuntimeError(f"WHAM failed to converge: residual={residual}")


def inefficiency(x):
    """Initial-positive-sequence FFT autocorrelation estimate (diagnostic)."""
    x = np.asarray(x, float)
    if len(x) < 4 or not np.all(np.isfinite(x)) or np.var(x) == 0:
        raise ValueError("At least four finite, varying observations required")
    x = x - x.mean()
    n = len(x)
    ft = np.fft.rfft(x, n=2 * n)
    ac = np.fft.irfft(ft * ft.conj())[:n] / np.arange(n, 0, -1)
    ac /= ac[0]
    total = 0.0
    for i in range(1, n - 1, 2):
        pair = ac[i] + ac[i + 1]
        if pair <= 0:
            break
        total += pair
    return max(1.0, 1.0 + 2.0 * total)


def align(pmf, reference_bin):
    p = np.asarray(pmf, float)
    if not np.isfinite(p[reference_bin]):
        raise ValueError("Reference bin not sampled")
    return p - p[reference_bin]


def rmse(a, b, mask):
    a, b = np.asarray(a)[mask], np.asarray(b)[mask]
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Assessment interval contains unsampled bins")
    return float(np.sqrt(np.mean((a - b) ** 2)))


def resample_blocks(series, block_length, rng):
    """Circular moving-block bootstrap, separately within each production series."""
    x = np.asarray(series)
    if not 1 <= block_length <= len(x):
        raise ValueError("Invalid block length")
    starts = rng.integers(len(x), size=int(np.ceil(len(x) / block_length)))
    indices = (starts[:, None] + np.arange(block_length)) % len(x)
    return x[indices.ravel()[:len(x)]]
