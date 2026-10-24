import numpy as np
import pytest
from umbrella.analysis import wham, overlap, inefficiency, align, rmse, resample_blocks


def test_wham_recovers_known_distribution_with_unequal_window_counts():
    x = np.linspace(-2, 2, 201)
    target = np.exp(-0.5*x*x)
    centers = np.linspace(-1.5, 1.5, 7)
    bias = 4*(x[None, :]-centers[:, None])**2
    distributions = target[None, :]*np.exp(-bias)
    distributions /= distributions.sum(axis=1, keepdims=True)
    counts = distributions*np.arange(1, 8)[:, None]*10000
    pmf, info = wham(counts, bias)
    np.testing.assert_allclose(align(pmf, 100), 0.5*x*x, atol=1e-7)
    assert info['residual_kbt'] < 1e-9


def test_wham_invariant_to_window_energy_constants():
    counts = np.array([[8, 4, 1], [1, 4, 8.]])
    bias = np.array([[0, 1, 2.], [2, 1, 0.]])
    a, _ = wham(counts, bias)
    b, _ = wham(counts, bias + np.array([[20], [-30]]))
    np.testing.assert_allclose(align(a, 1), align(b, 1), atol=1e-8)


def test_empty_bins_remain_unknown():
    pmf, _ = wham([[10, 0, 10]], [[0, 0, 0]])
    assert np.isnan(pmf[1])
    with pytest.raises(ValueError, match='Reference'):
        align(pmf, 1)


@pytest.mark.parametrize('counts,bias', [([[0, 0]], [[0, 0]]),
    ([[-1, 2]], [[0, 0]]), ([[1, np.nan]], [[0, 0]]),
    ([[1, 2]], [[0, np.inf]]), ([[1, 2]], [[0]])])
def test_invalid_wham_input(counts, bias):
    with pytest.raises(ValueError):
        wham(counts, bias)


def test_nonconvergence_raises():
    with pytest.raises(RuntimeError):
        wham([[10, 1], [1, 2]], [[0, 2], [2, 0]], max_iterations=1)


def test_overlap_extremes_and_partial_overlap():
    np.testing.assert_allclose(overlap([[10, 0], [5, 5], [0, 10], [0, 10]]), [.5, .5, 1])
    assert overlap([[1, 0], [0, 1]])[0] == 0
    with pytest.raises(ValueError):
        overlap([[0, 0], [1, 2]])


def test_inefficiency_detects_correlated_series():
    rng = np.random.default_rng(42)
    noise = rng.normal(size=50000)
    correlated = np.empty(len(noise))
    correlated[0] = noise[0]
    for i in range(1, len(noise)):
        correlated[i] = .9*correlated[i-1]+noise[i]
    assert 1 <= inefficiency(noise) < 1.2
    assert 14 < inefficiency(correlated) < 26  # AR(1): theoretical g = 19
    with pytest.raises(ValueError):
        inefficiency(np.ones(20))


def test_circular_block_bootstrap_and_seed():
    x = np.arange(20)
    a = resample_blocks(x, 5, np.random.default_rng(42))
    b = resample_blocks(x, 5, np.random.default_rng(42))
    np.testing.assert_array_equal(a, b)
    assert len(a) == len(x)
    for block in a.reshape(-1, 5):
        np.testing.assert_array_equal(np.diff(block) % len(x), np.ones(4))
    with pytest.raises(ValueError):
        resample_blocks(x, 21, np.random.default_rng(42))


def test_rmse_rejects_missing_assessment_bins():
    assert rmse([0, 2], [0, 0], [True, True]) == pytest.approx(np.sqrt(2))
    with pytest.raises(ValueError):
        rmse([np.nan], [0], [True])
