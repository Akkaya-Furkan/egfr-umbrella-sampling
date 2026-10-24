from pathlib import Path
import json
import numpy as np
from umbrella.analysis import R, align, histograms, wham

ROOT = Path(__file__).resolve().parents[1]


def test_saved_pmf_is_reconstruction_of_actual_coordinate_histograms():
    out = ROOT/'results/model'
    coordinates = np.load(out/'samples.npz')['x_nm']
    stored = np.load(out/'histograms.npz')
    meta = json.loads((out/'simulation.json').read_text())['config']
    counts = np.array([histograms(replica, stored['edges_nm']) for replica in coordinates])
    np.testing.assert_array_equal(counts, stored['counts'])
    assert np.isfinite(coordinates).all()
    assert counts.sum() == coordinates.size
    centers = np.array(meta['centers_nm'])
    bins = (stored['edges_nm'][1:]+stored['edges_nm'][:-1])/2
    bias = .5*meta['spring_kj_mol_nm2']*(bins[None, :]-centers[:, None])**2/(R*meta['temperature_K'])
    weighted = counts/stored['inefficiencies'][:, :, None]
    f, _ = wham(weighted.reshape(-1, len(bins)), np.tile(bias, (len(coordinates), 1)))
    f = align(f*R*meta['temperature_K'], np.argmin(abs(bins+.5)))
    saved = np.genfromtxt(out/'pmf.csv', delimiter=',', names=True)
    np.testing.assert_allclose(f, saved['pmf_kj_mol'], atol=1e-10, equal_nan=True)
