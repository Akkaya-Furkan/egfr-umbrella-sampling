import numpy as np
from umbrella.atomistic import assess


def test_missing_ensemble_cannot_be_reported_as_pmf():
    report, curve = assess([],np.arange(17),complete=False)
    assert not report['accepted']
    assert report['actual_windows']==0
    assert curve is None
    assert report['status']=='PMF yakınsamadı/kanıtlanmadı'


def test_numerical_wham_solution_does_not_override_short_duration_or_no_overlap():
    rng = np.random.default_rng(23)
    centers = np.linspace(1,2.6,17)
    xs = [rng.normal(c,.001,40) for c in centers]
    report,_ = assess(xs,centers,complete=True,bootstrap=0,
                      starts=[dict(min_heavy_nm=.25) for _ in centers])
    assert not report['accepted']
    assert not report['gates']['production_duration']
    assert not report['gates']['adjacent_overlap']
    assert not report['gates']['bootstrap']


def test_nonfinite_sample_is_rejected_before_histogramming():
    xs = [np.array([1.,2.,3.,np.nan]) for _ in range(17)]
    report,curve = assess(xs,np.linspace(1,2.6,17),complete=True)
    assert not report['gates']['finite_samples']
    assert not report['accepted']
    assert curve is None


def test_recorded_run_is_not_accepted():
    import json, pathlib
    path = pathlib.Path(__file__).parents[1]/'results/atomistic/analysis.json'
    report = json.loads(path.read_text())
    assert report['actual_windows']==17
    assert not report['accepted']
    assert not report['gates']['adjacent_overlap'] and not report['gates']['effective_samples']
    assert report['curve_artifact']=='UNACCEPTED_diagnostic.csv'
