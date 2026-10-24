#!/usr/bin/env python3
"""Re-read real saved samples; never publish an unaccepted reconstruction as PMF."""
import hashlib
import json
from pathlib import Path
import numpy as np
from umbrella.atomistic import assess
from umbrella.analysis import histograms, overlap


def main():
    out = Path('results/atomistic')
    sim_path = out/'simulation.json'
    sim = json.loads(sim_path.read_text()) if sim_path.exists() else {}
    config = sim.get('config',{})
    centers = config.get('centers_nm',list(range(17)))
    expected = [out/f'samples_{w:02d}.npz' for w in range(17)]
    existing = [p for p in expected if p.is_file()]
    windows = json.loads((out/'windows.json').read_text()) if (out/'windows.json').exists() else []
    report,curve = assess([np.load(p)['r_nm'] for p in existing],centers,
                           spring=config.get('spring_kj_mol_nm2',1000),
                           temperature=config.get('temperature_K',300),
                           interval_ps=config.get('sample_interval_ps',.1),
                           burn_ps=config.get('burn_ps',0),starts=windows,
                           complete=sim.get('status')=='production_complete')
    report.update(simulation_status=sim.get('status','not started'),
                  completed_window_ids=[int(p.stem.split('_')[-1]) for p in existing],
                  expected_series=17, verified_series=len(existing),
                  missing_window_ids=[i for i,p in enumerate(expected) if not p.exists()],
                  simulation_error=sim.get('error'),
                  simulation_elapsed_seconds=sim.get('elapsed_seconds'),
                  input_series={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in existing})
    report['partial_windows'] = {p.name: dict(samples=int(np.load(p).size), duration_ps=float(np.load(p).size*config.get('sample_interval_ps',.1))) for p in out.glob('partial_*.npy')}
    box = sim.get('box_nm')
    if box and existing:
        max_r = max(float(np.load(p)['r_nm'].max()) for p in existing)
        report['maximum_sampled_com_nm'] = max_r
        report['periodic_image_margin_nm'] = min(box)/2-max_r
        safe = report['periodic_image_margin_nm']>.1
        report['gates']['periodic_image_safety'] = safe
        if not safe:
            report['accepted'] = False
            report['status'] = 'PMF yakınsamadı/kanıtlanmadı'
            report['failures'].append('periodic_image_safety')
    for stale in ('pmf.csv','UNACCEPTED_diagnostic.csv'):
        (out/stale).unlink(missing_ok=True)
    if curve is not None:
        name = 'pmf.csv' if report['accepted'] else 'UNACCEPTED_diagnostic.csv'
        np.savetxt(out/name,curve,delimiter=',',header='com_distance_nm,conditional_radial_free_energy_kj_mol',comments='')
        report['curve_artifact'] = name
    if len(existing)>=2:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        series=[np.load(p)['r_nm'] for p in existing]
        edges=np.arange(np.floor(min(x.min() for x in series)*100)/100-.01,
                        np.ceil(max(x.max() for x in series)*100)/100+.0201,.01)
        counts=histograms(series,edges)
        np.savez_compressed(out/'histograms.npz',edges_nm=edges,counts=counts,window_ids=report['completed_window_ids'])
        fig,axes=plt.subplots(2,1,figsize=(10,7),layout='constrained')
        bins=(edges[:-1]+edges[1:])/2
        for w,h in zip(report['completed_window_ids'],counts):
            axes[0].plot(bins,h/h.sum(),label=str(w),linewidth=1)
        axes[0].set(xlabel='Protein–erlotinib COM distance (nm)',ylabel='Biased probability per bin',
                    title='EGFR pilot: biased histograms (not a PMF)')
        axes[0].legend(ncol=9,fontsize=7,title='Window ID')
        ov=overlap(counts)
        axes[1].bar(np.arange(len(ov)),ov)
        axes[1].axhline(.10,color='red',linestyle='--',label='Acceptance >=0.10')
        axes[1].set(xlabel='Adjacent pair (i, i+1)',ylabel='Histogram intersection',ylim=(0,1))
        axes[1].legend()
        fig.savefig(out/'window_histograms.png',dpi=160)
        plt.close(fig)
    (out/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))

if __name__ == '__main__': main()
