#!/usr/bin/env python3
"""Read back deliverables, hashes and numerical gates; fail closed."""
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    result = Path('results')
    model = json.loads((result/'model/analysis.json').read_text())
    egfr = json.loads((result/'egfr/preparation.json').read_text())
    expected = [result/p for p in [
        'data_manifest.json', 'environment.json', 'inhibitors.txt', 'model_run.log', 'test-results.xml',
        'analysis_run.log', 'model/samples.npz', 'model/simulation.json',
        'model/analysis.json', 'model/histograms.npz', 'model/pmf.csv',
        'model/pmf.png', 'model/window_histograms.png', 'model/convergence.png',
        'egfr/protein.pdb', 'egfr/ligand_AQ4.pdb', 'egfr/preparation.json']]
    expected += [result/f'egfr/window_{w:02d}.pdb' for w in range(17)]
    missing = [str(p) for p in expected if not p.is_file() or not p.stat().st_size]
    if missing:
        raise ValueError(f'Missing/empty deliverables: {missing}')
    files = {str(p):dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in expected}
    coords = np.load(result/'model/samples.npz')['x_nm']
    raw = Path('data/raw/1M17.pdb')
    manifest = json.loads((result/'data_manifest.json').read_text())
    # The PDB is redownloadable and ignored; verify its hash when present.
    pdb_verified = None if not raw.exists() else hashlib.sha256(raw.read_bytes()).hexdigest() == manifest['sha256']
    structural_gate = egfr['actual_windows']==egfr['expected_windows']==17 and egfr['max_projection_error_nm']<1e-4
    report = dict(expected_artifacts=len(expected), actual_artifacts=len(files),
                  unique_hashes=len({v['sha256'] for v in files.values()}),
                  coordinate_shape=list(coords.shape), finite_samples=int(np.isfinite(coords).sum()),
                  model_gates_passed=model['model_gates_passed'],
                  structural_gate_passed=structural_gate, pdb_checksum_verified=pdb_verified,
                  files=files)
    (result/'artifacts.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='files'}, indent=2))
    if not model['model_gates_passed'] or not structural_gate or pdb_verified is False:
        raise SystemExit('At least one acceptance gate FAILED; do not claim validation')


if __name__ == '__main__':
    main()
