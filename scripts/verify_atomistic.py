#!/usr/bin/env python3
"""Audit saved relaxed states independently of the live MD loop, then hash evidence."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import openmm as mm
from openmm import app, unit as u


def digest(path):
    return dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    out=Path(sys.argv[1] if len(sys.argv)>1 else 'results/atomistic')
    windows=json.loads((out/'windows.json').read_text()) if (out/'windows.json').exists() else []
    prm=app.AmberPrmtopFile('data/atomistic/solvated.prmtop')
    system=prm.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1*u.nanometer,constraints=app.HBonds)
    atoms=list(prm.topology.atoms())
    protein=[a.index for a in atoms if a.residue.name not in ('AQ4','WAT','HOH','Na+','Cl-','NA','CL') and a.element!=app.element.hydrogen]
    ligand=[a.index for a in atoms if a.residue.name=='AQ4' and a.element!=app.element.hydrogen]
    masses=np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())])
    verified=[]
    for record in windows:
        path=out/f"window_{record['window']:02d}.xml"
        if not path.exists():
            # A rejected clashing start can appear in windows.json without a saved state.
            if record['min_heavy_nm']>=.15: raise ValueError(f'Missing nonclashing state {path}')
            continue
        state=mm.XmlSerializer.deserialize(path.read_text())
        xyz=state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        box=np.diag(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
        if len(xyz)!=system.getNumParticles() or not np.isfinite(xyz).all():
            raise ValueError(f'Invalid state coordinates {path}')
        delta=np.average(xyz[ligand],axis=0,weights=masses[ligand])-np.average(xyz[protein],axis=0,weights=masses[protein])
        delta-=box*np.round(delta/box)
        com=float(np.linalg.norm(delta))
        minimum=float(cKDTree(xyz[protein]%box,boxsize=box).query(xyz[ligand]%box)[0].min())
        if abs(com-record['actual_com_nm'])>1e-10 or abs(minimum-record['min_heavy_nm'])>1e-10:
            raise ValueError(f'Start-state measurements differ on read-back: {path}')
        verified.append(dict(window=record['window'],com_nm=com,minimum_heavy_nm=minimum,**digest(path)))
    report=dict(expected_start_states=json.loads((out/'analysis.json').read_text())['expected_windows'],verified_start_states=len(verified),
                nonclashing_start_states=sum(r['minimum_heavy_nm']>=.15 for r in verified),
                states=verified)
    (out/'state_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    expected=['repair.json','parameterization.json','preparation_audit.json','simulation.json',
              'pulling.csv','analysis.json','state_verification.json']
    for optional in ('attempt.json','equilibration.json','equilibration_thermodynamics.csv',
                     'windows.json','thermodynamics.csv','pmf.csv','UNACCEPTED_diagnostic.csv',
                     'histograms.npz','window_histograms.png'):
        if (out/optional).exists(): expected.append(optional)
    analysis=json.loads((out/'analysis.json').read_text())
    expected += [f'samples_{i:02d}.npz' for i in analysis['completed_window_ids']]
    expected += list(analysis.get('partial_windows',{}))
    # Round 2 reuses round-1 preparation records stored in results/atomistic.
    shared=Path('results/atomistic')
    paths=[out/p if (out/p).exists() or out==shared else shared/p for p in expected]
    if any(not p.is_file() or p.stat().st_size==0 for p in paths): raise ValueError('Missing retained atomistic evidence')
    files={p.name:digest(p) for p in paths}
    for name,old in analysis['input_series'].items():
        if files[name]!=old: raise ValueError(f'Coordinate sample changed after analysis: {name}')
    manifest=dict(expected_artifacts=len(paths),verified_artifacts=len(files),
                  unique_hashes=len({v['sha256'] for v in files.values()}),files=files)
    (out/'artifacts.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='files'}))
    print(json.dumps({k:v for k,v in report.items() if k!='states'}))

if __name__=='__main__': main()
