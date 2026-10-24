#!/usr/bin/env python3
"""Read back input chemistry, package records, download hashes and retained files."""
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from collections import Counter
import parmed


def info(path):
    return dict(bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    base=Path('data/atomistic'); out=Path('results/atomistic')
    system=parmed.load_file(str(base/'solvated.prmtop'),str(base/'solvated.inpcrd'))
    leap_log=(base/'leap.log').read_text()
    exits=re.findall(r'Exiting LEaP: Errors = (\d+); Warnings = (\d+)',leap_log)
    if not exits or int(exits[-1][0])!=0: raise ValueError('LEaP did not finish without errors')
    ligand=[a for a in system.atoms if a.residue.name=='AQ4']
    counts=Counter(r.name for r in system.residues)
    sdf=(base/'AQ4.sdf').read_text().splitlines()
    crystal=[[float(l[i:i+10]) for i in (0,10,20)] for l in sdf[4:4+int(sdf[3][:3])]]
    import numpy as np
    dry=parmed.load_file(str(base/'dry.prmtop'),str(base/'dry.inpcrd'))
    heavy=np.array([[a.xx,a.xy,a.xz] for a in dry.atoms if a.residue.name=='AQ4' and a.atomic_number!=1])
    change=float(np.linalg.norm(heavy-np.array(crystal),axis=1).max())
    charge=sum(a.charge for a in system.atoms)
    ligand_charge=sum(a.charge for a in ligand)
    if len(ligand)!=52 or abs(ligand_charge)>1e-5 or abs(charge)>1e-5 or change>1e-5:
        raise ValueError('Topology ligand identity/charge/coordinate audit failed')
    if info(Path('data/raw/1M17.pdb'))['sha256'] != json.loads(Path('results/data_manifest.json').read_text())['sha256']:
        raise ValueError('Source PDB checksum changed')
    packages={}
    openmm_python=Path(os.environ.get('OPENMM_PYTHON',str(Path.home()/'.conda/envs/openmm_env/bin/python')))
    envs = {'ambertools':Path(sys.prefix),'openmm':openmm_python.parent.parent}
    for env, prefix in envs.items():
        records=[]
        for path in (prefix/'conda-meta').glob('*.json'):
            data=json.loads(path.read_text())
            if data['name'] in ('ambertools','openmm','pdbfixer','parmed','openbabel','openmm-hip'):
                records.append({k:data.get(k) for k in ('name','version','build','channel')})
        packages[env]=records
    live = subprocess.run([str(openmm_python),'-c',"import json,openmm,importlib.metadata as m; print(json.dumps(dict(openmm=openmm.__version__,pdbfixer=m.version('pdbfixer'))))"],check=True,capture_output=True,text=True)
    runtime_versions=json.loads(live.stdout)
    inputs={p.name:info(p) for p in [base/'solvated.prmtop',base/'solvated.inpcrd',base/'dry.prmtop',base/'dry.inpcrd',base/'ligand.mol2',base/'ligand.frcmod',base/'protein_amber.pdb']}
    remote_sources=[dict(url='https://files.rcsb.org/download/1M17.pdb',**info(Path('data/raw/1M17.pdb'))),dict(url='https://models.rcsb.org/v1/1m17/ligand?auth_asym_id=A&auth_seq_id=999&encoding=sdf',**info(base/'AQ4.sdf'))]
    report=dict(system_atoms=len(system.atoms),residue_counts=dict(counts),
                ligand_atoms=len(ligand),ligand_heavy_atoms=len(heavy),ligand_charge_e=ligand_charge,
                system_charge_e=charge,box_angstrom=system.box.tolist(),
                max_ligand_crystal_coordinate_change_angstrom=change,
                expected_prepared_inputs=7,verified_prepared_inputs=len(inputs),inputs=inputs,
                conda_installation_records=packages,runtime_versions=runtime_versions,remote_sources=remote_sources,
                preparation_errors_resolved=['Initial LEaP failed: NGLY terminal H lacked atom type; renamed H to H1.',
                                             'AM1-BCC text charges summed to -0.001e; uniform +0.000019230769e correction per atom.'],
                leap_final_errors=int(exits[-1][0]),leap_final_warnings=int(exits[-1][1]))
    (out/'preparation_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('inputs','packages','residue_counts')},indent=2))

if __name__=='__main__': main()
