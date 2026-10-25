#!/usr/bin/env python3
"""Round 2: independent re-sampling of the round-1 relaxed starts plus one window
added in the only overlap hole (pair 10-11). Gates and force constant unchanged."""
import argparse, json, time
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import app, unit as u

SRC = Path('results/atomistic'); OUT = Path('results/atomistic_round2')
DATA = Path('data/atomistic')


def save_json(name, value):
    t = (OUT/name).with_suffix('.tmp'); t.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n'); t.replace(OUT/name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--deadline-unix', type=float, required=True)
    ap.add_argument('--production-ps', type=float, default=500)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    started = time.time()
    r1 = json.loads((SRC/'simulation.json').read_text())
    old = r1['config']['centers_nm']
    extra = (old[10]+old[11])/2
    centers = old[:11]+[extra]+old[11:]
    source = list(range(11))+[10]+list(range(11, 17))   # round-1 start state per window
    config = dict(r1['config'], production_ps=args.production_ps, n_windows=len(centers),
                  centers_nm=centers, seed=r1['config']['seed']+1000, source_window=source,
                  round=2, note='extra window between round-1 windows 10 and 11, started from window 10 state')
    progress = dict(status='initializing', completed_windows=[], started_unix=started, config=config)
    save_json('simulation.json', progress)
    try:
        prm = app.AmberPrmtopFile(str(DATA/'solvated.prmtop'))
        system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer,
                                  constraints=app.HBonds, rigidWater=True)
        atoms = list(prm.topology.atoms())
        protein = [a.index for a in atoms if a.residue.name not in ('AQ4','WAT','HOH','Na+','Cl-','NA','CL') and a.element != app.element.hydrogen]
        ligand = [a.index for a in atoms if a.residue.name == 'AQ4' and a.element != app.element.hydrogen]
        ca = [a.index for a in atoms if a.name == 'CA' and a.index in protein]
        masses = np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())])
        wp, wl = masses[protein], masses[ligand]
        inp = app.AmberInpcrdFile(str(DATA/'solvated.inpcrd'))
        initial = np.array(inp.positions.value_in_unit(u.nanometer))
        bias = mm.CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-r0)^2')
        bias.addGlobalParameter('k', 0); bias.addGlobalParameter('r0', centers[0])
        bias.addGroup(protein, wp.tolist()); bias.addGroup(ligand, wl.tolist()); bias.addBond([0,1], [])
        bias.setUsesPeriodicBoundaryConditions(True); system.addForce(bias)
        rest = mm.CustomExternalForce('0.5*kca*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        rest.addGlobalParameter('kca', config['ca_restraint_kj_mol_nm2'])
        for n in ('x0','y0','z0'): rest.addPerParticleParameter(n)
        for i in ca: rest.addParticle(i, initial[i].tolist())
        system.addForce(rest)
        bar = mm.MonteCarloBarostat(1*u.bar, 300*u.kelvin, 25); bar.setFrequency(0); system.addForce(bar)
        integ = mm.LangevinMiddleIntegrator(300*u.kelvin, 1/u.picosecond, .002*u.picoseconds)
        sim = app.Simulation(prm.topology, system, integ, mm.Platform.getPlatformByName('HIP'), {'Precision':'mixed'})
        sim.loadState(str(SRC/'window_00.xml'))
        box = np.diag(sim.context.getState().getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
        progress.update(platform='HIP', openmm_version=mm.__version__, box_nm=box.tolist())

        def run_ps(ps):
            n = round(ps/.002)
            while n > 0:
                if time.time() > args.deadline_unix-30: raise TimeoutError('deadline reached; completed windows retained')
                k = min(1000, n); sim.step(k); n -= k

        def com():
            st = sim.context.getState(getPositions=True)
            pos = st.getPositions(asNumpy=True).value_in_unit(u.nanometer)
            L = np.diag(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
            d = np.average(pos[ligand],axis=0,weights=wl)-np.average(pos[protein],axis=0,weights=wp)
            d -= L*np.round(d/L); r = float(np.linalg.norm(d))
            if not np.isfinite(r) or r >= min(L)/2-.1: raise ValueError('Invalid radial sample')
            from scipy.spatial import cKDTree
            md = float(cKDTree(pos[protein] % L, boxsize=L).query(pos[ligand] % L)[0].min())
            return r, md

        records = []
        for w, (c, s) in enumerate(zip(centers, source)):
            sim.loadState(str(SRC/f'window_{s:02d}.xml'))
            sim.context.setParameter('k', config['spring_kj_mol_nm2']); sim.context.setParameter('r0', float(c))
            seed = config['seed']+w+1
            integ.setRandomNumberSeed(seed); sim.context.setVelocitiesToTemperature(300*u.kelvin, seed)
            progress.update(status='window_burn', active_window=w); save_json('simulation.json', progress)
            run_ps(config['burn_ps'])
            r, md = com()
            records.append(dict(window=w, center_nm=float(c), actual_com_nm=r, min_heavy_nm=md))
            save_json('windows.json', records)
            if md < .15: raise ValueError(f'Window {w} clash {md}')
            sim.saveState(str(OUT/f'window_{w:02d}.xml'))
            progress.update(status='window_production'); save_json('simulation.json', progress)
            xs = []
            for _ in range(round(args.production_ps/config['sample_interval_ps'])):
                run_ps(config['sample_interval_ps']); xs.append(com()[0])
            np.savez_compressed(OUT/f'samples_{w:02d}.npz', r_nm=np.array(xs))
            progress['completed_windows'].append(w)
            progress.update(status='window_complete', elapsed_seconds=time.time()-started); save_json('simulation.json', progress)
            print('window_complete', w, flush=True)
        progress['status'] = 'production_complete'
    except Exception as exc:
        progress.update(status='failed', error=f'{type(exc).__name__}: {exc}'); raise
    finally:
        progress['elapsed_seconds'] = time.time()-started; save_json('simulation.json', progress)

if __name__ == '__main__': main()
