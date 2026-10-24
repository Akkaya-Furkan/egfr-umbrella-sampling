#!/usr/bin/env python3
"""Bounded explicit-solvent EGFR pilot; checkpoint stages and scalar observations."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.spatial import cKDTree
import openmm as mm
from openmm import app, unit as u

OUT = Path('results/atomistic')
DATA = Path('data/atomistic')


def save_json(name, value):
    dest = OUT/name
    temp = dest.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(dest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--deadline-unix', type=float, required=True)
    parser.add_argument('--production-ps', type=float, default=100)
    parser.add_argument('--resume-equilibrated', action='store_true')
    args = parser.parse_args()
    started = time.time()
    config = dict(seed=23024, temperature_K=300, timestep_ps=.002,
                  friction_per_ps=1, spring_kj_mol_nm2=1000,
                  ca_restraint_kj_mol_nm2=100, burn_ps=20,
                  production_ps=args.production_ps, sample_interval_ps=.1,
                  n_windows=17, displacement_nm=1.6, pull_ps=400,
                  nvt_ps=20, npt_ps=80, added_salt_molar=0)
    progress = dict(status='initializing', completed_windows=[], started_unix=started,
                    config=config)
    save_json('simulation.json', progress)

    def stage(name, **extra):
        progress.update(status=name, elapsed_seconds=time.time()-started, **extra)
        save_json('simulation.json', progress)
        print(name, json.dumps(extra), flush=True)

    def budget():
        if time.time() > args.deadline_unix-30:
            raise TimeoutError('Task wall-clock budget reached; retained completed windows only')

    try:
        prm = app.AmberPrmtopFile(str(DATA/'solvated.prmtop'))
        inp = app.AmberInpcrdFile(str(DATA/'solvated.inpcrd'))
        system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer,
                                  constraints=app.HBonds, rigidWater=True)
        atoms = list(prm.topology.atoms())
        protein = [a.index for a in atoms if a.residue.name not in ('AQ4','WAT','HOH','Na+','Cl-','NA','CL') and a.element != app.element.hydrogen]
        ligand = [a.index for a in atoms if a.residue.name == 'AQ4' and a.element != app.element.hydrogen]
        ca = [a.index for a in atoms if a.name == 'CA' and a.index in protein]
        if len(ligand) != 29 or len(protein) < 2000:
            raise ValueError(f'Unexpected groups: protein={len(protein)}, ligand={len(ligand)}')
        masses = np.array([system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles())])
        wp, wl = masses[protein], masses[ligand]
        initial = np.array(inp.positions.value_in_unit(u.nanometer))
        r0 = float(np.linalg.norm(np.average(initial[ligand],axis=0,weights=wl)-np.average(initial[protein],axis=0,weights=wp)))
        centers = r0+np.linspace(0, config['displacement_nm'], config['n_windows'])
        config['centers_nm'] = centers.tolist()
        bias = mm.CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-r0)^2')
        bias.addGlobalParameter('k', 0)
        bias.addGlobalParameter('r0', r0)
        bias.addGroup(protein, wp.tolist()); bias.addGroup(ligand, wl.tolist())
        bias.addBond([0,1], [])
        bias.setUsesPeriodicBoundaryConditions(True)
        system.addForce(bias)
        restraint = mm.CustomExternalForce('0.5*kca*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('kca', config['ca_restraint_kj_mol_nm2'])
        for name in ('x0','y0','z0'): restraint.addPerParticleParameter(name)
        for i in ca: restraint.addParticle(i, initial[i].tolist())
        system.addForce(restraint)
        barostat = mm.MonteCarloBarostat(1*u.bar,300*u.kelvin,25)
        barostat.setRandomNumberSeed(config['seed'])
        barostat.setFrequency(0)
        system.addForce(barostat)
        integrator = mm.LangevinMiddleIntegrator(300*u.kelvin,1/u.picosecond,.002*u.picoseconds)
        integrator.setRandomNumberSeed(config['seed'])
        platform = mm.Platform.getPlatformByName('HIP')
        sim = app.Simulation(prm.topology,system,integrator,platform,{'Precision':'mixed'})
        sim.context.setPositions(inp.positions)
        sim.context.setPeriodicBoxVectors(*inp.boxVectors)
        progress.update(platform=platform.getName(), device=platform.getPropertyValue(sim.context,'DeviceName'),
                        system_atoms=system.getNumParticles(), protein_heavy_atoms=len(protein), ligand_heavy_atoms=len(ligand),
                        protein_ca_atoms=len(ca), openmm_version=mm.__version__)
        sim.reporters.append(app.StateDataReporter(str(OUT/'thermodynamics.log'),5000,step=True,time=True,
                                                   potentialEnergy=True,temperature=True,volume=True,speed=True))

        def run_ps(ps):
            steps = round(ps/.002)
            while steps > 0:
                budget(); n = min(1000,steps); sim.step(n); steps -= n

        def observe():
            state = sim.context.getState(getPositions=True,getEnergy=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
            box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer)
            lengths = np.diag(box)
            delta = np.average(pos[ligand],axis=0,weights=wl)-np.average(pos[protein],axis=0,weights=wp)
            delta -= lengths*np.round(delta/lengths)
            r = float(np.linalg.norm(delta))
            energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
            # Periodic minimum-image atom distances, even for unwrapped molecules.
            mind = float(cKDTree(pos[protein] % lengths,boxsize=lengths).query(pos[ligand] % lengths)[0].min())
            if not np.isfinite([r,energy,mind]).all(): raise ValueError('Nonfinite state')
            if r >= min(lengths)/2-.1: raise ValueError('Radial coordinate approaches periodic image limit')
            return r,energy,mind,lengths.tolist()

        if args.resume_equilibrated:
            sim.loadCheckpoint(str(OUT/'equilibrated.chk'))
            stage('resumed_equilibrated_checkpoint')
        else:
            stage('minimizing')
            sim.minimizeEnergy(maxIterations=2000)
            sim.context.setVelocitiesToTemperature(300*u.kelvin, config['seed'])
    
            stage('nvt'); run_ps(config['nvt_ps'])
            barostat.setFrequency(25); sim.context.reinitialize(preserveState=True)
            stage('npt'); run_ps(config['npt_ps'])
        # Freeze the equilibrated box for identical Hamiltonians in every window.
        barostat.setFrequency(0); sim.context.reinitialize(preserveState=True)
        sim.saveCheckpoint(str(OUT/'equilibrated.chk'))
        rr, ee, dd, box = observe()
        stage('equilibrated', equilibrated_distance_nm=rr, equilibrated_min_heavy_nm=dd, box_nm=box)
        sim.context.setParameter('k',config['spring_kj_mol_nm2'])
        sim.context.setParameter('r0',float(centers[0]))
        run_ps(20)
        # Continuous SMD, storing states only after physical dynamics (no translations).
        states = [sim.context.getState(getPositions=True,getVelocities=True,getParameters=True)]
        pull = []
        stage('pulling')
        for j in range(1,401):
            target = centers[0]+config['displacement_nm']*j/400
            sim.context.setParameter('r0',float(target)); run_ps(config['pull_ps']/400)
            r,e,d,box = observe(); pull.append([j*config['pull_ps']/400,target,r,d])
            if j % 25 == 0:
                states.append(sim.context.getState(getPositions=True,getVelocities=True,getParameters=True))
            np.savetxt(OUT/'pulling.csv',pull,delimiter=',',header='time_ps,target_nm,com_nm,min_heavy_nm',comments='')
        progress['pulled_states'] = len(states)
        if len(states) != config['n_windows']: raise ValueError('Incomplete SMD states')
        start_records = []
        for w,(center,state) in enumerate(zip(centers,states)):
            budget()
            sim.context.setState(state)
            sim.context.setParameter('r0',float(center))
            sim.context.setParameter('k',config['spring_kj_mol_nm2'])
            integrator.setRandomNumberSeed(config['seed']+w+1)
            sim.context.setVelocitiesToTemperature(300*u.kelvin,config['seed']+w+1)
            stage('window_burn', active_window=w)
            run_ps(config['burn_ps'])
            r,e,d,box = observe()
            start_records.append(dict(window=w,center_nm=float(center),actual_com_nm=r,min_heavy_nm=d))
            save_json('windows.json',start_records)
            if d < .15: raise ValueError(f'Window {w} severe clash after relaxation: {d} nm')
            sim.saveState(str(OUT/f'window_{w:02d}.xml'))
            stage('window_production', active_window=w)
            samples = []
            for sample in range(round(args.production_ps/config['sample_interval_ps'])):
                run_ps(config['sample_interval_ps'])
                state = sim.context.getState(getPositions=True)
                pos = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
                delta = np.average(pos[ligand],axis=0,weights=wl)-np.average(pos[protein],axis=0,weights=wp)
                lengths = np.diag(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.nanometer))
                delta -= lengths*np.round(delta/lengths)
                r = float(np.linalg.norm(delta))
                if not np.isfinite(r) or r >= min(lengths)/2-.1: raise ValueError('Invalid radial sample')
                samples.append(r)
                # Partial files are explicitly separate from complete production.
                if (sample+1)%100 == 0:
                    np.save(OUT/f'partial_{w:02d}.npy',np.array(samples))
            np.savez_compressed(OUT/f'samples_{w:02d}.npz',r_nm=np.array(samples))
            (OUT/f'partial_{w:02d}.npy').unlink(missing_ok=True)
            r,e,d,box = observe()
            progress['completed_windows'].append(w)
            stage('window_complete', active_window=w, last_min_heavy_nm=d)
        stage('production_complete')
    except Exception as exc:
        stage('failed',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        progress['finished_unix'] = time.time()
        progress['elapsed_seconds'] = time.time()-started
        save_json('simulation.json',progress)

if __name__ == '__main__': main()
