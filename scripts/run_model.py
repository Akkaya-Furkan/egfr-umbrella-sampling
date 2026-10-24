#!/usr/bin/env python3
"""Generate actual biased Langevin trajectories; run with systemd-inhibit."""
import json
from importlib import metadata as package_metadata
from pathlib import Path
import platform
import sys
import time
import numpy as np
import openmm as mm
from openmm import unit


def main():
    out = Path('results/model')
    out.mkdir(parents=True, exist_ok=True)
    environment = dict(python=sys.version, openmm_full_version=mm.version.full_version,
                       packages={name: package_metadata.version(name) for name in
                                 ['numpy', 'scipy', 'matplotlib', 'openmm', 'pytest']})
    (out.parent/'environment.json').write_text(json.dumps(environment, indent=2)+'\n')
    config = dict(temperature_K=300.0, mass_Da=10.0, friction_per_ps=10.0,
                  dt_ps=0.002, burn_steps=50000, production_steps=500000,
                  stride=50, spring_kj_mol_nm2=100.0,
                  centers_nm=np.linspace(-0.8, 0.8, 17).tolist(), replicas=2)
    samples = np.empty((2, 17, 10000))
    records = []
    start = time.monotonic()
    for rep in range(2):
        for w, center in enumerate(config['centers_nm']):
            seed = 230000 + rep * 100 + w
            system = mm.System()
            system.addParticle(config['mass_Da'])
            force = mm.CustomExternalForce('8*((x/0.5)^2-1)^2+0.5*k*(x-c)^2+500*(y*y+z*z)')
            force.addGlobalParameter('k', config['spring_kj_mol_nm2'])
            force.addGlobalParameter('c', center)
            force.addParticle(0, [])
            system.addForce(force)
            integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 10/unit.picosecond, 0.002*unit.picoseconds)
            integrator.setRandomNumberSeed(seed)
            context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
            context.setPositions([[center, 0, 0]])
            context.setVelocitiesToTemperature(300*unit.kelvin, seed + 1000)
            integrator.step(config['burn_steps'])
            for frame in range(samples.shape[-1]):
                integrator.step(config['stride'])
                samples[rep, w, frame] = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)[0, 0]
            energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            records.append(dict(replica=rep, window=w, seed=seed, samples=samples.shape[-1], final_energy_kj_mol=energy))
            print(f'replica={rep} window={w} samples={samples.shape[-1]} elapsed_s={time.monotonic()-start:.2f}', flush=True)
            del context, integrator
    np.savez_compressed(out/'samples.npz', x_nm=samples)
    metadata = dict(config=config, records=records, elapsed_seconds=time.monotonic()-start,
                    openmm_version=mm.__version__, platform='Reference', python=platform.python_version())
    (out/'simulation.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
