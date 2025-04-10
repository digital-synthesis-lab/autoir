import json
import os
import uuid

import click
from autoir.analyze import process_file
from autoir.timer import Timer
from autoir.create import gas_simulation
from autoir.openmm.simulator import DEFAULT_DIPOLES_FILE, OpenMMSimulator


@click.command("gas")
@click.argument("smiles")
@click.option("-o", "--output", default=None, help="Output directory")
@click.option("--time_step", default=2, help="Simulation timestep in fs (default: 2 fs)")
@click.option("--temperature", default=300, help="Simulation temperature in K (default: 300 K)")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--nvt_equi_steps", default=100_000, help="Number of equilibration steps")
@click.option("--nvt_prod_steps", default=1_000_000, help="Number of production steps")
@click.option(
    "--trj_freq", default=10_000, help="Number of steps for dumping the trajectory"
)
def liquid_sim(
    smiles,
    output,
    time_step,
    temperature,
    seed,
    nvt_equi_steps,
    nvt_prod_steps,
    trj_freq,
):
    """Run a liquid phase simulation for the given SMILES string."""
    sim_id = str(uuid.uuid4())

    if output is not None:
        sim_dir = output
    else:
        sim_dir = os.path.join(os.getcwd(), sim_id)

    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    click.echo("Creating gas box")

    # Generate LAMMPS data file and header
    topology, interchange = gas_simulation(smiles)

    sim = OpenMMSimulator(
        time_step=time_step,
        temperature=temperature,
        pressure=pressure,
        trj_freq=trj_freq,
    )

    with Timer() as t:
        sim.run(
            interchange,
            npt_equi_steps=0,
            nvt_equi_steps=nvt_equi_steps,
            nvt_prod_steps=nvt_prod_steps,
        )
    runtime = t.time

    click.echo("Processing and saving file")
    process_file(DEFAULT_DIPOLES_FILE, out_file="ir.csv")

    # Prepare parameters for reproducibility
    params = {
        "id": sim_id,
        "smiles": smiles,
        "n_mols": 1,
        "n_atoms": topology.n_atoms,
        "temperature": temperature,
        "pressure": 0.0,
        "seed": seed,
        "npt_equi_steps": 0,
        "nvt_equi_steps": nvt_equi_steps,
        "nvt_prod_steps": nvt_prod_steps,
        "trj_freq": trj_freq,
        "phase": "gas",
        "runtime": runtime,
    }

    with open("job.json", "w") as f:
        json.dump(params, f)
