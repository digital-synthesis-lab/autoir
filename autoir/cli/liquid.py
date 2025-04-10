import json
import os
import uuid

import click
from autoir.analyze import process_file
from autoir.timer import Timer
from autoir.create import DEFAULT_NUM_MOLS, liq_simulation
from autoir.openmm.simulator import DEFAULT_DIPOLES_FILE, OpenMMSimulator


@click.command("liquid")
@click.argument("smiles")
@click.option("-o", "--output", default=None, help="Output directory")
@click.option(
    "-n",
    "--n_mols",
    default=DEFAULT_NUM_MOLS,
    help="Number of molecules inside the box",
)
@click.option("--time_step", default=2, help="Simulation temperature in K")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--pressure", default=1.0, help="Simulation pressure in atm")
@click.option(
    "--target_density", default=0.5, help="Target density when packing the molecules"
)
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--npt_equi_steps", default=500_000, help="Number of equilibration steps")
@click.option("--nvt_equi_steps", default=50_000, help="Number of equilibration steps")
@click.option("--nvt_prod_steps", default=1_000_000, help="Number of production steps")
@click.option(
    "--trj_freq", default=10_000, help="Number of steps for dumping the trajectory"
)
def liquid_sim(
    smiles,
    output,
    n_mols,
    time_step,
    temperature,
    pressure,
    target_density,
    seed,
    npt_equi_steps,
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

    click.echo("Creating liquid box")

    # Generate LAMMPS data file and header
    topology, interchange = liq_simulation(
        smiles, n_mols=n_mols, target_density=target_density
    )

    sim = OpenMMSimulator(
        time_step=time_step,
        temperature=temperature,
        pressure=pressure,
        trj_freq=trj_freq,
    )

    with Timer() as t:
        sim.run(
            interchange,
            npt_equi_steps=npt_equi_steps,
            nvt_equi_steps=nvt_equi_steps,
            nvt_prod_steps=nvt_prod_steps,
        )
    runtime = t.time

    click.echo("Processing and saving files")
    # uses the last half of the production simulation
    prod = pd.read_csv("prod.csv")
    prod = prod.iloc[len(prod) // 2:]

    avg_E = prod['Potential Energy (kJ/mole)'].mean()
    avg_D = prod['Density (g/mL)'].mean()
    avg_V = prod['Box Volume (nm^3)'].mean()

    # computes and processes the infrared spectra
    df = process_file(DEFAULT_DIPOLES_FILE, out_file="ir.csv")
    wn, ir = smooth_ir(df)

    # Prepare parameters for reproducibility
    params = {
        "id": sim_id,
        "smiles": smiles,
        "n_mols": n_mols,
        "n_atoms": topology.n_atoms,
        "temperature": temperature,
        "pressure": pressure,
        "seed": seed,
        "npt_equi_steps": npt_equi_steps,
        "nvt_equi_steps": nvt_equi_steps,
        "nvt_prod_steps": nvt_prod_steps,
        "trj_freq": trj_freq,
        "phase": "liquid",
        "runtime": runtime,
        "avg_energy": avg_E,
        "avg_density": avg_D,
        "avg_volume": avg_V,
        "ir": ir.tolist()
    }

    with open("job.json", "w") as f:
        json.dump(params, f)
