import json
import os
import uuid
import pandas as pd

import click
from autoir.analyze import process_file, smooth_ir
from autoir.timer import Timer
from autoir.create import DEFAULT_NUM_MOLS, mix_simulation
from autoir.openmm.simulator import (
    DEFAULT_DIPOLES_FILE,
    DEFAULT_PROD_FILE,
    DEFAULT_EQUI_FILE,
    DEFAULT_AVGE_FILE,
    OpenMMSimulator,
)


@click.command("mixture")
@click.argument("smiles_1")
@click.argument("smiles_2")
@click.option("-o", "--output", default=None, help="Output directory")
@click.option(
    "-n",
    "--n_mols",
    default=DEFAULT_NUM_MOLS,
    help="Number of molecules inside the box",
)
@click.option("--ratio", default=0.5, help="Ratio between 1 and 2")
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
def mixture_sim(
    smiles_1,
    smiles_2,
    output,
    n_mols,
    ratio,
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
    """Run a liquid phase simulation of a mixture for the given SMILES strings."""
    sim_id = str(uuid.uuid4())

    if output is not None:
        sim_dir = output
    else:
        sim_dir = os.path.join(os.getcwd(), sim_id)

    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    click.echo("Creating mixture box")
    # Generate LAMMPS data file and header
    topology, interchange = mix_simulation(
        smiles_1, smiles_2, ratio=ratio, n_mols=n_mols, target_density=target_density
    )

    # Prepare parameters for the main input file
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

    # the volume is computed from the NVT simulation
    equi = pd.read_csv(DEFAULT_EQUI_FILE)
    avg_D = equi.iloc[-1]["Density (g/mL)"]
    avg_V = equi.iloc[-1]["Box Volume (nm^3)"]

    # uses the last half of the production simulation
    prod = pd.read_csv(DEFAULT_PROD_FILE)
    prod = prod.iloc[len(prod) // 2:]

    #avg_E = prod['Potential Energy (kJ/mole)'].mean()
    # average enegy from every time step
    avg_E = pd.read_csv(DEFAULT_AVGE_FILE)
    avg_E = avg_E.iloc[-1]["avgE (kJ/mol)"]

    # computes and processes the infrared spectra
    df = process_file(DEFAULT_DIPOLES_FILE, out_file="ir.csv")
    wn, ir = smooth_ir(df)

    params = {
        "id": sim_id,
        "smiles_1": smiles_1,
        "smiles_2": smiles_2,
        "n_mols": n_mols,
        "n_atoms": topology.n_atoms,
        "ratio": ratio,
        "temperature": temperature,
        "pressure": pressure,
        "seed": seed,
        "npt_equi_steps": npt_equi_steps,
        "nvt_equi_steps": nvt_equi_steps,
        "nvt_prod_steps": nvt_prod_steps,
        "trj_freq": trj_freq,
        "phase": "mixture",
        "runtime": runtime,
        "avg_energy": avg_E,
        "avg_density": avg_D,
        "avg_volume": avg_V,
        "ir": ir.tolist()
    }

    with open("job.json", "w") as f:
        json.dump(params, f)
