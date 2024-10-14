import json
import click
import uuid
import os
from autoir.create import mix_simulation, DEFAULT_NUM_MOLS
from autoir.render import render_input_file, write_input_file


@click.command("mixture")
@click.argument("smiles_1")
@click.argument("smiles_2")
@click.option("-o", "--output", default=None, help="Output directory")
@click.option("-n", "--n_mols", default=DEFAULT_NUM_MOLS, help="Number of molecules inside the box")
@click.option("--ratio", default=0.5, help="Ratio between 1 and 2")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--pressure", default=1.0, help="Simulation pressure in atm")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--equi_steps", default=50000, help="Number of equilibration steps")
@click.option("--prod_steps", default=300000, help="Number of production steps")
def mixture_sim(
    smiles_1, smiles_2, output, n_mols, ratio, temperature, pressure, seed, equi_steps, prod_steps
):
    """Run a liquid phase simulation of a binary molecular mixture for the given SMILES strings."""
    sim_id = str(uuid.uuid4())

    if output is not None:
        sim_dir = output
    else: 
        sim_dir = os.path.join(os.getcwd(), sim_id)

    # Change to the simulation directory
    os.chdir(sim_dir)

    # Generate LAMMPS data file and header
    mix_simulation(smiles_1, smiles_2, ratio=ratio, n_mols=n_mols)

    # Prepare parameters for the main input file
    params = {
        "id": sim_id,
        "smiles_1": smiles_1,
        "smiles_2": smiles_2,
        "n_mols": n_mols,
        "ratio": ratio,
        "temperature": temperature,
        "pressure": pressure,
        "seed": seed,
        "equi_steps": equi_steps,
        "prod_steps": prod_steps,
        "dump_freq": None,
        "phase": "mixture",
    }

    # Render and write the main input file
    write_input_file(params, "mixture.in")

    with open("job.json", "w") as f:
        json.dump(params, f)

    click.echo(f"Liquid phase simulation files created in directory: {sim_dir}")
