import json
import click
import uuid
import os
from autoir.create import mix_simulation
from autoir.render import render_input_file, write_input_file


@click.command("mixture")
@click.argument("smiles_1")
@click.argument("smiles_2")
@click.option("--ratio", default=0.5, help="Ratio between 1 and 2")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--pressure", default=1.0, help="Simulation pressure in atm")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--equi_steps", default=50000, help="Number of equilibration steps")
@click.option("--prod_steps", default=300000, help="Number of production steps")
def mixture_sim(
    smiles1, smiles2, ratio, temperature, pressure, seed, equi_steps, prod_steps
):
    """Run a liquid phase simulation for the given SMILES string."""
    sim_id = str(uuid.uuid4())
    sim_dir = os.path.join(os.getcwd(), sim_id)
    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    # Generate LAMMPS data file and header
    mix_simulation(smiles1, smiles2, ratio)

    # Prepare parameters for the main input file
    params = {
        "id": sim_id,
        "smiles1": smiles1,
        "smiles2": smiles2,
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
