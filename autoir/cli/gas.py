import json
import click
import uuid
import os
from autoir.create import gas_simulation
from autoir.render import render_input_file, write_input_file


@click.command("gas")
@click.argument("smiles")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--equi_steps", default=50000, help="Number of equilibration steps")
@click.option("--prod_steps", default=500000, help="Number of production steps")
def gas_sim(smiles, temperature, seed, equi_steps, prod_steps):
    """Run a gas phase simulation for the given SMILES string."""
    sim_id = str(uuid.uuid4())
    sim_dir = os.path.join(os.getcwd(), sim_id)
    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    # Generate LAMMPS data file and header
    gas_simulation(smiles)

    # Prepare parameters for the main input file
    params = {
        "smiles": smiles,
        "temperature": temperature,
        "seed": seed,
        "equi_steps": equi_steps,
        "prod_steps": prod_steps,
        "dump_freq": None,
        "phase": "gas",
    }

    # Render and write the main input file
    write_input_file(params, "gas.in")

    with open("job.json", "w") as f:
        json.dump(params, f)

    click.echo(f"Gas phase simulation files created in directory: {sim_dir}")
