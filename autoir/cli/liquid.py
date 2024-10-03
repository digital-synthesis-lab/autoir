import click
import uuid
import os
from autoir.create import liq_simulation
from autoir.render import render_input_file, write_input_file


@click.command("liquid")
@click.argument("smiles")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--pressure", default=1.0, help="Simulation pressure in atm")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--equi_steps", default=10000, help="Number of equilibration steps")
@click.option("--prod_steps", default=100000, help="Number of production steps")
def liquid_sim(smiles, temperature, pressure, seed, equi_steps, prod_steps):
    """Run a liquid phase simulation for the given SMILES string."""
    sim_id = str(uuid.uuid4())
    sim_dir = os.path.join(os.getcwd(), sim_id)
    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    # Generate LAMMPS data file and header
    liq_simulation(smiles)

    # Prepare parameters for the main input file
    params = {
        "temperature": temperature,
        "pressure": pressure,
        "seed": seed,
        "equi_steps": equi_steps,
        "prod_steps": prod_steps,
    }

    # Render and write the main input file
    write_input_file(params, "simulation.in")

    click.echo(f"Liquid phase simulation files created in directory: {sim_dir}")
