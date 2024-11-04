import json
import click
import uuid
import os
from autoir.create import liq_simulation, DEFAULT_NUM_MOLS
from autoir.render import render_input_file, write_input_file


@click.command("liquid")
@click.argument("smiles")
@click.option("-o", "--output", default=None, help="Output directory")
@click.option("-n", "--n_mols", default=DEFAULT_NUM_MOLS, help="Number of molecules inside the box")
@click.option("--temperature", default=300, help="Simulation temperature in K")
@click.option("--pressure", default=1.0, help="Simulation pressure in atm")
@click.option("--box_size", default=None, help="Size of the simulation box (optional)")
@click.option("--seed", default=12345, help="Random seed for the simulation")
@click.option("--equi_steps", default=50000, help="Number of equilibration steps")
@click.option("--prod_steps", default=100000, help="Number of production steps")
def liquid_sim(smiles, n_mols, output, temperature, pressure, box_size, seed, equi_steps, prod_steps):
    """Run a liquid phase simulation for the given SMILES string."""
    sim_id = str(uuid.uuid4())

    if output is not None:
        sim_dir = output
    else: 
        sim_dir = os.path.join(os.getcwd(), sim_id)
    os.makedirs(sim_dir, exist_ok=True)

    # Change to the simulation directory
    os.chdir(sim_dir)

    # Generate LAMMPS data file and header
    liq_simulation(smiles, n_mols=n_mols, box_size=box_size)

    # Prepare parameters for the main input file
    params = {
        "id": sim_id,
        "smiles": smiles,
        "n_mols": n_mols,
        "temperature": temperature,
        "pressure": pressure,
        "seed": seed,
        "equi_steps": equi_steps,
        "prod_steps": prod_steps,
        "dump_freq": None,
        "phase": "liquid",
    }

    # Render and write the main input file
    write_input_file(params, "liquid.in")

    with open("job.json", "w") as f:
        json.dump(params, f)

    click.echo(f"Liquid phase simulation files created in directory: {sim_dir}")
