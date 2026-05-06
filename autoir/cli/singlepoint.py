import click

from autoir.create import DEFAULT_FF
from autoir.singlepoint import single_point_energy


@click.command("singlepoint")
@click.argument("pdb_file", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-s",
    "--smiles",
    multiple=True,
    required=True,
    help="SMILES for each unique molecule in the PDB (repeat for multi-component systems).",
)
@click.option(
    "--frame",
    default=0,
    type=int,
    show_default=True,
    help="Zero-based MODEL index to evaluate.",
)
@click.option(
    "--force_field",
    default=DEFAULT_FF,
    show_default=True,
    help="OpenFF SMIRNOFF force-field filename.",
)
@click.option(
    "--platform",
    default=None,
    type=str,
    help="OpenMM platform (e.g. CUDA, CPU). Defaults to auto-selection.",
)
def singlepoint_cmd(pdb_file, smiles, frame, force_field, platform):
    """Compute the single-point potential energy (kJ/mol) for a PDB frame."""
    energy = single_point_energy(
        pdb_file=pdb_file,
        smiles=list(smiles),
        frame=frame,
        force_field=force_field,
        platform=platform,
    )
    click.echo(f"{energy:.6f}")
    return energy
