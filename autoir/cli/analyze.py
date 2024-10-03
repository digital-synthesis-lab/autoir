import click
import uuid
import os
from autoir.analyze import process_file


@click.command("analyze")
@click.option("-i", "--input", default="dipoles.csv", help="Input file")
@click.option("-o", "--output", default="ir.csv", help="Output file")
def analyze_file(input, output):
    """Compute the IR spectrum from the output LAMMPS simulation"""
    process_file(input, output)
