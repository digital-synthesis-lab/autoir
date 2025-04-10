import click
from autoir.analyze import process_file, DEFAULT_TRUNCATE


@click.command("analyze")
@click.option("-i", "--input", default="dipoles.csv", help="Input file")
@click.option("-o", "--output", default="ir.csv", help="Output file")
@click.option("-t", "--truncate_autocorr", default=DEFAULT_TRUNCATE, help="Number of steps to truncate the autocorrelation function")
def analyze_file(input, output, truncate_autocorr):
    """Compute the IR spectrum from the output LAMMPS simulation"""
    process_file(input, output, truncate_autocorr=truncate_autocorr)
