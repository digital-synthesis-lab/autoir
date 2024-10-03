import click

from autoir.cli.gas import gas_sim
from autoir.cli.liquid import liquid_sim
from autoir.cli.analyze import analyze_file


class AutoIRGroup(click.Group):
    pass


@click.command(cls=AutoIRGroup)
def autoir():
    """Command line interface for AutoIR"""


autoir.add_command(gas_sim)
autoir.add_command(liquid_sim)
autoir.add_command(analyze_file)
