import click

from autoir.cli.gas import gas_sim
from autoir.cli.liquid_sim import liquid_sim


class AutoIRGroup(click.Group):
    pass


@click.command(cls=AutoIRGroup)
def autoir():
    """Command line interface for AutoIR"""


autoir.add_command(gas_sim)
autoir.add_command(liquid_sim)
