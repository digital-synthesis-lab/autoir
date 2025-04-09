from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig


def write_lammps_file(interchange: Interchange):
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)
    return mdconfig
