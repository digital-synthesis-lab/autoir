from openff.toolkit import ForceField, Molecule, unit

from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig
from openff.interchange.components._packmol import (
    UNIT_CUBE,
    RHOMBIC_DODECAHEDRON,
    pack_box,
)


def prep_topology():
    #mol = Molecule.from_smiles("CCO")
    solv = Molecule.from_smiles("C(Cl)(Cl)(Cl)Cl")

    topology = pack_box(
        molecules=[solv],
        number_of_copies=[500],
        box_vectors=3.5 * UNIT_CUBE * unit.nanometer,
        # box_vectors=3.5 * RHOMBIC_DODECAHEDRON * unit.nanometer,
    )

    return topology


if __name__ == "__main__":
    topology = prep_topology()
    print("obtained topology")

    sage = ForceField("openff_unconstrained-2.0.0.offxml")

    interchange: Interchange = Interchange.from_smirnoff(
        force_field=sage, topology=topology
    )
    interchange.to_lammps("interchange.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="auto_generated.in", interchange=interchange)
