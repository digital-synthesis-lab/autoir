from openff.toolkit import ForceField, Molecule, unit, Topology

from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig
from openff.interchange.components._packmol import (
    UNIT_CUBE,
    RHOMBIC_DODECAHEDRON,
    pack_box,
)


def prep_topology():
    mol = Molecule.from_smiles("C(Cl)(Cl)(Cl)Cl")
    mol.generate_conformers(n_conformers=1, rms_cutoff=0.1 * unit.angstrom),
    topology = mol.to_topology()

    return topology


if __name__ == "__main__":
    topology = prep_topology()
    print("obtained topology")

    sage = ForceField("openff_unconstrained-2.0.0.offxml")

    interchange: Interchange = Interchange.from_smirnoff(
        force_field=sage, topology=topology
    )
    interchange.box = unit.Quantity([4, 4, 4], unit.nanometer)
    interchange.to_lammps("interchange.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="auto_generated.in", interchange=interchange)
