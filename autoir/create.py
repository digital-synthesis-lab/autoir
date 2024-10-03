from openff.toolkit import ForceField, Molecule, unit, Topology

from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig
from openff.interchange.components._packmol import (
    UNIT_CUBE,
    RHOMBIC_DODECAHEDRON,
    pack_box,
)


DEFAULT_FF = "openff_unconstrained-2.0.0.offxml"
DEFAULT_BOX_SIZE = 4.0  # nm for gas phase


def topology_gas(mol: Molecule):
    mol.generate_conformers(n_conformers=1, rms_cutoff=0.1 * unit.angstrom),
    topology = mol.to_topology()

    return topology


def topology_single_liquid(mol: Molecule, copies: int = 200, box_size: float = 2.6):
    topology = pack_box(
        molecules=[mol],
        number_of_copies=[copies],
        box_vectors=box_size * UNIT_CUBE * unit.nanometer,
    )

    return topology


def gas_simulation(smiles):
    mol = Molecule.from_smiles(smiles)
    topology = topology_gas(mol)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.box = unit.Quantity([DEFAULT_BOX_SIZE] * 3, unit.nanometer)
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)


def liq_simulation(smiles):
    mol = Molecule.from_smiles(smiles)

    # estimates the size of the box given the number of heavy atoms
    heavy_atoms = [
        at for at in mol.atoms
        if at.atomic_number > 1
    ]
    box_size = 2.7 * len(heavy_atoms) / 3

    topology = topology_single_liquid(mol, box_size=box_size)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)
