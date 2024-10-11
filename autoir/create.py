from typing import List
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


def estimate_box_size(mol: Molecule):
    heavy_atoms = [at for at in mol.atoms if at.atomic_number > 1]
    return 2.7 * len(heavy_atoms) / 3


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


def topology_binary(
    mol1: Molecule,
    mol2: Molecule,
    ratio: float = 0.5,
    copies: List[int] = 200,
    box_size: float = 2.6,
):
    topology = pack_box(
        molecules=[mol1, mol2],
        number_of_copies=[round(copies * ratio), round(copies * (1 - ratio))],
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

    box_size = estimate_box_size(mol)

    topology = topology_single_liquid(mol, box_size=box_size)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)


def mix_simulation(smiles1, smiles2, ratio: float = 0.5):
    mol1 = Molecule.from_smiles(smiles1)
    mol2 = Molecule.from_smiles(smiles2)

    # estimates the size of the box given the number of heavy atoms
    box_size = (estimate_box_size(mol1) + estimate_box_size(mol2)) / 2

    topology = topology_binary(mol1, mol2, ratio, box_size=box_size)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)
