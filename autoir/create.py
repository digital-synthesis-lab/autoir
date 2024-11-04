from typing import List
from rdkit.Chem import AllChem as Chem
from openff.toolkit import ForceField, Molecule, unit, Topology
from openff.interchange import Interchange
from openff.interchange.components.mdconfig import MDConfig
from openff.interchange.components._packmol import (
    UNIT_CUBE,
    RHOMBIC_DODECAHEDRON,
    pack_box,
)


DEFAULT_FF = "openff_unconstrained-2.0.0.offxml"
DEFAULT_BOX_SIZE_GAS = 4.0  # nm for gas phase
DEFAULT_BOX_SIZE_LIQ = 2.6  # nm for liq phase
DEFAULT_NUM_MOLS = 100


def estimate_box_size(smiles: str, n_mols: int = DEFAULT_NUM_MOLS, scaling: float = 2):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.EmbedMolecule(mol)
    vol = Chem.ComputeMolVolume(mol)
    total_vol = n_mols * scaling * vol
    size = total_vol ** (1 / 3)
    return size


def topology_gas(mol: Molecule):
    mol.generate_conformers(n_conformers=1, rms_cutoff=0.1 * unit.angstrom),
    topology = mol.to_topology()

    return topology


def topology_single_liquid(
    mol: Molecule,
    n_mols: int = DEFAULT_NUM_MOLS,
    box_size: float = DEFAULT_BOX_SIZE_LIQ,
):
    topology = pack_box(
        molecules=[mol],
        number_of_copies=[n_mols],
        box_vectors=box_size * UNIT_CUBE * unit.nanometer,
    )

    return topology


def topology_binary(
    mol1: Molecule,
    mol2: Molecule,
    ratio: float = 0.5,
    n_mols: List[int] = 200,
    box_size: float = DEFAULT_BOX_SIZE_LIQ,
):
    topology = pack_box(
        molecules=[mol1, mol2],
        number_of_copies=[round(n_mols * ratio), round(n_mols * (1 - ratio))],
        box_vectors=box_size * UNIT_CUBE * unit.nanometer,
    )

    return topology


def gas_simulation(smiles, box_size: float = DEFAULT_BOX_SIZE_GAS):
    mol = Molecule.from_smiles(smiles)
    topology = topology_gas(mol)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.box = unit.Quantity([box_size] * 3, unit.nanometer)
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)


def liq_simulation(smiles, n_mols: int = DEFAULT_NUM_MOLS, box_size: float = None):
    mol = Molecule.from_smiles(smiles)

    if box_size is None:
        box_size = estimate_box_size(smiles, n_mols=n_mols)

    topology = topology_single_liquid(mol, box_size=box_size, n_mols=n_mols)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)


def mix_simulation(
    smiles1, smiles2, ratio: float = 0.5, n_mols: int = DEFAULT_NUM_MOLS
):
    mol1 = Molecule.from_smiles(smiles1)
    mol2 = Molecule.from_smiles(smiles2)

    # estimates the size of the box given the number of heavy atoms
    box_size = (
        estimate_box_size(mol1, n_mols=n_mols) + estimate_box_size(mol2, n_mols=n_mols)
    ) / 2

    topology = topology_binary(mol1, mol2, ratio, box_size=box_size)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.to_lammps("out.lmp")
    mdconfig = MDConfig.from_interchange(interchange)
    mdconfig.write_lammps_input(input_file="header.in", interchange=interchange)
