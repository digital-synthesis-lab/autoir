from typing import List

from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box
from openff.toolkit import ForceField, Molecule, Topology, unit
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors

from . import const

DEFAULT_FF = "openff_unconstrained-2.0.0.offxml"
DEFAULT_BOX_SIZE_GAS = 4.0  # nm for gas phase
DEFAULT_DENSITY = 0.5  # g/cm3
DEFAULT_NUM_MOLS = 80


def estimate_box_size(
    smiles: str, n_mols: int = DEFAULT_NUM_MOLS, target_density: float = DEFAULT_DENSITY
):
    """Estimates the box size to obtain a targeted density (in g/cm3)"""
    mol = Chem.MolFromSmiles(smiles)
    mass = Descriptors.ExactMolWt(mol)  # g/mol
    box_mass = (mass / const.N_AVOGADRO) * n_mols  # g
    box_vol = box_mass / target_density  # cm^3
    box_vol = box_vol * 1e21  # nm^3
    size = box_vol ** (1 / 3)
    return size


def topology_gas(mol: Molecule) -> Topology:
    mol.generate_conformers(n_conformers=1, rms_cutoff=0.1 * unit.angstrom),
    topology = mol.to_topology()

    return topology


def topology_single_liquid(
    mol: Molecule,
    box_size: float,
    n_mols: int = DEFAULT_NUM_MOLS,
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
    box_size: float,
    ratio: float = 0.5,
    n_mols: List[int] = 200,
):
    topology = pack_box(
        molecules=[mol1, mol2],
        number_of_copies=[round(n_mols * ratio), round(n_mols * (1 - ratio))],
        box_vectors=box_size * UNIT_CUBE * unit.nanometer,
    )

    return topology


def topology_n_mols(
    mol_list: List[Molecule],
    ratio_list: List[float],
    box_size: float,
    n_mols: List[int] = 200,
):
    copies = [round(n_mols * ratio) for ratio in ratio_list]

    topology = pack_box(
        molecules=mol_list,
        number_of_copies=copies,
        box_vectors=box_size * UNIT_CUBE * unit.nanometer,
    )

    return topology


def gas_simulation(smiles, box_size: float = DEFAULT_BOX_SIZE_GAS):
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    topology = topology_gas(mol)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    interchange.box = unit.Quantity([box_size] * 3, unit.nanometer)
    return topology, interchange


def liq_simulation(
    smiles,
    n_mols: int = DEFAULT_NUM_MOLS,
    target_density: float = DEFAULT_DENSITY,
    **kwargs,
):
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    box_size = estimate_box_size(smiles, n_mols=n_mols, target_density=target_density)

    topology = topology_single_liquid(mol, box_size=box_size, n_mols=n_mols)
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    return topology, interchange


def mix_simulation(
    smiles1,
    smiles2,
    ratio: float = 0.5,
    n_mols: int = DEFAULT_NUM_MOLS,
    target_density: float = DEFAULT_DENSITY,
):
    mol1 = Molecule.from_smiles(smiles1, allow_undefined_stereo=True)
    mol2 = Molecule.from_smiles(smiles2, allow_undefined_stereo=True)

    # estimates the size of the box given the number of heavy atoms
    box1 = estimate_box_size(smiles1, n_mols=n_mols, target_density=target_density)
    box2 = estimate_box_size(smiles2, n_mols=n_mols, target_density=target_density)
    box_size = box1 * ratio + box2 * (1 - ratio)

    topology = topology_binary(
        mol1, mol2, ratio=ratio, n_mols=n_mols, box_size=box_size
    )
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    return topology, interchange


def n_mix_simulation(
    smiles_list: List[str],
    ratio_list: List[float] = None,
    n_mols: int = DEFAULT_NUM_MOLS,
    target_density: float = DEFAULT_DENSITY,
):
    n = len(smiles_list)

    if ratio_list is None:
        ratio_list = [1 / n for _ in range(n)]

    ratio_sum = sum(ratio_list)
    ratio_list = [x / ratio_sum for x in ratio_list]

    mol_list = [
        Molecule.from_smiles(smi, allow_undefined_stereo=True) for smi in smiles_list
    ]

    # estimates the size of the box given the target density
    box_size = sum([
        x * estimate_box_size(smi, n_mols=n_mols, target_density=target_density)
        for x, smi in zip(ratio_list, smiles_list)
    ])

    topology = topology_n_mols(
        mol_list, ratio_list=ratio_list, n_mols=n_mols, box_size=box_size
    )
    ff = ForceField(DEFAULT_FF)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )
    return topology, interchange
