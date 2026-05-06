"""Single-point energy evaluation from a PDB frame.

This module reproduces the parameterization pipeline used in :mod:`autoir.create`
(``Molecule.from_smiles`` -> ``Topology`` -> ``Interchange.from_smirnoff``) but
instead of running dynamics it simply evaluates the potential energy for the
coordinates stored in a PDB frame.

OpenFF cannot infer bond orders or formal charges from a bare PDB, so the
caller must supply one SMILES string per unique molecule in the PDB via the
``smiles`` argument to :func:`single_point_energy`.
"""

from typing import List, Optional, Union

from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology, unit

import openmm
from openmm import unit as openmm_unit
from openmm.app import PDBFile
from openmm.openmm import Platform

from .create import DEFAULT_BOX_SIZE_GAS, DEFAULT_FF
from .openmm.utils import get_default_platform_props


def _load_pdb_frame(pdb_file: str, frame: int):
    """Return (PDBFile, positions for ``frame``, box vectors or None).

    ``PDBFile`` loads every ``MODEL`` section of a multi-model PDB into its
    internal position list; ``getPositions(frame=i)`` selects frame ``i``.
    The box vectors come from the single ``CRYST1`` record and therefore apply
    to every frame.
    """
    pdb = PDBFile(pdb_file)

    n_frames = pdb.getNumFrames()
    if frame < 0 or frame >= n_frames:
        raise IndexError(
            f"Frame index {frame} out of range for PDB with {n_frames} frame(s)."
        )

    positions = pdb.getPositions(asNumpy=True, frame=frame)

    box_vectors = pdb.topology.getPeriodicBoxVectors()
    if box_vectors is not None:
        # Treat an all-zero box as "no box"
        diag = [box_vectors[i][i].value_in_unit(openmm_unit.nanometer) for i in range(3)]
        if all(abs(x) < 1e-8 for x in diag):
            box_vectors = None

    return pdb, positions, box_vectors


def single_point_energy(
    pdb_file: str,
    smiles: Union[str, List[str]],
    frame: int = 0,
    force_field: str = DEFAULT_FF,
    fallback_box_nm: float = DEFAULT_BOX_SIZE_GAS,
    platform: Optional[str] = None,
    platform_props: Optional[dict] = None,
) -> float:
    """Compute the OpenMM potential energy (kJ/mol) for a single PDB frame.

    Parameters
    ----------
    pdb_file:
        Path to a PDB file. Multi-``MODEL`` files are supported.
    smiles:
        SMILES string, or a list of SMILES strings, covering every unique
        molecule present in the PDB. Required because OpenFF cannot infer
        chemistry from PDB coordinates alone.
    frame:
        Zero-based index of the ``MODEL`` to evaluate. Defaults to the first
        frame.
    force_field:
        OpenFF SMIRNOFF force-field filename. Defaults to the same force field
        used by the rest of autoir (``openff_unconstrained-2.0.0.offxml``).
    fallback_box_nm:
        Cubic box edge used when the PDB does not define a periodic box. The
        default matches :func:`autoir.create.gas_simulation`.
    platform:
        Optional OpenMM platform name (``"CUDA"``, ``"CPU"``, ...). ``None``
        lets OpenMM pick the fastest available platform.
    platform_props:
        Optional overrides merged on top of
        :func:`autoir.openmm.utils.get_default_platform_props`.

    Returns
    -------
    float
        The potential energy in kJ/mol.
    """
    if isinstance(smiles, str):
        smiles_list = [smiles]
    else:
        smiles_list = list(smiles)

    if not smiles_list:
        raise ValueError("At least one SMILES string must be provided.")

    unique_molecules = [
        Molecule.from_smiles(smi, allow_undefined_stereo=True) for smi in smiles_list
    ]

    # Read the requested frame and (optional) box from the PDB.
    _, positions, box_vectors = _load_pdb_frame(pdb_file, frame)

    # Build an OpenFF topology with chemistry taken from the reference
    # molecules and connectivity/atom order taken from the PDB.
    topology = Topology.from_pdb(pdb_file, unique_molecules=unique_molecules)

    # Parameterize.
    ff = ForceField(force_field)
    interchange: Interchange = Interchange.from_smirnoff(
        force_field=ff, topology=topology
    )

    # Periodic box: prefer CRYST1 from the PDB, fall back to a cubic box.
    if box_vectors is not None:
        a = box_vectors[0][0].value_in_unit(openmm_unit.nanometer)
        b = box_vectors[1][1].value_in_unit(openmm_unit.nanometer)
        c = box_vectors[2][2].value_in_unit(openmm_unit.nanometer)
        interchange.box = unit.Quantity([a, b, c], unit.nanometer)
    else:
        interchange.box = unit.Quantity(
            [fallback_box_nm] * 3, unit.nanometer
        )

    # Build an OpenMM simulation solely to get a Context configured with the
    # force field; the integrator is never stepped.
    integrator = openmm.LangevinIntegrator(
        300 * openmm_unit.kelvin,
        1 / openmm_unit.picosecond,
        1 * openmm_unit.femtosecond,
    )

    if platform is None:
        simulation = interchange.to_openmm_simulation(integrator)
    else:
        if not isinstance(platform, str):
            raise ValueError(f"Argument platform {platform!r} is not a string")
        props = {**get_default_platform_props(platform), **(platform_props or {})}
        simulation = interchange.to_openmm_simulation(
            integrator,
            platformProperties=props,
            platform=Platform.getPlatformByName(platform),
        )

    # Override positions (and box, if any) with the requested PDB frame.
    simulation.context.setPositions(positions)
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)

    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(openmm_unit.kilojoules_per_mole)
