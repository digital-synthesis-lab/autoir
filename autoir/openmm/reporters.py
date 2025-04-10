import csv
import numpy as np
from openff.interchange import Interchange
from openmm import unit


class DipoleReporter:
    def __init__(self, file: str, report_interval: int, interchange: Interchange):
        self._out = open(file, "w", newline="")
        self._writer = csv.writer(self._out)
        self._writer.writerow(["Step", "DipoleX", "DipoleY", "DipoleZ"])
        self._report_interval = report_interval

        # Store the partial charges (in elementary charge units)
        charge_dict = interchange["Electrostatics"].charges
        num_atoms = interchange.topology.n_atoms
        self._charges = np.zeros(num_atoms)  # plain array of charges (in e)

        for key, value in charge_dict.items():
            atom_idx = key.atom_indices[0]  # TopologyKey holds a tuple like (i,)
            self._charges[atom_idx] = value.m

    def describeNextReport(self, simulation):
        steps = self._report_interval - simulation.currentStep % self._report_interval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)._value * 10  # Å
        charges_quantity = self._charges  # unit: elementary charge
        dipole = np.sum(
            charges_quantity[:, None] * positions, axis=0
        )  # shape (3,), units: e·Å

        self._writer.writerow(
            [
                simulation.context.getTime().value_in_unit(unit.femtosecond),
                dipole[0],
                dipole[1],
                dipole[2],
            ]
        )

    def __del__(self):
        self._out.close()
