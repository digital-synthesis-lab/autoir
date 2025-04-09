import time
import numpy as np
import csv
import pandas as pd
import openmm
from openmm import unit

from autoir.create import mix_simulation

_, interchange = mix_simulation("C1COCO1", "NCCN", n_mols=100, scaling=0.002)


class DipoleReporter:
    def __init__(self, file, reportInterval, interchange):
        self._out = open(file, "w", newline="")
        self._writer = csv.writer(self._out)
        self._writer.writerow(
            ["Step", "Dipole_X_Debye", "Dipole_Y_Debye", "Dipole_Z_Debye"]
        )
        self._reportInterval = reportInterval

        # Store the partial charges (in elementary charge units)
        charge_dict = interchange["Electrostatics"].charges
        num_atoms = interchange.topology.n_atoms
        self._charges = np.zeros(num_atoms)  # plain array of charges (in e)

        for key, value in charge_dict.items():
            atom_idx = key.atom_indices[0]  # TopologyKey holds a tuple like (i,)
            self._charges[atom_idx] = value.m

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)._value * 10  # Å
        # Compute dipole vector in e·nm
        # dipole = sum(q * pos for q, pos in zip(self._charges, positions))
        charges_quantity = self._charges  # unit: charge
        dipole = np.sum(
            charges_quantity[:, None] * positions, axis=0
        )  # shape (3,), units: e·nm

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


# Propagate the System with Langevin dynamics.
time_step = 2 * unit.femtoseconds  # simulation timestep
temperature = 300 * unit.kelvin  # simulation temperature
friction = 1 / unit.picosecond  # collision rate
integrator = openmm.LangevinIntegrator(temperature, friction, time_step)

# Define pressure and barostat frequency
pressure = 1.0 * unit.atmosphere
barostat_frequency = 100  # in steps


# Logging options.
trj_freq = 1000  # number of steps per written trajectory frame
data_freq = 100  # number of steps per written simulation statistics

# Set up an OpenMM simulation.
simulation = interchange.to_openmm_simulation(integrator)
# simulation = openmm.app.Simulation(interchange.topology.to_openmm(), system, integrator)
simulation.system.addForce(
    openmm.MonteCarloBarostat(pressure, temperature, barostat_frequency)
)
simulation.context.reinitialize(True)

# Randomize the velocities from a Boltzmann distribution at a given temperature.
simulation.context.setVelocitiesToTemperature(temperature)

# Configure the information in the output files.
state_data_reporter = openmm.app.StateDataReporter(
    "data.csv",
    data_freq,
    step=True,
    potentialEnergy=True,
    temperature=True,
    density=True,
    volume=True,
)
simulation.reporters.append(state_data_reporter)

# Length of the simulation.
num_steps = 20000  # number of integration steps to run


print("Initial equilibration")
start = time.process_time()

# Run the simulation
simulation.step(num_steps)

end = time.process_time()
print("Elapsed time %.2f seconds" % (end - start))
print("Done!")

# simulation.context.reinitialize(True)
for i, f in enumerate(simulation.system.getForces()):
    if isinstance(f, openmm.MonteCarloBarostat):
        f.setFrequency(0)
        simulation.system.removeForce(i)


# load the data
df = pd.read_csv("data.csv")
mean_vol = df.iloc[-100:]["Box Volume (nm^3)"].mean()
mean_a = np.power(mean_vol, 1 / 3)
lattice = np.eye(3) * mean_a
print(mean_vol, mean_a)
simulation.context.setPeriodicBoxVectors(*lattice)
simulation.context.reinitialize(True)

# Length of the simulation.
num_steps = 40000  # number of integration steps to run

print("NVT Equilibration")
start = time.process_time()

# Run the simulation
simulation.step(num_steps)

end = time.process_time()
print("Elapsed time %.2f seconds" % (end - start))


## Production
print("Production")
pdb_reporter = openmm.app.PDBReporter(
    "trajectory.pdb", trj_freq, enforcePeriodicBox=True
)
simulation.reporters.append(pdb_reporter)
dipole_reporter = DipoleReporter(
    "dipoles.csv", reportInterval=1, interchange=interchange
)
simulation.reporters.append(dipole_reporter)


num_steps = 400000  # number of integration steps to run
start = time.process_time()

# Run the simulation
simulation.step(num_steps)

end = time.process_time()
print("Elapsed time %.2f seconds" % (end - start))
print("Done!")
