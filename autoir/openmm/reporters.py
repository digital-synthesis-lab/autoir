import csv
import numpy as np
from openff.interchange import Interchange
import openmm
from openmm import unit


class DipoleReporter:
    def __init__(self, file: str, reportInterval: int, interchange: Interchange):
        self.output = file
        self._results = []
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
        charges_quantity = self._charges  # unit: elementary charge
        dipole = np.sum(
            charges_quantity[:, None] * positions, axis=0
        )  # shape (3,), units: e·Å

        self._results.append(
            [
                simulation.context.getTime().value_in_unit(unit.femtosecond),
                dipole[0],
                dipole[1],
                dipole[2],
            ]
        )

    def __del__(self):
        with open(self.output, "w", newline="") as f:
            self._writer = csv.writer(f)
            self._writer.writerow(["Step", "DipoleX", "DipoleY", "DipoleZ"])

            for line in self._results:
                self._writer.writerow(line)


class FastDataReporter:
    """StateDataReporter outputs information about a simulation, such as energy and temperature, to a file.

    To use it, create a StateDataReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """

    def __init__(
        self,
        file,
        reportInterval,
        step=False,
        time=False,
        potentialEnergy=False,
        kineticEnergy=False,
        totalEnergy=False,
        temperature=False,
        volume=False,
        density=False,
        progress=False,
        remainingTime=False,
        speed=False,
        elapsedTime=False,
        systemMass=None,
        totalSteps=None,
        append=False,
    ):
        """Create a StateDataReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to write frames
        step : bool=False
            Whether to write the current step index to the file
        time : bool=False
            Whether to write the current time to the file
        potentialEnergy : bool=False
            Whether to write the potential energy to the file
        kineticEnergy : bool=False
            Whether to write the kinetic energy to the file
        totalEnergy : bool=False
            Whether to write the total energy to the file
        temperature : bool=False
            Whether to write the instantaneous temperature to the file
        volume : bool=False
            Whether to write the periodic box volume to the file
        density : bool=False
            Whether to write the system density to the file
        progress : bool=False
            Whether to write current progress (percent completion) to the file.
            If this is True, you must also specify totalSteps.
        remainingTime : bool=False
            Whether to write an estimate of the remaining clock time until
            completion to the file.  If this is True, you must also specify
            totalSteps.
        speed : bool=False
            Whether to write an estimate of the simulation speed in ns/day to
            the file
        elapsedTime : bool=False
            Whether to write the elapsed time of the simulation in seconds to
            the file.
        separator : string=','
            The separator to use between columns in the file
        systemMass : mass=None
            The total mass to use for the system when reporting density.  If
            this is None (the default), the system mass is computed by summing
            the masses of all particles.  This parameter is useful when the
            particle masses do not reflect their actual physical mass, such as
            when some particles have had their masses set to 0 to immobilize
            them.
        totalSteps : int=None
            The total number of steps that will be included in the simulation.
            This is required if either progress or remainingTime is set to True,
            and defines how many steps will indicate 100% completion.
        append : bool=False
            If true, append to an existing file.  This has two effects.  First,
            the file is opened in append mode.  Second, the header line is not
            written, since there is assumed to already be a header line at the
            start of the file.
        """
        self._reportInterval = reportInterval
        self._results = []
        self.output = file

        if (progress or remainingTime) and totalSteps is None:
            raise ValueError(
                "Reporting progress or remaining time requires total steps to be specified"
            )

        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._separator = separator
        self._totalMass = systemMass
        self._totalSteps = totalSteps
        self._append = append
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = (
            potentialEnergy or kineticEnergy or totalEnergy or temperature
        )
        self._includes = ["energy"] if self._needEnergy else []

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        dict
            A dictionary describing the required information for the next report
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return {"steps": steps, "periodic": None, "include": self._includes}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            if not self._append:
                print(
                    '#"%s"' % ('"' + self._separator + '"').join(headers),
                    file=self._out,
                )
            self._results.append(headers)
            self._hasInitialized = True

        # Query for the values
        values = self._constructReportValues(simulation, state)
        self._results.append(values)

    def _constructReportValues(self, simulation, state):
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []

        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(state.getTime().value_in_unit(unit.picosecond))
        if self._potentialEnergy:
            values.append(
                state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            )
        if self._kineticEnergy:
            values.append(
                state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
            )
        if self._totalEnergy:
            values.append(
                (state.getKineticEnergy() + state.getPotentialEnergy()).value_in_unit(
                    unit.kilojoules_per_mole
                )
            )
        if self._temperature:
            integrator = simulation.context.getIntegrator()
            if hasattr(integrator, "computeSystemTemperature"):
                values.append(
                    integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
                )
            else:
                values.append(
                    (
                        2
                        * state.getKineticEnergy()
                        / (self._dof * unit.MOLAR_GAS_CONSTANT_R)
                    ).value_in_unit(unit.kelvin)
                )

        if self._volume or self._density:
            box = state.getPeriodicBoxVectors()
            volume = box[0][0] * box[1][1] * box[2][2]

        if self._volume:
            values.append(volume.value_in_unit(unit.nanometer**3))

        if self._density:
            values.append(
                (self._totalMass / volume).value_in_unit(
                    unit.gram / unit.item / unit.milliliter
                )
            )

        return values

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports

        Parameters
        - simulation (Simulation) The simulation to generate a report for
        """
        system = simulation.system
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0 * unit.dalton:
                    dof += 3
            for i in range(system.getNumConstraints()):
                p1, p2, distance = system.getConstraintParameters(i)
                if (
                    system.getParticleMass(p1) > 0 * unit.dalton
                    or system.getParticleMass(p2) > 0 * unit.dalton
                ):
                    dof -= 1
            if any(
                type(system.getForce(i)) == mm.CMMotionRemover
                for i in range(system.getNumForces())
            ):
                dof -= 3
            self._dof = dof
        if self._density:
            if self._totalMass is None:
                # Compute the total system mass.
                self._totalMass = 0 * unit.dalton
                for i in range(system.getNumParticles()):
                    self._totalMass += system.getParticleMass(i)
            elif not unit.is_quantity(self._totalMass):
                self._totalMass = self._totalMass * unit.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._step:
            headers.append("Step")
        if self._time:
            headers.append("Time (ps)")
        if self._potentialEnergy:
            headers.append("Potential Energy (kJ/mole)")
        if self._kineticEnergy:
            headers.append("Kinetic Energy (kJ/mole)")
        if self._totalEnergy:
            headers.append("Total Energy (kJ/mole)")
        if self._temperature:
            headers.append("Temperature (K)")
        if self._volume:
            headers.append("Box Volume (nm^3)")
        if self._density:
            headers.append("Density (g/mL)")
        return headers

    def __del__(self):
        with open(self.output, "w", newline="") as f:
            self._writer = csv.writer(f)
            for line in self._results:
                self._writer.writerow(line)
