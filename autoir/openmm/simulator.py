from openff.interchange import Interchange

import openmm
from openmm import unit

from .reporters import DipoleReporter
from .utils import deactivate_barostat, resize_box

DEFAULT_LOG_FILE = "data.csv"
DEFAULT_DIPOLES_FILE = "dipoles.csv"
DEFAULT_TRAJ_FILE = "trajectory.pdb"

DEFAULT_TIME_STEP = 2  # fs
DEFAULT_TEMPERATURE = 300  # K
DEFAULT_FRICTION = 1  # 1/ps
DEFAULT_PRESSURE = 1  # atm
DEFAULT_BAROSTAT_FREQ = 100  # steps
DEFAULT_LOG_FREQ = 20  # steps
DEFAULT_TRJ_FREQ = 1000  # steps


class OpenMMSimulator:
    def __init__(
        self,
        time_step: float = DEFAULT_TIME_STEP,
        temperature: float = DEFAULT_TEMPERATURE,
        friction: float = DEFAULT_FRICTION,
        pressure: float = DEFAULT_PRESSURE,
        barostat_freq: float = DEFAULT_BAROSTAT_FREQ,
        log_freq: float = DEFAULT_LOG_FREQ,
        trj_freq: float = DEFAULT_TRJ_FREQ,
        trj_file: str = DEFAULT_TRAJ_FILE,
        log_file: str = DEFAULT_LOG_FILE,
        dipole_file: str = DEFAULT_DIPOLES_FILE,
    ):
        self.time_step = time_step * unit.femtoseconds  # simulation timestep
        self.temperature = temperature * unit.kelvin  # simulation temperature
        self.friction = friction / unit.picosecond  # collision rate

        self.pressure = pressure * unit.atmosphere
        self.barostat_freq = barostat_freq
        self.trj_freq = trj_freq
        self.log_freq = log_freq

        self.trj_file = trj_file
        self.log_file = log_file
        self.dipole_file = dipole_file

        self.logger = self.get_logger()

    def get_logger(self):
        logger = logging.getLogger("AutoIR")
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "AutoIR - %(asctime)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)
        logger.handlers = []
        logger.addHandler(handler)
        return logger

    def create_simulation(self, interchange: Interchange):
        simulation = interchange.to_openmm_simulation(self.get_integrator())
        simulation.system.addForce(self.get_barostat())
        simulation.context.reinitialize(True)
        simulation.context.setVelocitiesToTemperature(self.temperature)
        simulation.reporters.append(self.get_data_reporter())
        return simulation

    def get_integrator(self):
        return openmm.LangevinIntegrator(
            self.temperature, self.friction, self.time_step
        )

    def get_barostat(self):
        return openmm.MonteCarloBarostat(
            self.pressure, self.temperature, self.barostat_freq
        )

    def get_data_reporter(self):
        return openmm.app.StateDataReporter(
            self.log_file,
            self.log_freq,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            volume=True,
        )

    def get_traj_reporter(self):
        return openmm.app.PDBReporter(
            self.trj_file, self.trj_freq, enforcePeriodicBox=True
        )

    def get_dipole_reporter(self, interchange: Interchange):
        return DipoleReporter(
            self.dipoles_file, reportInterval=1, interchange=interchange
        )

    def run(
        self,
        interchange: Interchange,
        npt_equi_steps=100_000,
        nvt_equi_steps=100_000,
        nvt_prod_steps=1_000_000,
        npt_equi_volume_steps=200,
    ):
        self.logger.info("Creating simulation")
        simulation = self.get_simulation(interchange)
        self.logger.info("NPT equilibration")
        simulation.step(npt_equi_steps)

        self.logger.info("NVT equilibration")
        deactivate_barostat(simulation)
        resize_box(simulation, log_file=self.log_file, last_n=npt_equi_volume_steps)
        simulation.step(nvt_equi_steps)

        self.logger.info("NVT Production")
        simulation.reporters.append(self.get_traj_reporter())
        simulation.reporters.append(self.get_dipole_reporter())
        simulation.step(nvt_prod_steps)
