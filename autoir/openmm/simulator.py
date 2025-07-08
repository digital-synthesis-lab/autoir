import logging
import os

from openff.interchange import Interchange

import openmm
from openmm import unit
from openmm.openmm import Platform

from .reporters import AverageEnergyReporter, DipoleReporter, FastDataReporter
from .utils import (
    deactivate_barostat,
    deactivate_data_reporters,
    get_default_platform_props,
    resize_box,
)

DEFAULT_PROD_FILE = "prod.csv"
DEFAULT_EQUI_FILE = "equi.csv"
DEFAULT_AVGE_FILE = "avgE.csv"
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
        equi_file: str = DEFAULT_EQUI_FILE,
        prod_file: str = DEFAULT_PROD_FILE,
        dipole_file: str = DEFAULT_DIPOLES_FILE,
        platform: str = None,
        platform_props: dict = None,
    ):
        self.time_step = time_step * unit.femtoseconds  # simulation timestep
        self.temperature = temperature * unit.kelvin  # simulation temperature
        self.friction = friction / unit.picosecond  # collision rate

        self.pressure = pressure * unit.atmosphere
        self.barostat_freq = barostat_freq
        self.trj_freq = trj_freq
        self.log_freq = log_freq

        self.trj_file = trj_file
        self.equi_file = equi_file
        self.prod_file = prod_file
        self.dipole_file = dipole_file

        self.platform = platform
        if type(platform_props) == dict:
            self.platform_props = {
                **get_default_platform_props(platform),
                **platform_props,
            }
        else:
            self.platform_props = get_default_platform_props(platform)

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

    def create_simulation(self, interchange: Interchange, barostat: bool = True):
        if self.platform is None:
            simulation = interchange.to_openmm_simulation(self.get_integrator())
        elif type(self.platform) != str:
            raise ValueError(f"Argument platform {self.platform} is not string")
        else:
            simulation = interchange.to_openmm_simulation(
                self.get_integrator(),
                platformProperties=self.platform_props,
                platform=Platform.getPlatformByName(self.platform),
            )

        if barostat:
            simulation.system.addForce(self.get_barostat())
            simulation.context.reinitialize(True)

        simulation.context.setVelocitiesToTemperature(self.temperature)
        return simulation

    def get_integrator(self):
        return openmm.LangevinIntegrator(
            self.temperature, self.friction, self.time_step
        )

    def get_barostat(self):
        return openmm.MonteCarloBarostat(
            self.pressure, self.temperature, self.barostat_freq
        )

    def get_traj_reporter(self):
        return openmm.app.PDBReporter(
            self.trj_file, self.trj_freq, enforcePeriodicBox=True
        )

    def get_dipole_reporter(self, interchange: Interchange):
        return DipoleReporter(
            self.dipole_file, reportInterval=1, interchange=interchange
        )

    def run(
        self,
        interchange: Interchange,
        npt_equi_steps=100_000,
        nvt_equi_steps=100_000,
        nvt_prod_steps=1_000_000,
        npt_equi_volume_steps=1000,
    ):
        self.logger.info("Creating simulation")
        has_npt = npt_equi_steps > 0
        simulation = self.create_simulation(interchange, barostat=has_npt)
        simulation.minimizeEnergy(maxIterations=100)
        simulation.reporters.append(
            FastDataReporter(
                self.equi_file,
                reportInterval=self.log_freq,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                density=True,
                volume=True,
            )
        )

        if has_npt:
            self.logger.info(f"NPT equilibration for {npt_equi_steps} steps")
            simulation.step(npt_equi_steps)
            del simulation.reporters[-1]
            deactivate_barostat(simulation)
            resize_box(
                simulation, log_file=self.equi_file, last_n=npt_equi_volume_steps
            )
            state = simulation.context.getState()
            simulation.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())

        self.logger.info(f"NVT equilibration for {nvt_equi_steps} steps")
        simulation.step(nvt_equi_steps)

        self.logger.info(f"NVT production for {nvt_prod_steps} steps")
        deactivate_data_reporters(simulation)
        simulation.reporters.append(self.get_traj_reporter())
        simulation.reporters.append(self.get_dipole_reporter(interchange))
        simulation.reporters.append(
            AverageEnergyReporter(
                "avgE.csv", reportInterval=1, startingStep=nvt_prod_steps // 2
            )
        )
        simulation.reporters.append(
            FastDataReporter(
                self.prod_file,
                reportInterval=self.log_freq,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=False,
                density=False,
                volume=False,
            )
        )
        simulation.step(nvt_prod_steps)

        self.logger.info(f"Production simulation done")
