import os
import numpy as np
import pandas as pd
import openmm

from .reporters import FastDataReporter


def resize_box(simulation, log_file: str = "data.csv", last_n: int = 2000):
    # load the data
    df = pd.read_csv(log_file)
    mean_vol = df.iloc[-last_n:]["Box Volume (nm^3)"].mean()
    mean_density = df.iloc[-last_n:]["Density (g/mL)"].mean()
    mean_a = np.power(mean_vol, 1 / 3)
    lattice = np.eye(3) * mean_a
    simulation.context.setPeriodicBoxVectors(*lattice)
    simulation.context.reinitialize(True)
    return simulation, mean_vol, mean_density


def deactivate_barostat(simulation):
    for i, f in enumerate(simulation.system.getForces()):
        if isinstance(f, openmm.MonteCarloBarostat):
            f.setFrequency(0)
            simulation.system.removeForce(i)

    return simulation


def deactivate_data_reporters(simulation):
    to_remove = []
    for i, f in enumerate(simulation.reporters):
        if isinstance(f, openmm.app.StateDataReporter):
            to_remove.append(i)

        elif isinstance(f, FastDataReporter):
            to_remove.append(i)

    for i in sorted(to_remove, reverse=True):
        simulation.reporters.pop(i)

    return simulation


def get_default_platform_props(platform: str) -> dict:
    if platform is None:
        return {}

    if platform.lower() == "cuda":
        devices = os.environ.get("CUDA_VISIBLE_DEVICES", None)
        if devices is None:
            return {}

        if "," in devices:
            return {
                "DeviceIndex": devices,
            }

        return {
            "DeviceIndex": int(devices),
        }
