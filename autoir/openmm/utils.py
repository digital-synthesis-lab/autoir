import numpy as np
import pandas as pd
import openmm


def resize_box(simulation, log_file: str = "data.csv", last_n: int = 2000):
    # load the data
    df = pd.read_csv(log_file)
    mean_vol = df.iloc[-last_n:]["Box Volume (nm^3)"].mean()
    mean_a = np.power(mean_vol, 1 / 3)
    lattice = np.eye(3) * mean_a
    simulation.context.setPeriodicBoxVectors(*lattice)
    simulation.context.reinitialize(True)
    return simulation


def deactivate_barostat(simulation):
    for i, f in enumerate(simulation.system.getForces()):
        if isinstance(f, openmm.MonteCarloBarostat):
            f.setFrequency(0)
            simulation.system.removeForce(i)

    return simulation
