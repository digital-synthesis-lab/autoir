import numpy as np
import pandas as pd
from scipy import fftpack, signal

from autoir.openmm.simulator import DEFAULT_DIPOLES_FILE

from . import const

DEFAULT_TRUNCATE = 20_000


def compute_autocorr(
    filename: str = DEFAULT_DIPOLES_FILE, truncate_autocorr: int = DEFAULT_TRUNCATE
):
    # Load data
    time, mu_x, mu_y, mu_z = np.loadtxt(
        filename, skiprows=1, unpack=True, delimiter=","
    )
    n = len(mu_x)

    # Shift the array if odd to perform the convolution
    shift = n % 2
    mu_x_shifted = np.zeros(n * 2 - shift)
    mu_y_shifted = np.zeros(n * 2 - shift)
    mu_z_shifted = np.zeros(n * 2 - shift)

    mu_x_shifted[n // 2 : n // 2 + n] = mu_x
    mu_y_shifted[n // 2 : n // 2 + n] = mu_y
    mu_z_shifted[n // 2 : n // 2 + n] = mu_z

    # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
    autocorr_x_full = signal.fftconvolve(mu_x_shifted, mu_x[::-1], mode="same")
    autocorr_y_full = signal.fftconvolve(mu_y_shifted, mu_y[::-1], mode="same")
    autocorr_z_full = signal.fftconvolve(mu_z_shifted, mu_z[::-1], mode="same")

    norm = np.arange(n, 0, -1)
    autocorr_x_full = autocorr_x_full[-n:] / norm
    autocorr_y_full = autocorr_y_full[-n:] / norm
    autocorr_z_full = autocorr_z_full[-n:] / norm

    autocorr_full = autocorr_x_full + autocorr_y_full + autocorr_z_full
    autocorr = autocorr_full[:truncate_autocorr]

    timestep = (
        time[1] - time[0]
    ) * 1.0e-15  # converts time from femtoseconds to seconds

    return autocorr, timestep


def compute_spectra(autocorr, timestep, T=300):
    # Calculate the FFTs of autocorrelation functions
    intensities = fftpack.dct(autocorr, type=1)[1:]
    frequencies = np.linspace(0, 0.5 / timestep, len(autocorr))[1:]
    wavenums = frequencies / (100.0 * const.LIGHTSPEED)

    # Calculate spectra
    correction = 1.0 - np.exp(-const.HBAR * frequencies / (const.kB * T))

    field_description = frequencies * correction
    quantum_correction = frequencies / correction
    spectra = intensities * field_description
    spectra_qm = spectra * quantum_correction
    return wavenums, spectra_qm


def get_spectra_from_file(
    inp_file: str = DEFAULT_DIPOLES_FILE, truncate_autocorr: int = DEFAULT_TRUNCATE
):
    autocorr, timestep = compute_autocorr(inp_file)
    wavenums, spectra = compute_spectra(autocorr, timestep)

    df = pd.DataFrame({"w": wavenums, "IR": spectra})
    return df


def process_file(
    inp_file: str = DEFAULT_DIPOLES_FILE,
    out_file: str = "ir.csv",
    truncate_autocorr: int = DEFAULT_TRUNCATE,
):
    df = get_spectra_from_file(inp_file)
    df.to_csv(out_file)
