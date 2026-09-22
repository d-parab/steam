"""Reading vibrational frequencies and converting them to Hz."""

import numpy as np
from ase import units


def convert_frequencies_to_hz(file_path):
    """Read frequency values in cm^-1 from a text file and convert them to Hz.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to a text file with one frequency per line, in cm^-1.

    Returns
    -------
    numpy.ndarray
        Frequencies in Hz.
    """
    c = units._c  # speed of light in m/s

    with open(file_path, 'r') as file:
        freq_cm_inv = np.array([float(line.strip()) for line in file])

    freq_Hz = c * 100 * freq_cm_inv

    return freq_Hz