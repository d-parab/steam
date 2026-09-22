"""Harmonic oscillator thermodynamics for surface intermediates."""

import numpy as np
from ase import units

from steam.extract_frequencies import convert_frequencies_to_hz


def harmonic_properties(file_path, temperature):
    """Zero-point energy, vibrational entropy, Helmholtz free energy,
    standard chemical potential and internal energy.

    The chemical potential contribution mu is identical to the Helmholtz
    free energy A; it is reported separately for convenience.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Text file with one frequency per line, in cm^-1.
    temperature : float
        Temperature in K.

    Returns
    -------
    dict
        ZPE, S, A, mu and U in J-based and eV-based units.
    """
    freq = convert_frequencies_to_hz(file_path)

    h = units._hplanck   # J s
    k = units._k         # J/K
    Na = units._Nav      # 1/mol
    e = units._e         # C
    R = k * Na           # J/(mol K)

    T = float(temperature)
    beta = 1.0 / (k * T)

    x = h * freq * beta
    exp_mx = np.exp(-x)

    zpe = np.sum(0.5 * h * freq) * Na
    S = R * np.sum(x * exp_mx / (1.0 - exp_mx) - np.log(1.0 - exp_mx))
    A = R * T * np.sum(0.5 * x + np.log(1.0 - exp_mx))
    U = A + T * S

    per_molecule = 1.0 / (e * Na)

    return {
        "temperature": T,
        "n_modes": len(freq),
        "ZPE_J_per_mol": zpe,
        "ZPE_eV": zpe * per_molecule,
        "S_J_per_mol_K": S,
        "S_eV_per_K": S * per_molecule,
        "A_J_per_mol": A,
        "A_eV": A * per_molecule,
        "mu_J_per_mol": A,
        "mu_eV": A * per_molecule,
        "U_J_per_mol": U,
        "U_eV": U * per_molecule,
    }


def format_report(results, file_path=None):
    """Return the results as a human-readable text block."""
    r = results
    lines = ["Harmonic oscillator thermodynamics", "=" * 64]
    if file_path is not None:
        lines.append(f"Frequency file : {file_path}")
    lines += [
        f"Temperature    : {r['temperature']:.2f} K",
        f"Modes used     : {r['n_modes']}",
        "",
        f"{'Zero-Point Energy (ZPE)':<38}"
        f"{r['ZPE_J_per_mol']:.3e} J/mol, {r['ZPE_eV']:.3e} eV/molecule",
        f"{'Total Vibrational Entropy (S)':<38}"
        f"{r['S_J_per_mol_K']:.3e} J/(mol K), {r['S_eV_per_K']:.3e} eV/(molecule K)",
        f"{'Total Helmholtz Free Energy (A)':<38}"
        f"{r['A_J_per_mol']:.3e} J/mol, {r['A_eV']:.3e} eV/molecule",
        f"{'Standard Chemical Potential (mu)':<38}"
        f"{r['mu_J_per_mol']:.3e} J/mol, {r['mu_eV']:.3e} eV/molecule",
        f"{'Total Internal Energy (U)':<38}"
        f"{r['U_J_per_mol']:.3e} J/mol, {r['U_eV']:.3e} eV/molecule",
    ]
    return "\n".join(lines)


def write_report(results, output_path, file_path=None):
    """Write the formatted results to a text file."""
    with open(output_path, "w") as f:
        f.write(format_report(results, file_path=file_path) + "\n")
    return output_path