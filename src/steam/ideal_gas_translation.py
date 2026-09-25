"""3D ideal gas translational thermodynamics."""

import numpy as np
from ase import units
from ase.io import read


def ideal_gas_translation_properties(contcar_file, temperature, pressure=1e5):
    """3D ideal gas translational contribution to the thermodynamic properties.

    Parameters
    ----------
    contcar_file : str
        VASP CONTCAR of the gas-phase molecule. The mass is the sum of
        all atoms in the file.
    temperature : float
        Temperature, K.
    pressure : float
        Standard-state pressure, Pa (default 1e5). Used for mu, S, A and U.

    Returns
    -------
    dict
    """
    k = units._k
    h = units._hplanck
    e = units._e
    Na = units._Nav
    R = k * Na

    T = float(temperature)
    P = float(pressure)

    atoms = read(contcar_file, format="vasp")
    m_amu = sum(atom.mass for atom in atoms)
    m = m_amu * units._amu

    # Standard chemical potential
    mu = (-1 * R * T) * (
        1.5 * np.log((2 * np.pi * m * k * T) / (h**2)) + np.log((k * T) / P)
    )

    # Entropy, Helmholtz free energy, internal energy
    q = np.power((2 * np.pi * m * k * T) / (h**2), 3 / 2) * (k * T / P)
    S = R * ((5 / 2) + np.log(q))
    A = -R * T * (1 + np.log(q))
    U = A + T * S

    per_molecule = 1.0 / (e * Na)

    return {
        "temperature": T,
        "pressure": P,
        "mass_amu": m_amu,
        "U_J_per_mol": U,
        "U_eV": U * per_molecule,
        "S_J_per_mol_K": S,
        "S_eV_per_K": S * per_molecule,
        "A_J_per_mol": A,
        "A_eV": A * per_molecule,
        "mu_J_per_mol": mu,
        "mu_eV": mu * per_molecule,
    }


def format_report(results, contcar_file=None):
    """Return the results as a human-readable text block."""
    r = results
    lines = ["Ideal gas translation thermodynamics (3D)", "=" * 64]
    if contcar_file is not None:
        lines.append(f"Structure      : {contcar_file}")
    lines += [
        f"Temperature    : {r['temperature']:.2f} K",
        f"Pressure       : {r['pressure']:.3e} Pa",
        f"Mass           : {r['mass_amu']:.4f} amu",
        "",
        f"{'Total Internal Energy (U)':<38}"
        f"{r['U_J_per_mol']:.3e} J/mol, {r['U_eV']:.3e} eV/molecule",
        f"{'Total Translational Entropy (S)':<38}"
        f"{r['S_J_per_mol_K']:.3e} J/(mol K), {r['S_eV_per_K']:.3e} eV/(molecule K)",
        f"{'Total Helmholtz Free Energy (A)':<38}"
        f"{r['A_J_per_mol']:.3e} J/mol, {r['A_eV']:.3e} eV/molecule",
        f"{'Standard Chemical Potential (mu)':<38}"
        f"{r['mu_J_per_mol']:.3e} J/mol, {r['mu_eV']:.3e} eV/molecule",
    ]
    return "\n".join(lines)


def write_report(results, output_path, contcar_file=None):
    """Write the formatted results to a text file."""
    with open(output_path, "w") as f:
        f.write(format_report(results, contcar_file=contcar_file) + "\n")
    return output_path