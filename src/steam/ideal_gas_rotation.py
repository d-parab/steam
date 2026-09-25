"""Ideal gas rotational thermodynamics (linear and non-linear molecules)."""

import numpy as np
from ase import units
from ase.io import read


def _constants():
    k = units._k
    h = units._hplanck
    Na = units._Nav
    e = units._e
    return k, h, Na, e, k * Na


def ideal_gas_rotation_nonlinear_properties(contcar_file, symmetry_number,
                                            temperature):
    """Rotational contribution for a non-linear ideal gas molecule.

    Parameters
    ----------
    contcar_file : str
        VASP CONTCAR of the gas-phase molecule.
    symmetry_number : float
        Rotational symmetry number.
    temperature : float
        Temperature, K.

    Returns
    -------
    dict
    """
    k, h, Na, e, R = _constants()

    T = float(temperature)
    sigma = float(symmetry_number)

    atoms = read(contcar_file, format="vasp")
    I_amu_A2 = atoms.get_moments_of_inertia()
    Ix, Iy, Iz = I_amu_A2 * units._amu * 1e-20   # kg m^2

    # Rotational partition function
    beta = 1 / (k * T)
    h_ = h / (2 * np.pi)
    q_rot = (np.sqrt(np.pi) / sigma) * np.sqrt((2 * Ix) / (beta * h_**2)) * \
            np.sqrt((2 * Iy) / (beta * h_**2)) * np.sqrt((2 * Iz) / (beta * h_**2))

    # Standard chemical potential
    mu = (-1 * R * T) * np.log(q_rot)

    # Entropy, Helmholtz free energy, internal energy
    S = R * ((3 / 2) + np.log(q_rot))
    A = -R * T * (np.log(q_rot))
    U = A + T * S

    per_molecule = 1.0 / (e * Na)

    return {
        "geometry": "non-linear",
        "temperature": T,
        "symmetry_number": sigma,
        "moments_of_inertia_amu_A2": I_amu_A2,
        "U_J_per_mol": U,
        "U_eV": U * per_molecule,
        "S_J_per_mol_K": S,
        "S_eV_per_K": S * per_molecule,
        "A_J_per_mol": A,
        "A_eV": A * per_molecule,
        "mu_J_per_mol": mu,
        "mu_eV": mu * per_molecule,
    }


def _linear_moment(I_amu_A2, tol):
    """Return the non-zero, equal moment of inertia (amu A^2) of a linear
    molecule. Raises ValueError if the moments are not those of a linear
    molecule within the relative tolerance tol."""
    I = np.asarray(I_amu_A2, dtype=float)
    I_max = I.max()
    for i in range(3):
        others = [I[j] for j in range(3) if j != i]
        if I[i] <= tol * I_max and np.isclose(others[0], others[1], rtol=tol):
            return others[0]
    raise ValueError(
        f"WARNING: Check moments of inertia {I.tolist()} amu A^2 - "
        "not a linear molecule within tolerance. Make sure the linear "
        "molecule is aligned along one of the unit cell axes (x, y or z)."
    )


def ideal_gas_rotation_linear_properties(contcar_file, symmetry_number,
                                         temperature, tol=1e-3):
    """Rotational contribution for a linear ideal gas molecule.

    Parameters
    ----------
    contcar_file : str
        VASP CONTCAR of the gas-phase molecule.
    symmetry_number : float
        Rotational symmetry number.
    temperature : float
        Temperature, K.
    tol : float
        Relative tolerance for identifying the zero moment and the two
        equal moments (default 1e-3).

    Returns
    -------
    dict
    """
    k, h, Na, e, R = _constants()

    T = float(temperature)
    sigma = float(symmetry_number)

    atoms = read(contcar_file, format="vasp")
    I_amu_A2 = atoms.get_moments_of_inertia()
    I_used_amu_A2 = _linear_moment(I_amu_A2, tol)
    I_ = I_used_amu_A2 * units._amu * 1e-20       # kg m^2

    # Rotational partition function
    beta = 1 / (k * T)
    h_ = h / (2 * np.pi)
    q_rot = (2 * I_) / (sigma * beta * h_**2)

    # Standard chemical potential
    mu = (-1 * R * T) * np.log(q_rot)

    # Entropy, Helmholtz free energy, internal energy
    S = R * (1 + np.log(q_rot))
    A = -R * T * (np.log(q_rot))
    U = A + T * S

    per_molecule = 1.0 / (e * Na)

    return {
        "geometry": "linear",
        "temperature": T,
        "symmetry_number": sigma,
        "moments_of_inertia_amu_A2": I_amu_A2,
        "moment_used_amu_A2": I_used_amu_A2,
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
    Ix, Iy, Iz = r["moments_of_inertia_amu_A2"]
    lines = [f"Ideal gas rotation thermodynamics ({r['geometry']})", "=" * 64]
    if contcar_file is not None:
        lines.append(f"Structure        : {contcar_file}")
    lines += [
        f"Temperature      : {r['temperature']:.2f} K",
        f"Symmetry number  : {r['symmetry_number']:g}",
        f"Moments (amu A^2): {Ix:.4f}, {Iy:.4f}, {Iz:.4f}",
    ]
    if r["geometry"] == "linear":
        lines.append(f"Moment used      : {r['moment_used_amu_A2']:.4f} amu A^2")
    lines += [
        "",
        f"{'Total Internal Energy (U)':<38}"
        f"{r['U_J_per_mol']:.3e} J/mol, {r['U_eV']:.3e} eV/molecule",
        f"{'Total Rotational Entropy (S)':<38}"
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