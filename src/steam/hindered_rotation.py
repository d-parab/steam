"""Hindered rotor thermodynamics for surface intermediates."""

import numpy as np
from ase import units
from ase.io import read
from scipy.special import iv

CONSTRAINT_WARNING = "WARNING: make sure all atoms are constrained except the adsorbate."


def _moment_of_inertia(atoms, atom_index):
    """In-plane moment of inertia (amu A^2) of the unconstrained atoms."""
    masses = atoms.get_masses()
    fixed = [i for c in atoms.constraints for i in c.index]
    x_c, y_c, _ = atoms[atom_index].position

    I_calc = 0.0
    for i, atom in enumerate(atoms):
        if i not in fixed:
            x, y, _ = atom.position
            I_calc += masses[i] * ((x - x_c) ** 2 + (y - y_c) ** 2)
    return I_calc


def barrier_to_frequency_rot(barrier_eV, contcar_file, atom_index, n_minima):
    """Convert a hindered rotation barrier (eV) to a frequency in cm^-1."""
    e = units._e
    c = units._c

    Wr = barrier_eV * e
    structure = read(contcar_file, format="vasp")
    I = _moment_of_inertia(structure, atom_index) * units._amu * 1e-20

    nu_r = (1 / (2 * np.pi)) * np.sqrt(n_minima**2 * Wr / (2 * I))

    return nu_r / (c * 100)


def hindered_rotation_properties(barrier_eV, contcar_file, atom_index,
                                 n_minima, symmetry_number, temperature):
    """Hindered rotation contribution to the thermodynamic properties.

    Parameters
    ----------
    barrier_eV : float
        Rotational barrier, eV.
    contcar_file : str
        VASP CONTCAR with all atoms constrained except the adsorbate.
    atom_index : int
        Index of the atom about which the molecule rotates.
    n_minima : int
        Number of equivalent minima in a full rotation.
    symmetry_number : int
        Symmetry number.
    temperature : float
        Temperature, K.

    Returns
    -------
    dict
    """
    h = units._hplanck
    k = units._k
    Na = units._Nav
    e = units._e
    c = units._c
    R = k * Na

    T = float(temperature)
    W = barrier_eV * e
    sigma = symmetry_number

    nu_cm = barrier_to_frequency_rot(barrier_eV, contcar_file, atom_index,
                                     n_minima)
    nu = nu_cm * (c * 100)

    structure = read(contcar_file, format="vasp")
    I_amu_A2 = _moment_of_inertia(structure, atom_index)

    ri = W / (h * nu)
    Ti = (k * T) / (h * nu)

    # Chemical potential
    A_const = np.pi * W / k
    B_const = ((h * nu) / k) * (0.5 * ri + 0.5 - 1 / (2 + 16 * ri))
    C_const = W / (2 * k)
    D_const = (h * nu) / k

    mu = -R * T * (
        0.5 * np.log(A_const / T)
        - B_const / T
        + np.log(iv(0, C_const / T))
        - np.log(sigma)
        - np.log(1 - np.exp(-D_const / T))
    )

    # Entropy
    bessel0 = iv(0, ri / (2 * Ti))
    bessel1 = iv(1, ri / (2 * Ti))

    S_HO = R * ((1 / Ti) / (np.exp(1 / Ti) - 1) - np.log(1 - np.exp(-1 / Ti)))
    delS = R * (
        -0.5
        - (ri * bessel1) / (2 * Ti * bessel0)
        + np.log(np.sqrt(np.pi * ri / Ti) * bessel0)
    )
    S = S_HO + delS

    # Helmholtz free energy
    A_HO = R * T * (1 / (2 * Ti) + np.log(1 - np.exp(-1 / Ti)))
    delA = R * T * (
        -1 / ((2 + 16 * ri) * Ti)
        + ri / (2 * Ti)
        - np.log(np.sqrt(np.pi * ri / Ti) * bessel0)
    )
    A = A_HO + delA

    # Internal energy
    U = A + T * S

    per_molecule = 1.0 / (e * Na)

    return {
        "temperature": T,
        "barrier_eV": barrier_eV,
        "atom_index": atom_index,
        "n_minima": n_minima,
        "symmetry_number": sigma,
        "moment_of_inertia_amu_A2": I_amu_A2,
        "nu_cm_inv": nu_cm,
        "mu_J_per_mol": mu,
        "mu_eV": mu * per_molecule,
        "S_HO_J_per_mol_K": S_HO,
        "S_HO_eV_per_K": S_HO * per_molecule,
        "delS_J_per_mol_K": delS,
        "delS_eV_per_K": delS * per_molecule,
        "S_J_per_mol_K": S,
        "S_eV_per_K": S * per_molecule,
        "A_HO_J_per_mol": A_HO,
        "A_HO_eV": A_HO * per_molecule,
        "delA_J_per_mol": delA,
        "delA_eV": delA * per_molecule,
        "A_J_per_mol": A,
        "A_eV": A * per_molecule,
        "U_J_per_mol": U,
        "U_eV": U * per_molecule,
    }


def format_report(r, contcar_file=None):
    """Return the results as a human-readable text block."""
    lines = ["Hindered rotation thermodynamics", "=" * 70,
             CONSTRAINT_WARNING, ""]
    if contcar_file is not None:
        lines.append(f"Structure        : {contcar_file}")
    lines += [
        f"Temperature      : {r['temperature']:.2f} K",
        f"Barrier          : {r['barrier_eV']:.4f} eV",
        f"Rotation axis    : atom {r['atom_index']}",
        f"Equivalent minima: {r['n_minima']}",
        f"Symmetry number  : {r['symmetry_number']}",
        f"Moment of inertia: {r['moment_of_inertia_amu_A2']:.4f} amu A^2",
        f"Frequency        : {r['nu_cm_inv']:.2f} cm^-1",
        "",
        f"{'Chemical potential (mu)':<32}"
        f"{r['mu_J_per_mol']:.3e} J/mol  |  {r['mu_eV']:.3e} eV/molecule",
        "",
        f"{'Vibrational entropy (S_HO)':<32}"
        f"{r['S_HO_J_per_mol_K']:.3f} J/(mol K)  |  {r['S_HO_eV_per_K']:.3e} eV/(molecule K)",
        f"{'Additional entropy (delS)':<32}"
        f"{r['delS_J_per_mol_K']:.3f} J/(mol K)  |  {r['delS_eV_per_K']:.3e} eV/(molecule K)",
        f"{'Total entropy (S)':<32}"
        f"{r['S_J_per_mol_K']:.3f} J/(mol K)  |  {r['S_eV_per_K']:.3e} eV/(molecule K)",
        "",
        f"{'Vibrational Helmholtz (A_HO)':<32}"
        f"{r['A_HO_J_per_mol']:.3f} J/mol  |  {r['A_HO_eV']:.3e} eV/molecule",
        f"{'Additional Helmholtz (delA)':<32}"
        f"{r['delA_J_per_mol']:.3f} J/mol  |  {r['delA_eV']:.3e} eV/molecule",
        f"{'Total Helmholtz (A)':<32}"
        f"{r['A_J_per_mol']:.3f} J/mol  |  {r['A_eV']:.3e} eV/molecule",
        "",
        f"{'Internal energy (U)':<32}"
        f"{r['U_J_per_mol']:.3f} J/mol  |  {r['U_eV']:.3e} eV/molecule",
    ]
    return "\n".join(lines)


def write_report(r, output_path, contcar_file=None):
    """Write the formatted results to a text file."""
    with open(output_path, "w") as f:
        f.write(format_report(r, contcar_file=contcar_file) + "\n")
    return output_path