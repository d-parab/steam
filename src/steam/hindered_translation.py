"""Hindered translator thermodynamics for surface intermediates."""

import numpy as np
from ase import units
from ase.io import read
from scipy.special import iv

CONSTRAINT_WARNING = "WARNING: make sure all atoms are constrained except the adsorbate."
GEOMETRY_WARNING = "WARNING: hard coded for (111) geometry."


def _adsorbate_mass(atoms):
    """Total mass (amu) of the unconstrained atoms."""
    masses = atoms.get_masses()
    fixed = [i for c in atoms.constraints for i in c.index]
    return sum(masses[i] for i in range(len(masses)) if i not in fixed)


def barrier_to_frequency_trans(barrier_eV, contcar_file, nn_distance):
    """Convert a hindered translation barrier (eV) to a frequency in cm^-1."""
    e = units._e
    c = units._c

    W_x = barrier_eV * e
    atoms = read(contcar_file, format="vasp")
    mass = _adsorbate_mass(atoms) * units._amu

    b = nn_distance * 1e-10
    as_factor = np.sqrt(3) * b**2 / 2

    nu_x = (W_x / (2 * mass * as_factor)) ** 0.5

    return nu_x / (c * 100)


def hindered_translation_properties(barrier_eV, contcar_file, nn_distance,
                                    temperature):
    """Hindered translation contribution to the thermodynamic properties.

    Parameters
    ----------
    barrier_eV : float
        Activation barrier for translation, eV.
    contcar_file : str
        VASP CONTCAR with all atoms constrained except the adsorbate.
    nn_distance : float
        Nearest neighbour distance, angstrom.
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

    nu_cm = barrier_to_frequency_trans(barrier_eV, contcar_file, nn_distance)
    nu = nu_cm * (c * 100)

    atoms = read(contcar_file, format="vasp")
    mass_amu = _adsorbate_mass(atoms)

    ri = W / (h * nu)
    Ti = (k * T) / (h * nu)

    # Chemical potential
    A_const = np.pi * W / k
    B_const = (ri + 1 - 2.0 / (2 + 16 * ri)) * (h * nu) / k
    C_const = 0.5 * W / k
    D_const = (h * nu) / k

    mu = -0.5 * R * T * (
        np.log(A_const / T)
        - B_const / T
        + np.log(iv(0, C_const / T) ** 2)
        - 2 * np.log(1 - np.exp(-D_const / T))
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
        "nn_distance": nn_distance,
        "adsorbate_mass_amu": mass_amu,
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
    lines = ["Hindered translation thermodynamics", "=" * 70,
             CONSTRAINT_WARNING, GEOMETRY_WARNING, ""]
    if contcar_file is not None:
        lines.append(f"Structure       : {contcar_file}")
    lines += [
        f"Temperature     : {r['temperature']:.2f} K",
        f"Barrier         : {r['barrier_eV']:.4f} eV",
        f"NN distance     : {r['nn_distance']:.4f} A",
        f"Adsorbate mass  : {r['adsorbate_mass_amu']:.4f} amu",
        f"Frequency       : {r['nu_cm_inv']:.2f} cm^-1",
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