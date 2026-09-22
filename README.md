# STEAM

Surface Thermodynamics and Entropy Analysis Module

## Installation

```bash
pip install git+https://github.com/d-parab/steam.git
```

## 'Zero point energy' and 'Harmonic oscillator' analysis

### Input format

A plain text file (eg. freq.txt) with one vibrational frequency per line, in cm^-1.

Frequencies are typically taken from a DFT vibrational frequency analysis.
Imaginary modes are not handled and should be removed before use.

### Command line

```bash
steam harmonic freq.txt -T 298.15 -o results.txt
```

Returns a dict with these keys:

| Key | Quantity | Units |
|---|---|---|
| `ZPE_J_per_mol`, `ZPE_eV` | Zero-point energy | J/mol, eV/molecule |
| `S_J_per_mol_K`, `S_eV_per_K` | Vibrational entropy | J/(mol K), eV/(molecule K) |
| `A_J_per_mol`, `A_eV` | Helmholtz free energy | J/mol, eV/molecule |
| `mu_J_per_mol`, `mu_eV` | Standard chemical potential | J/mol, eV/molecule |
| `U_J_per_mol`, `U_eV` | Internal energy | J/mol, eV/molecule |

All quantities are the vibrational contribution only, computed from
the harmonic oscillator partition function over the supplied modes. They
do not include the DFT electronic energy, which must be added separately.
For an adsorbate, the Helmholtz free
energy is then:

Helmholtz free energy = E_DFT + A
where `A` is the value reported here.

The contribution to standard chemical potential `mu` is identical to `A`
and is reported under both names for convenience.

## License

MIT

## Citation

If you use this code, please cite [add reference].
