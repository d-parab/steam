# STEAM

Surface Thermodynamics and Entropy Analysis Module

## Installation

```bash
pip install git+https://github.com/d-parab/steam.git
```
Requires Python 3.9+, NumPy, SciPy and ASE

## Development installation
```bash
git clone https://github.com/d-parab/steam.git
cd steam
pip install -e .
```

## Harmonic oscillator

### Command line

STEAM uses subcommands. Run `steam --help` or `steam <command> --help` for details.

```bash
steam harmonic freq.txt -T 298.15 -o harmonic.txt
```

| Option | Meaning |
|---|---|
| `freq_file` | text file of frequencies in cm^-1, one per line, imaginary nodes removed beforehand |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `-o`, `--output` | also write the report to this file |


## Hindered translation

### Command line

```bash
steam hindered-trans CONTCAR_Hind-trans -W 1.01 -b 2.83 -T 298.15 -o trans.txt
```

| Option | Meaning |
|---|---|
| `contcar` | VASP CONTCAR file |
| `-W`, `--barrier` | translation barrier in eV |
| `-b`, `--nn-distance` | nearest neighbour distance in A |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `-o`, `--output` | also write the report to this file |

All atoms must be constrained in the CONTCAR except the adsorbate. The
output is derived from the barrier assuming a **(111) surface
geometry** — modify the module for other geometries.

## Hindered rotation

### Command line

```bash
steam hindered-rot CONTCAR_Hind-rot -W 0.51 -a 45 -n 6 -s 2 -T 298.15 -o rot.txt
```

| Option | Meaning |
|---|---|
| `contcar` | VASP CONTCAR file |
| `-W`, `--barrier` | rotation barrier in eV |
| `-a`, `--atom-index` | index of the atom about which the molecule rotates |
| `-n`, `--n-minima` | number of equivalent minima in a full rotation |
| `-s`, `--symmetry-number` | symmetry number |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `-o`, `--output` | also write the report to this file |

All atoms must be constrained except the adsorbate.

## Ideal gas translation (3D)

### Command line

```bash
steam ig-trans CONTCAR-IG -T 298.15 -P 1e5 -o ig.txt
```

| Option | Meaning |
|---|---|
| `contcar` | VASP CONTCAR of the gas-phase molecule |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `-P`, `--pressure` | standard-state pressure in Pa (default 1e5) |
| `-o`, `--output` | also write the report to this file |

## Ideal gas rotation

### Non-linear molecules

```bash
steam ig-rot-nonlinear CONTCAR-IG -s 2 -T 298.15 -o igrot.txt
```

| Option | Meaning |
|---|---|
| `contcar` | VASP CONTCAR of the gas-phase molecule |
| `-s`, `--symmetry-number` | rotational symmetry number |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `-o`, `--output` | also write the report to this file |

### Linear molecules

```bash
steam ig-rot-linear CONTCAR-IG -s 2 -T 298.15 -o igrot.txt
```

| Option | Meaning |
|---|---|
| `contcar` | VASP CONTCAR of the gas-phase molecule |
| `-s`, `--symmetry-number` | rotational symmetry number |
| `-T`, `--temperature` | temperature in K (default 298.15) |
| `--tol` | relative tolerance for the linearity check (default 1e-3) |
| `-o`, `--output` | also write the report to this file |

The linear molecule must be aligned along one of the unit cell axes
(x, y or z) in the CONTCAR.

## Output

All quantities are the additional contribution, computed from
various assumptions for the vibrational modes. They do not include the
DFT electronic energy, which must be added separately.
For an adsorbate, the Helmholtz free
energy is then:

Helmholtz free energy = E_DFT + A
where `A` is the value reported here.

#### Reference article for formulae used in hindered translation and rotation codes
Lynza H. Sprowl, Charles T. Campbell, Líney Árnadóttir; Hindered Translator and Hindered Rotor Models for Adsorbates: Partition Functions and Entropies. J. Phys. Chem. C 12 May 2016; 120 (18): 9719–9731.

## Citation

If you use this code, please cite [add reference].
