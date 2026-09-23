# STEAM

Surface Thermodynamics and Entropy Analysis Module

## Installation

```bash
pip install git+https://github.com/d-parab/steam.git
```
Requires Python 3.9+, NumPy, SciPy and ASE

## Harmonic oscillator

### Command line

STEAM uses subcommands. Run `steam --help` or `steam <command> --help` for details.

```bash
steam harmonic freq.txt -T 298.15 -o results.txt
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

## Python API

```python
from steam import (
    harmonic_properties,
    hindered_translation_properties,
    hindered_rotation_properties,
)

r = harmonic_properties("freq.txt", temperature=298.15)
print(r["A_eV"], r["S_J_per_mol_K"])
```

Each returns a dict. The hindered modules also hold the separate
harmonic and anharmonic parts (`S_HO`, `delS`, `A_HO`, `delA`), which
are not printed in the report.

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


## License

MIT

## Citation

If you use this code, please cite [add reference].
