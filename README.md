# lammps_to_gro

**Convert LAMMPS data files to GROMACS format** - A Python tool for molecular dynamics simulation workflow integration. Transform LAMMPS data files with force field parameters into GROMACS topology (.itp/.top) and coordinate (.gro) files for seamless cross-platform MD simulations.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Keywords:** molecular dynamics, LAMMPS, GROMACS, file conversion, force field, topology, coordinates, simulation, computational chemistry, MD workflow

Convert a LAMMPS data file (with force field parameters in the `* Coeffs` sections) into GROMACS topology (`.itp`/`.top`) and coordinate (`.gro`) files.

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/askforarun/lammps_to_gro.git
cd lammps_to_gro

# Convert LAMMPS data to GROMACS format
# (Example: 90 = 10 ethanol molecules × 9 atoms each)
python lammps_to_gro.py \
  --infile your_data.lammps \
  --residue-sizes 90 \
  --residue-names MOL

# Run simulation in GROMACS
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```
---

## Table of Contents
- [What This Tool Does](#what-this-tool-does)
- [Force Field Mapping & Compatibility](#force-field-mapping--compatibility)
- [Platform Compatibility](#platform-compatibility)
- [What You Need](#what-you-need)
- [Command Reference & Options](#command-reference--options)
  - [All CLI Options](#all-cli-options)
  - [Core Options](#core-options)
  - [Output Control](#output-control)
  - [Residue Assignment](#residue-assignment)
  - [Function Types](#function-types-gromacs-funct-column)
  - [What It Generates](#what-it-generates)
  - [Unit Conversions](#unit-conversions)
- [Requirements & Installation](#requirements--installation)
  - [System Requirements](#system-requirements)
  - [Installation](#installation)
- [How It Works](#how-it-works)
- [Tutorial](#tutorial)
  - [CLI with GROMACS Execution](#cli-with-gromacs-execution)
  - [Python API Usage](#python-api-usage)
- [Validation](#validation)
- [Troubleshooting](#troubleshooting)
  - [Common Issues and Solutions](#common-issues-and-solutions)
  - [Exit Codes Reference](#exit-codes-reference)
- [Open Source](#open-source)
- [Citation](#citation)
- [License](#license)

---

## What This Tool Does
Run molecular dynamics simulations in GROMACS using systems prepared in LAMMPS.

**Typical workflow**
1) Prepare a LAMMPS data file (with force field parameters).  
2) Convert it to GROMACS format with this tool.  
3) Run simulation/analysis in GROMACS.

**Input → Output**
- LAMMPS data file → GROMACS topology and coordinate files.

**Generated files**
- `sorted_<input>` — sorted copy of LAMMPS data file (always generated)  
- `forcefield.itp` — `[defaults]`, `[atomtypes]` (LJ parameters)  
- `<molecule>.itp` — `[moleculetype]`, `[atoms]`, `[bonds]`, `[pairs]`, `[angles]`, `[dihedrals]`  
- `topol.top` — `#include`s + `[system]` + `[molecules]`  
- `system.gro` — coordinates + box vectors  
---

## What You Need

**Input file format**
- LAMMPS data file with force field parameters; assumes `atom_style full` (see data file format given below).  
- Recommended reading: <https://docs.lammps.org/atom_style.html>

**Required sections**
- Box dimensions (`xlo xhi`, `ylo yhi`, `zlo zhi`)
- `Masses`, `Atoms`, `Bonds`, `Angles`, `Dihedrals`, `Impropers` (if present)
- `* Coeffs` sections (`Pair`, `Bond`, `Angle`, `Dihedral`, `Improper`)
- Missing any required section aborts conversion (CLI exit code 4).
- Atom IDs in `Atoms` must be contiguous and 1-based (`1..N`).

**Atoms section format**
```
atom-ID  molecule-ID  atom-type  charge  x  y  z
```

**Example LAMMPS data file**
```
# LAMMPS data file for crosslinked polymer network
12.0 48.0 xlo xhi
12.0 48.0 ylo yhi
12.0 48.0 zlo zhi

Masses
1 12.0107
2 1.00794
3 15.9994

Atoms
1 1 1 -0.120 15.234 22.456 18.123
2 1 2  0.060 15.345 22.567 18.234
...

Bonds
1 1 1 2
2 1 2 3
...

Bond Coeffs
1 340.0 1.529
...

Pair Coeffs
1 0.0628 3.981
...
```

---

## Force Field Mapping & Compatibility

### Supported LAMMPS Force Field Styles
```bash
pair_style      lj/cut/coul/cut 9 9   # or lj/cut/coul/long 9 9 (any cutoff value)
bond_style      harmonic
angle_style      harmonic
dihedral_style  fourier     # periodic dihedrals
dihedral_style  opls        # OPLS dihedrals 
improper_style  cvff        # CVFF impropers 
```

### LAMMPS to GROMACS Mapping
Suggested reading: https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
| LAMMPS Style | GROMACS Function | Description |
|--------------|------------------|-------------|
| `pair_style lj/cut/coul/cut` | `nb-func 1` | Lennard-Jones + Coulomb cutoff |
| `bond_style harmonic` | `bond-funct 1` | Harmonic bonds: E = k(r - r₀)² |
| `angle_style harmonic` | `angle-funct 1` | Harmonic angles: E = k(θ - θ₀)² |
| `dihedral_style fourier` | `dihedral-funct 1` | Fourier dihedrals |
| `dihedral_style opls` | `dihedral-funct 3` | OPLS dihedrals converted to Ryckaert-Bellemans form |
| `improper_style cvff` | `improper-funct 4` | CVFF impropers mapped to periodic improper form (`phi_s`, `k`, `n`) |


### Combination Rules & 1-4 Scaling
Choose a preset with `--comb-rule`:

| Preset | `comb-rule` | `fudgeLJ` | `fudgeQQ` | Typical match |
|--------|-------------|-----------|-----------|---------------|
| `amber` (default) | `2` | `0.5` | `0.8333` | LAMMPS `pair_modify mix arithmetic` + `special_bonds ... coul 0.8333` |
| `opls` | `2` | `0.5` | `0.5` | LAMMPS `pair_modify mix arithmetic` + `special_bonds ... coul 0.5` |

Notes:
- Both presets use arithmetic LJ mixing (`comb-rule 2`): `σ_ij = (σ_i + σ_j)/2`, `ε_ij = sqrt(ε_i * ε_j)`.
- `fudgeLJ` is `0.5` for both presets (matches `special_bonds lj 0.0 0.0 0.5`).
- Only `fudgeQQ` changes between presets.

---

## Platform Compatibility
Validated on: Linux (Ubuntu, CentOS, Red Hat, Debian), macOS (Intel/Apple Silicon), Windows 10/11 (native Python or WSL2).

---

## Command Reference & Options

### All CLI Options
```bash
python lammps_to_gro.py --help
```

### Core Options
| Option | Description | Default |
|--------|-------------|---------|
| `--infile` | Input LAMMPS data file | Required |
| `--molecule` | Molecule name for `.itp` | `MOL` |
| `--system` | System name for topology | `SYSTEM` |
| `--comb-rule` | Nonbonded preset for `[ defaults ]` (`amber`→`2, 0.5, 0.8333`; `opls`→`2, 0.5, 0.5`) | `amber` |
| `--gro-out` | Output coordinate file | `system.gro` |
| `--forcefield-out` | Output forcefield file | `forcefield.itp` |
| `--itp-out` | Output molecule topology | `<molecule>.itp` |
| `--topol-out` | Output topology file | `topol.top` |

### Output Control
| Option | Description | Default |
|--------|-------------|---------|
| `--quiet` | Suppress info messages | `False` |

### Residue Assignment (Required)
| Option | Description | Example |
|--------|-------------|---------|
| `--residue-sizes` | Atoms per residue block | `73 73 31` |
| `--residue-names` | Residue names for blocks | `RES1 RES2 RES3` |
| `--repeat` | Repeat residue block pattern until all atoms are assigned | Flag (no value) |

**Residue Assignment Modes**
- **Non-repeat mode (default)**: Residue blocks are consumed once; `sizes` must sum to exactly `n_atoms`
- **Repeat mode (`--repeat`)**: Pattern repeats until all atoms are assigned; useful for periodic systems

**Examples**
- `--residue-sizes 73 73 31 --residue-names ASP BEN ETH` → Atoms 1-73 → ASP, 74-146 → BEN, 147-177 → ETH (non-repeat)
- `--residue-sizes 50 --residue-names MON --repeat` → Pattern (50 atoms per MON residue) repeats for all atoms

### What It Generates
- `sorted_<input>` — sorted copy of LAMMPS data file (always generated)
- `forcefield.itp` — `[defaults]`, `[atomtypes]` (LJ parameters)
- `<molecule>.itp` — `[moleculetype]`, `[atoms]`, `[bonds]`, `[pairs]`, `[angles]`, `[dihedrals]`
- `topol.top` — `#include`s + `[system]` + `[molecules]`
- `system.gro` — coordinates + box vectors

### Unit Conversions
Lengths and energies are converted from LAMMPS (real units) to GROMACS SI-based units:
- **Length**: Å → nm (× 0.1)
- **Energy**: kcal/mol → kJ/mol (× 4.184)
- **Bond force constant**: kcal/mol/Å² → kJ/mol/nm² (× 4.184 × 100 × 2)
- **Angle force constant**: kcal/mol/rad² → kJ/mol/rad² (× 4.184 × 2)

---

## Requirements & Installation

### System Requirements
- Python 3.8+
- Optional: `networkx` (for `--pairs-method networkx`)
- Optional: GROMACS (running converted files) and LAMMPS (creating inputs)

### Installation
```bash
git clone <repository-url>
cd lammps_to_gro
python lammps_to_gro.py --help
```

Optional dependency:
```bash
pip install networkx
```

---

## How It Works

1) **File Reading & Validation** — parse header box dimensions; validate section formats (Masses, Atoms, etc.); check coefficient consistency and atom type coverage; ensure every interaction type has corresponding parameters; creates `sorted_<input>` copy with key sections sorted by numeric IDs (always generated). Section parsing is header-based (so internal blank lines do not truncate sections), and comment lines (`# ...`) are ignored during data parsing.

2) **Data Parsing** — parse masses (type ID → mass), atoms (ID, mol ID, type, charge, coordinates), and all bonded interactions; convert units (Å→nm, kcal/mol→kJ/mol); validate total charge neutrality; scale force constants appropriately. Atom IDs are validated as contiguous and 1-based (`1..N`) before topology generation. Generate atom names using just the element symbol (e.g., `C`, `H`, `O`, `N`, `CL`). If element/type mapping is missing, atom naming falls back to the first character of the residue name, or `X` if unavailable. This simple format always fits within the GRO 5-character limit and works for systems of any size. GROMACS uses atom indices internally, so duplicate atom names are handled correctly.

3) **Residue Assignment** — assign atoms to residues using block-based patterns from required `--residue-sizes` and `--residue-names`; generate residue numbers and local atom indices; validate residue sizes/names arrays match in length; ensure total residue atoms equals system atoms (non-repeat mode) or pattern can repeat (repeat mode).  
   *Example*: `--residue-sizes 73 73 31 --residue-names ASP BEN ETH` → Atoms 1-73 → ASP, 74-146 → BEN, 147-177 → ETH (non-repeat).  
   *Example with repeat*: `--residue-sizes 50 --residue-names MON --repeat` → Pattern (50 atoms per MON residue) repeats for all atoms.  

4) **Force Field Conversion** — convert LAMMPS coefficients to GROMACS format; map atom types to GROMACS atomtype names; validate all atoms reference defined atom types; detect orphaned coefficient sections; warn about unused atom type definitions.

5) **1–4 Pair Generation** — generate 1-4 pairs using BFS (default) or NetworkX algorithms; identify atoms separated by exactly 3 bonds; validate all bond/angle/dihedral/improper references exist in atoms; check correct column counts for each bonded section.

6) **Output Generation** — write `forcefield.itp` with `[defaults]` and `[atomtypes]`; write `<molecule>.itp` with `[moleculetype]`, `[atoms]`, `[bonds]`, `[pairs]`, `[angles]`, `[dihedrals]`; write `topol.top` with include statements and system info; write `system.gro` with coordinates and box vectors.



---

## Tutorial

This tutorial demonstrates converting three LAMMPS data files to GROMACS format using `lammps_to_gro.py`:

1. **`multi_ethanol_data.lammps`** — 10 ethanol molecules
2. **`mixed_data.lammps`** — 1 aspirin, 2 benzene, and 5 ethanol molecules
3. **`mixed_abc_data_replicate.lammps`** — 1 aspirin, 1 benzene, and 1 ethanol replicated 8 times (illustrates the `--repeat` functionality, useful for large molecules with repeating residue patterns)

### CLI with GROMACS Execution

### First Example
```bash
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --molecule ETHANOL \
  --system POLYMER \
  --residue-sizes 90 \
  --residue-names ETH
```

### Second Example
Mixed system: 1 aspirin (21 atoms), 2 benzene (24 atoms), 5 ethanol (45 atoms).

```bash
python lammps_to_gro.py \
  --infile mixed_data.lammps \
  --molecule MIXED \
  --system SIMULATION \
  --residue-sizes 21 24 45 \
  --residue-names ASP BEN ETH
```

### Third Example
```bash
python lammps_to_gro.py \
  --infile mixed_abc_data_replicate.lammps \
  --molecule COMPLEX \
  --system SIMULATION \
  --residue-sizes 21 12 9 \
  --residue-names ASP BEN ETH \
  --repeat
```

**Run GROMACS**
```bash
# prepare
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr -maxwarn 4
# run
gmx mdrun -v -deffnm md
# analyze
gmx energy -f md.edr -o energy.xvg
```

Additional CLI examples:
```bash
# minimal conversion (multi_ethanol_data.lammps)
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --residue-sizes 90 \
  --residue-names ETH

# custom output filenames (mixed_data.lammps)
python lammps_to_gro.py \
  --infile mixed_data.lammps \
  --gro-out coords.gro \
  --forcefield-out ff.itp \
  --topol-out system.top \
  --molecule MIXED \
  --residue-sizes 21 24 45 \
  --residue-names ASP BEN ETH

# single residue type system (multi_ethanol_data.lammps)
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --residue-sizes 90 \
  --residue-names ETH

# with explicit GROMACS function types (all defaults shown)
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --molecule ETHANOL \
  --system POLYMER \
  --residue-sizes 90 \
  --residue-names ETH \
  --bond-funct 1 \
  --angle-funct 1 \
  --dihedral-funct 1 \
  --improper-funct 4 \
  --pairs-funct 1

# OPLS dihedrals with RB conversion (for OPLS-style dihedral coeffs: type k1 k2 k3 k4)
python lammps_to_gro.py \
  --infile data.benzene \
  --molecule benzene \
  --system benzene_system \
  --residue-sizes 12 \
  --residue-names BEN \
  --comb-rule opls \
  --dihedral-funct 3 \
  --improper-funct 4

# multiple identical residue blocks (mixed_data.lammps)
python lammps_to_gro.py \
  --infile mixed_data.lammps \
  --residue-sizes 21 24 45 \
  --residue-names ASP BEN ETH

# periodic polymer with repeat pattern (mixed_abc_data_replicate.lammps)
python lammps_to_gro.py \
  --infile mixed_abc_data_replicate.lammps \
  --residue-sizes 21 12 9 \
  --residue-names ASP BEN ETH \
  --repeat \
  --molecule COMPLEX

# mixed residue types with repeat (mixed_abc_data_replicate.lammps)
python lammps_to_gro.py \
  --infile mixed_abc_data_replicate.lammps \
  --residue-sizes 21 12 9 \
  --residue-names ASP BEN ETH \
  --repeat \
  --system COMPLEX

# using absolute paths (mixed_data.lammps)
python lammps_to_gro.py \
  --infile /Users/arunsrikanthsridhar/Downloads/lammps_to_gro/mixed_data.lammps \
  --molecule MIXED \
  --system SIMULATION \
  --residue-sizes 21 24 45 \
  --residue-names ASP BEN ETH \
  --gro-out /home/user/gromacs/coords.gro

# quiet mode (multi_ethanol_data.lammps)
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --residue-sizes 90 \
  --residue-names ETH \
  --quiet
```

### Python API Usage
```python
from lammps_to_gro import build_gromacs_itps_and_gro_from_lammps

# First example: multi_ethanol_data.lammps
build_gromacs_itps_and_gro_from_lammps(
    infile="multi_ethanol_data.lammps",
    molecule_name="ETHANOL",
    system_name="POLYMER",
    residue_block_sizes=(90,),
    residue_block_names=("ETH",),
)
```
Example with repeat pattern:
```python
from lammps_to_gro import build_gromacs_itps_and_gro_from_lammps

# Third example: mixed_abc_data_replicate.lammps with repeat
build_gromacs_itps_and_gro_from_lammps(
    infile="mixed_abc_data_replicate.lammps",
    molecule_name="COMPLEX",
    system_name="SIMULATION",
    residue_block_sizes=(21, 12, 9),
    residue_block_names=("ASP", "BEN", "ETH"),
    repeat_residue_blocks=True,
)
```
Additional API examples:
```python
# custom output paths (mixed_data.lammps)
build_gromacs_itps_and_gro_from_lammps(
    infile="mixed_data.lammps",
    molecule_name="MIXED",
    system_name="SIMULATION",
    residue_block_sizes=(21, 24, 45),
    residue_block_names=("ASP", "BEN", "ETH"),
    gro_outfile="/path/to/output/coordinates.gro",
    forcefield_outfile="/path/to/output/forcefield.itp",
    itp_outfile="/path/to/output/molecule.itp",
    topol_outfile="/path/to/output/topology.top",
)

# advanced conversion with custom function types and pairs (mixed_data.lammps)
build_gromacs_itps_and_gro_from_lammps(
    infile="mixed_data.lammps",
    molecule_name="MIXED",
    system_name="SIMULATION",
    residue_block_sizes=(21, 24, 45),
    residue_block_names=("ASP", "BEN", "ETH"),
    bond_funct=2,
    angle_funct=2,
    dihedral_funct=9,
    pairs_method="networkx",
)

# batch processing multiple systems
systems = [
    ("multi_ethanol_data.lammps", "ETHANOL", "POLYMER"),
    ("mixed_data.lammps", "MIXED", "SIMULATION"),
    ("mixed_abc_data_replicate.lammps", "COMPLEX", "SIMULATION"),
]

for data_file, mol_name, sys_name in systems:
    if data_file == "multi_ethanol_data.lammps":
        residue_sizes = (90,)
        residue_names = ("ETH",)
    elif data_file == "mixed_data.lammps":
        residue_sizes = (21, 24, 45)
        residue_names = ("ASP", "BEN", "ETH")
    else:  # mixed_abc_data_replicate.lammps
        residue_sizes = (21, 12, 9)
        residue_names = ("ASP", "BEN", "ETH")
    
    build_gromacs_itps_and_gro_from_lammps(
        infile=data_file,
        molecule_name=mol_name,
        system_name=sys_name,
        residue_block_sizes=residue_sizes,
        residue_block_names=residue_names,
    )
    print(f"Converted {data_file} to GROMACS format")

# batch processing with repeat patterns
polymers = [
    ("mixed_abc_data_replicate.lammps", "COMPLEX", "SIMULATION_A"),
    ("mixed_abc_data_replicate.lammps", "COMPLEX", "SIMULATION_B"),
]

for data_file, mol_name, sys_name in polymers:
    build_gromacs_itps_and_gro_from_lammps(
        infile=data_file,
        molecule_name=mol_name,
        system_name=sys_name,
        residue_block_sizes=(21, 12, 9),
        residue_block_names=("ASP", "BEN", "ETH"),
        repeat_residue_blocks=True,
    )
    print(f"Converted {data_file} with repeat pattern")
```

---

## Validation

The following validation uses the **multi-ethanol example** (`multi_ethanol_data.lammps`) from Example 1 above.

Energy comparison between LAMMPS (real units, kcal/mol) and GROMACS (kJ/mol converted to kcal/mol) for `multi_ethanol_data.lammps`.

**Commands used:**

```bash
# Convert LAMMPS data to GROMACS topology/coordinates
python lammps_to_gro.py \
  --infile multi_ethanol_data.lammps \
  --molecule ETHANOL \
  --system POLYMER \
  --residue-sizes 90 \
  --residue-names ETH

# LAMMPS single-point energy (input: example_lammps_input.lmp)
lmp < example_lammps_input.lmp

# GROMACS single-point energy decomposition
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -v -deffnm md
echo "1 2 3 4 5 6 7 8" | gmx energy -f md.edr -o energy.xvg
```

| Interaction | LAMMPS (kcal/mol) | GROMACS (kJ/mol) | GROMACS (kcal/mol) | Match |
|-------------|-------------------|------------------|-------------------|-------|
| Bond | 1.113 | 4.656 | 1.113 | ✅ |
| Angle | 2.475 | 10.356 | 2.475 | ✅ |
| Dihedral | 0.015 | 0.064 | 0.015 | ✅ |
| LJ | 1.523 | 6.374 | 1.524 | ✅ |
| Coulomb | -38.499 | -161.08 | -38.50 | ✅ |
| **Total** | -33.37 | -139.63 | -33.38 | ✅ |

Note: This comparison assumes `pair_modify mix arithmetic` in LAMMPS (GROMACS `comb-rule 2`). Small differences are expected due to numerical precision and floating-point handling; the total potential energy matches within ~0.1%.

---

### OPLS/Ryckaert-Bellemans Validation

Energy comparison for `data.benzene` with OPLS dihedrals converted to Ryckaert-Bellemans (`--dihedral-funct 3`, `--comb-rule opls`; see additional CLI examples).

**Commands used:**

```bash
# Convert LAMMPS data to GROMACS topology/coordinates
python lammps_to_gro.py \
  --infile data.benzene \
  --molecule benzene \
  --system benzene_system \
  --residue-sizes 12 \
  --residue-names BEN \
  --comb-rule opls \
  --dihedral-funct 3 \
  --improper-funct 4

# Keep generated system.gro for GROMACS validation.
# LAMMPS single-point energy (input: in.benzene, special_bonds ... coul 0.5)
lmp < in.benzene

# GROMACS single-point energy decomposition
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -v -deffnm md
echo "1 2 3 4 5 6 7 8 9" | gmx energy -f md.edr -o energy_benzene_opls_exact.xvg
```

| Interaction | LAMMPS (kcal/mol) | GROMACS (kJ/mol) | GROMACS (kcal/mol) | Match |
|-------------|-------------------|------------------|-------------------|-------|
| Bond | 0.2255 | 0.9436 | 0.2255 | ✅ |
| Angle | 0.001318 | 0.005521 | 0.001320 | ✅ |
| RB Dihedral | 0.0003246 | 0.001345 | 0.0003215 | ✅ |
| Improper | 1.61e-05 | 6.73e-05 | 1.61e-05 | ✅ |
| LJ (14+SR) | 5.3283 | 22.2938 | 5.3284 | ✅ |
| Coulomb (14+SR) | 2.7207 | 11.3834 | 2.7208 | ✅ |
| **Total** | 8.2762 | 34.6277 | 8.2762 | ✅ |

Note: This comparison assumes `special_bonds ... coul 0.5` in LAMMPS (and `--comb-rule opls` in conversion). Small differences are expected due to coordinate precision and floating-point handling; the total potential energy matches within ~0.2%.

**OPLS → RB conversion formulas used:**
```
C0 = K2 + (K1 + K3)/2
C1 = (-K1 + 3*K3)/2
C2 = -K2 + 4*K4
C3 = -2*K3
C4 = -4*K4
C5 = 0
```

Where K1-K4 are the LAMMPS OPLS coefficients (converted to kJ/mol).

## Troubleshooting

**Charge neutrality error** — ensure charges sum to zero.  
**Missing required sections** — add any absent `Masses`/`Atoms`/`Bonds`/`Angles`/`Dihedrals` and all `* Coeffs`.  
**Missing coefficient coverage** — provide coeffs for every type referenced.  
**Atom ID validation error** — ensure `Atoms` IDs are contiguous and 1-based (`1..N`) with no gaps.  
**Residue assignment mismatch** — sizes must sum to atom count (non-repeat) or be positive (repeat); names length must match sizes.  
**CLI argument errors** — re-run with `--help`.  
**File not found** — verify `--infile` path.

### Exit Codes Reference
| Exit | Meaning | Common cause |
|------|---------|--------------|
| 0 | Success | Normal completion |
| 1 | Unexpected error | Internal exception |
| 2 | File not found | Missing input |
| 3 | Invalid input | CLI/format/validation error |
| 4 | Conversion failed | Runtime errors; missing sections |

---

## Open Source

MIT-licensed project. Repo: <https://github.com/askforarun/lammps_to_gro>  
Contributions welcome via issues/PRs.

## Contributing

We welcome contributions! Here's how you can help:

- **Report Bugs**: Open an issue with detailed description and error messages
- **Request Features**: Suggest enhancements or new force field support
- **Submit Pull Requests**: 
  - Fork the repository
  - Create a feature branch
  - Make your changes with tests
  - Submit a pull request
- **Improve Documentation**: Help us improve the README and code comments
- **Share Your Use Cases**: Let us know how you're using the tool in your research

### Development Setup

```bash
git clone https://github.com/askforarun/lammps_to_gro.git
cd lammps_to_gro
# Test with example files
python lammps_to_gro.py --infile multi_ethanol_data.lammps --residue-sizes 90 --residue-names ETH
```

---

## Citation

**BibTeX**
```bibtex
@software{lammps_to_gro,
  title = {lammps_to_gro},
  author = {Sridhar, Arun Srikanth},
  type = {software},
  doi = {10.5281/zenodo.18329896},
  url = {https://github.com/askforarun/lammps_to_gro}
}
```
**Plain text** — Sridhar, Arun Srikanth. *lammps_to_gro*. Zenodo, 2024. doi:10.5281/zenodo.18329896

---

## License
MIT License (see `LICENSE`).
