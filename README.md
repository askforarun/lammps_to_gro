# lammps_to_gro

Convert a LAMMPS data file (with force field parameters in the `* Coeffs` sections) into GROMACS topology (`.itp`/`.top`) and coordinate (`.gro`) files.

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
pair_modify     mix arithmetic tail no
bond_style      harmonic
angle_style      harmonic
dihedral_style  fourier
```

### LAMMPS to GROMACS Mapping
Suggested reading: https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
| LAMMPS Style | GROMACS Function | Description |
|--------------|------------------|-------------|
| `pair_style lj/cut/coul/cut` | `nb-func 1` | Lennard-Jones + Coulomb cutoff |
| `bond_style harmonic` | `bond-funct 1` | Harmonic bonds: E = k(r - r₀)² |
| `angle_style harmonic` | `angle-funct 1` | Harmonic angles: E = k(θ - θ₀)² |
| `dihedral_style fourier` | `dihedral-funct 1` | Fourier dihedrals |
| `improper_funct 4` | `improper-funct 4` | Periodic impropers (default) |

### Combination Rules & 1-4 Scaling
| Parameter | Default | Description |
|-----------|---------|-------------|
| `comb-rule` | `2` | Arithmetic mixing: σ_ij = (σ_i + σ_j)/2, ε_ij = √(ε_i × ε_j). Matches LAMMPS `pair_modify mix arithmetic`. |
| `fudgeLJ` | `0.5` | 1-4 LJ scaling factor. Matches LAMMPS `special_bonds lj 0.0 0.0 0.5`. |
| `fudgeQQ` | `0.8333` | 1-4 Coulomb scaling factor. Matches LAMMPS `special_bonds coul 0.0 0.0 0.8333`. |

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

1) **File Reading & Validation** — parse header box dimensions; validate section formats (Masses, Atoms, etc.); check coefficient consistency and atom type coverage; ensure every interaction type has corresponding parameters; creates `sorted_<input>` copy with key sections sorted by numeric IDs (always generated).

2) **Data Parsing** — parse masses (type ID → mass), atoms (ID, mol ID, type, charge, coordinates), and all bonded interactions; convert units (Å→nm, kcal/mol→kJ/mol); validate total charge neutrality; scale force constants appropriately. Generate atom names using pattern `element_residue_first_letter + local_index` (e.g., residue "ASP", carbon atom at index 5 → "C_A5"). The output preserves original LAMMPS atom IDs for global numbering (first column in `.gro` file), while local indices in atom names reset to 1 for each new residue, ensuring compatibility with bonded interactions.

3) **Residue Assignment** — assign atoms to residues using block-based patterns (or default single residue); generate residue numbers and local atom indices; validate residue sizes/names arrays match in length; ensure total residue atoms equals system atoms (non-repeat mode) or pattern can repeat (repeat mode).  
   *Example*: `--residue-sizes 73 73 31 --residue-names ASP BEN ETH` → Atoms 1-73 → ASP, 74-146 → BEN, 147-177 → ETH (non-repeat).  
   *Example with repeat*: `--residue-sizes 50 --residue-names MON --repeat` → Pattern (50 atoms per MON residue) repeats for all atoms.  
   *Default*: All atoms → single residue (first 3 chars of molecule name).

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
```bash
python lammps_to_gro.py \
  --infile mixed_data.lammps \
  --molecule MIXED \
  --system SIMULATION \
  --residue-sizes 21 24 45 \ # 21 atoms (1 aspirin), 24 atoms (2 benzene × 12 atoms each), 45 atoms (5 ethanol × 9 atoms each)
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

### Validation
The following validation uses the **multi-ethanol example** (`multi_ethanol_data.lammps`) from Example 1 above.

**LAMMPS:**
```bash
lmp < example_lammps_input.lmp
```

**GROMACS:**
```bash
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -v -deffnm md
echo "1 2 3 4 5 6 7 8" | gmx energy -f md.edr -o energy.xvg
```

Energy comparison between LAMMPS (real units, kcal/mol) and GROMACS (kJ/mol converted to kcal/mol) for `multi_ethanol_data.lammps`:

| Interaction | LAMMPS (kcal/mol) | GROMACS (kJ/mol) | GROMACS (kcal/mol) | Match |
|-------------|-------------------|------------------|-------------------|-------|
| Bond | 1.113 | 4.656 | 1.113 | ✅ |
| Angle | 2.475 | 10.356 | 2.475 | ✅ |
| Dihedral | 0.015 | 0.064 | 0.015 | ✅ |
| LJ | 1.523 | 6.374 | 1.524 | ✅ |
| Coulomb | -38.499 | -161.08 | -38.50 | ✅ |
| **Total** | -33.37 | -139.63 | -33.38 | ✅ |

Note: Requires `pair_modify mix arithmetic` in LAMMPS to match GROMACS `comb-rule 2`.

---

## Troubleshooting

**Charge neutrality error** — ensure charges sum to zero.  
**Missing required sections** — add any absent `Masses`/`Atoms`/`Bonds`/`Angles`/`Dihedrals` and all `* Coeffs`.  
**Missing coefficient coverage** — provide coeffs for every type referenced.  
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
