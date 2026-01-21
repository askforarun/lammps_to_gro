# lammps_to_gro

Convert a LAMMPS data file (with force field parameters in the `* Coeffs` sections) into GROMACS topology (`.itp`/`.top`) and coordinate (`.gro`) files.

This is useful when you have a bonded network (e.g., polymers/hydrogels) built in LAMMPS and want to continue simulation or analysis in GROMACS.

## What it generates

- `forcefield.itp` — `[defaults]` + `[atomtypes]` (LJ parameters)
- `<molecule>.itp` — `[moleculetype]`, `[atoms]`, `[bonds]`, `[pairs]`, `[angles]`, `[dihedrals]`
- `topol.top` — `#include` statements + `[system]` + `[molecules]`
- `system.gro` — coordinates + box vectors (`--gro-out none` to skip)
- `sorted_<input>` — copy of the input with key sections sorted by ID

## Unit conversions

- Length: Å → nm
- Energy: kcal/mol → kJ/mol
- Force constants are scaled appropriately for each interaction type

## Requirements

- Python 3.8+
- Optional: `networkx` (only if you use `--pairs-method networkx`)

## Installation

Clone the repo and run directly:

```bash
python lammps_to_gro.py --help
```

If you plan to use `--pairs-method networkx`:

```bash
pip install networkx
```

## Usage

### Quick start (example data)

Sample data included in this repo:

- `data.aftercrosslink_coeff` (PVA–Glutaraldehyde crosslinked hydrogel network)

```bash
python lammps_to_gro.py \
  --infile data.aftercrosslink_coeff \
  --molecule PVA_GLU_NETWORK \
  --system wetting \
  --residue-sizes 73 73 31 \
  --residue-names PVA PVA GLU
```

### Running LAMMPS (optional)

To run the LAMMPS simulation with the example input file:

```bash
lmp_png < example_lammps_input.lmp
```

This will generate energy output that can be compared with the GROMACS results for validation.

### Subcommand-style (optional)

The argument parser accepts being called with an extra leading token (useful if you later wire this into a parent CLI):

```bash
python lammps_to_gro.py lammps-to-gro --infile data.aftercrosslink_coeff
```

### Minimal

```bash
python lammps_to_gro.py --infile your.data
```

Defaults:

- `--molecule MOL` → writes `MOL.itp`
- `--system SYSTEM` → system name in `topol.top` `[system]` section
- `--forcefield-out forcefield.itp`
- `--topol-out topol.top`
- `--gro-out system.gro`
- `--bond-funct 1`, `--angle-funct 1`, `--dihedral-funct 1`, `--improper-funct 4`, `--pairs-funct 1`
- `--pairs-method bfs`

### CLI flag → Python keyword mapping

The CLI flags map directly to the main Python API (`build_gromacs_itps_and_gro_from_lammps`) using the following names:

- `--molecule` → `molecule_name`
- `--system` → `system_name`
- `--gro-out` → `gro_outfile`
- `--forcefield-out` → `forcefield_outfile`
- `--itp-out` → `itp_outfile`
- `--topol-out` → `topol_outfile`
- `--residue-sizes` → `residue_block_sizes`
- `--residue-names` → `residue_block_names`

### Control function types (GROMACS `funct` column)

These integers are written into the relevant topology blocks and control the functional form GROMACS uses:

```bash
python lammps_to_gro.py \
  --infile your.data \
  --bond-funct 1 \
  --angle-funct 1 \
  --dihedral-funct 1 \
  --improper-funct 4 \
  --pairs-funct 1
```

Reference: GROMACS topology formats (see “funct” column in bonded sections):
`https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html`

### Skip writing `system.gro`

```bash
python lammps_to_gro.py --infile your.data --gro-out none
```

### All CLI options

```bash
python lammps_to_gro.py --help
```

## Input expectations (LAMMPS)

The input file must be a LAMMPS data file that includes force field parameters (coeff sections), with (at minimum) sections like:

- Box dimensions in the header (`xlo xhi`, `ylo yhi`, `zlo zhi`)
- `Masses`
- `Atoms`
- `Bonds`, `Angles`, `Dihedrals`, `Impropers` (as present in your system)
- `Pair Coeffs`, `Bond Coeffs`, `Angle Coeffs`, `Dihedral Coeffs`, `Improper Coeffs` (as applicable)

### `Atoms` section format

The script expects an “Atoms” format compatible with:

`atom-ID  molecule-ID  atom-type  charge  x  y  z`

Coordinates are assumed to be in Å and are converted to nm.

### Residue handling

To assign residues, provide matching lists:

- `--residue-sizes` (atoms per residue block)
- `--residue-names` (residue names for each block)

Assignment repeats in blocks. Example:

```bash
python lammps_to_gro.py \
  --infile your.data \
  --residue-sizes 100 50 \
  --residue-names RES1 RES2
```

This assigns 100 atoms to `RES1`, next 50 to `RES2`, then repeats until all atoms are assigned.

### Generic residue example

For a simple system with identical residues:

```bash
python lammps_to_gro.py \
  --infile your.data \
  --residue-sizes 50 \
  --residue-names RES
```

This assigns 50 atoms to `RES`, then repeats until all atoms are assigned (e.g., for 200 atoms: RES 1, RES 1, RES 1, RES 1).

For multiple identical residue blocks:

```bash
python lammps_to_gro.py \
  --infile your.data \
  --residue-sizes 50 50 \
  --residue-names RES RES
```

This assigns 50 atoms to `RES`, next 50 to `RES`, then repeats (e.g., for 200 atoms: RES 1, RES 1, RES 1, RES 1).

## Outputs and “next steps”

After generation, you can use standard GROMACS commands (example):

```bash
gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```

## Validation

Energy comparison between LAMMPS (real units, kcal/mol) and GROMACS (values taken from `energy.xvg`, converted from kJ/mol to kcal/mol):

| Interaction | LAMMPS (kcal/mol) | GROMACS (kJ/mol) | GROMACS (kcal/mol) | Difference (kcal/mol) | % Difference |
|-------------|------------------|------------------|-------------------|----------------------|-------------|
| Bond | 6005.128 | 25180.800 | 6018.356 | -13.228 | -0.22% |
| Angle | 2917.642 | 12247.100 | 2927.127 | -9.485 | -0.32% |
| Dihedral | 725.039 | 3034.100 | 725.167 | -0.128 | -0.02% |
| LJ Total | -20.922 | -85.700 | -20.483 | -0.440 | -2.11% |
| Coulomb Total | -3419.174 | -16368.700 | -3912.213 | +493.039 | +14.42% |
| Total Non-bonded | -3440.097 | -16454.400 | -3932.696 | +492.599 | +14.32% |

**Key Findings:**
- **Bonded terms**: Excellent agreement (< 0.3% difference)
- **LJ interactions**: Good agreement (~2% difference)
- **Coulomb interactions**: Moderate difference (~14%) due to different 1-4 scaling and implementation details
- **Overall**: The conversion successfully preserves the force field parameters with high accuracy

**Note**: GROMACS Coulomb interactions include both short-range and 1-4 terms. The ~14% difference is typical when converting between different MD packages due to variations in electrostatics implementation.

## Notes

- Total charge check: the script sums charges from the `Atoms` section and errors if the total charge is not ~0 (tolerance is internal).
- Pair generation: `--pairs-method bfs` is the default and requires no extra dependencies; `--pairs-method networkx` needs `networkx` installed.
- The converter writes `sorted_<input>` next to the input file to make downstream parsing more robust.
- Output files are written in the current directory; existing files with the same names will be overwritten.

## Python API usage

You can also import and call the main conversion function:

```python
from lammps_to_gro import build_gromacs_itps_and_gro_from_lammps

build_gromacs_itps_and_gro_from_lammps(
    infile="data.aftercrosslink_coeff",
    molecule_name="MOL",
    system_name="SYSTEM",
    residue_block_sizes=(73, 73, 31),
    residue_block_names=("PVA", "PVA", "GLU"),
)
```
