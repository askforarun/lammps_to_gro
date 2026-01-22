#!/usr/bin/env python3
#
# LAMMPS -> GROMACS converter
#
# Input: LAMMPS data file with force field parameters (coeff sections).
# Output: GROMACS topology and coordinate files (forcefield.itp, <molecule>.itp, topol.top, system.gro).
# See README.md for detailed usage/examples.
#

import argparse
import logging
import os
import sys
from collections import deque

# ============================================================================
# MODULE-LEVEL CONSTANTS
# ============================================================================

# Atomic numbers for element symbols
ATOMIC_NUMBERS = {
    "H": 1, "C": 6, "N": 7, "O": 8,
    "F": 9, "P": 15, "S": 16,
    "CL": 17, "BR": 35, "I": 53
}

# Default LAMMPS atom type ID -> GROMACS atom type name mapping
DEFAULT_TYPE_NAME_MAP = {
    1: "c3", 2: "oh", 3: "hc", 4: "h1", 5: "ho", 6: "c6",
    7: "os", 8: "h2", 9: "c3", 10: "c6", 11: "c3", 12: "c6",
}

# Default GROMACS atom type name -> element symbol mapping
DEFAULT_ELEMENT_MAP = {
    "c3": "C", "c6": "C",
    "hc": "H", "h1": "H", "h2": "H",
    "oh": "O", "ho": "H", "os": "O",
}

# Unit conversion factors
ANGSTROM_TO_NM = 0.1
KCAL_TO_KJ = 4.184

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# ============================================================================
# HELPER FUNCTIONS (MODULE LEVEL)
# ============================================================================

def get_atomic_number(atom_name, element_map):
    # Resolve a GROMACS atomtype name (e.g., "c3") to an atomic number via element_map.
    if atom_name not in element_map:
        raise ValueError(f"Missing element_map entry for atomtype '{atom_name}'")
    elem = element_map[atom_name].upper()
    if elem not in ATOMIC_NUMBERS:
        raise ValueError(f"Unsupported element '{elem}'")
    return ATOMIC_NUMBERS[elem]


def read_section(lines, key):
    # Return the non-blank content lines for a named LAMMPS section (e.g., "Atoms", "Bonds").
    for i, line in enumerate(lines):
        if line.strip().startswith(key):
            out = []
            j = i + 1
            # Skip blank lines after header
            while j < len(lines) and not lines[j].strip():
                j += 1
            # Collect non-blank lines
            while j < len(lines) and lines[j].strip():
                out.append(lines[j])
                j += 1
            return out
    return []


def parse_box_from_header(lines):
    # Read x/y/z box bounds from the header and return box lengths in nm.
    xlo = xhi = ylo = yhi = zlo = zhi = None
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Masses"):
            break
        t = s.split()
        if len(t) == 4:
            if t[2] == "xlo" and t[3] == "xhi":
                xlo, xhi = float(t[0]), float(t[1])
            elif t[2] == "ylo" and t[3] == "yhi":
                ylo, yhi = float(t[0]), float(t[1])
            elif t[2] == "zlo" and t[3] == "zhi":
                zlo, zhi = float(t[0]), float(t[1])

    if None in (xlo, xhi, ylo, yhi, zlo, zhi):
        raise RuntimeError("Box bounds (xlo/xhi/ylo/yhi/zlo/zhi) not found in data file")

    # Convert Angstrom to nm
    return (
        (xhi - xlo) * ANGSTROM_TO_NM,
        (yhi - ylo) * ANGSTROM_TO_NM,
        (zhi - zlo) * ANGSTROM_TO_NM
    )


def normalize_blocks(sizes, names):
    # Normalize residue block definitions to tuples and validate consistent lengths.
    if isinstance(sizes, int):
        sizes = (sizes,)
    if isinstance(names, str):
        names = (names,)
    if len(sizes) != len(names):
        raise ValueError("residue_block_sizes and residue_block_names must have same length")
    return tuple(map(int, sizes)), tuple(map(str, names))


def build_residue_assignment(n_atoms, sizes, names):
    # Assign each atom ID to (resnr, resname, local_index) using repeating residue blocks.
    assignment = {}
    aid, resnr, block_idx = 1, 1, 0
    while aid <= n_atoms:
        size = sizes[block_idx % len(sizes)]
        name = names[block_idx % len(names)]
        for local in range(1, size + 1):
            if aid > n_atoms:
                break
            assignment[aid] = (resnr, name, local)
            aid += 1
        resnr += 1
        block_idx += 1
    return assignment


def find_1_4_pairs_bfs(n_atoms, edges):
    # Compute 1-4 pairs: atom pairs separated by exactly 3 bonds, using a BFS per atom.
    # Build adjacency list
    adj = [[] for _ in range(n_atoms + 1)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)

    pairs = []
    for start in range(1, n_atoms + 1):
        dist = [-1] * (n_atoms + 1)
        dist[start] = 0
        queue = deque([start])

        while queue:
            u = queue.popleft()
            if dist[u] == 3:
                continue
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = dist[u] + 1
                    queue.append(v)

        # Collect pairs where distance is exactly 3
        for v in range(start + 1, n_atoms + 1):
            if dist[v] == 3:
                pairs.append((start, v))

    return pairs


def find_1_4_pairs_networkx(n_atoms, edges):
    # Alternative 1-4 pair finder using NetworkX shortest path lengths (requires `networkx`).
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(1, n_atoms + 1))
    G.add_edges_from(edges)

    pairs = set()
    for u in range(1, n_atoms + 1):
        sp = nx.single_source_shortest_path_length(G, u, cutoff=3)
        for v, d in sp.items():
            if d == 3 and u < v:
                pairs.add((u, v))

    return sorted(pairs)


def validate_atom_types(masses, type_name_map):
    # Ensure every LAMMPS atom type ID (from Masses) maps to a GROMACS atomtype name.
    missing = set(masses.keys()) - set(type_name_map.keys())
    if missing:
        raise ValueError(
            f"type_name_map missing atom types: {sorted(missing)}\n"
            f"Add mappings for these LAMMPS type IDs to type_name_map."
        )


def _sorted_copy_path(infile):
    # Place `sorted_<basename>` alongside the input file.
    in_dir = os.path.dirname(os.path.abspath(infile))
    base = os.path.basename(infile)
    return os.path.join(in_dir, f"sorted_{base}")


def _find_section_bounds(lines, key):
    # Return (start_idx, end_idx) of the content region for a section key, excluding header/blank lines.
    for i, line in enumerate(lines):
        if line.strip().startswith(key):
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            start = j
            while j < len(lines) and lines[j].strip():
                j += 1
            end = j
            return start, end
    return None


def write_sorted_lammps_data_copy(infile, preview_lines=3):
    # Write `sorted_<infile>` where key sections are sorted by their leading numeric IDs.
    with open(infile) as f:
        lines = f.readlines()

    section_keys = [
        "Masses",
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "Atoms",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
    ]

    preview = {}
    for key in section_keys:
        bounds = _find_section_bounds(lines, key)
        if bounds is None:
            continue
        start, end = bounds
        segment = lines[start:end]

        keep = []
        sortable = []
        for line in segment:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                keep.append(line)
                continue
            first = stripped.split()[0]
            try:
                idx = int(first)
            except ValueError:
                keep.append(line)
                continue
            sortable.append((idx, line))

        sortable.sort(key=lambda x: x[0])
        sorted_segment = keep + [line for _, line in sortable]
        lines[start:end] = sorted_segment

        if preview_lines and sorted_segment:
            preview[key] = [l.rstrip("\n") for l in sorted_segment[:preview_lines]]

    outpath = _sorted_copy_path(infile)
    with open(outpath, "w") as f:
        f.write("".join(lines))

    logger.info(f"Wrote {outpath}")
    for key, plines in preview.items():
        logger.info(f"Preview ({key}):")
        for pl in plines:
            logger.info(f"  {pl}")

    return outpath


# ============================================================================
# MAIN CONVERSION FUNCTION
# ============================================================================

def build_gromacs_itps_and_gro_from_lammps(
    infile,
    # Type mappings (use defaults if None)
    type_name_map=None,
    element_map=None,
    # Output files
    forcefield_outfile="forcefield.itp",
    itp_outfile=None,
    topol_outfile="topol.top",
    gro_outfile="system.gro",
    # System parameters
    system_name="SYSTEM",
    molecule_name="MOL",
    nrexcl=3,
    defaults=(1, 2, "yes", 0.5, 0.8333),
    # Residue layout (None = single residue for all atoms)
    residue_block_sizes=None,
    residue_block_names=None,
    # Function types for bonded interactions (GROMACS funct column)
    bond_funct=1,      # 1=harmonic, 2=G96, 3=Morse, etc.
    angle_funct=1,     # 1=harmonic, 5=Urey-Bradley, etc.
    dihedral_funct=1,  # 1=proper, 3=Ryckaert-Bellemans, 9=multi, etc.
    improper_funct=4,  # 2=harmonic, 4=periodic
    pairs_funct=1,     # 1=LJ-14, 2=excluded
    pairs_method="bfs",
    # Validation
    total_charge_tolerance=1e-6,
    # Side outputs
    write_sorted_infile_copy=True,
):
    # End-to-end conversion pipeline:
    # 1) Optionally write `sorted_<infile>` (sections sorted by numeric IDs)
    # 2) Parse header + sections (Masses/Atoms + bonded sections + coeff sections)
    # 3) Convert units (Å->nm, kcal/mol->kJ/mol) and validate total charge
    # 4) Generate bonded terms + 1-4 pairs
    # 5) Write GROMACS outputs: forcefield.itp, <molecule>.itp, topol.top, system.gro
    if itp_outfile is None:
        itp_outfile = f"{molecule_name}.itp"

    # Use defaults if not provided
    if type_name_map is None:
        type_name_map = DEFAULT_TYPE_NAME_MAP.copy()
    if element_map is None:
        element_map = DEFAULT_ELEMENT_MAP.copy()

    logger.info(f"Reading LAMMPS data file: {infile}")

    if pairs_method not in ("bfs", "networkx"):
        raise ValueError("pairs_method must be 'bfs' or 'networkx'")

    if write_sorted_infile_copy:
        write_sorted_lammps_data_copy(infile)

    # -------------------- Read file --------------------
    with open(infile) as f:
        lines = f.readlines()

    box_nm = parse_box_from_header(lines)
    logger.info(f"Box dimensions (nm): {box_nm[0]:.3f} x {box_nm[1]:.3f} x {box_nm[2]:.3f}")

    # Read sections
    masses_l = read_section(lines, "Masses")
    atoms_l = read_section(lines, "Atoms")
    bonds_l = read_section(lines, "Bonds")
    angles_l = read_section(lines, "Angles")
    diheds_l = read_section(lines, "Dihedrals")
    impropers_l = read_section(lines, "Impropers")
    pairc_l = read_section(lines, "Pair Coeffs")
    bondc_l = read_section(lines, "Bond Coeffs")
    anglec_l = read_section(lines, "Angle Coeffs")
    dihc_l = read_section(lines, "Dihedral Coeffs")
    improperc_l = read_section(lines, "Improper Coeffs")

    if not atoms_l or not masses_l:
        raise RuntimeError("Atoms or Masses section missing from data file")

    # -------------------- Parse masses --------------------
    masses = {int(l.split()[0]): float(l.split()[1]) for l in masses_l}
    logger.info(f"Found {len(masses)} atom types")

    # Validate atom types
    validate_atom_types(masses, type_name_map)

    # -------------------- Parse atoms --------------------
    atoms = []
    qtot = 0.0
    for l in atoms_l:
        t = l.split()
        aid, mol, atype = int(t[0]), int(t[1]), int(t[2])
        q = float(t[3])
        x, y, z = float(t[4]) * ANGSTROM_TO_NM, float(t[5]) * ANGSTROM_TO_NM, float(t[6]) * ANGSTROM_TO_NM
        atoms.append((aid, mol, atype, q, x, y, z))
        qtot += q

    if abs(qtot) > total_charge_tolerance:
        raise RuntimeError(f"Total charge is not zero: {qtot:.6f}")

    atoms.sort()
    n_atoms = len(atoms)
    logger.info(f"Found {n_atoms} atoms, total charge: {qtot:.6e}")

    # -------------------- Residue assignment --------------------
    if (residue_block_sizes is None) != (residue_block_names is None):
        raise ValueError("residue_block_sizes and residue_block_names must be provided together")

    if residue_block_sizes is None and residue_block_names is None:
        residue_block_sizes = (n_atoms,)
        residue_block_names = (molecule_name[:3],)

    sizes, names = normalize_blocks(residue_block_sizes, residue_block_names)
    res_assign = build_residue_assignment(n_atoms, sizes, names)

    # -------------------- Parse coefficients --------------------
    # Pair coeffs: type -> (sigma_nm, epsilon_kJ)
    paircoeff = {
        int(l.split()[0]): (
            float(l.split()[2]) * ANGSTROM_TO_NM,
            float(l.split()[1]) * KCAL_TO_KJ
        )
        for l in pairc_l
    }

    # Bond coeffs: type -> (r0_nm, k_kJ_nm2)
    bondcoeff = {
        int(l.split()[0]): (
            float(l.split()[2]) * ANGSTROM_TO_NM,
            float(l.split()[1]) * KCAL_TO_KJ * 100 * 2  # kcal/mol/A^2 -> kJ/mol/nm^2
        )
        for l in bondc_l
    }

    # Angle coeffs: type -> (theta0_deg, k_kJ_rad2)
    anglecoeff = {
        int(l.split()[0]): (
            float(l.split()[2]),
            float(l.split()[1]) * KCAL_TO_KJ * 2
        )
        for l in anglec_l
    }

    # Dihedral coeffs: type -> [(phi, k, mult), ...]
    dihcoeff = {}
    for l in dihc_l:
        t = l.split()
        dtype = int(t[0])
        nterms = int(t[1])
        terms = []
        i = 2
        for _ in range(nterms):
            k = float(t[i]) * KCAL_TO_KJ
            mult = int(t[i + 1])
            phi = float(t[i + 2])
            terms.append((phi, k, mult))
            i += 3
        dihcoeff[dtype] = terms

    # Improper coeffs: type -> (k, d, n) for cvff style or (k, chi0) for harmonic
    # Assuming cvff style: K d n (d=1 or -1, n=multiplicity)
    impropercoeff = {}
    for l in improperc_l:
        t = l.split()
        itype = int(t[0])
        k = float(t[1]) * KCAL_TO_KJ
        d = int(t[2])  # +1 or -1
        n = int(t[3])  # multiplicity
        # For GROMACS periodic improper (funct=4): phi_s = 0 or 180 based on d
        phi_s = 0.0 if d > 0 else 180.0
        impropercoeff[itype] = (phi_s, k, n)

    # -------------------- Parse bonds --------------------
    bonds, edges = [], []
    for l in bonds_l:
        t = l.split()
        bt, i, j = int(t[1]), int(t[2]), int(t[3])
        bonds.append((i, j, bt))
        edges.append((i, j))

    logger.info(f"Found {len(bonds)} bonds, {len(angles_l)} angles, {len(diheds_l)} dihedrals, {len(impropers_l)} impropers")

    # -------------------- Find 1-4 pairs --------------------
    logger.info(f"Finding 1-4 pairs using {pairs_method} method...")
    if pairs_method == "bfs":
        pairs = find_1_4_pairs_bfs(n_atoms, edges)
    else:
        try:
            pairs = find_1_4_pairs_networkx(n_atoms, edges)
        except ModuleNotFoundError as e:
            raise RuntimeError(
                "pairs_method='networkx' requires the 'networkx' package to be installed."
            ) from e
    logger.info(f"Found {len(pairs)} 1-4 pairs")

    # ==========================================================
    # WRITE OUTPUT FILES
    # ==========================================================

    # -------------------- forcefield.itp --------------------
    ff_lines = [
        "[ defaults ]\n",
        "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ\n",
        " ".join(map(str, defaults)) + "\n\n",
        "[ atomtypes ]\n",
        "; name at.num mass charge ptype sigma epsilon\n"
    ]

    seen_by_name = {}
    for tid in sorted(masses):
        name = type_name_map[tid]
        atnum = get_atomic_number(name, element_map)
        sig, eps = paircoeff[tid]
        key = (atnum, masses[tid], sig, eps)

        if name in seen_by_name:
            if seen_by_name[name] != key:
                raise RuntimeError(
                    f"Inconsistent parameters for duplicate atomtype '{name}'.\n"
                    f"LAMMPS types mapping to same name must have identical parameters."
                )
            continue

        seen_by_name[name] = key
        ff_lines.append(f"{name:6s} {atnum:5d} {masses[tid]:10.5f} 0.0  A "
                       f"{sig:10.6f} {eps:10.6f}\n")

    with open(forcefield_outfile, "w") as f:
        f.write("".join(ff_lines))
    logger.info(f"Wrote {forcefield_outfile}")

    # -------------------- molecule.itp --------------------
    mol_lines = [
        "[ moleculetype ]\n",
        "; name nrexcl\n",
        f"{molecule_name} {nrexcl}\n\n",
        "[ atoms ]\n",
        "; nr type resnr residue atom cgnr charge mass\n"
    ]

    for aid, molid, atype, q, *_ in atoms:
        resnr, resname, local = res_assign[aid]
        aname = f"{resname[0]}{local}"
        mol_lines.append(
            f"{aid:5d} {type_name_map[atype]:6s} {resnr:5d} "
            f"{resname:6s} {aname:6s} {aid:5d} {q: .6f} {masses[atype]:.5f}\n"
        )

    mol_lines.append("\n[ bonds ]\n; ai aj funct r k\n")
    for i, j, bt in bonds:
        r0, k = bondcoeff[bt]
        mol_lines.append(f"{i:5d} {j:5d} {bond_funct} {r0:.6f} {k:.6f} ; type {bt}\n")

    mol_lines.append("\n[ pairs ]\n; ai aj funct\n")
    for i, j in pairs:
        mol_lines.append(f"{i:5d} {j:5d} {pairs_funct}\n")

    mol_lines.append("\n[ angles ]\n; ai aj ak funct theta k\n")
    for l in angles_l:
        t = l.split()
        at, i, j, k = int(t[1]), int(t[2]), int(t[3]), int(t[4])
        th, ka = anglecoeff[at]
        mol_lines.append(f"{i:5d} {j:5d} {k:5d} {angle_funct} {th:.5f} {ka:.6f} ; type {at}\n")

    mol_lines.append("\n[ dihedrals ]\n; ai aj ak al funct phi k mult\n")
    for l in diheds_l:
        t = l.split()
        dt, i, j, k, m = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])
        for phi, kd, mult in dihcoeff[dt]:
            mol_lines.append(
                f"{i:5d} {j:5d} {k:5d} {m:5d} {dihedral_funct} {phi:.6f} {kd:.6f} {mult} ; type {dt}\n"
            )

    # Impropers (second [dihedrals] section with improper function type)
    if impropers_l and impropercoeff:
        mol_lines.append("\n[ dihedrals ] ; impropers\n; ai aj ak al funct phi k mult\n")
        for l in impropers_l:
            t = l.split()
            it, i, j, k, m = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])
            phi_s, ki, mult = impropercoeff[it]
            mol_lines.append(
                f"{i:5d} {j:5d} {k:5d} {m:5d} {improper_funct} {phi_s:.6f} {ki:.6f} {mult} ; improper type {it}\n"
            )
    else:
        logger.warning("No Impropers section found in LAMMPS data file - skipping impropers in ITP")

    with open(itp_outfile, "w") as f:
        f.write("".join(mol_lines))
    logger.info(f"Wrote {itp_outfile}")

    # -------------------- topol.top --------------------
    with open(topol_outfile, "w") as f:
        f.write(f'; Topology generated by lammps_to_gro.py\n\n')
        f.write(f'#include "{forcefield_outfile}"\n')
        f.write(f'#include "{itp_outfile}"\n\n')
        f.write("[ system ]\n")
        f.write(f"{system_name}\n\n")
        f.write("[ molecules ]\n")
        f.write(f"{molecule_name} 1\n")
    logger.info(f"Wrote {topol_outfile}")

    # -------------------- system.gro --------------------
    lx, ly, lz = box_nm
    with open(gro_outfile, "w") as f:
        f.write(f"{molecule_name}\n")
        f.write(f"{n_atoms:5d}\n")
        for aid, molid, atype, q, x, y, z in atoms:
            resnr, resname, local = res_assign[aid]
            aname = f"{resname[0]}{local}"
            f.write(
                f"{resnr % 100000:5d}{resname:>5s}{aname:>5s}"
                f"{aid % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )
        f.write(f"{lx:10.5f}{ly:10.5f}{lz:10.5f}\n")
    logger.info(f"Wrote {gro_outfile}")

    logger.info("Conversion complete!")


# ============================================================================
# CLI INTERFACE
# ============================================================================

def parse_args(argv=None):
    # CLI argument parsing for direct script usage.
    if argv is None:
        argv = sys.argv[1:]

    # Allow subcommand-style invocation (e.g., "lammps-to-gro ...") by dropping
    # the subcommand token if present.
    if argv and argv[0] in {"lammps-to-gro", "lammps_to_gro"}:
        argv = argv[1:]

    parser = argparse.ArgumentParser(
        description="Convert LAMMPS data file to GROMACS topology and coordinate files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Direct script usage
  %(prog)s --infile data.aftercrosslink_coeff --molecule PVA_GLU_NETWORK --system wetting \\
           --residue-sizes 73 73 31 --residue-names PVA PVA GLU --pairs-method bfs

  # Subcommand-style (if wired into a parent CLI)
  lammps-to-gro --infile data.aftercrosslink_coeff --molecule PVA_GLU_NETWORK --system wetting \\
                --residue-sizes 73 73 31 --residue-names PVA PVA GLU --pairs-method bfs
        """
    )

    # Required arguments
    parser.add_argument(
        "--infile", "-i",
        required=True,
        help="Input LAMMPS data file with force field parameters (coeff sections)"
    )

    # Output options
    parser.add_argument(
        "--molecule", "-m",
        default="MOL",
        help="Molecule name for topology (default: MOL)"
    )
    parser.add_argument(
        "--system", "-s",
        default="SYSTEM",
        help="System name (default: SYSTEM)"
    )
    parser.add_argument(
        "--forcefield-out",
        default="forcefield.itp",
        help="Output forcefield file (default: forcefield.itp)"
    )
    parser.add_argument(
        "--itp-out",
        default=None,
        help="Output molecule topology file (default: <molecule>.itp)"
    )
    parser.add_argument(
        "--topol-out",
        default="topol.top",
        help="Output system topology file (default: topol.top)"
    )
    parser.add_argument(
        "--gro-out",
        default="system.gro",
        help="Output coordinate file (default: system.gro)"
    )

    # Residue options
    parser.add_argument(
        "--residue-sizes",
        type=int,
        nargs="+",
        default=None,
        help="Atoms per residue type (e.g., 146 31)"
    )
    parser.add_argument(
        "--residue-names",
        nargs="+",
        default=None,
        help="Residue type names (e.g., RES1 RES2)"
    )

    # Function types for bonded interactions
    parser.add_argument(
        "--bond-funct",
        type=int,
        default=1,
        help="GROMACS bond function type: 1=harmonic, 2=G96, 3=Morse (default: 1)"
    )
    parser.add_argument(
        "--angle-funct",
        type=int,
        default=1,
        help="GROMACS angle function type: 1=harmonic, 5=Urey-Bradley (default: 1)"
    )
    parser.add_argument(
        "--dihedral-funct",
        type=int,
        default=1,
        help="GROMACS dihedral function type: 1=proper, 3=RB, 9=multi (default: 1)"
    )
    parser.add_argument(
        "--improper-funct",
        type=int,
        default=4,
        help="GROMACS improper function type: 2=harmonic, 4=periodic (default: 4)"
    )
    parser.add_argument(
        "--pairs-funct",
        type=int,
        default=1,
        help="GROMACS pairs function type: 1=LJ-14 (default: 1)"
    )

    # Algorithm options
    parser.add_argument(
        "--pairs-method",
        choices=["bfs", "networkx"],
        default="bfs",
        help="Algorithm for finding 1-4 pairs (default: bfs)"
    )

    # Verbosity
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress info messages"
    )

    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    if args.quiet:
        logger.setLevel(logging.WARNING)

    # Handle itp output filename
    itp_out = args.itp_out if args.itp_out else f"{args.molecule}.itp"
    
    # Handle gro output
    gro_out = args.gro_out

    # Handle residue blocks
    residue_sizes = tuple(args.residue_sizes) if args.residue_sizes else None
    residue_names = tuple(args.residue_names) if args.residue_names else None

    try:
        # Run conversion
        build_gromacs_itps_and_gro_from_lammps(
            infile=args.infile,
            forcefield_outfile=args.forcefield_out,
            itp_outfile=itp_out,
            topol_outfile=args.topol_out,
            gro_outfile=gro_out,
            molecule_name=args.molecule,
            system_name=args.system,
            residue_block_sizes=residue_sizes,
            residue_block_names=residue_names,
            bond_funct=args.bond_funct,
            angle_funct=args.angle_funct,
            dihedral_funct=args.dihedral_funct,
            improper_funct=args.improper_funct,
            pairs_funct=args.pairs_funct,
            pairs_method=args.pairs_method,
        )

        # Print summary
        print("\nGenerated files:")
        print(f"  {_sorted_copy_path(args.infile)}")
        print(f"  {args.forcefield_out}")
        print(f"  {itp_out}")
        print(f"  {args.topol_out}")
        print(f"  {gro_out}")

        print("\nNext steps (GROMACS commands):")
        print(f"  gmx grompp -f md.mdp -c {gro_out} -p {args.topol_out} -o md.tpr")
        print("  gmx mdrun -v -deffnm md")

    except FileNotFoundError as e:
        print(f"\nERROR: Input file not found\n  {e}")
    except ValueError as e:
        print(f"\nERROR: Invalid input\n  {e}")
    except RuntimeError as e:
        print(f"\nERROR: Conversion failed\n  {e}")
    except Exception as e:
        print(f"\nERROR: Unexpected error\n  {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
