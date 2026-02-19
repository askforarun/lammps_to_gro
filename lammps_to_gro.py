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

# Mass ranges for element detection (atomic mass units)
MASS_TO_ELEMENT = [
    (1.0, 1.1, "H"),      # Hydrogen ~1.008
    (12.0, 12.1, "C"),    # Carbon ~12.011
    (14.0, 14.1, "N"),    # Nitrogen ~14.007
    (15.9, 16.1, "O"),    # Oxygen ~15.999
    (19.0, 19.1, "F"),    # Fluorine ~18.998
    (30.9, 31.1, "P"),    # Phosphorus ~30.974
    (32.0, 32.1, "S"),    # Sulfur ~32.065
    (35.4, 35.5, "CL"),   # Chlorine ~35.453
    (79.9, 80.0, "BR"),   # Bromine ~79.904
    (126.9, 127.0, "I"),  # Iodine ~126.904
]


def detect_element_from_mass(mass):
    """Detect element symbol from atomic mass."""
    for low, high, element in MASS_TO_ELEMENT:
        if low <= mass <= high:
            return element
    raise ValueError(f"Cannot detect element from mass {mass}. Add entry to MASS_TO_ELEMENT.")


def auto_generate_type_maps(masses):
    """Auto-generate type_name_map and element_map from masses.
    
    Creates unique atomtype names like C1, C2, H1, H2 based on element and type ID.
    """
    type_name_map = {}
    element_map = {}
    element_counts = {}
    
    for tid in sorted(masses.keys()):
        mass = masses[tid]
        element = detect_element_from_mass(mass)
        
        # Count occurrences of each element to create unique names
        if element not in element_counts:
            element_counts[element] = 0
        element_counts[element] += 1
        
        # Create atomtype name like C1, C2, H1, H2
        atomtype_name = f"{element.lower()}{element_counts[element]}"
        type_name_map[tid] = atomtype_name
        element_map[atomtype_name] = element
    
    return type_name_map, element_map

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

# Known LAMMPS section headers for robust section scanning
SECTION_HEADERS = (
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
    "Velocities",
)

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


def get_element_symbol(atom_type, type_name_map, element_map):
    # Get element symbol for a given atom type
    atomtype_name = type_name_map.get(atom_type, None)
    if atomtype_name is None:
        raise ValueError(f"Missing type_name_map entry for atom type {atom_type}")
    return element_map[atomtype_name].upper()


def generate_atom_name(resname, local_index, atom_type, type_name_map, element_map):
    """Generate atom name using just the element symbol.

    Format: element only (H, C, O, N, CL, etc.)
    Simple and always fits within 5-char GRO limit.
    GROMACS uses atom indices internally, so duplicate names are fine.
    """
    try:
        element = get_element_symbol(atom_type, type_name_map, element_map)
    except (ValueError, KeyError):
        element = resname[0] if resname else "X"
    
    return element


def format_gro_name(name, width=5, keep_tail=False):
    # Normalize names for strict fixed-width GRO fields.
    # Keep underscores so atom names stay consistent between .itp and .gro.
    s = "".join(ch for ch in str(name) if ch.isalnum() or ch == "_") or "X"
    if len(s) > width:
        s = s[-width:] if keep_tail else s[:width]
    return f"{s:>{width}s}"


def is_section_header_line(line):
    # Return True if a stripped line looks like a LAMMPS section header.
    s = line.strip()
    return any(s == h or s.startswith(h + " ") for h in SECTION_HEADERS)


def read_section(lines, key):
    # Return the non-blank content lines for a named LAMMPS section (e.g., "Atoms", "Bonds").
    for i, line in enumerate(lines):
        if line.strip().startswith(key):
            out = []
            j = i + 1
            # Skip blank lines after header
            while j < len(lines) and not lines[j].strip():
                j += 1
            # Collect lines until the next section header; allow internal blank lines.
            while j < len(lines):
                s = lines[j].strip()
                if is_section_header_line(s):
                    break
                if s:
                    out.append(lines[j])
                j += 1
            return out
    return []


def iter_data_lines(section_lines):
    # Yield only non-empty, non-comment lines from a section.
    for line in section_lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        yield s


def has_data_lines(section_lines):
    # True when a section contains at least one non-comment data line.
    return any(True for _ in iter_data_lines(section_lines))


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


def build_residue_assignment(n_atoms, sizes, names, repeat=False):
    """
    Assign each atom ID to (resnr, resname, local_index).
    If repeat is True, the residue block pattern is repeated until all atoms are assigned.
    Otherwise, the blocks are consumed once; assumes sum(sizes) == n_atoms.
    """
    assignment = {}
    aid, resnr, block_idx = 1, 1, 0

    if repeat:
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
    else:
        for size, name in zip(sizes, names):
            for local in range(1, size + 1):
                assignment[aid] = (resnr, name, local)
                aid += 1
            resnr += 1

    return assignment


def find_1_4_pairs_bfs(n_atoms, edges):
    # Compute 1-4 pairs: atom pairs separated by exactly 3 bonds, using a BFS per atom.
    # Optimized to avoid O(N) scans per atom by collecting distance-3 nodes on the fly.
    adj = [[] for _ in range(n_atoms + 1)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)

    pairs = []
    seen = [0] * (n_atoms + 1)
    dist = [0] * (n_atoms + 1)
    epoch = 0

    for start in range(1, n_atoms + 1):
        epoch += 1
        queue = deque([start])
        seen[start] = epoch
        dist[start] = 0

        while queue:
            u = queue.popleft()
            if dist[u] == 3:
                continue
            for v in adj[u]:
                if seen[v] == epoch:
                    continue
                seen[v] = epoch
                dist[v] = dist[u] + 1
                if dist[v] == 3:
                    if v > start:
                        pairs.append((start, v))
                else:
                    queue.append(v)

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


def validate_header_box(lines):
    """Validate that box dimensions are present and properly formatted."""
    box_found = {"x": False, "y": False, "z": False}
    
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Masses"):
            break
        
        t = s.split()
        if len(t) == 4:
            if t[2] == "xlo" and t[3] == "xhi":
                try:
                    xlo, xhi = float(t[0]), float(t[1])
                    if xhi <= xlo:
                        raise ValueError(f"Invalid box dimensions: xhi ({xhi}) <= xlo ({xlo})")
                    box_found["x"] = True
                except ValueError:
                    raise ValueError(f"Invalid box dimension format: {s}")
            elif t[2] == "ylo" and t[3] == "yhi":
                try:
                    ylo, yhi = float(t[0]), float(t[1])
                    if yhi <= ylo:
                        raise ValueError(f"Invalid box dimensions: yhi ({yhi}) <= ylo ({ylo})")
                    box_found["y"] = True
                except ValueError:
                    raise ValueError(f"Invalid box dimension format: {s}")
            elif t[2] == "zlo" and t[3] == "zhi":
                try:
                    zlo, zhi = float(t[0]), float(t[1])
                    if zhi <= zlo:
                        raise ValueError(f"Invalid box dimensions: zhi ({zhi}) <= zlo ({zlo})")
                    box_found["z"] = True
                except ValueError:
                    raise ValueError(f"Invalid box dimension format: {s}")
    
    missing = [dim for dim, found in box_found.items() if not found]
    if missing:
        raise ValueError(f"Missing box dimensions for: {', '.join(missing)}")
    
    return True


def validate_section_format(lines, section_name, expected_columns=None):
    """Validate that a section has the correct format and expected number of columns."""
    section_lines = read_section(lines, section_name)
    
    if not section_lines:
        return f"Section '{section_name}' not found"
    
    errors = []
    for i, line in enumerate(section_lines, 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
            
        parts = stripped.split()
        if expected_columns and len(parts) < expected_columns:
            errors.append(f"{section_name} line {i}: expected >= {expected_columns} columns, got {len(parts)}: '{stripped}'")
            
        # Validate numeric columns (skip first column which is usually an ID)
        if section_name == "Masses":
            try:
                int(parts[0])  # atom type ID
                float(parts[1])  # mass
            except (ValueError, IndexError):
                errors.append(f"{section_name} line {i}: invalid mass format: '{stripped}' (hint: id mass)")
        
        elif section_name == "Atoms":
            try:
                int(parts[0])  # atom ID
                int(parts[1])  # molecule ID
                int(parts[2])  # atom type
                float(parts[3])  # charge
                float(parts[4])  # x
                float(parts[5])  # y
                float(parts[6])  # z
            except (ValueError, IndexError):
                errors.append(f"{section_name} line {i}: invalid atom format in first 7 columns: '{stripped}' (hint: id mol type q x y z)")
    
    return errors if errors else True


def validate_pair_coeffs(section_lines):
    """Validate Pair Coeffs lines: id epsilon sigma (>=3 columns)."""
    errors = []
    for i, line in enumerate(section_lines, 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 3:
            errors.append(f"Pair Coeffs line {i}: expected >= 3 columns (id epsilon sigma), got {len(parts)}: '{s}'")
            continue
        try:
            int(parts[0])
            float(parts[1])
            float(parts[2])
        except ValueError:
            errors.append(f"Pair Coeffs line {i}: non-numeric epsilon/sigma: '{s}' (hint: id epsilon sigma)")
    return errors if errors else True


def validate_bond_coeffs(section_lines):
    """Validate Bond Coeffs lines: id k r0 (>=3 columns)."""
    errors = []
    for i, line in enumerate(section_lines, 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 3:
            errors.append(f"Bond Coeffs line {i}: expected >= 3 columns (id k r0), got {len(parts)}: '{s}'")
            continue
        try:
            int(parts[0])
            float(parts[1])
            float(parts[2])
        except ValueError:
            errors.append(f"Bond Coeffs line {i}: non-numeric k or r0: '{s}' (hint: id k r0)")
    return errors if errors else True


def validate_angle_coeffs(section_lines):
    """Validate Angle Coeffs lines: id k theta0 (>=3 columns)."""
    errors = []
    for i, line in enumerate(section_lines, 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 3:
            errors.append(f"Angle Coeffs line {i}: expected >= 3 columns (id k theta0), got {len(parts)}: '{s}'")
            continue
        try:
            int(parts[0])
            float(parts[1])
            float(parts[2])
        except ValueError:
            errors.append(f"Angle Coeffs line {i}: non-numeric k or theta0: '{s}' (hint: id k theta0)")
    return errors if errors else True


def convert_opls_to_rb(k1, k2, k3, k4):
    """Convert OPLS dihedral coefficients (K1-K4) to Ryckaert-Bellemans format (C0-C5).
    
    LAMMPS OPLS format uses K values that are directly equivalent to GROMACS F values.
    (The 0.5 factor in the LAMMPS potential formula is already accounted for in K.)
    
    RB format: V = C0 + C1*cos(psi) + C2*cos^2(psi) + C3*cos^3(psi) + C4*cos^4(psi) + C5*cos^5(psi)
    where psi = phi - 180 (RB uses trans=0 convention)
    
    Conversion formulas (from GROMACS manual):
    C0 = F2 + (F1 + F3)/2
    C1 = (-F1 + 3*F3)/2
    C2 = -F2 + 4*F4
    C3 = -2*F3
    C4 = -4*F4
    C5 = 0
    
    Where F1=K1, F2=K2, F3=K3, F4=K4
    """
    f1, f2, f3, f4 = k1, k2, k3, k4
    
    c0 = f2 + (f1 + f3) / 2
    c1 = (-f1 + 3 * f3) / 2
    c2 = -f2 + 4 * f4
    c3 = -2 * f3
    c4 = -4 * f4
    c5 = 0.0
    
    return c0, c1, c2, c3, c4, c5


def validate_dihedral_coeffs(section_lines, dihedral_funct=1):
    """Validate Dihedral Coeffs lines based on dihedral_funct parameter.
    
    If dihedral_funct=3: expects OPLS format (id k1 k2 k3 k4)
    Otherwise: expects periodic format (id nterms (k mult phi)*)
    """
    if not section_lines:
        return True
    
    errors = []
    if dihedral_funct == 3:
        # OPLS style: id k1 k2 k3 k4 (5 columns total)
        for i, line in enumerate(section_lines, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) != 5:
                errors.append(f"Dihedral Coeffs line {i}: OPLS style requires 5 columns (id k1 k2 k3 k4), got {len(parts)}: '{s}'")
                continue
            try:
                int(parts[0])  # dihedral type ID
                for j in range(1, 5):
                    float(parts[j])  # k1-k4
            except ValueError:
                errors.append(f"Dihedral Coeffs line {i}: non-numeric coefficients: '{s}' (hint: id k1 k2 k3 k4)")
    else:
        # Periodic style: id nterms (k mult phi)*
        for i, line in enumerate(section_lines, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                errors.append(f"Dihedral Coeffs line {i}: expected >= 2 columns (id nterms), got {len(parts)}: '{s}'")
                continue
            try:
                dtype = int(parts[0])
                nterms = int(parts[1])
            except ValueError:
                errors.append(f"Dihedral Coeffs line {i}: non-numeric id or nterms: '{s}'")
                continue

            expected = 2 + 3 * nterms
            if len(parts) < expected:
                errors.append(
                    f"Dihedral Coeffs line {i}: expected {expected} columns for {nterms} terms, got {len(parts)}: '{s}' (hint: id nterms k mult phi ...)"
                )
                continue

            idx = 2
            for term in range(nterms):
                try:
                    float(parts[idx])      # k
                    int(parts[idx + 1])    # mult
                    float(parts[idx + 2])  # phi
                except ValueError:
                    errors.append(
                        f"Dihedral Coeffs line {i}: non-numeric term {term + 1} (k mult phi) starting at column {idx + 1}: '{s}'"
                    )
                    break
                idx += 3
    return errors if errors else True


def validate_improper_coeffs(section_lines):
    """Validate Improper Coeffs lines: id k d n (>=4 columns)."""
    errors = []
    for i, line in enumerate(section_lines, 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 4:
            errors.append(f"Improper Coeffs line {i}: expected >= 4 columns (id k d n), got {len(parts)}: '{s}'")
            continue
        try:
            int(parts[0])
            float(parts[1])
            int(parts[2])
            int(parts[3])
        except ValueError:
            errors.append(f"Improper Coeffs line {i}: non-numeric k/d/n: '{s}' (hint: id k d n)")
    return errors if errors else True


def validate_coefficient_consistency(lines):
    """Validate that each bonded section has its matching coefficient section and vice versa."""
    mapping = {
        "Bonds": "Bond Coeffs",
        "Angles": "Angle Coeffs",
        "Dihedrals": "Dihedral Coeffs",
        "Impropers": "Improper Coeffs",
    }
    errors = []
    
    for bonded_section, coeff_section in mapping.items():
        bonded_present = has_data_lines(read_section(lines, bonded_section))
        coeff_present = has_data_lines(read_section(lines, coeff_section))
        
        if bonded_present and not coeff_present:
            errors.append(f"Found {bonded_section} section but missing {coeff_section} section")
        if coeff_present and not bonded_present:
            errors.append(f"Found {coeff_section} section but no corresponding {bonded_section} section")
    
    return errors if errors else True


def validate_residue_assignment(n_atoms, sizes, names, repeat=False):
    """Validate residue assignment parameters."""
    if sizes is None and names is None:
        return True  # Default assignment
    
    if len(sizes) != len(names):
        raise ValueError(f"residue_sizes ({len(sizes)}) and residue_names ({len(names)}) must have same length")
    
    total_atoms = sum(sizes)
    if total_atoms > n_atoms:
        raise ValueError(f"Residue assignment requires {total_atoms} atoms but only {n_atoms} atoms found")
    if total_atoms < n_atoms and not repeat:
        raise ValueError(
            "Residue assignment covers only "
            f"{total_atoms} atoms but {n_atoms} atoms found. "
            "Residue sizes must sum to exactly n_atoms, or use --repeat to extend pattern"
        )
    if total_atoms <= 0:
        raise ValueError("Residue sizes must sum to a positive number")
    
    # Check for invalid sizes
    for i, size in enumerate(sizes):
        if size <= 0:
            raise ValueError(f"Residue size {i+1} is invalid: {size} (must be > 0)")
    
    # Check for empty names
    for i, name in enumerate(names):
        if not name or not name.strip():
            raise ValueError(f"Residue name {i+1} is empty")
    
    # Warn if --repeat is redundant
    if repeat and total_atoms == n_atoms:
        logger.warning("--repeat is redundant: residue pattern size equals total atoms. Pattern will execute exactly once.")
    
    return True


def validate_atom_type_consistency(atoms, masses):
    """Validate that all atoms reference valid atom types."""
    atom_types = set()
    for atom_line in iter_data_lines(atoms):
        parts = atom_line.split()
        if len(parts) >= 3:
            try:
                atom_types.add(int(parts[2]))  # atom type column
            except ValueError:
                raise ValueError(f"Invalid atom type in line: '{atom_line}'")
    
    mass_types = set(masses.keys())
    
    # Check for atom types without mass definitions
    missing_masses = atom_types - mass_types
    if missing_masses:
        raise ValueError(f"Atoms reference undefined atom types: {sorted(missing_masses)}")
    
    # Check for mass definitions not used by any atoms
    unused_masses = mass_types - atom_types
    if unused_masses:
        logger.warning(f"Unused atom type definitions: {sorted(unused_masses)}")
    
    return True


def validate_bonded_connectivity(atoms, bonds, angles, dihedrals, impropers):
    """Validate that bonded interactions reference valid atom IDs."""
    if not atoms:
        return True
    
    atom_ids = set()
    for atom_line in iter_data_lines(atoms):
        parts = atom_line.split()
        if len(parts) >= 1:
            try:
                atom_ids.add(int(parts[0]))
            except ValueError:
                raise ValueError(f"Invalid atom ID in line: '{atom_line}'")
    
    def check_section(section_lines, section_name, num_atoms):
        if not section_lines:
            return True
        
        for i, line in enumerate(iter_data_lines(section_lines), 1):
            parts = line.split()
            # Format: interaction_id type_id atom1 atom2 [atom3] [atom4]
            expected_cols = 2 + num_atoms  # interaction_id + type_id + atom IDs
            if len(parts) < expected_cols:
                raise ValueError(f"{section_name} line {i}: Expected at least {expected_cols} columns, got {len(parts)}")
            
            # Check atom IDs (start at index 2, after interaction_id and type_id)
            for j in range(2, 2 + num_atoms):
                try:
                    atom_id = int(parts[j])
                except ValueError:
                    raise ValueError(f"{section_name} line {i}: Invalid atom ID '{parts[j]}'")
                if atom_id not in atom_ids:
                    raise ValueError(f"{section_name} line {i}: Atom ID {atom_id} not found in Atoms section")
        
        return True
    
    # Validate each bonded section
    check_section(bonds, "Bonds", 2)      # bond: 2 atom IDs
    check_section(angles, "Angles", 3)    # angle: 3 atom IDs
    check_section(dihedrals, "Dihedrals", 4)  # dihedral: 4 atom IDs
    check_section(impropers, "Impropers", 4)   # improper: 4 atom IDs
    
    return True


def validate_atom_types(masses, type_name_map):
    # Ensure every LAMMPS atom type ID (from Masses) maps to a GROMACS atomtype name.
    missing = set(masses.keys()) - set(type_name_map.keys())
    if missing:
        raise ValueError(
            f"type_name_map missing atom types: {sorted(missing)}\n"
            f"Add mappings for these LAMMPS type IDs to type_name_map."
        )


def validate_contiguous_atom_ids(atom_records):
    # Enforce contiguous 1..N atom IDs required by residue assignment and BFS pair generation.
    ids = sorted(aid for aid, *_ in atom_records)
    n_atoms = len(ids)
    expected = list(range(1, n_atoms + 1))
    if ids != expected:
        missing = sorted(set(expected) - set(ids))
        extras = sorted(set(ids) - set(expected))
        raise ValueError(
            "Atom IDs must be contiguous and 1-based (expected 1..N). "
            f"Found range {ids[0] if ids else 'NA'}..{ids[-1] if ids else 'NA'} "
            f"for N={n_atoms}. Missing IDs: {missing[:10]}"
            f"{'...' if len(missing) > 10 else ''}. "
            f"Out-of-range IDs: {extras[:10]}{'...' if len(extras) > 10 else ''}."
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

        pending = []
        sortable = []
        for line in segment:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                pending.append(line)
                continue
            first = stripped.split()[0]
            try:
                idx = int(first)
            except ValueError:
                pending.append(line)
                continue
            sortable.append((idx, pending + [line]))
            pending = []

        sortable.sort(key=lambda x: x[0])
        sorted_segment = []
        for _, block in sortable:
            sorted_segment.extend(block)
        sorted_segment.extend(pending)
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
    defaults=(1, 2, "yes", 0.5, 0.8333),  # comb-rule 2 = arithmetic mixing (matches LAMMPS mix arithmetic)
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
    repeat_residue_blocks=False,
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

    logger.info(f"Reading LAMMPS data file: {infile}")

    if pairs_method not in ("bfs", "networkx"):
        raise ValueError("pairs_method must be 'bfs' or 'networkx'")

    if write_sorted_infile_copy:
        sorted_infile = write_sorted_lammps_data_copy(infile)
    else:
        sorted_infile = infile

    # -------------------- Read file --------------------
    with open(sorted_infile) as f:
        lines = f.readlines()

    box_nm = parse_box_from_header(lines)
    logger.info(f"Box dimensions (nm): {box_nm[0]:.3f} x {box_nm[1]:.3f} x {box_nm[2]:.3f}")

    # -------------------- Early validation --------------------
    logger.info("Performing early validation...")
    
    # Validate header and box dimensions
    validate_header_box(lines)
    
    # Validate section formats
    masses_validation = validate_section_format(lines, "Masses", 2)
    if masses_validation is not True:
        raise ValueError(f"Masses section validation failed:\n" + "\n".join(masses_validation))
    
    atoms_validation = validate_section_format(lines, "Atoms", 7)
    if atoms_validation is not True:
        raise ValueError(f"Atoms section validation failed:\n" + "\n".join(atoms_validation))

    # Validate coefficient sections
    pairc_validation = validate_pair_coeffs(read_section(lines, "Pair Coeffs"))
    if pairc_validation is not True:
        raise ValueError(f"Pair Coeffs validation failed:\n" + "\n".join(pairc_validation))

    bondc_validation = validate_bond_coeffs(read_section(lines, "Bond Coeffs"))
    if bondc_validation is not True:
        raise ValueError(f"Bond Coeffs validation failed:\n" + "\n".join(bondc_validation))

    anglec_validation = validate_angle_coeffs(read_section(lines, "Angle Coeffs"))
    if anglec_validation is not True:
        raise ValueError(f"Angle Coeffs validation failed:\n" + "\n".join(anglec_validation))

    dihc_validation = validate_dihedral_coeffs(read_section(lines, "Dihedral Coeffs"), dihedral_funct)
    if dihc_validation is not True:
        raise ValueError(f"Dihedral Coeffs validation failed:\n" + "\n".join(dihc_validation))

    improperc_validation = validate_improper_coeffs(read_section(lines, "Improper Coeffs"))
    if improperc_validation is not True:
        raise ValueError(f"Improper Coeffs validation failed:\n" + "\n".join(improperc_validation))
    
    # Validate coefficient consistency
    coeff_validation = validate_coefficient_consistency(lines)
    if coeff_validation is not True:
        raise ValueError(f"Coefficient consistency validation failed:\n" + "\n".join(coeff_validation))
    
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

    # Required vs optional sections
    required_sections = {
        "Masses": has_data_lines(masses_l),
        "Atoms": has_data_lines(atoms_l),
        "Bonds": has_data_lines(bonds_l),
        "Angles": has_data_lines(angles_l),
        "Dihedrals": has_data_lines(diheds_l),
        "Pair Coeffs": has_data_lines(pairc_l),
        "Bond Coeffs": has_data_lines(bondc_l),
        "Angle Coeffs": has_data_lines(anglec_l),
        "Dihedral Coeffs": has_data_lines(dihc_l),
    }
    missing_required = [name for name, present in required_sections.items() if not present]
    if missing_required:
        raise RuntimeError("Required sections missing or empty: " + ", ".join(missing_required))

    optional_sections = {
        "Impropers": has_data_lines(impropers_l),
        "Improper Coeffs": has_data_lines(improperc_l),
    }
    missing_optional = [name for name, present in optional_sections.items() if not present]
    if missing_optional:
        logger.warning("Missing or empty optional sections: " + ", ".join(missing_optional))

    # -------------------- Parse masses --------------------
    masses = {}
    for s in iter_data_lines(masses_l):
        t = s.split()
        masses[int(t[0])] = float(t[1])
    logger.info(f"Found {len(masses)} atom types")

    # Auto-generate type maps from masses if not provided
    if type_name_map is None:
        type_name_map, auto_element_map = auto_generate_type_maps(masses)
        if element_map is None:
            element_map = auto_element_map
        logger.info(f"Auto-generated atomtype names: {type_name_map}")
    elif element_map is None:
        element_map = DEFAULT_ELEMENT_MAP.copy()

    # Validate atom type consistency
    validate_atom_type_consistency(atoms_l, masses)
    
    # Validate bonded connectivity
    validate_bonded_connectivity(atoms_l, bonds_l, angles_l, diheds_l, impropers_l)

    # Validate atom types
    validate_atom_types(masses, type_name_map)

    # -------------------- Parse atoms --------------------
    atoms = []
    qtot = 0.0
    for s in iter_data_lines(atoms_l):
        t = s.split()
        aid, mol, atype = int(t[0]), int(t[1]), int(t[2])
        q = float(t[3])
        x, y, z = float(t[4]) * ANGSTROM_TO_NM, float(t[5]) * ANGSTROM_TO_NM, float(t[6]) * ANGSTROM_TO_NM
        atoms.append((aid, mol, atype, q, x, y, z))
        qtot += q

    if abs(qtot) > total_charge_tolerance:
        raise RuntimeError(f"Total charge is not zero: {qtot:.6f}")

    atoms.sort()
    n_atoms = len(atoms)
    validate_contiguous_atom_ids(atoms)
    logger.info(f"Found {n_atoms} atoms, total charge: {qtot:.6e}")

    # -------------------- Residue assignment --------------------
    if residue_block_sizes is None or residue_block_names is None:
        raise ValueError("residue_block_sizes and residue_block_names are required")

    sizes, names = normalize_blocks(residue_block_sizes, residue_block_names)
    
    # Validate residue assignment
    validate_residue_assignment(n_atoms, sizes, names, repeat=repeat_residue_blocks)
    
    res_assign = build_residue_assignment(n_atoms, sizes, names, repeat=repeat_residue_blocks)

    # -------------------- Parse coefficients --------------------
    # Pair coeffs: type -> (sigma_nm, epsilon_kJ)
    paircoeff = {}
    for s in iter_data_lines(pairc_l):
        t = s.split()
        paircoeff[int(t[0])] = (
            float(t[2]) * ANGSTROM_TO_NM,
            float(t[1]) * KCAL_TO_KJ
        )

    # Bond coeffs: type -> (r0_nm, k_kJ_nm2)
    bondcoeff = {}
    for s in iter_data_lines(bondc_l):
        t = s.split()
        bondcoeff[int(t[0])] = (
            float(t[2]) * ANGSTROM_TO_NM,
            float(t[1]) * KCAL_TO_KJ * 100 * 2  # kcal/mol/A^2 -> kJ/mol/nm^2
        )

    # Angle coeffs: type -> (theta0_deg, k_kJ_rad2)
    anglecoeff = {}
    for s in iter_data_lines(anglec_l):
        t = s.split()
        anglecoeff[int(t[0])] = (
            float(t[2]),
            float(t[1]) * KCAL_TO_KJ * 2
        )

    # Dihedral coeffs: handle based on dihedral_funct parameter
    dihcoeff = {}
    
    if dihedral_funct == 3:
        # RB function: expect OPLS format (type k1 k2 k3 k4) and convert to C0-C4
        logger.info("Using RB dihedral function - converting OPLS coefficients to Ryckaert-Bellemans")
        for s in iter_data_lines(dihc_l):
            t = s.split()
            dtype = int(t[0])
            k1 = float(t[1]) * KCAL_TO_KJ
            k2 = float(t[2]) * KCAL_TO_KJ
            k3 = float(t[3]) * KCAL_TO_KJ
            k4 = float(t[4]) * KCAL_TO_KJ
            # Convert to RB format
            c0, c1, c2, c3, c4, c5 = convert_opls_to_rb(k1, k2, k3, k4)
            dihcoeff[dtype] = (c0, c1, c2, c3, c4, c5)
    else:
        # Default periodic function: expect periodic format (id nterms (k mult phi)*)
        logger.info("Using periodic dihedral function")
        for s in iter_data_lines(dihc_l):
            t = s.split()
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
    for s in iter_data_lines(improperc_l):
        t = s.split()
        itype = int(t[0])
        k = float(t[1]) * KCAL_TO_KJ
        d = int(t[2])  # +1 or -1
        n = int(t[3])  # multiplicity
        # For GROMACS periodic improper (funct=4): phi_s = 0 or 180 based on d
        phi_s = 0.0 if d > 0 else 180.0
        impropercoeff[itype] = (phi_s, k, n)

    # Validate coefficient coverage
    missing_paircoeff = set(masses.keys()) - set(paircoeff.keys())
    if missing_paircoeff:
        raise RuntimeError(f"Missing Pair Coeffs for atom types: {sorted(missing_paircoeff)}")

    bond_types_in_use = {int(s.split()[1]) for s in iter_data_lines(bonds_l)}
    missing_bondcoeff = bond_types_in_use - set(bondcoeff.keys())
    if missing_bondcoeff:
        raise RuntimeError(f"Missing Bond Coeffs for bond types: {sorted(missing_bondcoeff)}")

    angle_types_in_use = {int(s.split()[1]) for s in iter_data_lines(angles_l)}
    missing_anglecoeff = angle_types_in_use - set(anglecoeff.keys())
    if missing_anglecoeff:
        raise RuntimeError(f"Missing Angle Coeffs for angle types: {sorted(missing_anglecoeff)}")

    dihedral_types_in_use = {int(s.split()[1]) for s in iter_data_lines(diheds_l)}
    missing_dihedralcoeff = dihedral_types_in_use - set(dihcoeff.keys())
    if missing_dihedralcoeff:
        raise RuntimeError(f"Missing Dihedral Coeffs for dihedral types: {sorted(missing_dihedralcoeff)}")

    if impropers_l:
        improper_types_in_use = {int(s.split()[1]) for s in iter_data_lines(impropers_l)}
        missing_impropercoeff = improper_types_in_use - set(impropercoeff.keys())
        if missing_impropercoeff:
            raise RuntimeError(f"Missing Improper Coeffs for improper types: {sorted(missing_impropercoeff)}")

    # -------------------- Parse bonds --------------------
    bonds, edges = [], []
    for s in iter_data_lines(bonds_l):
        t = s.split()
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
        aname = generate_atom_name(resname, local, atype, type_name_map, element_map)
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
    for s in iter_data_lines(angles_l):
        t = s.split()
        at, i, j, k = int(t[1]), int(t[2]), int(t[3]), int(t[4])
        th, ka = anglecoeff[at]
        mol_lines.append(f"{i:5d} {j:5d} {k:5d} {angle_funct} {th:.5f} {ka:.6f} ; type {at}\n")

    mol_lines.append("\n[ dihedrals ]\n")
    if dihedral_funct == 3:
        # RB format: ai aj ak al funct c0 c1 c2 c3 c4 c5
        mol_lines.append("; ai aj ak al funct c0 c1 c2 c3 c4 c5\n")
        for s in iter_data_lines(diheds_l):
            t = s.split()
            dt, i, j, k, m = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])
            c0, c1, c2, c3, c4, c5 = dihcoeff[dt]
            mol_lines.append(
                f"{i:5d} {j:5d} {k:5d} {m:5d} {dihedral_funct} {c0:.6f} {c1:.6f} {c2:.6f} {c3:.6f} {c4:.6f} {c5:.6f} ; type {dt}\n"
            )
    else:
        # Periodic format: ai aj ak al funct phi k mult
        mol_lines.append("; ai aj ak al funct phi k mult\n")
        for s in iter_data_lines(diheds_l):
            t = s.split()
            dt, i, j, k, m = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])
            for phi, kd, mult in dihcoeff[dt]:
                mol_lines.append(
                    f"{i:5d} {j:5d} {k:5d} {m:5d} {dihedral_funct} {phi:.6f} {kd:.6f} {mult} ; type {dt}\n"
                )

    # Impropers (second [dihedrals] section with improper function type)
    if impropers_l and impropercoeff:
        mol_lines.append("\n[ dihedrals ] ; impropers\n; ai aj ak al funct phi k mult\n")
        for s in iter_data_lines(impropers_l):
            t = s.split()
            it, i, j, k, m = int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5])
            phi_s, ki, mult = impropercoeff[it]
            mol_lines.append(
                f"{i:5d} {j:5d} {k:5d} {m:5d} {improper_funct} {phi_s:.6f} {ki:.6f} {mult} ; improper type {it}\n"
            )

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
            # Use the same element-tagged naming scheme as the ITP (fallback to legacy if unresolved)
            aname = generate_atom_name(
                resname,
                local,
                atype,
                type_name_map,
                element_map,
            )
            resname5 = format_gro_name(resname, width=5, keep_tail=False)
            aname5 = format_gro_name(aname, width=5, keep_tail=True)
            f.write(
                f"{resnr % 100000:5d}{resname5}{aname5}"
                f"{aid % 100000:5d}{x:8.4f}{y:8.4f}{z:8.4f}\n"
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

    # Residue assignment (required)
    parser.add_argument(
        "--residue-sizes",
        type=int,
        nargs="+",
        required=True,
        help="Atoms per residue type (e.g., 146 31)"
    )
    parser.add_argument(
        "--residue-names",
        type=str,
        nargs="+",
        required=True,
        help="Residue names for each type (e.g., PVA GLU)"
    )
    parser.add_argument(
        "--repeat",
        "-R",
        action="store_true",
        help="Repeat the residue block pattern until all atoms are assigned"
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
        help="GROMACS dihedral function type: 1=proper periodic, 3=OPLS/RB, 9=multi (default: 1)"
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
    parser.add_argument(
        "--comb-rule",
        type=str,
        choices=["amber", "opls"],
        default="amber",
        help=(
            "Nonbonded preset for GROMACS [defaults]: "
            "amber->comb-rule 2, fudgeLJ 0.5, fudgeQQ 0.8333 (default); "
            "opls->comb-rule 2, fudgeLJ 0.5, fudgeQQ 0.5"
        ),
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

    try:
        return parser.parse_args(argv)
    except SystemExit as e:
        # argparse exits with code 2 for argument errors
        # Map to code 3 for "invalid input" to avoid conflict with file not found (2)
        if e.code == 2:
            sys.exit(3)
        else:
            raise


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
    nonbonded_presets = {
        "amber": (2, 0.5, 0.8333),
        "opls": (2, 0.5, 0.5),
    }
    comb_rule_value, fudge_lj, fudge_qq = nonbonded_presets[args.comb_rule]

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
            defaults=(1, comb_rule_value, "yes", fudge_lj, fudge_qq),
            residue_block_sizes=residue_sizes,
            residue_block_names=residue_names,
            repeat_residue_blocks=args.repeat,
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
        
        sys.exit(0)  # Successful completion

    except FileNotFoundError as e:
        print(f"\nERROR: Input file not found\n  {e}")
        sys.exit(2)
    except ValueError as e:
        print(f"\nERROR: Invalid input\n  {e}")
        sys.exit(3)  # Also used by argparse errors
    except RuntimeError as e:
        print(f"\nERROR: Conversion failed\n  {e}")
        sys.exit(4)
    except Exception as e:
        print(f"\nERROR: Unexpected error\n  {type(e).__name__}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
