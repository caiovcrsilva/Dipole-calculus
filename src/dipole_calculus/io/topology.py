"""Topology and ITP charge parsing utilities for GROMACS files."""

import os
import numpy as np


def extract_charges_from_file(file_path, charges_dict, base_dir):
    """Recursively parse a .top or .itp file to extract atomic charges per molecule type."""
    reading_molecule = False
    reading_atoms = False
    current_molecule_name = ""
    temp_charges = []

    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                # Recursively process #include directives
                if line.startswith("#include"):
                    include_path = line.split('"')[1]
                    full_path = os.path.join(base_dir, include_path)
                    extract_charges_from_file(full_path, charges_dict, base_dir)
                    continue

                # Handle section headers starting with '['
                if line.startswith("["):
                    if reading_atoms:
                        reading_atoms = False
                        if current_molecule_name and temp_charges:
                            charges_dict[current_molecule_name] = temp_charges

                    if "moleculetype" in line:
                        reading_molecule = True
                        continue

                    if "atoms" in line:
                        reading_atoms = True
                        temp_charges = []
                        continue

                    continue

                # Read molecule name under [ moleculetype ]
                if reading_molecule:
                    current_molecule_name = line.split()[0]
                    reading_molecule = False
                    continue

                # Read atomic charge from column 7 under [ atoms ]
                if reading_atoms:
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            temp_charges.append(float(parts[6]))
                        except ValueError:
                            pass

            # Save the last molecule parsed in the file
            if reading_atoms and current_molecule_name and temp_charges:
                charges_dict[current_molecule_name] = temp_charges

    except FileNotFoundError:
        print(f"\t  WARNING: File not found -> {file_path}", flush=True)

    return charges_dict


def build_charges_from_top(top_path):
    """Orchestrate .top parsing and expand charges according to the [ molecules ] directive."""
    print(f"--- Reading Topology: {top_path} ---", flush=True)
    charges_dict = {}
    base_dir = os.path.dirname(top_path) or "."

    charges_dict = extract_charges_from_file(top_path, charges_dict, base_dir)

    total_charges = []
    reading_molecules = False

    # Read [ molecules ] section at the end of the .top file
    with open(top_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue

            if line.startswith("[") and "molecules" in line:
                reading_molecules = True
                continue

            if reading_molecules and not line.startswith("["):
                parts = line.split()
                if len(parts) >= 2:
                    mol_name = parts[0]
                    count = int(parts[1])

                    if mol_name in charges_dict:
                        print(f"\t  -> {mol_name}: {len(charges_dict[mol_name])} atoms ({count}x)", flush=True)
                        total_charges.extend(charges_dict[mol_name] * count)
                    else:
                        print(f"\t  ERROR: Topology of molecule '{mol_name}' not found!", flush=True)

    return np.array(total_charges)


def parse_charges_from_itp(itp_path):
    """Extract atomic charges directly from a standalone .itp file."""
    charges = []
    in_atoms = False
    try:
        with open(itp_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue
                if line.startswith("[") and "atoms" in line:
                    in_atoms = True
                    continue
                if in_atoms and line.startswith("["):
                    in_atoms = False
                    break
                if in_atoms:
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            charges.append(float(parts[6]))
                        except ValueError:
                            pass
    except OSError:
        return []
    return charges


def build_topology_charges(config_list):
    """Build full charge array from a list of (itp_path, molecule_count) tuples."""
    full_charges = []
    print("--- Reading charges ---", flush=True)
    for itp, count in config_list:
        q = parse_charges_from_itp(itp)
        print(f"  -> {itp}: {len(q)} atoms. ({count}x)", flush=True)
        for _ in range(count):
            full_charges.extend(q)
    return np.array(full_charges)
