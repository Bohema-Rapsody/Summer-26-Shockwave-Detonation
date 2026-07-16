"""
read_data.py

Reads a LAMMPS write_data file (atom_style full)
into a LAMMPSData object.
"""

from LAMMPS_data_classes import Atom, Velocity, Bond, Box, LAMMPSData


# --------------------------------------------------
# Helper Functions
# --------------------------------------------------

def find_section(lines, section_name):
    """
    Finds the line number of a section heading.

    Returns None if the section doesn't exist.
    """

    for i, line in enumerate(lines):

        if line.strip().startswith(section_name):
            return i

    return None

def find_next_section(lines, start):

    section_names = [
        "Masses",
        "Pair Coeffs",
        "Bond Coeffs",
        "Atoms",
        "Velocities",
        "Bonds"
    ]

    for i in range(start + 1, len(lines)):

        text = lines[i].strip()

        for section in section_names:
            if text.startswith(section):
                return i

    return len(lines)


# --------------------------------------------------
# Main reader
# --------------------------------------------------

def read_data(filename):

    system = LAMMPSData()

    with open(filename, "r") as file:
        lines = file.readlines()

    # --------------------------------------------------
    # Header
    # --------------------------------------------------

    for line in lines:

        words = line.split()

        if len(words) == 4:

            if words[2] == "xlo":
                system.box = Box(
                    float(words[0]),
                    float(words[1]),
                    0,
                    0,
                    0,
                    0
                )

            elif words[2] == "ylo":

                system.box.ylo = float(words[0])
                system.box.yhi = float(words[1])

            elif words[2] == "zlo":

                system.box.zlo = float(words[0])
                system.box.zhi = float(words[1])

    # --------------------------------------------------
    # Masses
    # --------------------------------------------------

    start = find_section(lines, "Masses")

    if start is not None:

        end = find_next_section(lines, start)

        for line in lines[start + 2:end]:

            words = line.split()

            # Skip blank lines
            if not words:
                continue

            # A valid mass line has exactly two entries:
            # atom_type  mass
            if len(words) != 2:
                continue

            system.masses[int(words[0])] = float(words[1])

    # --------------------------------------------------
    # Pair Coeffs
    # --------------------------------------------------

    start = find_section(lines, "Pair Coeffs")

    if start is not None:

        end = find_next_section(lines, start)

        for line in lines[start + 2:end]:

            if line.strip():

                system.pair_coeffs.append(line.rstrip())

    # --------------------------------------------------
    # Bond Coeffs
    # --------------------------------------------------

    start = find_section(lines, "Bond Coeffs")

    if start is not None:

        end = find_next_section(lines, start)

        for line in lines[start + 2:end]:

            if line.strip():

                system.bond_coeffs.append(line.rstrip())

    # --------------------------------------------------
    # Atoms
    # --------------------------------------------------

    start = find_section(lines, "Atoms")

    end = find_next_section(lines, start)

    for line in lines[start + 2:end]:

        if line.strip() == "":
            continue

        words = line.split()

        atom = Atom(

            id=int(words[0]),
            molecule=int(words[1]),
            atom_type=int(words[2]),

            charge=float(words[3]),

            x=float(words[4]),
            y=float(words[5]),
            z=float(words[6]),

            ix=int(words[7]),
            iy=int(words[8]),
            iz=int(words[9])
        )

        system.atoms.append(atom)

    # --------------------------------------------------
    # Velocities
    # --------------------------------------------------

    start = find_section(lines, "Velocities")

    end = find_next_section(lines, start)

    for line in lines[start + 2:end]:

        if line.strip() == "":
            continue

        words = line.split()

        velocity = Velocity(

            atom_id=int(words[0]),

            vx=float(words[1]),
            vy=float(words[2]),
            vz=float(words[3])
        )

        system.velocities.append(velocity)

    # --------------------------------------------------
    # Bonds
    # --------------------------------------------------

    start = find_section(lines, "Bonds")

    end = find_next_section(lines, start)

    for line in lines[start + 2:end]:

        if line.strip() == "":
            continue

        words = line.split()

        bond = Bond(

            id=int(words[0]),

            bond_type=int(words[1]),

            atom1=int(words[2]),
            atom2=int(words[3])
        )

        system.bonds.append(bond)

    return system