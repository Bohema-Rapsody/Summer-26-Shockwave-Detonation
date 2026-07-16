"""
transform.py

Functions for modifying the coordinates of a LAMMPSData object.
"""

from LAMMPS_data_classes import LAMMPSData
from read_data_file import read_data


def translate(system: LAMMPSData,
              dx: float = 0.0,
              dy: float = 0.0,
              dz: float = 0.0):
    """
    Translate every atom in the system.

    Parameters
    ----------
    system : LAMMPSData
        The system to translate.

    dx, dy, dz : float
        Distance (Å) to translate the system.
    """

    for atom in system.atoms:

        atom.x += dx
        atom.y += dy
        atom.z += dz


def renumber(system: LAMMPSData,
             atom_offset: int = 0,
             molecule_offset: int = 0,
             bond_offset: int = 0):
    """
    Renumber all IDs in the system.

    Parameters
    ----------
    atom_offset : int
        Added to every atom ID.

    molecule_offset : int
        Added to every molecule ID.

    bond_offset : int
        Added to every bond ID.
    """

    # -------------------------
    # Atoms
    # -------------------------

    for atom in system.atoms:

        atom.id += atom_offset
        atom.molecule += molecule_offset

    # -------------------------
    # Velocities
    # -------------------------

    for velocity in system.velocities:

        velocity.atom_id += atom_offset

    # -------------------------
    # Bonds
    # -------------------------

    for bond in system.bonds:

        bond.id += bond_offset

        bond.atom1 += atom_offset
        bond.atom2 += atom_offset


#Main
shock = read_data("N2_shock/data/shock_formed_1.data")
equilib = read_data("N2_shock/data/N2_equilib.data")

print(equilib.atoms[0])

translate(equilib, dx=8000)

print(equilib.atoms[0])

renumber(equilib,
          atom_offset=177212,
          molecule_offset=88606,
          bond_offset=88606)

print(equilib.atoms[0])