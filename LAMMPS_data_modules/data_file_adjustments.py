"""
transform.py

Functions for modifying the coordinates of a LAMMPSData object.
"""

from LAMMPS_data_modules.LAMMPS_data_classes import LAMMPSData
from LAMMPS_data_modules.read_data_file import read_data
from LAMMPS_data_modules.write_data_file import write_data
from copy import deepcopy
from itertools import islice


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

        if atom.molecule == 0:
            continue
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


def combine(system1: LAMMPSData,
            system2: LAMMPSData,
            box_param):

    combined = LAMMPSData()

    # Copy metadata from first system
    combined.masses = deepcopy(system1.masses)
    combined.pair_coeffs = deepcopy(system1.pair_coeffs)
    combined.bond_coeffs = deepcopy(system1.bond_coeffs)

    # Combine atoms
    combined.atoms = system1.atoms + system2.atoms

    # Combine velocities
    combined.velocities = system1.velocities + system2.velocities

    # Combine bonds
    combined.bonds = system1.bonds + system2.bonds

    # Copy box
    
    #combined.box = deepcopy(system1.box)

    combined.box = box_param

    return combined


def unwrap_x(system: LAMMPSData):
    """
    Reconstruct periodic molecules onto the low-x side
    of the simulation box.

    Parameters
    ----------
    system : LAMMPSData
        Periodic equilibrium gas system.
    """

    Lx = system.box.xhi - system.box.xlo

    molecules = system.build_molecules()
    mol_to_be_del = []

    for molecule in molecules.values():

        if len(molecule.atoms) != 2:
            raise ValueError(
                f"Molecule {molecule.id} does not contain 2 atoms."
            )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]


        if atom1.ix == atom2.ix:
            continue

        # Move both atoms into their true positions first
        #
        # x_true = x + ix*Lx
        #
        for atom in molecule.atoms:

            atom.x += atom.ix * Lx

            # Remove x-image information
            atom.ix = 0


        # Now check if molecule is positioned too far right
        #
        # We want the molecule close to the low-x side.
        #

        while (atom1.x < system.box.xlo) or (atom2.x < system.box.xlo):

            atom1.x += Lx
            atom2.x += Lx

        while (atom1.x >= system.box.xhi) or (atom2.x >= system.box.xhi):

            atom1.x -= Lx
            atom2.x -= Lx

        print ("Corrected molecule ",molecule.id,", Atom 1: ",atom1.id,atom1.x," Atom 2: ",atom2.id,atom2.x)

        if (atom1.x + atom2.x)/2 < system.box.xlo:
            mol_to_be_del.append(molecule.id)
        
        
    system.delete_molecule(mol_to_be_del)


def cut_box(system: LAMMPSData,cx,cy,cz,length,keep_inside=True):

    molecules = system.build_molecules()
    mol_to_be_del = []
    #print(dict(islice(molecules.items(), 5)))
    #print(next(iter(molecules)))
    for mol in molecules.values():
        #print(mol)
        if (
            abs(mol.x_com - cx) > length/2 and
            abs(mol.y_com - cy) > length/2 and
            abs(mol.z_com - cz) > length/2
            ):
            if keep_inside:
                mol_to_be_del.append(mol.id)

        else:
            if not keep_inside:
                mol_to_be_del.append(mol.id)

    system.delete_molecule(mol_to_be_del)

