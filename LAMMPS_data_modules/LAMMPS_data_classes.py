"""
objects.py

Defines the Python classes used to store data from a
LAMMPS write_data file (atom_style full).
"""

from dataclasses import dataclass, field
from collections import defaultdict
from copy import deepcopy


# -------------------------------------------------
# Atom
# -------------------------------------------------

@dataclass
class Atom:
    """
    Stores one atom from the Atoms section.

    Format (atom_style full):

    id mol type charge x y z ix iy iz
    """

    id: int
    molecule: int
    atom_type: int

    charge: float

    x: float
    y: float
    z: float

    ix: int
    iy: int
    iz: int


# -------------------------------------------------
# Velocity
# -------------------------------------------------

@dataclass
class Velocity:
    """
    Stores one velocity entry.

    Format:

    atomID vx vy vz
    """

    atom_id: int

    vx: float
    vy: float
    vz: float


# -------------------------------------------------
# Bond
# -------------------------------------------------

@dataclass
class Bond:
    """
    Stores one bond.

    Format:

    bondID bondType atom1 atom2
    """

    id: int

    bond_type: int

    atom1: int
    atom2: int

# -------------------------------------------------
# Molecule
# -------------------------------------------------

@dataclass
class Molecule:
    """
    Stores one molecule.

    At present this simply contains the atoms belonging to the molecule,
    but additional properties (centre of mass, etc.) can easily be added.
    """

    id: int

    atoms: list[Atom] = field(default_factory=list)

    bond: Bond | None = None

    # Calculated properties
    x_com: float = 0.0
    y_com: float = 0.0
    z_com: float = 0.0

    def calculate_com(self, masses):
        """
        Calculate the centre of mass.

        Parameters
        ----------
        masses : dict
            Dictionary of atom_type -> mass
        """

        total_mass = 0.0

        x = 0.0
        y = 0.0
        z = 0.0

        for atom in self.atoms:

            m = masses[atom.atom_type]

            total_mass += m

            x += m * atom.x
            y += m * atom.y
            z += m * atom.z

        self.x_com = x / total_mass
        self.y_com = y / total_mass
        self.z_com = z / total_mass

# -------------------------------------------------
# Simulation Box
# -------------------------------------------------

@dataclass
class Box:

    xlo: float
    xhi: float

    ylo: float
    yhi: float

    zlo: float
    zhi: float


# -------------------------------------------------
# Complete LAMMPS system
# -------------------------------------------------

@dataclass
class LAMMPSData:

    atoms: list[Atom] = field(default_factory=list)

    velocities: list[Velocity] = field(default_factory=list)

    bonds: list[Bond] = field(default_factory=list)

    masses: dict = field(default_factory=dict)

    pair_coeffs: list[str] = field(default_factory=list)

    bond_coeffs: list[str] = field(default_factory=list)

    box: Box | None = None

    def sort_atm_ID(self):
        """Sort atoms, velocities and bonds by their IDs."""

        self.atoms.sort(key=lambda atom: atom.id)

        self.velocities.sort(key=lambda velocity: velocity.atom_id)

        self.bonds.sort(key=lambda bond: bond.id)


    def build_molecules(self):
        """
        Build a dictionary of Molecule objects.

        Returns
        -------
        dict[int, Molecule]
        """

        molecules = {}

        for atom in self.atoms:

            mol_id = atom.molecule

            if mol_id == 0:
                continue

            if mol_id not in molecules:
                molecules[mol_id] = Molecule(id=mol_id)

            molecules[mol_id].atoms.append(atom)

        for molecule in molecules.values():
            molecule.calculate_com(self.masses)

        #print(next(iter(molecules)))
        
        return molecules
    
    def delete_molecule(self, molecule_ids):
    
        """
        Remove an entire molecule from the system.

        Deletes:
        - atoms
        - velocities
        - bonds

        Parameters
        ----------
        molecule_id : int
            Molecule ID to remove
        """


        molecule_ids = set(molecule_ids)

        atom_ids = {
            atom.id
            for atom in self.atoms
            if atom.molecule in molecule_ids
        }

        self.atoms = [
            atom for atom in self.atoms
            if atom.id not in atom_ids
        ]

        self.velocities = [
            vel for vel in self.velocities
            if vel.atom_id not in atom_ids
        ]

        self.bonds = [
            bond for bond in self.bonds
            if bond.atom1 not in atom_ids
            and bond.atom2 not in atom_ids
        ]

        print("Deleted molecules: ",molecule_ids)


    def convert_real_to_metal(self):
        """
        Convert a LAMMPS data structure from REAL units to METAL units.

        Only quantities stored in the data file are converted.

        REAL:
            velocity = Å/fs

        METAL:
            velocity = Å/ps

        Therefore:
            v_metal = 1000 * v_real
        """

        print("Converting REAL units -> METAL units")

        # -------------------------
        # Velocities
        # -------------------------

        for vel in self.velocities:
            vel.vx *= 1000.0
            vel.vy *= 1000.0
            vel.vz *= 1000.0

        # -------------------------
        # Store unit system
        # -------------------------

        self.units = "metal"

        print("Done.")

    def swap_types(self, type_map):
        """
        Swap atom types according to a mapping.

        Example:
            {1:2, 2:1}
        """

        # Swap atom types
        for atom in self.atoms:
            if atom.atom_type in type_map:
                atom.atom_type = type_map[atom.atom_type]

        # Swap masses
        old_masses = self.masses.copy()

        for old_type, new_type in type_map.items():
            self.masses[new_type] = old_masses[old_type]