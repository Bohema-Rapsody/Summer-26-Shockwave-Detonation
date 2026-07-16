"""
objects.py

Defines the Python classes used to store data from a
LAMMPS write_data file (atom_style full).
"""

from dataclasses import dataclass, field


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