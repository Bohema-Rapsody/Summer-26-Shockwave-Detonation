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

    system.box.xhi += dx
    system.box.xlo += dx

    system.box.yhi += dy
    system.box.ylo += dy

    system.box.zhi += dz
    system.box.zlo += dz

    print("translated by:",dx,dy,dz)
    


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
    system.consecutive_atm_ID()
    # -------------------------
    # Atoms
    # -------------------------

    for atom in system.atoms:

        atom.id += atom_offset

        if atom.molecule != 0:
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

    print("Offset system by: Atoms:",atom_offset,"Molecules:",molecule_offset,"Bonds:",bond_offset)
    #system.rebuild_mol_id()


def combine(system1: LAMMPSData,
            system2: LAMMPSData,
            box_param, offset_x=0.0,
            offset_y=0.0, offset_z=0.0):

    pre_combine_renum(system1,system2)

    combined = LAMMPSData()

    # Copy metadata from first system
    combined.masses = {**deepcopy(system1.masses)}
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
    combined.box.xhi += offset_x
    combined.box.yhi += offset_y
    combined.box.zhi += offset_z

    combined.consecutive_atm_ID()
    print('Combine complete, box params',combined.box)
    
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
                f"Molecule {molecule.id} does not contain 2 atoms. {molecule.atoms}"
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

def unwrap_y(system: LAMMPSData):
    """
    Reconstruct periodic molecules onto the low-y side
    of the simulation box.

    Parameters
    ----------
    system : LAMMPSData
        Periodic equilibrium gas system.
    """

    Ly = system.box.yhi - system.box.ylo

    molecules = system.build_molecules()
    mol_to_be_del = []

    for molecule in molecules.values():

        if len(molecule.atoms) != 2:
            raise ValueError(
                f"Molecule {molecule.id} does not contain 2 atoms."
            )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]


        if atom1.iy == atom2.iy:
            continue

        # Move both atoms into their true positions first
        #
        # y_true = y + iy*Ly
        #
        for atom in molecule.atoms:

            atom.y += atom.iy * Ly

            # Remove y-image information
            atom.iy = 0


        # Now check if molecule is positioned too far right
        #
        # We want the molecule close to the low-y side.
        #

        while (atom1.y < system.box.ylo) or (atom2.y < system.box.ylo):

            atom1.y += Ly
            atom2.y += Ly

        while (atom1.y >= system.box.yhi) or (atom2.y >= system.box.yhi):

            atom1.y -= Ly
            atom2.y -= Ly

        print ("Corrected molecule ",molecule.id,", Atom 1: ",atom1.id,atom1.y," Atom 2: ",atom2.id,atom2.y)

        if (atom1.y + atom2.y)/2 < system.box.ylo:
            mol_to_be_del.append(molecule.id)
        
        
    system.delete_molecule(mol_to_be_del)


def unwrap_z(system: LAMMPSData):
    """
    Reconstruct periodic molecules onto the low-z side
    of the simulation box.

    Parameters
    ----------
    system : LAMMPSData
        Periodic equilibrium gas system.
    """

    Lz = system.box.zhi - system.box.zlo

    molecules = system.build_molecules()
    mol_to_be_del = []

    for molecule in molecules.values():

        if len(molecule.atoms) != 2:
            raise ValueError(
                f"Molecule {molecule.id} does not contain 2 atoms."
            )

        atom1 = molecule.atoms[0]
        atom2 = molecule.atoms[1]


        if atom1.iz == atom2.iz:
            continue

        # Move both atoms into their true positions first
        #
        # z_true = z + iz*Lz
        #
        for atom in molecule.atoms:

            atom.z += atom.iz * Lz

            # Remove z-image information
            atom.iz = 0


        # Now check if molecule is positioned too far right
        #
        # We want the molecule close to the low-z side.
        #

        while (atom1.z < system.box.zlo) or (atom2.z < system.box.zlo):

            atom1.z += Lz
            atom2.z += Lz

        while (atom1.z >= system.box.zhi) or (atom2.z >= system.box.zhi):

            atom1.z -= Lz
            atom2.z -= Lz

        print ("Corrected molecule ",molecule.id,", Atom 1: ",atom1.id,atom1.z," Atom 2: ",atom2.id,atom2.z)

        if (atom1.z + atom2.z)/2 < system.box.zlo:
            mol_to_be_del.append(molecule.id)
        
        
    system.delete_molecule(mol_to_be_del)

def cut_box(system: LAMMPSData,cx,cy,cz,length,keep_inside):

    molecules = system.build_molecules()
    mol_to_be_del = []
    atom_to_be_del = []
    #print(dict(islice(molecules.items(), 5)))
    #print(next(iter(molecules)))
    for mol in molecules.values():
        #print(mol)
        for atom in mol.atoms:
            if (
                abs(atom.x - cx) < length/2 and
                abs(atom.y - cy) < length/2 and
                abs(atom.z - cz) < length/2
                ):
                if not keep_inside:
                    mol_to_be_del.append(mol.id)

            elif keep_inside:
                    mol_to_be_del.append(mol.id)

    system.delete_molecule(mol_to_be_del)

    for atom in system.atoms:
        if atom.atom_type == 2:
            #atom_to_be_del.append(atom.id)
            
            if (
                abs(atom.x - cx) < (length/2) and
                abs(atom.y - cy) < (length/2) and
                abs(atom.z - cz) < (length/2)
                ):
                if not keep_inside:
                    atom_to_be_del.append(atom.id)

            elif keep_inside:
                 atom_to_be_del.append(atom.id)
            

    system.delete_atoms(atom_to_be_del)

def trim(system:LAMMPSData,clo,chi,dir = 'x',offset=0):
    molecules = system.build_molecules()
    mol_to_be_del = []
    atom_to_be_del = []

    for mol in molecules.values():


        for atom in mol.atoms:

            if dir =='x':
                check_dim = atom.x

            if dir =='y':
                check_dim = atom.y

            if dir =='z':
                check_dim = atom.z

            if clo+offset < check_dim < chi-offset:
                continue
            else:
                mol_to_be_del.append(mol.id)

    system.delete_molecule(mol_to_be_del)

    for atom in system.atoms:
        if atom.atom_type == 2:
                        
            if dir =='x':
                check_dim = atom.x

            if dir =='y':
                check_dim = atom.y

            if dir =='z':
                check_dim = atom.z
            
            if clo+offset < check_dim < chi-offset:
                continue
            else:
                atom_to_be_del.append(atom.id)
            

    system.delete_atoms(atom_to_be_del)

    if dir =='x':
        system.box.xlo = clo
        system.box.xhi = chi

    if dir =='y':
        system.box.ylo = clo
        system.box.yhi = chi

    if dir =='z':
        system.box.zlo = clo
        system.box.zhi = chi
    
def pre_combine_renum(base_system:LAMMPSData,append_system:LAMMPSData):
    base_system.consecutive_atm_ID()

    atom_num = len(base_system.atoms)
    mol_num = len(base_system.build_molecules())
    bond_num = len(base_system.bonds)

    renumber(append_system,
            atom_offset=atom_num,
            molecule_offset=mol_num,
            bond_offset=bond_num)

def insert_section(base_system:LAMMPSData,insert_system:LAMMPSData,insert):
    #insert shows the position at which the insert_system will be inserted into base-system
    base_low = deepcopy(base_system)
    print('copied')
    #base_high = deepcopy(base_system)
    insert_width = insert_system.box.xhi-insert_system.box.xlo
    translate(insert_system,insert-insert_system.box.xlo,0,0)

    trim(base_low,base_low.box.xlo,insert,dir='x')
    trim(base_system,insert,base_system.box.xhi,dir='x')
    translate(base_system,insert_width,0,0)

    base_low = combine(base_low,insert_system,box_param=base_low.box,offset_x=insert_width)
    

    return(combine(base_low,base_system,box_param=base_low.box,offset_x=(base_system.box.xhi-base_system.box.xlo)))