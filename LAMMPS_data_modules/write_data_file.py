"""
write_data.py

Writes a LAMMPSData object to a write_data file.
"""

from LAMMPS_data_modules.LAMMPS_data_classes import LAMMPSData


def write_data(system: LAMMPSData,
               filename):
        
    #if xhi is None:
    #    xhi = system.box.xhi

    with open(filename, "w") as f:

        #Header
        f.write("LAMMPS data file\n\n")

        f.write(f"{len(system.atoms)} atoms\n")
        f.write("1 atom types\n\n")

        f.write(f"{len(system.bonds)} bonds\n")
        f.write("1 bond types\n\n")

        #f.write(f"{system.box.xlo} {xhi} xlo xhi\n")
        f.write(f"{system.box.xlo} {system.box.xhi} xlo xhi\n")
        f.write(f"{system.box.ylo} {system.box.yhi} ylo yhi\n")
        f.write(f"{system.box.zlo} {system.box.zhi} zlo zhi\n\n")

        #Masses
        f.write("Masses\n\n")

        for atom_type, mass in system.masses.items():
            f.write(f"{atom_type} {mass}\n")

        f.write("\n")

        #Pair Coeffs
        f.write("Pair Coeffs\n\n")

        for line in system.pair_coeffs:
            f.write(line + "\n")

        f.write("\n")

        #Bond Coeffs
        f.write("Bond Coeffs\n\n")

        for line in system.bond_coeffs:
            f.write(line + "\n")

        f.write("\n")

        #Atoms
        f.write("Atoms # full\n\n")

        for atom in system.atoms:

            f.write(
                f"{atom.id} "
                f"{atom.molecule} "
                f"{atom.atom_type} "
                f"{atom.charge} "
                f"{atom.x:.8f} "
                f"{atom.y:.8f} "
                f"{atom.z:.8f} "
                f"{atom.ix} "
                f"{atom.iy} "
                f"{atom.iz}\n"
            )

        f.write("\n")

        #Velocities
        f.write("Velocities\n\n")

        for vel in system.velocities:

            f.write(
                f"{vel.atom_id} "
                f"{vel.vx:.12f} "
                f"{vel.vy:.12f} "
                f"{vel.vz:.12f}\n"
            )

        f.write("\n")

        #Bonds
        f.write("Bonds\n\n")

        for bond in system.bonds:

            f.write(
                f"{bond.id} "
                f"{bond.bond_type} "
                f"{bond.atom1} "
                f"{bond.atom2}\n"
            )

