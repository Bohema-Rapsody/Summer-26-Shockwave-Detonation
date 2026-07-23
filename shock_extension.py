from LAMMPS_data_modules.data_file_adjustments import *

#Main
shock = read_data("N2_shock/data/shock_formed_3.moldata")
particle_equilib = read_data("Cu-N2/data/Cu_N2_combined.moldata")

shock.convert_real_to_metal()
#print(equilib.atoms[0])

#unwrap_x(particle_equilib)

translate(particle_equilib, dx=shock.box.xhi)

#print(equilib.atoms[0])

#Renumbering of atoms in the Cu-N2 sim
atom_num = len(shock.atoms)
mol_num = len(shock.build_molecules())
bond_num = len(shock.bonds)

renumber(particle_equilib,
          atom_offset=atom_num,
          molecule_offset=mol_num,
          bond_offset=bond_num)

#print(equilib.atoms[0])

combined_sim = combine(shock, particle_equilib, box_param=shock.box, offset_x=particle_equilib.box.xhi)

#print(combined_sim.atoms[0])

combined_sim.consecutive_atm_ID()

#print(combined_sim.atoms[0])

write_data(combined_sim, "N2_shock/data/shock_Cu_1.moldata")