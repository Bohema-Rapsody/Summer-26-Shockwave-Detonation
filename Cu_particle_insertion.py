from LAMMPS_data_modules.data_file_adjustments import *

#Main
particle = read_data("Ti-N2/data/Ti_equilib_50.moldata")
N2_shock = read_data("N2_shock/data/shock_Cu_2.moldata")


#Conversion modules to ensure compatibility
particle.swap_types({1:2, 2:1})
N2_shock.convert_real_to_metal()
#unwrap_x(N2_equilib)

#deleting uneeded molecules
cut_box(particle,0.0,0.0,0.0,80,keep_inside=True)
cut_box(N2_shock,17500.0,0.0,0.0,82,keep_inside=False)

translate(particle,17500.0,0,0)

#Renumbering of atoms in the Cu-N2 sim
atom_num = len(N2_shock.atoms)
mol_num = len(N2_shock.build_molecules())
bond_num = len(N2_shock.bonds)

renumber(particle,
          atom_offset=atom_num,
          molecule_offset=mol_num,
          bond_offset=bond_num)

combined_sim_Cu_N2 = combine(particle,N2_shock,box_param=N2_shock.box)

combined_sim_Cu_N2.consecutive_atm_ID()

write_data(combined_sim_Cu_N2, "Ti-N2/data/Ti_N2_shock_ready.moldata")