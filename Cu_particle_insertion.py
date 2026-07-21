from LAMMPS_data_modules.data_file_adjustments import *

#Main
particle = read_data("Cu-N2/data/Cu_equilib.data")
N2_equilib = read_data("N2_shock/data/N2_equilib.data")


#Conversion modules to ensure compatibility
particle.swap_types({1:2, 2:1})
N2_equilib.convert_real_to_metal()
#unwrap_x(N2_equilib)

#deleting uneeded molecules
cut_box(particle,0.0,0.0,0.0,80,keep_inside=True)
cut_box(N2_equilib,1500.0,0.0,0.0,82,keep_inside=False)

translate(particle,1500.0,0,0)

#Renumbering of atoms in the Cu-N2 sim
atom_num = len(N2_equilib.atoms)
mol_num = len(N2_equilib.build_molecules())
bond_num = len(N2_equilib.bonds)

renumber(particle,
          atom_offset=atom_num,
          molecule_offset=mol_num,
          bond_offset=bond_num)

combined_sim_Cu_N2 = combine(particle,N2_equilib,box_param=N2_equilib.box)

combined_sim_Cu_N2.consecutive_atm_ID()

write_data(combined_sim_Cu_N2, "Cu-N2/data/Cu_N2_combined.moldata")