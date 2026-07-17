from LAMMPS_data_modules.data_file_adjustments import *

#Main
shock = read_data("N2_shock/data/shock_formed_1.data")
equilib = read_data("N2_shock/data/N2_equilib.data")

print(equilib.atoms[0])

unwrap_x(equilib)

translate(equilib, dx=8000)

print(equilib.atoms[0])

renumber(equilib,
          atom_offset=177212,
          molecule_offset=88606,
          bond_offset=88606)

print(equilib.atoms[0])

combined_sim = combine(shock, equilib)

print(combined_sim.atoms[0])

combined_sim.sort_atm_ID()

print(combined_sim.atoms[0])

write_data(combined_sim, "N2_shock/data/shock_combined_1.data")