from LAMMPS_data_modules.data_file_adjustments import *
from math import ceil
from copy import deepcopy
#Main


N2_shock = read_data("N2_Shock/data/shock_Ti_300_100ps.moldata")
N2_equilib = read_data("N2_Shock/data/large_equilib_900x900x1000.moldata")


N2_shock.rebuild_mol_id()
N2_shock.consecutive_atm_ID()

trim(N2_equilib,0,500,dir='x')
translate(N2_equilib,22000,0,0)


#unwrap_x(N2_shock)
#unwrap_y(N2_shock)
#unwrap_z(N2_shock)


combined_sim = combine(N2_shock,N2_equilib,box_param=N2_shock.box,offset_x=500)
combined_sim.rebuild_mol_id()
combined_sim.consecutive_atm_ID()
print("Dims: ",combined_sim.box)

write_data(combined_sim,'N2_Shock/data/shock_ready_Ti_300_100ps.moldata')