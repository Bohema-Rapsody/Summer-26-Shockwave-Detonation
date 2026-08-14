from LAMMPS_data_modules.data_file_adjustments import *
from math import ceil
from copy import deepcopy
#Main


N2_shock = read_data("N2_Shock/data/shock_ready_Ti_300_200ps.moldata")

translate(N2_shock,-15631.0,0,0)

write_data(N2_shock,'N2_Shock/data/shock_ready_Ti_300_200ps.moldata')