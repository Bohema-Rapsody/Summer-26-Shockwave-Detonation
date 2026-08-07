from LAMMPS_data_modules.data_file_adjustments import *
N = read_data("N2_Shock/data/shock_formed_300.moldata") 
print('read')
write_data(N, "N2_Shock/data/shock_formed_300.moldata")