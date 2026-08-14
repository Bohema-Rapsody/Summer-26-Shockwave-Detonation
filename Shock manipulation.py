from LAMMPS_data_modules.data_file_adjustments import *
from math import ceil
from copy import deepcopy
#Main


N2_shock = read_data("N2_Shock/data/shock_Ti_300_200ps.moldata")
N2_equilib = read_data("N2_Shock/data/large_equilib_900x900x1000.moldata")
N2_post_shock = read_data("N2_Shock/data/large_post-shock_900x900x2000.moldata")

post_shock_width = N2_post_shock.box.xhi - N2_post_shock.box.xlo
equilib_width = N2_equilib.box.xhi - N2_equilib.box.xlo

N2_shock.rebuild_mol_id()
N2_shock.consecutive_atm_ID()


post_shock_req = 1000
post_shock_insert = 18000
post_shock_req_i = ceil(post_shock_req/post_shock_width)
print('Copy post-shock:',post_shock_req_i)

pre_shock_req = 3400
pre_shock_req_i =  ceil(pre_shock_req/equilib_width)
print('Copy pre-shock:',pre_shock_req_i)

#dim_req = 306
#dim_req_i = ceil(dim_req/(N2_shock.box.yhi-N2_shock.box.ylo))
#print('Copy dims:',dim_req_i)

#manufacture post-shock
N2_post_shock_copy = deepcopy(N2_post_shock)
for i in range(1,post_shock_req_i):

    append_shock = deepcopy(N2_post_shock)
    translate(append_shock,i*post_shock_width,0,0)

    N2_post_shock_copy=combine(N2_post_shock_copy,N2_post_shock,box_param=N2_post_shock.box,offset_x=post_shock_width)

trim(N2_post_shock_copy,N2_post_shock_copy.box.xlo,N2_post_shock_copy.box.xlo+post_shock_req,dir='x',offset=2)
#translate(N2_post_shock,post_shock_insert,0,0)
print(N2_post_shock_copy.box)

N2_shock_copy = insert_section(N2_shock,N2_post_shock_copy,post_shock_insert)
print(N2_shock_copy.box)
del N2_post_shock,N2_post_shock_copy,N2_shock

#manufacture pre-shock
N2_pre_shock = deepcopy(N2_equilib)
for i in range(1,pre_shock_req_i):

    append_shock = deepcopy(N2_equilib)
    translate(append_shock,i*equilib_width,0,0)

    N2_pre_shock=combine(N2_pre_shock,append_shock,box_param=N2_equilib.box,offset_x=equilib_width)

trim(N2_pre_shock,N2_pre_shock.box.xlo,N2_pre_shock.box.xlo+pre_shock_req,dir='x',offset=2)
translate(N2_pre_shock,N2_shock_copy.box.xhi,0,0)
print(N2_pre_shock.box)

final_shock = combine(N2_shock_copy,N2_pre_shock,box_param=N2_shock_copy.bsox,offset_x=pre_shock_req)

#unwrap_x(N2_shock)
#unwrap_y(N2_shock)
#unwrap_z(N2_shock)

print("Dims: ",final_shock.box)

translate(final_shock,-15631.0,0,0)

write_data(final_shock,'N2_Shock/data/shock_ready_Ti_300_200ps.moldata')