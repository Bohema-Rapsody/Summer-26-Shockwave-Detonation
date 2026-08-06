from LAMMPS_data_modules.data_file_adjustments import *
from math import ceil
from copy import deepcopy
#Main

piston_pos = 12789
post_shock_start = 13700
post_shock_end = 15700
width = post_shock_end-post_shock_start
particle_pos = 17500

N2_shock = read_data("N2_Shock/data/shock_Cu_2.moldata")
N2_equilib = read_data("N2_Shock/data/N2_equilib.moldata")
N2_shock.rebuild_mol_id()
N2_equilib.rebuild_mol_id()
N2_equilib.consecutive_atm_ID()
N2_shock.consecutive_atm_ID()

#deleting uneeded molecules (Cu particle)
cut_box(N2_shock,particle_pos,0.0,0.0,80,keep_inside=False)


unwrap_x(N2_equilib)
unwrap_y(N2_equilib)
unwrap_z(N2_equilib)
unwrap_x(N2_shock)
unwrap_y(N2_shock)
unwrap_z(N2_shock)



post_shock_req = 5500
post_shock_req = ceil(post_shock_req/width) + 1
print('Copy post-shock:',post_shock_req)

pre_shock_req = 0
pre_shock_req =  ceil(pre_shock_req/N2_equilib.box.xhi)
print('Copy pre-shock:',pre_shock_req)

dim_req = 900
dim_req = ceil(dim_req/(N2_shock.box.xhi-N2_shock.box.xlo))
print('Copy dims:',dim_req)




#Separating into shocks
Shock_low = deepcopy(N2_shock)
trim(Shock_low,N2_shock.box.xlo,post_shock_start)

Post_shock_section = deepcopy(N2_shock)
trim(Post_shock_section,post_shock_start+2,post_shock_end-2)

Shock_high = deepcopy(N2_shock)
trim(Shock_high,post_shock_end,N2_shock.box.xhi)


#Add post_shock section recursively
for i in range(0,post_shock_req):
    append_shock = deepcopy(Post_shock_section)
    translate(append_shock,i*width)
    Shock_low = combine(Shock_low,Post_shock_section,box_param=Shock_low.box,offset_x=width)
    print('Current box max x:',Shock_low.box.xhi)


#combine for the full length
translate(Shock_high,(post_shock_req-1)*width)
combined_sim = combine(Shock_low,Shock_high,box_param=Shock_low.box,offset_x=(Shock_high.box.xhi-Shock_high.box.xlo))
combined_sim.consecutive_atm_ID()

print('Final box max x',combined_sim.box.xhi)

write_data(combined_sim, "N2_Shock/data/shock_formed_300.moldata")
