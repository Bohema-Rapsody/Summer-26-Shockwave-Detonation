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
#trim end - minimise pre-shock length to 100ps
#trim(N2_shock,N2_shock.box.xlo,20000,dir='x')


unwrap_x(N2_equilib)
unwrap_y(N2_equilib)
unwrap_z(N2_equilib)
unwrap_x(N2_shock)
unwrap_y(N2_shock)
unwrap_z(N2_shock)



post_shock_req = 1500
post_shock_req = ceil(post_shock_req/width) + 1
print('Copy post-shock:',post_shock_req)

pre_shock_req = 0
pre_shock_req =  ceil(pre_shock_req/N2_equilib.box.xhi)
print('Copy pre-shock:',pre_shock_req)

dim_req = 306
dim_req_i = ceil(dim_req/(N2_shock.box.yhi-N2_shock.box.ylo))
print('Copy dims:',dim_req_i)




#Separating into shocks
Shock_low = deepcopy(N2_shock)
trim(Shock_low,N2_shock.box.xlo,post_shock_start,dir='x')
Shock_low.consecutive_atm_ID()

Post_shock_section = deepcopy(N2_shock)
trim(Post_shock_section,post_shock_start,post_shock_end,dir='x',offset=2)
Post_shock_section.consecutive_atm_ID()

Shock_high = deepcopy(N2_shock)
trim(Shock_high,post_shock_end,N2_shock.box.xhi,dir='x')
Shock_high.consecutive_atm_ID()


#Add post_shock section recursively
for i in range(0,post_shock_req):
    append_shock = deepcopy(Post_shock_section)
    translate(append_shock,i*width)

    Shock_low = combine(Shock_low,append_shock,box_param=Shock_low.box,offset_x=width)
    Shock_low.consecutive_atm_ID()
    print('Current box max x:',Shock_low.box.xhi)


#combine for the full length
translate(Shock_high,(post_shock_req-1)*width)
combined_sim = combine(Shock_low,Shock_high,box_param=Shock_low.box,offset_x=(Shock_high.box.xhi-Shock_high.box.xlo))
combined_sim.consecutive_atm_ID()

print('Final box max x',combined_sim.box.xhi)

#write_data(combined_sim, "N2_Shock/data/shock_formed_300.moldata")

y_stack = deepcopy(combined_sim)
width = y_stack.box.yhi - y_stack.box.ylo
for i in range(1,dim_req_i):
    append_shock = deepcopy(combined_sim)
    translate(append_shock,dy=i*width)

    y_stack = combine(y_stack,append_shock,box_param=y_stack.box,offset_y=width)
    y_stack.consecutive_atm_ID()
    print('Current box max y:',y_stack.box.yhi)


z_stack = deepcopy(y_stack)
width = z_stack.box.zhi - z_stack.box.zlo
for i in range(1,dim_req_i):
    append_shock = deepcopy(y_stack)
    translate(append_shock,dz=i*width)

    z_stack = combine(z_stack,append_shock,box_param=z_stack.box,offset_z=width)
    z_stack.consecutive_atm_ID()
    print('Current box max z:',z_stack.box.zhi)

z_stack.consecutive_atm_ID()

#Trim to size

cut_size = (z_stack.box.yhi - z_stack.box.ylo - dim_req)/2
trim(z_stack,z_stack.box.ylo+cut_size,z_stack.box.yhi-cut_size, dir='y',offset=2)
trim(z_stack,z_stack.box.zlo+cut_size,z_stack.box.zhi-cut_size,dir='z',offset=2)
print('complete trim to size:',dim_req,'x',dim_req)

write_data(z_stack, "N2_Shock/data/shock_formed_100.moldata")

#write_data(N, "N2_Shock/data/shock_formed_300.moldata")