from LAMMPS_data_modules.data_file_adjustments import *

#Main
#particle = read_data("Ti-N2/data/Ti_equilib_50_3.moldata")
particle = read_data("Ti-N2/data/Ti_equilib_300.moldata")

N2_shock = read_data("N2_Shock/data/shock_formed_300.moldata")


#N = 3599
#N = 5648
#N = 1197215
#N = 241628
N = 781539
CONV = 1.0364269e-4 
k_B = 8.617e-5
mass = 47.867
#mass = 14.007
#mass = 63.546

KE = 0.0
velocity_dict = particle.build_velocities()
for atom in particle.atoms:
    if atom.atom_type == 1:
        vel = velocity_dict[atom.id]
        KE += 0.5 * mass * (vel.vx**2 + vel.vy**2 + vel.vz**2) * CONV

T = 2*KE/(3*(N-1)*k_B)
print(T)
#exit()
#Conversion modules to ensure compatibility
particle.swap_types({1:2, 2:1})
#N2_shock.convert_real_to_metal()
#unwrap_x(N2_equilib)

#deleting uneeded molecules
#Shock dims - particle at x = 17000 + 2000 = 19000
# -128.0 lo 896.0 hi = 384 mid
#box size = 400 - ensures clean copy (box side length)

centre = (N2_shock.box.ylo+N2_shock.box.yhi)/2
print("centre at:",centre)
particle_loc = 19000
particle_size = 400 #cube side length

cut_box(particle,0.0,0.0,0.0,particle_size,keep_inside=True)
cut_box(N2_shock,particle_loc,centre,centre,particle_size+2,keep_inside=False)  
N2_shock.consecutive_atm_ID()
particle.consecutive_atm_ID()
print("cut box complete")



translate(particle,particle_loc,centre,centre)

#write_data(particle, "Ti-N2/data/Ti_N2_shock_ready.moldata")

print("done translate")

combined_sim = combine(particle,N2_shock,box_param=N2_shock.box)

combined_sim.consecutive_atm_ID()

print("done combining")
'''
KE = 0.0
velocity_dict = combined_sim_Cu_N2.build_velocities()
for atom in combined_sim_Cu_N2.atoms:
    if atom.atom_type == 2:
        vel = velocity_dict[atom.id]
        KE += 0.5 * mass * (vel.vx**2 + vel.vy**2 + vel.vz**2) * CONV

T = 2*KE/(3*(N-1)*k_B)
print(T)

min_dist = 1e9

for ti in combined_sim_Cu_N2.atoms:
    for n in combined_sim_Cu_N2.atoms:
        d = abs((ti.x - n.x)**2 +(ti.y - n.y)**2 +(ti.z - n.z)**2)**0.5
        min_dist = min(min_dist, d)

    print(min_dist)

print(min_dist)
'''

write_data(combined_sim, "Ti-N2/data/Ti_shock_ready_300.moldata")
#write_data(N2_shock, "Ti-N2/data/Ti_N2_shock_ready.moldata")