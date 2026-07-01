import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from log_reading import read_lammps_log

ideal_gas_data = read_lammps_log("N2/ideal_law_log.LAMMPS")

mean_values = []

for data in ideal_gas_data:
    mean_values.append({
        "Temp": data["Temp"].mean(),
        "Dens": data["Density"].mean(),
        "Pres": data["Press"].mean()
    })

mean_values = pd.DataFrame(mean_values).drop_duplicates()

M = 14.007*2
k_B = 1.381e-23
N_A = 6.022e+23
mean_values["DensTemp"] = mean_values["Dens"]*mean_values["Temp"]
#corrected to acount for LAMMPS real units
mean_values["IdealPres"] = (mean_values["DensTemp"]*1e+06*k_B*N_A/M)*9.86923e-6

mean_values = np.log10(mean_values)
#mean_values = mean_values.groupby(["Temp", "Dens"], as_index=False).mean()
print (mean_values)


# Create pressure matrices
P_sim = mean_values.pivot_table(index="Temp",
                          columns="Dens",
                          values="Pres",
                          aggfunc="mean").values

P_ideal = mean_values.pivot_table(index="Temp",
                            columns="Dens",
                            values="IdealPres",
                            aggfunc="mean").values


# Create 2D grids
dens = np.array(mean_values["Dens"])
temp = np.array(mean_values["Temp"])
P_sim = np.array(mean_values["Pres"])
P_ideal = np.array(mean_values["IdealPres"])

#DENS, TEMP = np.meshgrid(dens, temp)

#print(dens,temp)
#print(np.array(mean_values["Dens"]),np.array(mean_values["Temp"]))
#print(DENS, TEMP)
#print(DENS, TEMP[:,0])

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

#scatter, plot_wireframe, plot_surface, plot_trisurf

# Simulated pressure
ax.plot_trisurf(
    dens,
    temp,
    P_sim,
    color='b',      # color map
    edgecolor='none',    # remove triangle edges
    linewidth=0.2,       # edge width
    antialiased=True,
    alpha=0.9 
 

)

# Ideal gas pressure
ax.plot_trisurf(
    dens,
    temp,
    P_ideal,
    color='r',      # color map
    edgecolor='k',    # remove triangle edges
    linewidth=0.5,       # edge width
    antialiased=True,
    alpha=0.9
   
)


ax.set_xlabel("Density")
ax.set_ylabel("Temperature")
ax.set_zlabel("Pressure")

plt.show()
