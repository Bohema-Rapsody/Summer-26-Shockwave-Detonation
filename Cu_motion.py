import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

#profiles = {}

with open("N2_shock/profiles/Cu_particle.profile") as f:
    data = []

    while True:

        line = f.readline()

        if not line:
            break

        if line.startswith("#"):
            continue

        words = line.split()

        if len(words) == 11:

        
            data.append(
                list(map(float, f.readline().split()))
            )


    data = pd.DataFrame(
        data,
        columns=[
            "timestep",
            "comx",
            "comy",
            "comz",
            "vx",
            "vy",
            "vz",
            "temp",
            "ke",
            "pe",
            "r_gyr"
        ],
    )


plt.rcParams["figure.figsize"] = (17, 8)
fig, axs = plt.subplots(2,2)
 


axs[0,0].plot(
    data["timestep"],
    data["comx"]
)

axs[0,0].set_xlim(0, 400000)
axs[0,0].set_ylim(17400, 19000)
axs[0,0].set_ylabel("CoM position")
#axs[0,0].set_xlabel("Timestep (fs)")
axs[0,0].grid(True)


axs[1,0].plot(
    data["timestep"],
    data["vx"]
)

axs[1,0].set_xlim(0, 400000)
axs[1,0].set_ylim(-1, 10)
axs[1,0].set_ylabel("CoM velocity")
axs[1,0].set_xlabel("Timestep (fs)")
axs[1,0].grid(True)


axs[0,1].plot(
    data["timestep"],
    data["temp"]
)

axs[0,1].set_xlim(0, 400000)
axs[0,1].set_ylim(0, 1300)
axs[0,1].set_ylabel("Temperature (CoM corrected K)")
#axs[0,1].set_xlabel("Timestep (fs)")
axs[0,1].grid(True)


axs[1,1].plot(
    data["timestep"],
    data["pe"]
)

axs[1,1].set_xlim(0, 400000)
axs[1,1].set_ylim(-19000, -17000)
axs[1,1].set_ylabel("Potential Energy")
axs[1,1].set_xlabel("Timestep (fs)")
axs[1,1].grid(True)

plt.show()
