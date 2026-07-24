import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from math import trunc

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


    cu_data = pd.DataFrame(
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


profiles = {}

with open("N2_shock/profiles/N2_gas.profile") as f:
    while True:

        line = f.readline()

        if not line:
            break

        if line.startswith("#"):
            continue

        words = line.split()

        if len(words) == 3:

            timestep = int(words[0])
            nchunks = int(words[1])

            # skip header
            #f.readline()

            data = []

            for i in range(nchunks):
                data.append(
                    list(map(float, f.readline().split()))
                )

            profiles[timestep] = pd.DataFrame(
                data,
                columns=[
                    "chunk",
                    "x",
                    "count",
                    "vx",
                    "temp",
                    "ndensity",
                    "mdensity"
                ],
            )


def data_plot():
    fig, axs = plt.subplots(2,2)

    axs[0,0].plot(
        cu_data["timestep"],
        cu_data["comx"]
    )

    axs[0,0].set_xlim(0, 400000)
    axs[0,0].set_ylim(17400, 19000)
    axs[0,0].set_ylabel("CoM position")
    #axs[0,0].set_xlabel("Timestep (fs)")
    axs[0,0].grid(True)


    axs[1,0].plot(
        cu_data["timestep"],
        cu_data["vx"]
    )

    axs[1,0].set_xlim(0, 400000)
    axs[1,0].set_ylim(-1, 10)
    axs[1,0].set_ylabel("CoM velocity")
    axs[1,0].set_xlabel("Timestep (fs)")
    axs[1,0].grid(True)


    axs[0,1].plot(
        cu_data["timestep"],
        cu_data["temp"]
    )

    axs[0,1].set_xlim(0, 400000)
    axs[0,1].set_ylim(0, 1300)
    axs[0,1].set_ylabel("Temperature (CoM corrected K)")
    #axs[0,1].set_xlabel("Timestep (fs)")
    axs[0,1].grid(True)


    axs[1,1].plot(
        cu_data["timestep"],
        cu_data["pe"]
    )

    axs[1,1].set_xlim(0, 400000)
    axs[1,1].set_ylim(-19000, -17000)
    axs[1,1].set_ylabel("Potential Energy")
    axs[1,1].set_xlabel("Timestep (fs)")
    axs[1,1].grid(True)

    plt.show()


def vel_plot():

    hit_time = 83000

    #predicted vel

    U = 14.21 #Ang/ps
    M = 5.77e-22
    rho = 43.35
    d = 5e-9
    A = np.pi*(d**2)/4
    mu = 5.33e-5
    C_d = 3
    C_cunn = 1.8

    tau_stokes = M/(3*mu*np.pi*d)
    tau_stokes_corr = tau_stokes*C_cunn
    tau_CD = 2*M*(1e15)/(C_d*rho*A*U*100)

    time_list = [(t-hit_time) for t in cu_data["timestep"]]

    vel_pred_stokes = []
    for t in time_list:
        vel_pred_stokes.append(U*(1-np.e**(-(t)*1e-15/tau_stokes)))

    vel_pred_stokes_corr = []
    for t in time_list:
        vel_pred_stokes_corr.append(U*(1-np.e**(-(t)*1e-15/tau_stokes_corr)))

    vel_pred_CD = []
    for t in time_list:
        vel_pred_CD.append(U*(t/(t+tau_CD)))

    

    vel_pred_adjust = [U - v for v in vel_pred_stokes]
    vel_measured_adjust = [U - v for v in cu_data["vx"]]


    vel_loglin(time_list, vel_pred_CD, cu_data["vx"], vel_pred_stokes, vel_pred_stokes_corr)
    #vel_loglog(time_list, vel_pred_stokes)
    #vel_loglin(time_list, vel_pred_adjust, vel_measured_adjust)


    plt.xlim(0, 350)
    plt.ylim(0, 1100)
    plt.ylabel("CoM velocity (m/s)")
    plt.xlabel("Time (ps)")
    plt.grid(True)

    #plt.legend()
    plt.show()

def vel_loglog(time_list,vel_pred):

    plt.plot(
        time_list,
        cu_data["vx"],
        label = "measured"
    )


    plt.plot(
        time_list,
        vel_pred,
        label = "predicted"
    )

    plt.xscale('log')
    plt.yscale('log')


def vel_loglin(time_list, vel_pred_adjust, vel_measured_adjust, vel_pred_2, vel_pred_3):

    plt.plot(
        [t/1000 for t in time_list],
        [vel*100 for vel in vel_measured_adjust],
        label = "measured",
        #marker = "o"
        
    )

    #plt.scatter([t/1000 for t in time_list][::10],[vel*100 for vel in vel_measured_adjust][::10])


    plt.plot(
        [t/1000 for t in time_list],
        [vel*100 for vel in vel_pred_2],
        label = "Stokes",
        linestyle = "dashed"
    )


    plt.plot(
        [t/1000 for t in time_list],
         [vel*100 for vel in vel_pred_adjust],
        label = "CD const",
        linestyle = "dashed"
    )

    
    plt.plot(
        [t/1000 for t in time_list],
        [vel*100 for vel in vel_pred_3],
        label = "Stokes Corrected (Cunningham)",
        linestyle = "dashed"
    )


    #plt.yscale('log')


def temp_gas():

    temp_N2 = []
    for i in range(0,len(cu_data["timestep"])):
        pos = trunc(cu_data["comx"][i]/10)*10 + 5
        pos_temp = list(profiles[cu_data["timestep"][i]]["x"]).index(pos)
        temp_N2.append(profiles[cu_data["timestep"][i]]["temp"][pos_temp])


    plt.plot(
            cu_data["timestep"],
            cu_data["temp"],
            label = "measured"
        )
    
    plt.plot(
        cu_data["timestep"],
        temp_N2,
        label = "N2 gas"
    )

    plt.xlim(0, 400000)
    plt.ylim(0,1800)
    plt.ylabel("Temperature (CoM corrected K)")
    plt.xlabel("Timestep (fs)")
    plt.grid(True)

    plt.legend()
    plt.show()

plt.rcParams["figure.figsize"] = (12, 6)


#data_plot()
vel_plot()
#temp_gas()