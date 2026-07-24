import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

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


def determine_conds(profile,check_val,init_tol=0.8,tol=0.8):
    # pre-shock number density 0.000338, post-shock number density 0.00211  -recall two atoms per molecule
    # pre-shock mass density 0.00786, post-shock mass density 0.0491
    # pre-shock temp 300, post-shock temp 1350


    #check condition, first valid index (start), averaged post-shock value, final valid index (end)
    init_table = [check_val,'','','']

    #calulating shock threshold values
    for i in range(0,len(profile)):
        val = profile[i]

        if init_table[1] != '':
            #checks whether following values fits expected average
            if (init_table[2]*tol < val < init_table[2]/tol) or (i-init_table[1] < 6):
                fac = 1/(i+1-init_table[1])**1.1
                init_table[2] = init_table[2]*(1-fac) + val*fac
            else:
                init_table[3] = i
                break

        elif init_table[0]*init_tol < val < init_table[0]/init_tol:
            #sets initial valid value
            init_table[1] = i
            init_table[2] = val
        
    return init_table

def plot_tangents(step,init_table,x_profile,v_profile):
    global shock_data
    if not init_table[0][3] == '':
        plt.hlines(init_table[0][2],x_profile[init_table[0][1]],x_profile[init_table[0][3]])

    num_bins = len(x_profile)-1
    if not init_table[1][3] == '':
        plt.hlines(init_table[1][2],x_profile[num_bins-init_table[1][3]],x_profile[num_bins-init_table[1][1]])


    if not (init_table[0][3] == '' or init_table[1][3] == ''):
        start = init_table[0][3] + 3
        end = num_bins-init_table[1][3] - 3

        if end-start >5:
            coef = np.polyfit(x_profile[start:end], v_profile[start:end], 1)
            func = np.poly1d(coef)

            plt.plot(x_profile, func(x_profile))

            if coef[0] != 0:
                x_post = (init_table[0][2] - coef[1])/coef[0]
                x_pre = (init_table[1][2] - coef[1])/coef[0]
                print("Time step: ",step,"post-shock front at: ",x_post," Ang, pre-shock front at: ",x_pre," Ang\n" \
                "Total shock width: ",abs(x_pre-x_post)," Ang, calculated with ",end-start," values\n")
                shock_data.append([init_table[0][2],init_table[1][2],abs(x_pre-x_post)])
                return

    print("Time step:",step,"Shock wave undefined\n")

plt.rcParams["figure.figsize"] = (17, 8)
fig, axs = plt.subplots(2)
#shock animation
shock_data = []
for step in sorted(profiles):

    #Determining thickness from data
    #init_ndens = [determine_conds(profiles[step]["ndensity"],0.00211),determine_conds(list(reversed(profiles[step]["ndensity"])),0.000338,tol=0.6)]
    #init_temp = [determine_conds(profiles[step]["temp"],1350),determine_conds(list(reversed(profiles[step]["temp"])),300,tol=0.6)]

    #plot_tangents(step,init_ndens,profiles[step]["x"],profiles[step]["ndensity"])
    #plot_tangents(step,init_temp,profiles[step]["x"],profiles[step]["temp"])#-----------------------------------------------------------------------CHANGE
    #Note, shock size should be x = 268.0 Ang

    

    
    axs[0].plot(
        profiles[step]["x"],
        profiles[step]["ndensity"]
    )

    axs[1].plot(
        profiles[step]["x"],
        profiles[step]["temp"]
    )

    # Set x-limits on both plots
    axs[0].set_xlim(10000, 24000)
    axs[1].set_xlim(10000, 24000)

    # Set y-limits
    axs[0].set_ylim(0, 0.003)
    axs[1].set_ylim(0, 1800)

    # Labels
    axs[0].set_ylabel("Number density")
    axs[1].set_ylabel("Temperature (CoM corrected K)")
    axs[1].set_xlabel("x (Å)")   # Bottom subplot only

    # Grid
    axs[0].grid(True)
    axs[1].grid(True)

    plt.pause(0.1)

    axs[0].cla()
    axs[1].cla()
    
shock_data = np.array(shock_data)
#print("Average calculated post-shock conditions: ",f"{stats.trim_mean(shock_data[:,0], 0.2):.3g}"," units, Predicted: ",init_temp[0][0])#-----------CHANGE
#print("Average calculated pre-shock conditions: ",f"{stats.trim_mean(shock_data[:,1], 0.2):.3g}"," units, Predicted: ",init_temp[1][0])
#print("Average calculated shock width: ",f"{stats.trim_mean(shock_data[:,2], 0.2):.3g}"," Ang")



#single profile

#print(next(reversed(profiles.keys())))
profile = profiles[next(reversed(profiles.keys()))]

plt.figure(figsize=(8,4))

plt.plot(
    profile["x"],
    #profile["ndensity"],
    profile["temp"],
)

plt.xlim(10000,24000)
plt.xlabel("x (Å)")
#plt.ylim(0,0.003)
#plt.ylabel("Number density")
plt.ylim(0,1800)
plt.ylabel("Temperature (CoM corrected K)")
plt.grid()

plt.show()
