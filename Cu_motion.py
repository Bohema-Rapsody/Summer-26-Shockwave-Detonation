import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from math import trunc

#profiles = {}

#constants:

U = 14.21 #Ang/ps
M = 5.77e-22
N = 5473
rho = 43.35
d = 5e-9
A = np.pi*(d**2)/4
mu = 5.33e-5
Cp_N2 = 1245.1
Cp_Cu = 460
k_N2 = 0.0875

Re = rho*(U*100)*d/mu
Pr = Cp_N2*mu/k_N2

C_d = 3
C_cunn = 1.8
k_B = 8.617e-5
q = 1.6e-19

#time constants
tau_stokes = M/(3*mu*np.pi*d)
tau_stokes_corr = tau_stokes*C_cunn
tau_CD = 2*M*(1e15)/(C_d*rho*A*U*100)

hit_time = 83000 #potential to calculate completely


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
                list(map(float, line.split()))
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


# Load only the columns we need
# Columns:
# 0 Temperature
# 2 Density
# 8 Cp
# 11 Viscosity
# 12 Thermal conductivity

N2_data = np.loadtxt(
    "N2_transport_data.data",
    skiprows=1,
    usecols=(0, 2, 8, 11, 12)
)

#temperature = data[:, 0]
#density = data[:, 1]
#cp = data[:, 2]
#viscosity = data[:, 3]
#conductivity = data[:, 4]



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



    #predicted vel
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

    temp_N2 = current_gas_temp()

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


def current_gas_temp():
    temp_N2 = []
    for i in range(0,len(cu_data["timestep"])):
        pos = trunc(cu_data["comx"][i]/10)*10 + 5
        pos_temp = list(profiles[cu_data["timestep"][i]]["x"]).index(pos)
        temp_N2.append(profiles[cu_data["timestep"][i]]["temp"][pos_temp])

    return temp_N2


def energy_stagnation():

    time_list = [(t-hit_time) for t in cu_data["timestep"]]

    ke_change = [e-cu_data["ke"][hit_time/1000] for e in cu_data["ke"]]
    pe_change = [e-cu_data["pe"][hit_time/1000] for e in cu_data["pe"]]
    total_e_change_KE_PE = [k+p for k,p in zip(ke_change,pe_change)]

    '''
    plt.plot(
        time_list,
        ke_change,
        label = "KE change"
    )

    
    plt.plot(
        time_list,
        pe_change,
        label = "PE change"
    )
    
    
    plt.plot(
        time_list,
        total_e_change_KE_PE,
        label = "E Total change (KE+PE)"
    )
    '''
    E_velocity = kinetic_energy_breakdown(time_list)
    E_stokes_model = energy_transfer_model(time_list)
    
    plt.plot(
        time_list,
        [(T-E)*q/(M*Cp_Cu) for T,E in zip(total_e_change_KE_PE,E_velocity)],
        label = "E total & Stokes model difference (velocity term)"
    )
    

    plt.xlim(0, 350000)
    #plt.ylim(0,1800)
    #plt.ylabel("Energy (eV)")
    plt.ylabel("Temperature (K)")
    plt.xlabel("Timestep (fs)")
    plt.grid(True)

    plt.legend()
    plt.show()


def kinetic_energy_breakdown(time_list):

    v_ke_change = [M*(v*100)**2/(2*q) for v in cu_data["vx"]]
    T_ke_change = [(3*N-3)*k_B*(T-cu_data["temp"][hit_time/1000])/2 for T in cu_data["temp"]]
    total_e_change_v_T = [v+T for v,T in zip(v_ke_change,T_ke_change)]

    
    return v_ke_change
    plt.plot(
        time_list,
        v_ke_change,
        label = "KE change (velocity)"
    )

    
    
    plt.plot(
        time_list,
        T_ke_change,
        label = "KE change (temperature)"
    )
    
    
    plt.plot(
        time_list,
        total_e_change_v_T,
        label = "KE Total change (vel+temp)"
    )
    
    

    

def energy_transfer_model(time_list):


    t_list = time_list[int(hit_time/1000):]
    vel_data = cu_data["vx"][int(hit_time/1000):]

    #Stagnation model
    stagnation_delta_E = [((U-v)*100)**3*rho*np.pi*d**2/(8*q) for v in vel_data]
    stagnation_E_intgr = integrator(stagnation_delta_E,1e-12)

    stagnation_E_model = [Re*C_cunn*M*(U*100)**2/(72*q) *(1-np.e**(-3*t*1e-15/tau_stokes_corr)) for t in time_list]

    '''
    plt.plot(
            t_list,
            stagnation_E_intgr,
            label = "Stagnation Energy Integrated"
        )
    
    
    
    plt.plot(
            time_list,
            stagnation_E_model,
            label = "Stagnation Energy Model (Stokes corrected)"
        )
    '''

    #Conduction model
    conduction_delta_E = []
    conduction_E_intgr = []

    
    for i in range(int(hit_time/1000),len(time_list)):
        #Nu = 3.06 at max flow
        h_Nu, T_g = update_flow_parameters(i)
        #conduction_delta_E.append([h*np.pi*d**2*(T_g - cu_data["temp"][i])/q for h in h_Nu])
        conduction_delta_E.append([h*np.pi*d**2*(T_g - cu_data["temp"][i])/(M*Cp_Cu) for h in h_Nu])

    #print(len(conduction_delta_E), conduction_delta_E[0])
    for i in range(0,len(conduction_delta_E[0])):
        conduction_E_intgr.append(integrator([dE[i] for dE in conduction_delta_E],1e-12))

    plt.plot(
            time_list,
            [T-cu_data["temp"][int(hit_time/1000)] for T in cu_data["temp"]],
            label = "measured"
        )


    plt.plot(
            t_list,
            conduction_E_intgr[0],
            label = "Conduction Energy (Ranz-Marshall)"
        )

    plt.plot(
            t_list,
            conduction_E_intgr[1],
            label = "Conduction Energy (Faeth)"
        )

    plt.plot(
            t_list,
            conduction_E_intgr[2],
            label = "Conduction Energy (Whitaker)"
        )

   
    return stagnation_E_model
    
def update_flow_parameters(step):

    #from NIST nitrgoen data: #see N2_transport_data.data
    #at 1400K, rho = 44.83, mu = 5.216e-5, k = 0.0854, Cp = 1240
    #at 1600K, rho = 39.45, mu = 5.679e-5, k = 0.0938, Cp = 1260
    #step = int(time/1000)
    gas_temp = current_gas_temp()[step]
    current_vel = cu_data["vx"][step]

    surf_temp = cu_data["temp"][step]

    ave_temp = (gas_temp*0.5 + surf_temp*0.5)

    temp_eval = ave_temp

    rho_v = linear_approx(N2_data[:,0],N2_data[:,1],temp_eval)
    Cp_v = linear_approx(N2_data[:,0],N2_data[:,2],temp_eval)*1000
    mu_v = linear_approx(N2_data[:,0],N2_data[:,3],temp_eval)
    mu_inf = linear_approx(N2_data[:,0],N2_data[:,3],gas_temp)
    mu_surf = linear_approx(N2_data[:,0],N2_data[:,3],surf_temp)
    k_v = linear_approx(N2_data[:,0],N2_data[:,4],temp_eval)


    Re_v = rho_v*((U-current_vel)*100)*d/mu_v
    Pr_v = Cp_v*mu_v/k_v

    #Simpler model
    #Ranz-Marshall Correlation
    Nu_v_RM = 2 + 0.6*Re_v**(1/2)*Pr_v**(1/3)

    # Test correlation obtained from Faeth https://www.sciencedirect.com/science/article/pii/0360128577900120?via%3Dihub
    Nu_v_Fa = 2 + 0.555*Re_v**(1/2)*Pr_v**(1/3)/(1+1.232/(Re_v*Pr_v**(4/3)))**(1/2) 
    
    #Whitaker Correlation
    Nu_v_Wh = 2 + (0.4*Re_v**(1/2) + 0.06*Re_v**(2/3))*Pr_v**(0.4)*(mu_inf/mu_surf)**(1/4)

    Nu_list = [Nu_v_RM, Nu_v_Fa, Nu_v_Wh]
    #Nu_list = [Nu_v_RM]
    h_Nu_v = [Nu_v*k_v/d for Nu_v in Nu_list]

    print(step, temp_eval, current_vel, rho_v, mu_v, k_v, Cp_v, Re_v, Pr_v, Nu_list, h_Nu_v)
    return h_Nu_v,ave_temp

def linear_approx(temp_list, prop_list,T):
    for i in range(len(temp_list)-1):
        if temp_list[i] <= T <= temp_list[i+1]:

            frac = (T - temp_list[i]) / (temp_list[i+1] - temp_list[i])

            return prop_list[i] + frac * (prop_list[i+1] - prop_list[i])

    raise ValueError("Temperature outside table.")

def integrator(delta_list,step):

    integrated_list = [0]
    for i in range(0,len(delta_list)-1):
        integrated_list.append(step*delta_list[i]+integrated_list[-1])

    return integrated_list


plt.rcParams["figure.figsize"] = (12, 6)
#data_plot()
#vel_plot()
#temp_gas()
energy_stagnation()