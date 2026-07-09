import pandas as pd
import matplotlib.pyplot as plt

profiles = {}

with open("N2_Shock\dens.profile") as f:
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
                    "density",
                ],
            )


#shock animation
for step in sorted(profiles):

    plt.clf()

    plt.plot(
        profiles[step]["x"],
        profiles[step]["density"]/2
        
        #profiles[step]["temp"]
    )


    plt.xlim(0,4500)
    plt.xlabel("x (Å)")
    plt.ylim(0,0.003)
    plt.ylabel("Number density")
    #plt.ylim(0,1800)
    #plt.ylabel("Temperature (CoM corrected K)")
    plt.grid()

    plt.pause(0.01)



#single profile
#print(next(reversed(profiles.keys())))
profile = profiles[next(reversed(profiles.keys()))]

plt.figure(figsize=(8,4))

plt.plot(
    profile["x"],
    #profile["density"]/2,
    profile["temp"],
)

plt.xlim(0,4500)
plt.xlabel("x (Å)")
#plt.ylim(0,0.003)
#plt.ylabel("Number density")
plt.ylim(0,1800)
plt.ylabel("Temperature (CoM corrected K)")
plt.grid()

plt.show()