import pandas as pd
import matplotlib.pyplot as plt

profiles = {}

with open("dens.profile") as f:
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
            f.readline()

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
        profiles[step]["density"]
    )

    plt.xlim(0,1500)
    plt.ylim(0,0.003)

    plt.pause(0.01)


#single profile
profile = profiles[200000]

plt.figure(figsize=(8,4))

plt.plot(
    profile["x"],
    profile["density"],
)

plt.xlabel("x (Å)")
plt.ylabel("Number density")
plt.grid()

plt.show()