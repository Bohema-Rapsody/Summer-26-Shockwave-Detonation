import pandas as pd
import os
import tempfile


def xyz_dump_clean(filename):

    timestep_positions = {}

    # -------------------------
    # First pass: find timesteps
    # -------------------------
    with open(filename, "r") as f:

        while True:

            position = f.tell()

            line = f.readline()

            if not line:
                break

            line = line.strip()

            if not line:
                continue

            natoms = int(line)

            header = f.readline()

            timestep = int(
                header.split("Timestep:")[1].strip()
            )

            # Last occurrence wins
            timestep_positions[timestep] = position

            # Skip coordinates
            for _ in range(natoms):
                f.readline()

    print(f"Found {len(timestep_positions)} unique timesteps.")

    # --------------------------------
    # Create temporary file in same dir
    # --------------------------------
    directory = os.path.dirname(os.path.abspath(filename))

    fd, temp_filename = tempfile.mkstemp(
        suffix=".tmp",
        dir=directory
    )

    os.close(fd)

    try:

        # -------------------------
        # Second pass: write clean
        # -------------------------
        with open(filename, "r") as f, open(temp_filename, "w") as out:

            for timestep in sorted(timestep_positions):

                f.seek(timestep_positions[timestep])

                natoms_line = f.readline()
                header = f.readline()

                out.write(natoms_line)
                out.write(header)

                natoms = int(natoms_line.strip())

                for _ in range(natoms):
                    out.write(f.readline())

        # -------------------------
        # Replace original
        # -------------------------
        os.replace(temp_filename, filename)

        print(f"Cleaned file: {filename}")

    except Exception:

        # If something goes wrong, remove temporary file
        if os.path.exists(temp_filename):
            os.remove(temp_filename)

        raise







def multi_line_profile_clean(file_path):
    #Removes duplicate timesteps for large, mulit-line database

    profiles = {}
    header = []

    with open(file_path, "r") as f:

        while True:

            line = f.readline()

            if not line:
                break

            if line.startswith("#"):
                header.append(line)
                continue


            words = line.split()

            if len(words) == 3:

                timestep = int(words[0])
                nchunks = int(words[1])
                total_count = float(words[2])

                data = []

                for i in range(nchunks):
                    data.append(
                        list(map(float, f.readline().split()))
                    )

                profiles[timestep] = {
                    "nchunks": nchunks,
                    "total_count": total_count,
                    "data": data
                }


    # Sort by timestep
    profiles = dict(sorted(profiles.items()))


    # ------------------------------------------------------------
    # Write cleaned file
    # ------------------------------------------------------------

    with open(file_path, "w") as f:

        for line in header:
            f.write(line)

        for timestep, profile in profiles.items():

            data = profile["data"]

            f.write(
                f"{timestep} {len(data)} {profile['total_count']}\n"
            )

            for row in data:

                f.write(
                    " ".join(
                        f"{value:.10g}" for value in row
                    )
                    + "\n"
                )


    print(f"Cleaned file:  {file_path}")
    print(f"Unique timesteps: {len(profiles)}")









def single_line_profile_clean(file_path):
    #Removes duplicate timesteps for simple, single-line database

        data = {}
        header = []

        with open(file_path, "r") as f:
            for line in f:
                # Preserve all comment/header lines
                if line.startswith("#"):
                    header.append(line)
                    continue

                words = line.split()

                if not words:
                    continue

                # First word is always the timestep
                timestep = int(words[0])

                # Later occurrences overwrite earlier ones
                data[timestep] = line

        # Write sorted by timestep
        with open(file_path, "w") as f:
            for line in header:
                f.write(line)

            for timestep in sorted(data):
                f.write(data[timestep])

        print(f"Cleaned file: {file_path}")
        print(f"Unique timesteps: {len(data)}")



#Dump clean used for .xyz files
xyz_dump_clean("N2_shock/dump/Ti_dump_300.xyz")
xyz_dump_clean("N2_shock/dump/all_dump_Ti_300.xyz")

#Multi used for Particle Grid, N2 gas chunks, narrow N2 gas chunks
multi_line_profile_clean("N2_shock/profiles/N2_gas_Ti_300.profile")
multi_line_profile_clean("N2_shock/profiles/N2_gas_narrow_Ti_300.profile")
multi_line_profile_clean("N2_shock/profiles/Ti_300_grid_properties.profile")

#Single used for whole particle data, gas shell data
single_line_profile_clean("N2_shock/profiles/Ti_particle_300.profile")
single_line_profile_clean("N2_shock/profiles/gas_shell_temperature_Ti_300.profile")
