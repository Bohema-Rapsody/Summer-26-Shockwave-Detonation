import pandas as pd

def read_lammps_log(filename):
    tables = []

    with open(filename) as f:

        reading = False
        columns = None
        data = []

        for line in f:

            line = line.strip()

            # Start of a new thermo table
            if line.startswith("Step"):

                columns = line.split()
                data = []
                reading = True
                continue

            if reading:

                # End of current table
                cutoff = 40
                if line.startswith("Loop time") and len(data) > cutoff:
                    tables.append(pd.DataFrame(data[cutoff:], columns=columns)) #truncates the first 10 data point for stability
                    reading = False
                    continue

                parts = line.split()

                if len(parts) != len(columns):
                    continue

                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    pass

    return tables

log_data = read_lammps_log("Ar/ideal_law_log.LAMMPS")

print(log_data[0].head())
#print(df["Press"])