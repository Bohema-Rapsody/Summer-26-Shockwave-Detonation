import csv

from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

# Read trajectory
pipeline = import_file("N2_Shock/dump/Cu_dump.xyz")

# Perform CNA
pipeline.modifiers.append(CommonNeighborAnalysisModifier())

# Create output file
with open("Cu-N2/data/Cu_structure_analysis.csv", "w", newline="") as f:

    writer = csv.writer(f)

    # Header
    writer.writerow([
        "Frame",
        "Other",
        "FCC",
        "HCP",
        "BCC",
        "ICO"  
    ])

    # Loop through every frame
    for frame in range(pipeline.source.num_frames):

        data = pipeline.compute(frame)

        table = data.tables["structures"]

        # Extract the counts
        counts = table["Count"]

        writer.writerow([
            frame,
            int(counts[0]),   # Other
            int(counts[1]),   # FCC
            int(counts[2]),   # HCP
            int(counts[3]),   # BCC
            int(counts[4])    # ICO
        ])

print("Finished.")