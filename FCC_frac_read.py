import csv

from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier
from ovito.modifiers import PolyhedralTemplateMatchingModifier

# Read trajectory
pipeline = import_file("N2_Shock/dump/Ti_dump_300.xyz")

# Perform CNA
pipeline.modifiers.append(CommonNeighborAnalysisModifier())
#pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())

save_file = "Ti-N2/data/Ti_300_structure_analysis_CNA.csv"
# Create output file
with open(save_file, "w", newline="") as f:

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