import pandas as pd
from collections import Counter
import time

# File paths
time_open = time.time()
filepath = r"C:\Users\caleb\OneDrive\Brown\ProjectG\processed_data.h5"
time_close = time.time()
df = pd.read_excel(r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_gap2.xlsx")
fp = pd.read_hdf(filepath, key="data")
# Single test gene
ENSG00000228794 = {"Chromosome": "chr1", "Start": 825138, "Stop": 868835, "Strand": "-"}

start_time = time.time()# Initialize a counter for the number of times each chromosome appears
protein_counts = Counter()
num_lines = 0
for i in range(len(fp)):
    if fp.loc[i, "Chromosome"] == ENSG00000228794["Chromosome"]:
        if fp.loc[i, "Strand"] == ENSG00000228794["Strand"]:
            # Check for overlap between the start and stop indices
            if (
                fp.loc[i, "Start"] >= ENSG00000228794["Start"] and fp.loc[i, "Stop"] <= ENSG00000228794["Stop"]
            ) or (
                fp.loc[i, "Start"] <= ENSG00000228794["Start"] and fp.loc[i, "Stop"] >= ENSG00000228794["Start"]
            ) or (
                fp.loc[i, "Stop"] >= ENSG00000228794["Start"] and fp.loc[i, "Start"] <= ENSG00000228794["Stop"]
            ):
                protein_counts[fp.loc[i, "Protein"]] += 1
        num_lines += 1
end_time = time.time()

# Print the results
print(protein_counts)
print("Ran through all", num_lines, "lines in", end_time-start_time, "seconds")
print("Opened file in", time_close-time_open, "seconds")
