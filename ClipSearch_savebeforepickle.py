import pandas as pd
from collections import Counter
import time

# File paths
fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.txt"
df = pd.read_excel(r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_gap2.xlsx")

# The following 5 lines are not necessary but seem to improve efficiency
# Convert the relevant columns to lists for faster access
ENSG = df["ENSG"].tolist()  # Convert the column to a list
Chromosome = df["Chromosome"].tolist()  # Convert the column to a list
Start = df["Start"].tolist()  # Convert the column to a list
Stop = df["End"].tolist()  # Convert the column to a list
Strand = df["Strand"].tolist()  # Convert the column to a list

# Initalize lists (this may have to be modified to increase efficiency)
chromosomes = []
start_index = []
stop_index = []
strands = []
proteins = []
# Read the file and count the occurrences of each line
# ENSG00000228794 = {"Chromosome": "chr1", "Start": 825138, "Stop": 868835, "Strand": "-"}
key = {"ENSG": "ENSG00000272145","Chromosome": "chr1", "Start": 40669089, "Stop": 40692086, "Strand": "+"}
time_open = time.time()
# Open file and turn into lists
with open(fp) as f: # The total number of lines in this dataset is 38,361,545
    count = 0
    for line in f:
        field = line.strip().split("\t")
        # Extract Relevant Fields
        chromosome = field[0]
        start = int(field[1])
        stop = int(field[2])
        strand = field[4]
        protein = field[5]
        

        # Append to lists
        chromosomes.append(chromosome)
        start_index.append(start)
        stop_index.append(stop)
        strands.append(strand)
        proteins.append(protein)
        # print(f"Chromosome: {chromosome}, Start: {start}, Stop: {stop}, Strand: {strand}")
        count += 1
        # if count >= 1000000: # Only took 4 seconds to read 1 million lines and append to lists
        #     break
# key = ENSG00000228794
print('Key:', key)
time_close = time.time()
start_time = time.time()# Initialize a counter for the number of times each chromosome appears
protein_counts = Counter()
for i in range(len(chromosomes)):
    if chromosomes[i] == key["Chromosome"]:
        # print("Chromosomes match")
        if strands[i] == key["Strand"]:
            # print("Strands match")
            # Next I will check if there is any overlap between the start and stop indices
            if start_index[i] >= key["Start"] and stop_index[i] <= key["Stop"]:
                # print("Indices match")                
                protein_counts[proteins[i]] += 1
            elif start_index[i] <= key["Start"] and stop_index[i] >= key["Start"]:
                print("Indices match")
                protein_counts[proteins[i]] += 1
            elif start_index[i] <= key["Stop"] and stop_index[i] >= key["Stop"]:
                print("Indices match")
                protein_counts[proteins[i]] += 1
            elif start_index[i] >= key["Stop"]:
                break
end_time = time.time()
print(protein_counts)
print("Ran through all", len(chromosomes), "lines in", end_time-start_time, "seconds")
print("Opened file in", time_close-time_open, "seconds")

#         # This next if statement should check for overlap between the start and stop indices
#           If everything matches, then this is a protein binding partner. It needs to be stored in a counter
#           Syntax will look something like this:
#           protein_counts = Counter()
#           protein_counts[chromosome] += 1 # This will count the number of times each chromosome appears
#   
# else:
#     print("This is not a protein binding partner")
