# Pre-process data
import pandas as pd

# File path to the large dataset
fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.txt"

# Initialize lists to store data
chromosomes = []
start_index = []
stop_index = []
strands = []
proteins = []
print("Opening file...")
print("lists initialized")
# Open the file and process it
with open(fp) as f:
    for line in f:
        field = line.strip().split("\t")
        # Extract relevant fields
        chromosomes.append(field[0])
        start_index.append(int(field[1]))
        stop_index.append(int(field[2]))
        strands.append(field[4])
        proteins.append(field[5])
print("File opened and data extracted")
# Create a DataFrame from the lists
df = pd.DataFrame({
    "Chromosome": chromosomes,
    "Start": start_index,
    "Stop": stop_index,
    "Strand": strands,
    "Protein": proteins
})
print("DataFrame created")

# Save the DataFrame to an HDF5 file
output_file = r"C:\Users\caleb\OneDrive\Brown\ProjectG\processed_data.h5"
df.to_hdf(output_file, key="data", mode="w", format="table")
print(f"Data saved to {output_file} in HDF5 format")