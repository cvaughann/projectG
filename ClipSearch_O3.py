### ClipSearch_O3
import pandas as pd
import time
from collections import Counter

# HDF5 file path.
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"

# Define the query key.
# key = {"ENSG": "ENSG00000228794","Chromosome": "chr1", "Start": 825138, "Stop": 868835, "Strand": "-"}
# key = {"ENSG": "ENSG00000272145","Chromosome": "chr1", "Start": 40669089, "Stop": 40692086, "Strand": "+"}
key = {"ENSG": "ENSG00000176728", "Chromosome": 'chrY', 'Start': 18772706, 'Stop': 19077555, 'Strand': '-'}

def query_data(hdf5_fp, key):
    # Use HDF5's partial I/O to load only rows matching our Chromosome and Strand.
    query_clause = f'(Chromosome == "{key["Chromosome"]}") & (Strand == "{key["Strand"]}")'
    load_start = time.time()
    df_subset = pd.read_hdf(hdf5_fp, key="data", where=query_clause)
    load_end = time.time()
    print(f"Opened HDF5 partial data in {load_end - load_start:.2f} seconds.")

    # Display the key being used.
    print("Key:", key)

    # Use vectorized operations to test for overlaps.
    overlap_condition = (df_subset['Start'] <= key["Stop"]) & (df_subset['Stop'] >= key["Start"])
    query_start = time.time()
    df_overlap = df_subset[overlap_condition]
    protein_counts = Counter(df_overlap['Protein'])
    query_end = time.time()

    print(protein_counts)
    print(f"Ran through {len(df_subset)} filtered rows in {query_end - query_start:.2f} seconds.")

if __name__ == "__main__":
    query_data(hdf5_fp, key)
