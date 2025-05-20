import pandas as pd
import time
from collections import Counter

# File paths
excel_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_gap2.xlsx"
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"

# Read the Excel file.
# Expected columns: 'ENSG', 'Chromosome', 'Start', 'End', and 'Strand2'
df = pd.read_excel(excel_fp)

# Adjust the Chromosome values so they match the HDF5 format (i.e., 'chr1', 'chr2', etc.)
df['Chromosome'] = df['Chromosome'].astype(str).apply(lambda x: x if x.lower().startswith('chr') else "chr" + x)

# Build a dictionary of gene query keys.
gene_dict = {}
# For testing purposes, use only the first 10 genes; remove .head(10) to process all keys.
for idx, row in df.iterrows():
    gene_dict[row["ENSG"]] = {
        "ENSG": row["ENSG"],
        "Chromosome": row["Chromosome"],
        "Start": row["Start"],
        "Stop": row["End"],      # Map 'End' in Excel to 'Stop'
        "Strand": row["Strand2"]   # Map 'Strand2' in Excel to 'Strand'
    }

# print("Gene dictionary for processing:")
# print(gene_dict)

# Load the entire HDF5 dataset once into memory.
print("Loading entire HDF5 dataset into memory...")
load_start = time.time()
df_full = pd.read_hdf(hdf5_fp, key="data")
load_end = time.time()
print(f"Loaded full dataset in {load_end - load_start:.2f} seconds, total rows: {len(df_full)}")

# Initialize a global protein counter
global_protein_counter = Counter()
results = {}

# Process each gene query against the full dataset.
query_total_start = time.time()
for gene_id, key in gene_dict.items():
    # Filter by Chromosome and Strand from the full dataset.
    df_subset = df_full[(df_full['Chromosome'] == key["Chromosome"]) & (df_full['Strand'] == key["Strand"])]
    
    # Use vectorized filtering to check for overlaps:
    # An overlap exists if: (df_subset.Start <= key.Stop) and (df_subset.Stop >= key.Start)
    overlap_condition = (df_subset['Start'] <= key["Stop"]) & (df_subset['Stop'] >= key["Start"])
    df_overlap = df_subset[overlap_condition]
    
    # Count overlapping proteins.
    protein_counts = Counter(df_overlap['Protein'])
    results[gene_id] = protein_counts
    
    # Update the global counter with the counts from this gene.
    global_protein_counter.update(protein_counts)
    
    # print(f"\nResults for gene {gene_id}:")
    # print(protein_counts)
query_total_end = time.time()

print("\nGlobal protein counter (summed over all queries):")
print(global_protein_counter)
print(f"\nProcessed {len(gene_dict)} gene queries in {query_total_end - query_total_start:.2f} seconds.")
print(f"Total time taken: {query_total_end - load_start:.2f} seconds.")
print(f"Total number of unique proteins found: {len(global_protein_counter)}")
print(f"Total number of unique genes processed: {len(gene_dict)}")
