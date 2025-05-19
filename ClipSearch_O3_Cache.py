import pandas as pd
import time
from collections import Counter
time_start = time.time()
def query_data_cached(hdf5_fp, key, cache):
    """
    For the given gene key, load (or retrieve from cache) the HDF5 data corresponding to 
    its (Chromosome, Strand) group and then perform an overlap query.
    Returns a Counter of proteins that overlap the query interval.
    """
    group_key = (key["Chromosome"], key["Strand"])
    
    # If this group hasn't been loaded yet, read and cache it.
    if group_key not in cache:
        print(f"Loading data for group {group_key} from disk...")
        query_clause = f'(Chromosome == "{group_key[0]}") & (Strand == "{group_key[1]}")'
        cache[group_key] = pd.read_hdf(hdf5_fp, key="data", where=query_clause)
    df_subset = cache[group_key]
    
    # Determine overlap: intervals overlap if (df_subset.Start <= key.Stop) and (df_subset.Stop >= key.Start)
    overlap_condition = (df_subset['Start'] <= key["Stop"]) & (df_subset['Stop'] >= key["Start"])
    df_overlap = df_subset[overlap_condition]
    
    return Counter(df_overlap['Protein'])

# File paths
excel_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_no_gap2.xlsx"
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"

# Read the Excel file, which we expect to have columns: 'ENSG', 'Chromosome', 'Start', 'End', 'Strand2'.
df = pd.read_excel(excel_fp)

# Adjust the Chromosome column to include the 'chr' prefix if missing.
df['Chromosome'] = df['Chromosome'].astype(str).apply(lambda x: x if x.lower().startswith('chr') else "chr" + x)

# Build a dictionary of gene keys.
gene_dict = {}

for idx, row in df.iterrows():
    gene_dict[row["ENSG"]] = {
        "ENSG": row["ENSG"],
        "Chromosome": row["Chromosome"],
        "Start": row["Start"],
        "Stop": row["End"],      # Map 'End' in Excel to 'Stop'
        "Strand": row["Strand2"]   # Map 'Strand2' (Excel) to 'Strand'
    }

# print("Gene dictionary for processing:")
# print(gene_dict)

# Set up a cache dictionary to hold HDF5 subsets keyed by (Chromosome, Strand).
cache = {}

# Global protein counter to sum all protein counts across every gene.
global_protein_counter = Counter()

# Process each gene and update the global counter.
for gene_id, key in gene_dict.items():
    protein_counts = query_data_cached(hdf5_fp, key, cache)
    global_protein_counter.update(protein_counts)
    
    # print(f"\nResults for gene {gene_id}:")
    # print(protein_counts)

# After processing, print out the overall global protein counts.
print("\nGlobal protein counter (summed over all queries):")
print(global_protein_counter)
time_stop = time.time()
print("Total time taken:", time_stop - time_start, "seconds")
print("Total number of genes processed:", len(gene_dict))
print("Total number of proteins found:", len(global_protein_counter))
print("Total number of unique proteins found:", len(set(global_protein_counter.elements())))

