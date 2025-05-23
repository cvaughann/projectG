import pandas as pd
import time
from collections import Counter
import numpy as np
time_start = time.time()
### Make user confirm files are closed before running
user_answer = input("Please confirm that the Excel file is closed before running this script (y/n): ")
if user_answer.lower() != 'y':
    print("Please close the Excel file and run the script again.")
    exit()
''' Old Function
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
    overlap_condition = (
        ((df_subset['Start'] <= key["Hstart"]) & (df_subset['Stop'] >= key["Hstart"])) |
        ((df_subset['Start'] <= key["Wstart"]) & (df_subset['Stop'] >= key["Wstart"])) |
        ((df_subset['Start'] <= key["Cstart"]) & (df_subset['Stop'] >= key["Cstart"]))
    )    
    # overlap_conditionW = (df_subset['Start'] <= key["Wstart"]) & (df_subset['Stop'] >= key["Wstart"])
    # overlap_conditionC = (df_subset['Start'] <= key["Cstart"]) & (df_subset['Stop'] >= key["Cstart"])
    df_overlap = df_subset[overlap_condition]
    
    return Counter(df_overlap['Protein'])
'''
def query_data_cached(hdf5_fp, key, cache, *, 
                      return_detail=False,     # opt-in switch for the caller
                      inclusive=True):         # set False if you treat intervals as half-open
    """
    For the given gene key, load (or retrieve from cache) the HDF5 data corresponding to 
    its (Chromosome, Strand) group and then perform an overlap query.
    Returns a Counter of proteins that overlap the query interval.
    """
    group_key = (key["Chromosome"], key["Strand"])
        # ---------------------------------------------------------------
    # Skip records whose (Chromosome, Strand) are both NaN
    # ---------------------------------------------------------------
    # works for real NaNs *and* for the string "nan"
    if (pd.isna(group_key[0]) or str(group_key[0]).lower() == "nan") and \
       (pd.isna(group_key[1]) or str(group_key[1]).lower() == "nan"):

        # Nothing to do → return empty results in the same format
        return (Counter(), pd.DataFrame()) if return_detail else Counter()
    # If this group hasn't been loaded yet, read and cache it.
    if group_key not in cache:
        print(f"Loading data for group {group_key} from disk...")
        query_clause = f'(Chromosome == "{group_key[0]}") & (Strand == "{group_key[1]}")'
        cache[group_key] = pd.read_hdf(hdf5_fp, key="data", where=query_clause)
    df_subset = cache[group_key]
        # ---------- build masks for the three strand-specific intervals
    intervals = [
        ("H", key["Hstart"], key["Hstop"]),
        ("W", key["Wstart"], key["Wstop"]),
        ("C", key["Cstart"], key["Cstop"]),
    ]

    # masks = [
    #     (df_subset["Start"] <= hi) & (df_subset["Stop"] >= lo)
    #     for _tag, lo, hi in intervals]

    # overlap_any = masks[0] | masks[1] | masks[2]
    # if not overlap_any.any():
    #     return (Counter(), pd.DataFrame()) if return_detail else Counter()

    # ---------- for every strand that overlaps, compute the span ----
    detail_frames = []
    for tag, lo_raw, hi_raw in intervals:
        lo = pd.to_numeric(lo_raw, errors="coerce")
        hi = pd.to_numeric(hi_raw, errors="coerce")

        if pd.isna(lo) or pd.isna(hi):
            continue  # skip if either end is NaN

        lo, hi = int(lo), int(hi)  # convert to integers
        if lo > hi:
            lo, hi = hi, lo
        mask = (df_subset["Start"] <= hi) & (df_subset["Stop"] >= lo)
        if not mask.any():        # skip if this strand doesn’t hit anything
            continue

        df_tag = df_subset.loc[mask].copy()

        df_tag["StrandTag"]    = tag
        df_tag["OverlapStart"] = np.maximum(df_tag["Start"], lo)
        df_tag["OverlapStop"]  = np.minimum(df_tag["Stop"],  hi)
        df_tag["OverlapLen"]   = (
            df_tag["OverlapStop"] - df_tag["OverlapStart"] +
            (1 if inclusive else 0)
        )
        # df_tag["ENSG"] = key["ENSG"]     # keep gene ID for later grouping
        detail_frames.append(df_tag)

    detail_df = (
        pd.concat(detail_frames, ignore_index=True)
        if detail_frames else pd.DataFrame()
    )
    # ---------- new guard: bail out if nothing to work with -------------
    if detail_df.empty:
        return Counter(), detail_df if return_detail else Counter()
    # --------------------------------------------------------------------
    protein_counts = Counter(detail_df["Protein"])

    if return_detail:
        cols = ["Protein", "StrandTag",
                "Start", "Stop",
                "OverlapStart", "OverlapStop", "OverlapLen"]
        return protein_counts, detail_df[cols]
    else:
        return protein_counts

# File paths
excel_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\gap_triplex_genomic_mappings.xlsx"
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"

# Read the Excel file, which we expect to have columns: 'ENST', 'Chromosome', 'H_genomic_start', 'H-genomic_end', 'C_genomic_start', 'C-genomic_end','W_genomic_start', 'W-genomic_end','Strand'.
df = pd.read_excel(excel_fp)
'''
# Adjust the Chromosome column to include the 'chr' prefix if missing.
df['Chromosome'] = df['Chromosome'].astype(str).apply(lambda x: x if x.lower().startswith('chr') else "chr" + x)
'''
# Build a dictionary of gene keys.
gene_dict = {}

for idx, row in df.iterrows():
    gene_dict[row["ENST"]] = {
        "ENST": row["ENST"],
        "Chromosome": row["Chromosome"],
        "Hstart": row["H_genomic_start"],  # Map 'H Start Index' in Excel to 'Start'
        "Cstart": row["C_genomic_start"],  # Map 'C Start Index' in Excel to 'Start'    
        "Wstart": row["W_genomic_start"],  # Map 'W Start Index' in Excel to 'Start'
        "Hstop": row["H_genomic_end"],      # Map 'End' in Excel to 'Stop'
        "Cstop": row["C_genomic_end"],      # Map 'End' in Excel to 'Stop'
        "Wstop": row["W_genomic_end"],      # Map 'End' in Excel to 'Stop'
        "Strand": row["Strand"]   # Map 'Strand' (Excel) to 'Strand'
    }

# print("Gene dictionary for processing:")
# print(gene_dict)

# Set up a cache dictionary to hold HDF5 subsets keyed by (Chromosome, Strand).
cache = {}

# Global protein counter to sum all protein counts across every gene.
global_protein_counter = Counter()

all_details = []
# Process each gene and update the global counter.
for transcript_id, key in gene_dict.items():
    protein_counts, detail_df = query_data_cached(hdf5_fp, key, cache, return_detail=True, inclusive=True)
    global_protein_counter.update(protein_counts)
    
    if detail_df.empty:
        continue
    detail_df["ENST"] = transcript_id         # <-- add it right here
    all_details.append(detail_df)  # Collect all detail DataFrames
overlaps = pd.concat(all_details, ignore_index=True)
    # print(f"\nResults for gene {transcript_id}:")
    # print(protein_counts)
bases_by_protein = overlaps.groupby("Protein")["OverlapLen"].sum()
bases_by_protein.to_csv("bases_by_protein_gap_triplex.csv")
print(bases_by_protein)
# After processing, print out the overall global protein counts.
print("\nGlobal protein counter (summed over all queries):")
print(global_protein_counter)
time_stop = time.time()
print("Total time taken:", time_stop - time_start, "seconds")
print("Total number of genes processed:", len(gene_dict))
print("Total number of proteins found:", len(global_protein_counter))
print("Total number of unique proteins found:", len(set(global_protein_counter.elements())))

