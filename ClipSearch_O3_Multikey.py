import pandas as pd
import sys
import os
from ClipSearch_O3 import query_data  # Use relative import if ClipSearch_O3 is in the same directory
# Add the directory containing ClipSearch_O3 to the Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
sys.path.append(parent_dir)
excel_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_gap2.xlsx"
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"
# Read the Excel file. We expect columns: 'ENSG', 'Chromosome', 'Start', 'End', and 'Strand2'.
df = pd.read_excel(excel_fp)
df['Chromosome'] = df['Chromosome'].astype(str).apply(lambda x: x if x.lower().startswith('chr') else "chr" + x) # Change to match h5

gene_dict = {}
for idx, row in df.head(10).iterrows():
    global_counter = ()
    gene_dict[row["ENSG"]] = {
        "ENSG": row["ENSG"],
        "Chromosome": row["Chromosome"],
        "Start": row["Start"],
        "Stop": row["End"],     # Map 'End' to 'Stop'
        "Strand": row["Strand2"]    # Map 'Strand2' to 'Strand'
    }
print(gene_dict)
for key in gene_dict:
    query_data(hdf5_fp, gene_dict[key])