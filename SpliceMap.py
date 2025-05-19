#!/usr/bin/env python3
"""
Script: map_triplex_to_genome.py

Reads a triplex motif file (Excel) with cDNA coordinates for Hoogsteen, Crick, and Watson strands,
merges with transcript metadata (ENST IDs, display names), and maps each
cDNA interval back to genomic coordinates via the Ensembl REST API.

Outputs a CSV with original motifs and their genomic mappings.
"""
import pandas as pd
import requests
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# === CONFIGURATION ===
TRIPLEX_FILE  = r"C:\Users\caleb\OneDrive\Brown\ProjectG\with gap lncRNA.xlsx"
METADATA_FILE = r"C:\Users\caleb\OneDrive\Brown\ProjectG\lncRNA_transcript_lookup.xlsx"
OUTPUT_FILE   = "triplex_genomic_mappings.csv"

# === Read input files ===
df_trip = pd.read_excel(TRIPLEX_FILE)
df_meta = pd.read_excel(METADATA_FILE)

# Merge on ENST
df = df_trip.merge(df_meta[['ENST', 'display_name']], on='ENST', how='left')
missing = df['ENST'].isna().sum()
if missing > 0:
    print(f"⚠️ Warning: {missing} triplex rows have no matching ENST ID in metadata.")

motifs = [
    ('H', 'H Start Index', 'H End Index'),
    ('C', 'C Start Index', 'C End Index'),
    ('W', 'W Start Index', 'W End Index'),
]
for prefix, _, _ in motifs:
    df[f"{prefix}_chrom"] = None
    df[f"{prefix}_genomic_start"] = None
    df[f"{prefix}_genomic_end"] = None
    df[f"{prefix}_strand"] = None

SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}
session = requests.Session()
session.headers.update(HEADERS)

lock = threading.Lock()
motif_cache = {}
tasks = []

def map_motif(row_idx, enst, prefix, cdna_start, cdna_end):
    motif_key = (enst, cdna_start, cdna_end)
    with lock:
        if motif_key in motif_cache:
            mappings = motif_cache[motif_key]
        else:
            uri = f"/map/cdna/{enst}/{cdna_start}..{cdna_end}"
            resp = session.get(SERVER + uri)
            if not resp.ok:
                print(f"❌ Error mapping {enst} {prefix} {cdna_start}-{cdna_end}: {resp.status_code}")
                return
            mappings = resp.json().get("mappings", [])
            motif_cache[motif_key] = mappings
    if not mappings:
        print(f"⚠️ No mapping returned for {enst} {prefix} {cdna_start}-{cdna_end}")
        return

    block = mappings[0]
    chrom  = block['seq_region_name']
    gstart = block['start']
    gend   = block['end']
    strand = block['strand']
    gstr   = '+' if strand == 1 else '-'

    with lock:
        df.at[row_idx, f"{prefix}_chrom"]          = f"chr{chrom}"
        df.at[row_idx, f"{prefix}_genomic_start"] = gstart
        df.at[row_idx, f"{prefix}_genomic_end"]   = gend
        df.at[row_idx, f"{prefix}_strand"]        = gstr

# === Parallel mapping ===
overall_start = time.time()
with ThreadPoolExecutor(max_workers=10) as executor:
    for idx, row in df.iterrows():
        enst = row['ENST']
        if pd.isna(enst):
            continue
        for prefix, start_col, end_col in motifs:
            cdna_start = row.get(start_col)
            cdna_end   = row.get(end_col)
            if pd.isna(cdna_start) or pd.isna(cdna_end):
                continue
            cdna_start = int(cdna_start)
            cdna_end   = int(cdna_end)
            tasks.append(executor.submit(map_motif, idx, enst, prefix, cdna_start, cdna_end))

    for i, future in enumerate(as_completed(tasks), start=1):
        if i % 100 == 0:
            elapsed = time.time() - overall_start
            print(f"⏱ Completed {i} motif mappings in {elapsed:.1f}s")

# Save
df.to_csv(OUTPUT_FILE, index=False)
print(f"✅ Triplex-to-genome mapping saved to {OUTPUT_FILE} ({len(df)} rows) | Total time: {time.time() - overall_start:.1f}s")
