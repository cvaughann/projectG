#!/usr/bin/env python3
"""
fetch_transcript_loc_example.py

Fetches and prints chromosome, start/end coordinates, and strand
for the hard-coded transcript ID ENST00000445118.7 when run.
"""

import requests
import sys

API_BASE = "https://rest.ensembl.org"

def get_transcript_location(transcript_id: str):
    """
    Query Ensembl REST API /lookup endpoint for the given transcript ID.
    Returns a tuple: (chromosome, start, end, strand)
    Raises requests.HTTPError on HTTP errors.
    """
    url = f"{API_BASE}/lookup/id/{transcript_id}"
    headers = {"Accept": "application/json"}
    resp = requests.get(url, headers=headers)
    resp.raise_for_status()
    data = resp.json()
    return (
        data["seq_region_name"],
        data["start"],
        data["end"],
        data["strand"]
    )

if __name__ == "__main__":
    # --- Example run on ENST00000445118.7 ---
    tid = "ENST00000445118.7"
    try:
        chrom, start, end, strand = get_transcript_location(tid)
    except requests.HTTPError as e:
        print(f"ERROR: Could not fetch '{tid}': {e}", file=sys.stderr)
        sys.exit(1)

    strand_sym = "+" if strand == 1 else "-"
    print(f"Transcript: {tid}")
    print(f"Chromosome: {chrom}")
    print(f"Coordinates: {chrom}:{start}–{end}")
    print(f"Strand: {strand_sym} ({strand})")
