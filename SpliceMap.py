#!/usr/bin/env python3
import os
import json
import pandas as pd
import httpx
import asyncio
import time

# === CONFIGURATION ===
TRIPLEX_FILE  = r"C:\Users\caleb\OneDrive\Brown\ProjectG\no gap lncRNA.xlsx"
METADATA_FILE = r"C:\Users\caleb\OneDrive\Brown\ProjectG\nogap_lncRNA_transcript_lookup.xlsx"
OUTPUT_FILE   = "nogap_triplex_genomic_mappings.xlsx"
CACHE_FILE    = "motif_cache.json"

SERVER = "https://rest.ensembl.org"
HEADERS = {
    "Accept": "application/json",
    "User-Agent": "your.email@example.com"
}

# Concurrency / rate-limiting settings
CONCURRENCY_LIMIT = 5       # max simultaneous requests
MAX_RETRIES       = 5       # for exponential backoff
INITIAL_BACKOFF   = 1.0     # seconds

motifs = [
    ('H', 'H Start Index', 'H End Index'),
    ('C', 'C Start Index', 'C End Index'),
    ('W', 'W Start Index', 'W End Index'),
]

# === Read input files ===
print("🔄 Reading input files...")
df_trip = pd.read_excel(TRIPLEX_FILE)
df_meta = pd.read_excel(METADATA_FILE)

df = df_trip.merge(df_meta[['ENST', 'display_name']], on='ENST', how='left')
missing = df['ENST'].isna().sum()
if missing > 0:
    print(f"⚠️ Warning: {missing} triplex rows lack ENST metadata.")

# Initialize output columns
for prefix, _, _ in motifs:
    df[f"{prefix}_chrom"] = None
    df[f"{prefix}_genomic_start"] = None
    df[f"{prefix}_genomic_end"] = None
    df[f"{prefix}_strand"] = None

# === Step 1: Deduplicate unique motifs ===
print("🔍 Deduplicating motifs...")
unique_motifs = set()
index_map = {}
for idx, row in df.iterrows():
    enst = row['ENST']
    if pd.isna(enst): continue
    for prefix, start_col, end_col in motifs:
        start = row.get(start_col)
        end   = row.get(end_col)
        if pd.notna(start) and pd.notna(end):
            key = (enst, int(start), int(end))
            unique_motifs.add(key)
            index_map.setdefault(key, []).append((idx, prefix))

# === Load or initialize cache ===
print("🗄️ Loading cache...")
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'r') as f:
        cache = json.load(f)
else:
    cache = {}

# Separate cached vs to-fetch
print("↔️ Splitting cached vs to-fetch motifs...")
cached_results = []
to_fetch = []
for key in unique_motifs:
    enst, start, end = key
    key_str = f"{enst}:{start}:{end}"
    if key_str in cache:
        cached_results.append((enst, start, end, cache[key_str]))
    else:
        to_fetch.append(key)
print(f"✅ {len(cached_results)} cached; {len(to_fetch)} to fetch.")

# === Step 2: Async mapping with concurrency limit & backoff ===
semaphore = asyncio.Semaphore(CONCURRENCY_LIMIT)

async def async_map_motif(client, enst, start, end):
    uri = f"/map/cdna/{enst}/{start}..{end}"
    backoff = INITIAL_BACKOFF
    for attempt in range(1, MAX_RETRIES + 1):
        async with semaphore:
            r = await client.get(SERVER + uri)
        if r.status_code == 200:
            return (enst, start, end, r.json().get("mappings", []))
        if r.status_code == 429:
            retry_after = int(r.headers.get("Retry-After", backoff))
            print(f"⚠️ 429 rate-limit on {enst}:{start}-{end}; sleeping {retry_after}s (attempt {attempt})...")
            await asyncio.sleep(retry_after)
            backoff *= 2
            continue
        print(f"❌ HTTP {r.status_code} for {enst}:{start}-{end}")
        break
    return (enst, start, end, [])

async def fetch_all(motif_list):
    async with httpx.AsyncClient(headers=HEADERS, limits=httpx.Limits(max_connections=CONCURRENCY_LIMIT)) as client:
        tasks = [async_map_motif(client, *key) for key in motif_list]
        print(f"🚀 Launching {len(tasks)} tasks with concurrency={CONCURRENCY_LIMIT}...")
        return await asyncio.gather(*tasks)

# Fetch missing motifs
fetched_results = []
if to_fetch:
    start = time.time()
    fetched_results = asyncio.run(fetch_all(to_fetch))
    print(f"⏱ Fetched {len(fetched_results)} in {time.time()-start:.1f}s")

# Update cache
def save_cache():
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)
for enst, start, end, mappings in fetched_results:
    cache[f"{enst}:{start}:{end}"] = mappings
save_cache()

# === Step 3: Apply mappings ===
print("📝 Applying mappings to DataFrame...")
all_results = cached_results + fetched_results
for enst, start, end, mappings in all_results:
    if not mappings: continue
    block = mappings[0]
    chrom = f"chr{block['seq_region_name']}"
    gstart, gend = block['start'], block['end']
    strand = '+' if block['strand']==1 else '-'
    for idx, prefix in index_map[(enst, start, end)]:
        df.at[idx, f"{prefix}_chrom"] = chrom
        df.at[idx, f"{prefix}_genomic_start"] = gstart
        df.at[idx, f"{prefix}_genomic_end"] = gend
        df.at[idx, f"{prefix}_strand"] = strand

# === Save output ===
print("💾 Saving results...")
df.to_excel(OUTPUT_FILE, index=False)
print(f"✅ Done. Output -> {OUTPUT_FILE}")