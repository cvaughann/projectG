import requests
import pandas as pd

# Read your input Excel (replace filename/path as needed)
df = pd.read_excel(r"C:\Users\caleb\OneDrive\Brown\ProjectG\no gap lncRNA.xlsx")
enst_ids = df['ENST'].dropna().astype(str).tolist()

# Function to split a list into chunks
def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

# Ensembl REST API endpoint
url = "http://rest.ensembl.org/lookup/id"
headers = {"Content-Type": "application/json"}

# Prepare lists to collect results
Processed_ENST   = []
ENSG             = []
DisplayName      = []
Chromosome       = []   # ← new list for chr-prefixed seq_region_name
Start            = []
End              = []
Strand           = []   # ← will hold "+" or "−"

# Batch up to 990 at a time (API limit is 1000 per request)
for chunk in chunk_list(enst_ids, 990):
    payload = {"ids": chunk}
    resp = requests.post(url, headers=headers, json=payload)
    
    if resp.status_code != 200:
        print(f"Error: {resp.status_code} – {resp.text}")
        continue
    
    results = resp.json()
    for tid in chunk:
        tr = results.get(tid)
        if tr:
            Processed_ENST.append(tid)
            ENSG.append(tr.get("Parent"))
            DisplayName.append(tr.get("display_name"))
            # prefix "chr"
            Chromosome.append("chr" + str(tr.get("seq_region_name")))
            Start.append(tr.get("start"))
            End.append(tr.get("end"))
            # map strand value to symbol
            raw_strand = tr.get("strand")
            if raw_strand == 1:
                Strand.append("+")
            elif raw_strand == -1:
                Strand.append("−")
            else:
                Strand.append(str(raw_strand))
        else:
            print(f"⚠️ No data for {tid}")

# Build output DataFrame, now including chromosome & formatted strand
df_out = pd.DataFrame({
    "ENST":           Processed_ENST,
    "ENSG":           ENSG,
    "display_name":   DisplayName,
    "chromosome":     Chromosome,
    "start":          Start,
    "end":            End,
    "strand":         Strand
})

# Write to Excel
output_file = r"C:\Users\caleb\OneDrive\Brown\ProjectG\lncRNA_transcript_lookup.xlsx"
df_out.to_excel(output_file, index=False)

print(f"✅ Saved {len(df_out)} records to {output_file}")
