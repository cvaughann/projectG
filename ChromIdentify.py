#"C:\Users\caleb\OneDrive\Brown\ProjectG\no gap lncRNA.xlsx"
## With Start and Stop Codons
## Chromosome Identifier

import pandas as pd
import requests
import time

# Load the Excel file
df = pd.read_excel('no gap lncRNA.xlsx')

# Extract the "ENSG" column (update column name if necessary)
ensg_column = df["ENSG"]  # Ensure your Excel file has this column

# Dictionary for caching API results
chromosome_cache = {}

def get_chromosome_ensembl(ensg_id):
    """Fetches the chromosome number, start index, and stop index for a given canonical transcript ID using the Ensembl API."""
    if ensg_id in chromosome_cache:
        return chromosome_cache[ensg_id]  # Return cached value if available

    server = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{ensg_id}?content-type=application/json"

    response = requests.get(server + endpoint, headers={"Content-Type": "application/json"})

    if response.status_code == 200:
        data = response.json()
        chromosome = data.get("seq_region_name")  # Extract chromosome number
        start = data.get("start")  # Extract start index
        stop = data.get("end")  # Extract stop index
        chromosome_cache[ensg_id] = (chromosome, start, stop)  # Store in cache
        return chromosome, start, stop
    else:
        return f"Error: {response.status_code}", None, None

def assign_chromosomes(df):
    processed_count = 0  # Counter for processed rows
    start_time = time.time()

    for index, row in df.iterrows():
        ensg_id = row["ENSG"]
        chromosome, start, stop = get_chromosome_ensembl(ensg_id)
        df.at[index, "Chromosome"] = chromosome
        df.at[index, "Gene Start"] = start
        df.at[index, "Gene Stop"] = stop
        processed_count += 1

        # Display elapsed time every 100 rows processed
        if processed_count % 100 == 0:
            elapsed_time = time.time() - start_time
            print(f"Processed {processed_count} rows in {elapsed_time:.2f} seconds...")

# Start timing the execution
start_time = time.time()

# Assign chromosome numbers
assign_chromosomes(df)

# Save the updated Excel file
df.to_excel('lnc_st_ch_nogapfirst2k.xlsx', index=False)

# End timing the execution
end_time = time.time()
execution_time = end_time - start_time

print(f"Chromosome numbers added and saved to 'updated_with_chromosome.xlsx'")
print(f"Total execution time: {execution_time:.2f} seconds")
