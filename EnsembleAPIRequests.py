import requests
import json
import pandas as pd

# Get data from excel file
df = pd.read_excel("with gap lncRNA.xlsx")
ensg_ids = df['ENSG'].tolist()  # Convert the column to a list

# Function to split a list into chunks
def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

# API details
url = "http://rest.ensembl.org/lookup/id"
headers = {"Content-Type": "application/json"}

# Initialize lists to store results
Chromosome = []
Start = []
End = []
Strand = []
Processed_ENSG = []

# Process IDs in chunks of 990 to avoid exceeding the API limit of 1000
for chunk in chunk_list(ensg_ids, 990):
    data = {"ids": chunk}
    response = requests.post(url, headers=headers, json=data)
    
    if response.status_code == 200:
        result = response.json()
        for id in chunk:
            gene_data = result.get(id)
            if gene_data:
                Processed_ENSG.append(id)
                Chromosome.append(gene_data.get("seq_region_name"))
                Start.append(gene_data.get("start"))
                End.append(gene_data.get("end"))
                Strand.append(gene_data.get("strand"))
            else:
                print(f"Data for ID {id} not found in the response.")
    else:
        print(f"Error: {response.status_code}, {response.text}")

# Create a DataFrame from the extracted data
df_result = pd.DataFrame({
    "Transcript": 
    "ENSG": Processed_ENSG,
    "Chromosome": Chromosome,
    "Start": Start,
    "End": End,
    "Strand": Strand
})

# Save the DataFrame to an Excel file
df_result = df_result.drop_duplicates() # Remove duplicates
output_file = "lnc_data_complete_gap2.xlsx"
df_result.to_excel(output_file, index=False)
print(f"Data saved to {output_file}")