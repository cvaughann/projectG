# RBPMatch
import pandas as pd
# def RBPMatch()
df = pd.read_excel(r"C:\Users\caleb\OneDrive\Brown\ProjectG\Good Excel Files for CLIP comparison\lnc_data_complete_gap2.xlsx")
for _, row in df.iterrows():
    ENSG_data = {
        "ENSG": row["ENSG"],
        "Chromosome": row["Chromosome"],
        "Start": row["Start"],
        "Stop": row["End"],
        "Strand": row["Strand2"]
    }
print(ENSG_data["Stop"][1])