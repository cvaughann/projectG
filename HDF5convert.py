### HDF5convert

import pandas as pd
import time
import os

# File paths: adjust these as needed.
text_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.txt"
hdf5_fp = r"C:\Users\caleb\OneDrive\Brown\ProjectG\HumanClipDB\human.h5"

def convert_text_to_hdf5(text_fp, hdf5_fp):
    """Convert the large text file into an HDF5 file for fast partial I/O."""
    start_time = time.time()
    print("Reading text file...")
    # Read only the necessary columns.
    df = pd.read_csv(
        text_fp,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 4, 5],
        names=["Chromosome", "Start", "Stop", "Strand", "Protein"],
        dtype={
            "Chromosome": str,
            "Start": "int32",
            "Stop": "int32",
            "Strand": str,
            "Protein": str
        }
    )
    read_time = time.time()
    print(f"Finished reading text file in {read_time - start_time:.2f} seconds.")

    # Save the DataFrame to an HDF5 file using table format.
    print("Saving to HDF5 file...")
    df.to_hdf(hdf5_fp, key="data", mode="w", format="table", data_columns=["Chromosome", "Strand"])
    end_time = time.time()
    print(f"Saved HDF5 file in {end_time - read_time:.2f} seconds.")

if __name__ == "__main__":
    if os.path.exists(hdf5_fp):
        print("HDF5 file already exists. Delete it if you want to reconvert.")
    else:
        convert_text_to_hdf5(text_fp, hdf5_fp)
