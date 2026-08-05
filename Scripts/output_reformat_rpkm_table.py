import os
import sys
import pandas as pd

# This code flattens the ECs so that they'll duplicate the row info if there's more than 1 EC per gene
if __name__ == "__main__":

    # Check for correct number of arguments
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_rpkm_file> <output_rpkm_file>")
        sys.exit(1)

    input_rpkm = sys.argv[1]
    output_rpkm = sys.argv[2]

    # --- Performance Improvement 1: Use Pandas for fast I/O ---
    try:
        # Read the entire file into a DataFrame
        # Assuming the input is tab-separated and has a header row
        df = pd.read_csv(input_rpkm, sep="\t", engine="c")
    except Exception as e:
        print(f"Error reading input RPKM file: {e}")
        sys.exit(1)

    # --- Performance Improvement 2: Use highly optimized Pandas string methods ---
    # The EC column (assumed to be the 4th column, index 3, after reading)
    # The column name 'EC' is assumed to be the 4th header name in the file.
    # We must find the EC column name dynamically.
    if len(df.columns) < 4:
        print("Input file must have at least 4 columns (GeneID, ..., EC).")
        sys.exit(1)
        
    # Assuming the 4th column (index 3) is the one with the '|' separated ECs
    ec_col_name = df.columns[3] 

    # 'str.split()' and 'explode()' are the most efficient way to flatten the data.
    # 1. Split the EC column by the '|' delimiter, turning the strings into lists.
    df[ec_col_name] = df[ec_col_name].str.split("|")

    # 2. 'explode' the DataFrame based on the EC list. This duplicates all
    #    other row information for each EC in the list, achieving the desired
    #    "flattening" and duplicating behavior.
    df_flattened = df.explode(ec_col_name)

    # Clean up empty strings that might result from splitting (e.g., "A|B|" -> ["A", "B", ""])
    # This step is optional but makes the output cleaner.
    df_flattened = df_flattened[df_flattened[ec_col_name] != '']
    
    # --- Performance Improvement 3: Fast CSV Writing ---
    # Write the resulting DataFrame back to the output file
    df_flattened.to_csv(output_rpkm, sep="\t", index=False)