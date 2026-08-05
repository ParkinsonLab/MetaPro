import pandas as pd
import sys
import os

# does a simple groupby and sums up the taxa found
if __name__ == "__main__":
    # Check for correct number of arguments
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_taxa_file> <output_summary_file>")
        sys.exit(1)

    taxa_class_file = sys.argv[1]
    output_file = sys.argv[2]

    # --- Performance Improvement 1: Faster CSV Reading and Direct Selection ---
    # Use 'usecols' to read only the necessary columns (index 2 for 'taxa')
    # and use 'names' to set the column name directly, skipping the intermediate
    # step of reading all columns and then renaming and dropping 'class'/'read'.
    # Note: The original code inferred 'taxa' was the 3rd column (index 2)
    # after assigning names "class", "read", "taxa" to the 3 columns it read.
    # Assuming the input file has three columns and we only need the third one.
    try:
        # Use 'engine="c"' for faster reading, and 'dtype' to specify column type
        taxa_df = pd.read_csv(
            taxa_class_file,
            sep="\t",
            on_bad_lines='skip',
            quoting=3,
            engine="c",
            usecols=[2], # Assuming 'taxa' is the 3rd column (0-indexed)
            names=["taxa"],
            dtype={"taxa": str} # Enforcing string type for the grouping column
        )
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)


    # --- Performance Improvement 2: More Efficient Grouping/Counting ---
    # The 'value_counts()' method is often faster than:
    # 1. creating a new column 'count'=1
    # 2. performing a 'groupby().sum()'
    # It directly counts the occurrences of unique values in the 'taxa' column.
    taxa_counts_series = taxa_df["taxa"].value_counts()

    # Convert the resulting Series back to a DataFrame for writing to CSV
    # 'rename_axis' resets the index name (which is 'taxa') to a column
    # 'reset_names' renames the value column (which is the count)
    taxa_df_summed = taxa_counts_series.rename_axis('taxa').reset_index(name='count')

    # --- Output ---
    taxa_df_summed.to_csv(output_file, sep="\t", index=False)
    # print(taxa_df_summed)