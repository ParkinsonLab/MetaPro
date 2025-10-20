#!/usr/bin/env python
#This file is just supposed to repopulate the mRNA with the right number of duplicates
#but in order to do so, it must find out what's mRNA.  
#There's better ways to do this, but we're told that it's such a minor thing, we should leave it be.
#It also sorta ruins the modularity aspect, since we now have to include more than the original 3 files, but.... yeah...
#We need a clever solution for this.
#Also, this needs to check for file existence.  If there's missing files in the arg, don't run
#also this thing is totally broken
import sys
import pandas as pd
import numpy as np
import os
from itertools import islice

# --- Performance Improvement 1: Optimized FASTQ Parsing ---
def fastq_to_dataframe(filename):
    """Reads a FASTQ file and returns a DataFrame of records (ID, seq, junk, quality)."""
    # Use standard Python I/O and itertools for fast reading
    records = []
    
    try:
        with open(filename, 'r') as f:
            # Use islice to read lines in chunks of 4 (FASTQ record)
            while True:
                # Read 4 lines at a time
                lines = list(islice(f, 4))
                if not lines:
                    break
                
                # Clean ID line: get the first word and strip whitespace
                id_line = lines[0].strip()
                # Split only once to get the first part (the ID, e.g., @seq1)
                record_id = id_line.split(" ")[0]
                
                # Check for malformed record (e.g., incomplete FASTQ)
                if len(lines) != 4:
                     # This should not happen if the file is a clean multiple of 4 lines
                    raise ValueError(f"FASTQ file {filename} is truncated or malformed.")
                    
                records.append([
                    record_id,
                    lines[1].strip(),
                    lines[2].strip(),
                    lines[3].strip()
                ])
                
        # Create DataFrame from the list of records
        return pd.DataFrame(records, columns=["ID", "seq", "junk", "quality"])
    except FileNotFoundError:
        print(f"Error: Input file not found: {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading or parsing FASTQ file {filename}: {e}")
        sys.exit(1)


def repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename):
    
    # Use the faster parsing function
    ref_df = fastq_to_dataframe(ref_filename)
    mRNA_df = fastq_to_dataframe(mRNA_filename)

    # --- Performance Improvement 2: Optimized Clustering File Parsing ---
    cluster_map = {}
    with open(cluster_filename, "r") as clustr_read:
        rep = ""
        for line in clustr_read:
            line = line.strip()
            if not line or line.startswith(">"):
                continue

            # Optimized ID extraction based on the known format: ID is between '>' and '...'
            start_index = line.find(">") + 1
            end_index = line.find("...")
            
            # Use 'try/except' block to avoid unexpected crashes on malformed lines
            try:
                if start_index > 0 and end_index != -1 and end_index > start_index:
                    seq_id = line[start_index:end_index]
                    
                    if line.startswith("0"):
                        # Representative sequence: '0\t100%\t>seq1...'
                        rep = seq_id
                        cluster_map[rep] = [seq_id]
                    elif rep:
                        # Clustered sequence: '\t\t>seq2...'
                        cluster_map.setdefault(rep, []).append(seq_id)
            except:
                # Silently skip lines that fail the parsing logic
                continue

    # --- Optimized Repopulation Logic ---
    reduplicated_ids = set()

    # The mRNA_df IDs have a leading '@' from the FASTQ format (e.g., '@seq1').
    # The cluster map IDs do *not* have the leading '@'.
    # We use a vectorized operation to remove the '@' for matching.
    mRNA_ids_cleaned = mRNA_df["ID"].str[1:].tolist()

    for sequence in mRNA_ids_cleaned:
        if sequence in cluster_map:
            # If the sequence is a representative sequence, add all original IDs (with '@')
            # The original logic includes the representative itself if it has other members, 
            # and only itself if it was a singleton, which is covered by iterating the cluster_map list.
            for seq_id in cluster_map[sequence]:
                reduplicated_ids.add("@" + seq_id)
        else:
            # The sequence from mRNA_df was not a cluster representative, 
            # add it back as a unique read that passed filtering.
            reduplicated_ids.add("@" + sequence)

    # --- Faster Export ---
    # Filter the reference DataFrame using the set of desired IDs (fast with isin)
    data_df = ref_df[ref_df.ID.isin(reduplicated_ids)]
    
    # Convert the four columns into a single flat list of strings, which is the FASTQ format (4 lines per record)
    fastq_lines = data_df[['ID', 'seq', 'junk', 'quality']].values.flatten().tolist()

    # Write the lines to the file using efficient string joining
    with open(output_filename, "w") as out_file:
        # Use '\n'.join() to write the entire file content in one operation
        out_file.write('\n'.join(fastq_lines) + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} <ref_file> <mRNA_file> <cluster_file> <output_file>")
        sys.exit(1)
        
    ref_filename = sys.argv[1]      #in: the file that the duplicates will come from
    mRNA_filename = sys.argv[2]     #in: the file that will contain the reads needing to be duplicated
    cluster_filename = sys.argv[3]  #in: the cluster file showing what was in-fact duplicated
    output_filename = sys.argv[4]   #out: the final output
    
    # Check for file existence/size before running the main logic
    if not os.path.exists(mRNA_filename):
        print(f"Input mRNA file not found: {mRNA_filename}")
        sys.exit(1)

    if (os.path.getsize(mRNA_filename) == 0):
        with open(output_filename, "w") as out_file:
            print("input file empty. skipping")
    else:
        # Check for other necessary files
        if not all(os.path.exists(f) for f in [ref_filename, cluster_filename]):
            print("Missing ref_filename or cluster_filename. Aborting.")
            sys.exit(1)
            
        repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename)