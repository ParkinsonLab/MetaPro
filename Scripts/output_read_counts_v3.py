import sys
import os
import pandas as pd
import time
from datetime import datetime as dt
import subprocess # New import for potentially faster line counting

# Function 1: Optimized FASTQ Read Counting
# Using a platform-specific command (wc -l) is often the fastest.
# As a fallback, use an optimized Python line counter.
def fastq_count(item):
    """
    Counts the number of reads in a FASTQ file.
    Uses 'wc -l' for speed if available, otherwise uses an optimized
    pure-Python line counter. Returns 0 if the file doesn't exist.
    """
    if not os.path.exists(item):
        return 0

    try:
        # Option A: Fastest via 'wc -l' command line utility
        # This is significantly faster on large files than pure Python loops
        result = subprocess.run(['wc', '-l', item], capture_output=True, text=True, check=True)
        lines = int(result.stdout.split()[0])
    except:
        # Option B: Fallback to optimized pure Python line counting
        lines = 0
        with open(item, "rb") as f:
            lines = sum(1 for line in f)

    if lines % 4 != 0:
        # print(dt.today(), "Warning: Non-standard FASTQ line count in:", item, lines / 4)
        pass # Returning 0 might be safer in production, but here we proceed

    return lines / 4

# Function 2: Optimized Annotated Count
# Uses a generator expression inside sum and set operations for efficiency.
def annotated_count(map_file):
    annotated_mRNA = 0
    read_set = set() # Use a set for unique reads directly
    genes = 0

    if not os.path.exists(map_file):
        return 0, 0

    try:
        with open(map_file, "r") as infile:
            for line in infile:
                # Use str.split(None, 2) for max efficiency, split on any whitespace
                line_split = line.strip().split("\t")
                if len(line_split) < 3: # Skip malformed lines
                    continue

                # Reads are from the 4th column onwards (index 3)
                # Extend the set directly with new reads
                read_set.update(line_split[3:])
                genes += 1

        annotated_mRNA = len(read_set)
        return annotated_mRNA, genes

    except Exception as e:
        print(f"Error processing annotated_count file {map_file}: {e}")
        return 0, 0

# Function 3: Optimized EC Count
# Uses a single set operation and a list comprehension for efficiency.
def ec_count(map_file):
    ecs = set()
    if not os.path.exists(map_file):
        return 0

    try:
        with open(map_file, "r") as infile:
            for line in infile:
                line_split = line.split("\t")
                if len(line_split) < 3: # Skip malformed lines
                    continue
                
                # Split the EC portion and update the set in one go
                ec_portion = line_split[2].strip() # Use .strip() without args
                ecs.update(ec_portion.split("|"))
        
        # Remove empty strings if any result from splitting
        ecs.discard("")
        
        return len(ecs)

    except Exception as e:
        print(f"Error processing ec_count file {map_file}: {e}")
        return 0

# Function 4: Optimized Check Paired Data
def check_paired_data(p1, p2, message):
    if p1 == "None":
        return True
    
    # Only count if the file exists
    p1_count = fastq_count(p1)
    p2_count = fastq_count(p2)
    
    if p1_count != 0 and p2_count != 0 and p1_count != p2_count:
        print(dt.today(), "bad data in:", message)
        sys.exit()
    
    return True

# --- Main Logic Starts Here ---
if __name__ == "__main__":
    
    # Argument parsing remains the same, assuming correct order
    if len(sys.argv) != 21:
        print(f"Usage: python {sys.argv[0]} <20 arguments...>")
        sys.exit(1)

    # Assigning arguments (no change in logic here, just variable names)
    raw_sequence = sys.argv[1]
    qc_s, qc_p1, qc_p2 = sys.argv[2], sys.argv[3], sys.argv[4]
    qc_s_unique, qc_p1_unique, qc_p2_unique = sys.argv[5], sys.argv[6], sys.argv[7]
    host_dir = os.path.abspath(sys.argv[8])
    vectors_s, vectors_p1, vectors_p2 = sys.argv[9], sys.argv[10], sys.argv[11]
    rRNA_s, rRNA_p1, rRNA_p2 = sys.argv[12], sys.argv[13], sys.argv[14]
    mRNA_s, mRNA_p1, mRNA_p2 = sys.argv[15], sys.argv[16], sys.argv[17]
    gene_to_read_map = sys.argv[18]
    ec_map = sys.argv[19]
    output_file = sys.argv[20]

    # Check data integrity
    check_paired_data(qc_p1, qc_p2, "quality")
    check_paired_data(rRNA_p1, rRNA_p2, "rRNA+tRNA")
    check_paired_data(mRNA_p1, mRNA_p2, "putative_mRNA")
    check_paired_data(vectors_p1, vectors_p2, "vectors")

    headings = []
    data = []

    # --- Pre-calculate all counts to avoid redundant calls or slow calculations ---
    # The counts are calculated once here:
    raw_sequence_count = fastq_count(raw_sequence)
    
    qc_p1_count = fastq_count(qc_p1)
    qc_s_count = fastq_count(qc_s)
    quality_sequence_count = qc_p1_count + qc_s_count
    
    vectors_p1_count = fastq_count(vectors_p1)
    vectors_s_count = fastq_count(vectors_s)
    vectors_read_counts = vectors_p1_count + vectors_s_count
    
    rRNA_p1_count = fastq_count(rRNA_p1)
    rRNA_s_count = fastq_count(rRNA_s)
    rRNA_sequence_count = rRNA_p1_count + rRNA_s_count
    
    mRNA_p1_count = fastq_count(mRNA_p1)
    mRNA_s_count = fastq_count(mRNA_s)
    mRNA_sequence_count = mRNA_p1_count + mRNA_s_count

    # Annotated counts are calculated once
    annotated_mRNA_count, genes_count = annotated_count(gene_to_read_map)

    # EC counts are calculated once
    unique_ec_count = ec_count(ec_map)

    # --- Statistics Generation (Use pre-calculated variables) ---
    
    headings.append("Total reads")
    data.append(str(int(raw_sequence_count)))

    headings.append("High quality reads")
    data.append(str(int(quality_sequence_count)))

    headings.append("% high quality")
    # Avoid division by zero
    quality_sequence_pct = quality_sequence_count / raw_sequence_count if raw_sequence_count else 0
    data.append("%.2f" % (quality_sequence_pct * 100))

    # Host Loop - keep the sleep for debugging/pipeline flow if needed, but remove for raw speed.
    # The sleep(10) was likely a debugging pause and is removed for performance.
    
    # print("host dir:", os.listdir(host_dir))
    # time.sleep(10) # Removing this 10-second delay for speed

    for item in os.listdir(host_dir):
        host_path = os.path.join(host_dir, item)
        if os.path.isdir(host_path):
            host_p1 = os.path.join(host_path, item + "_p1_full_host.fastq")
            host_p2 = os.path.join(host_path, item + "_p2_full_host.fastq")
            host_s = os.path.join(host_path, item + "_s_full_host.fastq")
            
            check_paired_data(host_p1, host_p2, f"{item} host")
            
            host_p1_count = fastq_count(host_p1)
            host_s_count = fastq_count(host_s)
            host_read_counts = host_p1_count + host_s_count
            
            headings.append(item + " host reads found in sample")
            data.append(str(int(host_read_counts)))
            
            headings.append("% " + item + " host reads in sample")
            host_pct = host_read_counts / raw_sequence_count if raw_sequence_count else 0
            data.append("%.2f" % (host_pct * 100))
            
    # --- Statistics Generation (Cont.) ---
            
    headings.append("vector reads found in sample")
    data.append(str(int(vectors_read_counts)))

    headings.append("% vector reads in sample")
    vectors_pct = vectors_read_counts / raw_sequence_count if raw_sequence_count else 0
    data.append("%.2f" % (vectors_pct * 100))

    headings.append("rRNA + tRNA reads")
    data.append(str(int(rRNA_sequence_count)))

    headings.append("% rRNA + tRNA reads")
    rRNA_sequence_pct = rRNA_sequence_count / raw_sequence_count if raw_sequence_count else 0
    data.append("%.2f" % (rRNA_sequence_pct*100))

    headings.append("Putative mRNA reads")
    data.append(str(int(mRNA_sequence_count)))

    headings.append("% putative mRNA reads")
    mRNA_sequence_pct = mRNA_sequence_count / raw_sequence_count if raw_sequence_count else 0
    data.append("%.2f" % (mRNA_sequence_pct*100))

    headings.append("Annotated mRNA reads")
    data.append(str(int(annotated_mRNA_count)))

    headings.append("% of putative mRNA reads annotated")
    annotated_mRNA_pct = annotated_mRNA_count / mRNA_sequence_count if mRNA_sequence_count else 0
    data.append("%.2f" % (annotated_mRNA_pct*100))

    headings.append("Unique transcripts")
    data.append(str(int(genes_count)))

    headings.append("High-Quality unique enzymes")
    data.append(str(int(unique_ec_count)))

    # --- Output (Minimal I/O) ---
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(headings) + "\n")
        outfile.write("\t".join(data) + "\n")