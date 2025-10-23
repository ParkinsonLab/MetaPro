import os
import sys
from datetime import datetime as dt

def parse_fastq_streaming(filepath):
    """
    Generator that yields FASTQ records one at a time without loading entire file into memory.
    Returns (id, full_record) tuples.
    Fixed version that handles incomplete records at file end and preserves FASTQ formatting.
    """
    with open(filepath, 'r') as f:
        while True:
            # Read all 4 lines (includes trailing '\n')
            header = f.readline()
            if not header:
                break
            
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            
            # Check if we have a complete record
            if not seq or not plus or not qual:
                print("Warning: Incomplete FASTQ record at end of file %s" % filepath)
                print("Header: %r" % header)
                print("Seq: %r" % seq)
                print("Plus: %r" % plus)
                print("Qual: %r" % qual)
                break
            
            # --- Create stripped versions for validation and ID extraction ONLY ---
            stripped_header = header.strip()
            stripped_seq = seq.strip()
            stripped_plus = plus.strip()
            stripped_qual = qual.strip()
            
            # Additional validation
            if not stripped_header.startswith('@'):
                print("Warning: Invalid header line: %s" % stripped_header)
                continue
            if not stripped_plus.startswith('+'):
                print("Warning: Invalid plus line: %s" % stripped_plus)
                continue
            if len(stripped_seq) != len(stripped_qual):
                print("Warning: Sequence and quality lengths don't match: %d vs %d" % (len(stripped_seq), len(stripped_qual)))
                continue
            
            # Extract ID (everything before first space)
            read_id = stripped_header.split()[0] if ' ' in stripped_header else stripped_header
            
            # --- Reconstruct the full record using the original, un-stripped lines ---
            # The original lines already contain the necessary '\n' for correct FASTQ format
            full_record = header + seq + plus + qual
            
            yield (read_id, full_record)

def get_read_ids(filepath):
    """
    First pass: collect all read IDs from a FASTQ file.
    Memory efficient - only stores IDs, not full records.
    """
    ids = set()
    print("Scanning IDs from %s" % filepath)
    
    count = 0
    for read_id, _ in parse_fastq_streaming(filepath):
        ids.add(read_id)
        count += 1
        if count % 1000000 == 0:  # Progress indicator for large files
            print("  Processed %d reads..." % count)
    
    print("Found %d unique IDs in %s" % (len(ids), filepath))
    return ids

def filter_for_orphans(p0_path_i, p1_path_i, orphans_path_i, p0_path_o, p1_path_o, orphans_path_o):
    """
    Memory-efficient version that processes files in streaming fashion.
    Uses two passes to minimize memory usage with large files.
    """
    
    # Handle existing orphans first (FIXED: streams content to be memory efficient)
    if os.path.exists(orphans_path_i) and os.path.getsize(orphans_path_i) > 0:
        size = os.path.getsize(orphans_path_i)
        print("%s Found existing orphans file with %d bytes" % (dt.today(), size))
        # Stream content instead of reading all to memory
        with open(orphans_path_i, 'r') as infile, open(orphans_path_o, 'w') as outfile:
            for line in infile:
                outfile.write(line)
        print("Copied existing orphans from %s to %s" % (orphans_path_i, orphans_path_o))
    else:
        if os.path.exists(orphans_path_i):
            print("%s empty singletons file exists" % dt.today())
        else:
            print("%s doesn't exist, starting with empty orphans file" % orphans_path_i)
        open(orphans_path_o, 'w').close()
    
    print("Processing %s and %s" % (p0_path_i, p1_path_i))
    
    # Pass 1: Get all read IDs from both files (memory efficient)
    print("Pass 1: Collecting read IDs...")
    start_time = dt.now()
    
    ids_0 = get_read_ids(p0_path_i)
    ids_1 = get_read_ids(p1_path_i)
    
    # Find intersection - reads present in both files
    common_ids = ids_0.intersection(ids_1)
    print("Found %d matching pairs" % len(common_ids))
    print("File 0 orphans: %d" % (len(ids_0) - len(common_ids)))
    print("File 1 orphans: %d" % (len(ids_1) - len(common_ids)))
    print("ID collection time: %s" % (dt.now() - start_time))
    
    # Pass 2: Stream through files again and write output based on ID sets
    print("Pass 2: Writing output files...")
    start_time = dt.now()
    
    # Process first file
    matched_count_0 = 0
    orphan_count_0 = 0
    
    print("Processing file 0...")
    with open(p0_path_o, 'w') as matched_out, open(orphans_path_o, 'a') as orphan_out:
        for read_id, full_record in parse_fastq_streaming(p0_path_i):
            if read_id in common_ids:
                matched_out.write(full_record)
                matched_count_0 += 1
            else:
                orphan_out.write(full_record)
                orphan_count_0 += 1
            
            if (matched_count_0 + orphan_count_0) % 1000000 == 0:
                print("  Processed %d reads from file 0..." % (matched_count_0 + orphan_count_0))
    
    # Process second file
    matched_count_1 = 0
    orphan_count_1 = 0
    
    print("Processing file 1...")
    with open(p1_path_o, 'w') as matched_out, open(orphans_path_o, 'a') as orphan_out:
        for read_id, full_record in parse_fastq_streaming(p1_path_i):
            if read_id in common_ids:
                matched_out.write(full_record)
                matched_count_1 += 1
            else:
                orphan_out.write(full_record)
                orphan_count_1 += 1
            
            if (matched_count_1 + orphan_count_1) % 1000000 == 0:
                print("  Processed %d reads from file 1..." % (matched_count_1 + orphan_count_1))
    
    print("File processing time: %s" % (dt.now() - start_time))
    print("Saved %d matching pairs to %s" % (matched_count_0, p0_path_o))
    print("Saved %d matching pairs to %s" % (matched_count_1, p1_path_o))
    print("Appended %d new orphans to %s" % (orphan_count_0 + orphan_count_1, orphans_path_o))
    
    # Show final orphan file size
    if os.path.exists(orphans_path_o):
        print("Total orphans file size: %d bytes" % os.path.getsize(orphans_path_o))
    
    # Clean up memory
    del ids_0, ids_1, common_ids


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py <p0_in> <p1_in> <orphans_in> <p0_out> <p1_out> <orphans_out>")
        sys.exit(1)
    
    p0_path_i = os.path.abspath(sys.argv[1])
    p1_path_i = os.path.abspath(sys.argv[2])
    orphans_path_i = os.path.abspath(sys.argv[3])
    
    p0_path_o = os.path.abspath(sys.argv[4])
    p1_path_o = os.path.abspath(sys.argv[5])
    orphans_path_o = os.path.abspath(sys.argv[6])
    
    print("Input files:")
    print("  p0 in: %s" % p0_path_i)
    print("  p1 in: %s" % p1_path_i)
    print("  orphans in: %s" % orphans_path_i)
    print("Output files:")
    print("  p0 out: %s" % p0_path_o)
    print("  p1 out: %s" % p1_path_o)
    print("  orphans out: %s" % orphans_path_o)
    
    start_time = dt.now()
    filter_for_orphans(p0_path_i, p1_path_i, orphans_path_i, p0_path_o, p1_path_o, orphans_path_o)
    end_time = dt.now()
    print("Total processing time: %s" % (end_time - start_time))