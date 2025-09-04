import os
import sys
from datetime import datetime as dt

def parse_fastq_streaming(filepath):
    """
    Generator that yields FASTQ records one at a time without loading entire file into memory.
    Returns (id, full_record) tuples.
    Fixed version that handles incomplete records at file end.
    """
    with open(filepath, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            
            # Check if we have a complete record (all 4 lines present and non-empty)
            if not seq or not plus or not qual:
                print("Warning: Incomplete FASTQ record at end of file", filepath)
                print("Header:", repr(header))
                print("Seq:", repr(seq))
                print("Plus:", repr(plus))
                print("Qual:", repr(qual))
                break
            
            header = header.strip()
            seq = seq.strip()
            plus = plus.strip()
            qual = qual.strip()
            
            # Additional validation
            if not header.startswith('@'):
                print("Warning: Invalid header line:", header)
                continue
            if not plus.startswith('+'):
                print("Warning: Invalid plus line:", plus)
                continue
            if len(seq) != len(qual):
                print("Warning: Sequence and quality lengths don't match:", len(seq), "vs", len(qual))
                continue
            
            # Extract ID (everything before first space)
            read_id = header.split()[0] if ' ' in header else header
            # Fixed: Ensure proper FASTQ format - exactly 4 lines with proper newlines
            full_record = header + '\n' + seq + '\n' + plus + '\n' + qual + '\n'
            
            yield (read_id, full_record)

def get_read_ids(filepath):
    """
    First pass: collect all read IDs from a FASTQ file.
    Memory efficient - only stores IDs, not full records.
    """
    ids = set()
    print("Scanning IDs from", filepath)
    
    count = 0
    for read_id, _ in parse_fastq_streaming(filepath):
        ids.add(read_id)
        count += 1
        if count % 1000000 == 0:  # Progress indicator for large files
            print("  Processed", count, "reads...")
    
    print("Found", len(ids), "unique IDs in", filepath)
    return ids

def filter_for_orphans(p0_path_i, p1_path_i, orphans_path_i, p0_path_o, p1_path_o, orphans_path_o):
    """
    Memory-efficient version that processes files in streaming fashion.
    Uses two passes to minimize memory usage with large files.
    """
    
    # Handle existing orphans first - FIXED VERSION
    if os.path.exists(orphans_path_i) and os.path.getsize(orphans_path_i) > 0:
        print(dt.today(), "Found existing orphans file with", os.path.getsize(orphans_path_i), "bytes")
        
        # Read existing orphans and ensure proper formatting
        with open(orphans_path_i, 'r') as infile:
            content = infile.read()
            
            # Remove any trailing whitespace and ensure single final newline
            content = content.rstrip()
            if content:
                content += '\n'
        
        # Write cleaned content
        with open(orphans_path_o, 'w') as outfile:
            outfile.write(content)
            
        print("Copied and cleaned existing orphans from", orphans_path_i, "to", orphans_path_o)
        
        # Validate the copied file has correct line count
        with open(orphans_path_o, 'r') as f:
            line_count = sum(1 for _ in f)
        
        if line_count % 4 != 0:
            print("WARNING: Existing orphans file has", line_count, "lines, not divisible by 4!")
            print("Attempting to fix by removing incomplete records...")
            
            # Fix by removing incomplete records at the end
            lines_to_keep = (line_count // 4) * 4
            with open(orphans_path_o, 'r') as f:
                lines = f.readlines()
            
            with open(orphans_path_o, 'w') as f:
                f.writelines(lines[:lines_to_keep])
            
            print("Fixed: kept", lines_to_keep, "lines (removed", line_count - lines_to_keep, "incomplete lines)")
            
    else:
        if os.path.exists(orphans_path_i):
            print(dt.today(), "empty singletons file exists")
        else:
            print(orphans_path_i, "doesn't exist, starting with empty orphans file")
        open(orphans_path_o, 'w').close()
    
    print("Processing", p0_path_i, "and", p1_path_i)
    
    # Pass 1: Get all read IDs from both files (memory efficient)
    print("Pass 1: Collecting read IDs...")
    start_time = dt.now()
    
    ids_0 = get_read_ids(p0_path_i)
    ids_1 = get_read_ids(p1_path_i)
    
    # Find intersection - reads present in both files
    common_ids = ids_0.intersection(ids_1)
    print("Found", len(common_ids), "matching pairs")
    print("File 0 orphans:", len(ids_0) - len(common_ids))
    print("File 1 orphans:", len(ids_1) - len(common_ids))
    print("ID collection time:", dt.now() - start_time)
    
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
                print("  Processed", matched_count_0 + orphan_count_0, "reads from file 0...")
    
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
                print("  Processed", matched_count_1 + orphan_count_1, "reads from file 1...")
    
    print("File processing time:", dt.now() - start_time)
    print("Saved", matched_count_0, "matching pairs to", p0_path_o)
    print("Saved", matched_count_1, "matching pairs to", p1_path_o)
    print("Appended", orphan_count_0 + orphan_count_1, "new orphans to", orphans_path_o)
    
    # FINAL VALIDATION: Check that output orphans file is properly formatted
    if os.path.exists(orphans_path_o):
        with open(orphans_path_o, 'r') as f:
            final_line_count = sum(1 for _ in f)
        
        print("Total orphans file size:", os.path.getsize(orphans_path_o), "bytes")
        print("Total orphans file lines:", final_line_count)
        
        if final_line_count % 4 != 0:
            print("ERROR: Final orphans file has", final_line_count, "lines, not divisible by 4!")
            print("This indicates a bug - the file is malformed.")
        else:
            print("SUCCESS: Orphans file is properly formatted (", final_line_count//4, "complete FASTQ records)")
    
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
    print("  p0 in:", p0_path_i)
    print("  p1 in:", p1_path_i)
    print("  orphans in:", orphans_path_i)
    print("Output files:")
    print("  p0 out:", p0_path_o)
    print("  p1 out:", p1_path_o)
    print("  orphans out:", orphans_path_o)
    
    start_time = dt.now()
    filter_for_orphans(p0_path_i, p1_path_i, orphans_path_i, p0_path_o, p1_path_o, orphans_path_o)
    end_time = dt.now()
    print("Total processing time:", end_time - start_time)