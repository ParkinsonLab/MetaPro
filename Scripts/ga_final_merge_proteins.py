#!/usr/bin/env python3
"""
Merge all .faa files from two directories into one output file.
Remove duplicates based on sequence content.
"""

import os
import sys
from datetime import datetime as dt


def merge_faa_files(dir1, dir2, output_file):
    """
    Merge all .faa files from two directories, removing duplicates.
    """
    print(dt.today(), "Starting merge")
    
    # Find all .faa files
    faa_files = []
    for directory in [dir1, dir2]:
        if os.path.exists(directory):
            for filename in os.listdir(directory):
                if filename.endswith('.faa'):
                    faa_files.append(os.path.join(directory, filename))
    
    print(dt.today(), "Found", len(faa_files), ".faa files")
    
    # Track sequences and IDs to avoid duplicates
    seen_sequences = set()
    seen_ids = set()
    written_count = 0
    
    with open(output_file, 'w') as out:
        for faa_file in faa_files:
            print(dt.today(), "Processing", faa_file)
            
            current_header = ""
            current_seq = []
            
            with open(faa_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('>'):
                        # Write previous sequence if not duplicate
                        if current_header and current_seq:
                            sequence = ''.join(current_seq)
                            if sequence not in seen_sequences and current_header not in seen_ids:
                                out.write(current_header + "\n")
                                out.write(sequence + "\n")
                                seen_sequences.add(sequence)
                                seen_ids.add(current_header)
                                written_count += 1
                        
                        # Start new sequence
                        current_header = line
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Handle last sequence in file
                if current_header and current_seq:
                    sequence = ''.join(current_seq)
                    if sequence not in seen_sequences and current_header not in seen_ids:
                        out.write(current_header + "\n")
                        out.write(sequence + "\n")
                        seen_sequences.add(sequence)
                        seen_ids.add(current_header)
                        written_count += 1
    
    print(dt.today(), "Merged", written_count, "unique sequences to", output_file)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <dir1> <dir2> <output_file>")
        sys.exit(1)
    
    dir1, dir2, output_file = sys.argv[1:4]
    merge_faa_files(dir1, dir2, output_file)