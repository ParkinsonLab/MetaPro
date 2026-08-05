#!/usr/bin/env python3
"""
Merge all nucleotide FASTA files in a directory, removing duplicates.
Usage: python merge_fasta_files.py <input_directory> <output_file>
"""

import os
import sys


def main():
    if len(sys.argv) != 3:
        print("Usage: python merge_fasta_files.py <input_directory> <output_file>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # Get all files in directory
    files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) 
             if os.path.isfile(os.path.join(input_dir, f))]
    
    unique_sequences = set()
    
    with open(output_file, 'w') as out_f:
        for file_path in files:
            with open(file_path, 'r') as f:
                header = None
                sequence = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Process previous sequence
                        if header and sequence:
                            seq_str = ''.join(sequence).upper()
                            if seq_str not in unique_sequences:
                                unique_sequences.add(seq_str)
                                out_f.write(header + '\n')
                                for i in range(0, len(seq_str), 80):
                                    out_f.write(seq_str[i:i+80] + '\n')
                        
                        header = line
                        sequence = []
                    elif line:
                        sequence.append(line)
                
                # Process last sequence in file
                if header and sequence:
                    seq_str = ''.join(sequence).upper()
                    if seq_str not in unique_sequences:
                        unique_sequences.add(seq_str)
                        out_f.write(header + '\n')
                        for i in range(0, len(seq_str), 80):
                            out_f.write(seq_str[i:i+80] + '\n')


if __name__ == "__main__":
    main()