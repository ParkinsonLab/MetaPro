import pandas as pd
import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import math

"""
This code reads in a FASTA, and does stats on the reads.
The only function here is to find N50
"""

class contig_stats:
    def __init__(self, input_file):
        import pandas as pd
import numpy as np

class contig_stats:
    def __init__(self, input_file):
        # 1. Read the file into a list of lines (Standard Python reading)
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # 2. Process the lines to create two lists: headers and sequences
        headers = []
        sequences = []
        current_seq = []
        
        for line in lines:
            if line.startswith('>'):
                # Found a new header, save the previous sequence (if one exists)
                if current_seq:
                    sequences.append("".join(current_seq))
                
                # Start the new sequence and header
                headers.append(line)
                current_seq = []
            else:
                # Accumulate sequence data
                current_seq.append(line)
        
        # 3. Save the last sequence
        if current_seq:
            sequences.append("".join(current_seq))

        # 4. Create the final DataFrame
        self.df_orig = pd.DataFrame({
            "names": headers,
            "seq": sequences
        })
        
        # At this point, self.df_orig has two clean columns: 'names' and 'seq'
        # The remaining logic in your original script (pivot, apply, drop) is now unnecessary.
        # You can continue with your statistics calculation on this clean DataFrame.
        # Example: print(self.df_orig.head())
        
    def contig_stats(self):
        
        #get N50 from the contig
        #N50 is data pulled from contig lengths assembled end-to-end
        df_0 = pd.DataFrame(self.df_orig["seq"])
        df_0["len"] = df_0["seq"].apply(lambda x: len(x))
        df_0.sort_values(by=["len"], ascending = False, inplace = True)
        total_contig_length = df_0["len"].sum()
        #print(df_0)
        #print("total contig length:", total_contig_length)
        halfway_point = total_contig_length / 2
        #print("halfway:", halfway_point)
        cur_sum = total_contig_length
        n50_val = 0
        l50_val = 0
        count = 0
        for index, item in df_0.iterrows():
            prev_cur_sum = cur_sum
            cur_sum -= item["len"]
            l50_val += 1
            #print(l50_val, prev_cur_sum, " -> " , cur_sum , "(", item["len"], ")")
            if(cur_sum < halfway_point):
                #print("n50 found:", item["len"])
                n50_val = item["len"]
                break
        return [n50_val, l50_val]
        
            
if __name__ == "__main__":
    contig_file = sys.argv[1]
    output_file = sys.argv[2]
    stat_obj = contig_stats(contig_file)
    stat_report = stat_obj.contig_stats()
    
    stat_file = open(output_file, "w")
    n50_line = "N50: " + str(stat_report[0]) + " \n"
    l50_line = "L50: " + str(stat_report[1]) + " \n"
    stat_file.write(n50_line)
    stat_file.write(l50_line)
    stat_file.close()