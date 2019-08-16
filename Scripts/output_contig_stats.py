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
        self.df_orig = pd.read_csv(input_file, error_bad_lines=False, header=None, sep="\n")  # import the fasta
        self.df_orig.columns = ["row"]
        #There's apparently a possibility for NaNs to be introduced in the raw fasta.  We have to strip it before we process (from DIAMOND proteins.faa)
        self.df_orig.dropna(inplace=True)
        new_df = pd.DataFrame(self.df_orig.loc[self.df_orig.row.str.contains('>')])  # grab all the IDs
        new_df.columns = ["names"]
        new_data_df = self.df_orig.loc[~self.df_orig.row.str.contains('>')]  # grab the data
        new_data_df.columns = ["data"]
        self.df_orig = new_df.join(new_data_df, how='outer')  # join them into a 2-col DF
        self.df_orig["names"] = self.df_orig.fillna(method='ffill')  # fill in the blank spaces in the name section
        self.df_orig.dropna(inplace=True)  # remove all rows with no sequences
        self.df_orig.index = self.df_orig.groupby('names').cumcount()  # index it for transform
        temp_columns = self.df_orig.index  # save the index names for later
        self.df_orig = self.df_orig.pivot(values='data', columns='names')  # pivot
        self.df_orig = self.df_orig.T  # transpose
        self.df_orig["seq"] = self.df_orig[self.df_orig.columns[:]].apply(lambda x: "".join(x.dropna()), axis=1)  # consolidate all cols into a single sequence
        self.df_orig.drop(temp_columns, axis=1, inplace=True)
        #not really needed
        #self.df_orig["names"] = self.df_orig.index
        #self.df_orig.index = range(self.df_orig.shape[0])
        # At this point, we've already got the number of reads.
        
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