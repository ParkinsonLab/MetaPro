import pandas as pd
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import math




class read_quality_metrics:
    def __init__(self, input_file):
        self.input_file = input_file
        self.input_suffix = os.path.split(self.input_file)[1].split(".")[1]
        
        if(self.input_suffix == "fastq"):
<<<<<<< HEAD
            self.df_file = pd.read_csv(self.input_file, header = None, names=None, sep='\n', skip_blank_lines=False)
=======
            self.df_file = pd.read_csv(self.input_file, header = None, names=None, sep='\n', skip_blank_lines=False, quoting = 3)
>>>>>>> 4d5286c... committing final-ish code.
            self.df_orig = pd.DataFrame(self.df_file.values.reshape(int(len(self.df_file)/4), 4))
            self.df_orig.columns = ["ID", "seq", "junk", "quality"]
            self.df_orig.drop(columns =["ID", "seq", "junk"], inplace = True)
            
            
        elif(self.input_suffix == "fasta"):
            self.df_orig = pd.read_csv(self.input_file, error_bad_lines=False, header=None, sep="\n")  # import the fasta
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
            
        
        
    def string_to_ascii_array(self, line):
        new_line = ""
        for item in line:
            new_line += str(ord(item)) + ","
        return new_line[:-1]
    
    def avg_ascii(self, line):
        sum = 0
        for item in line:
            sum += ord(item)
        return sum / len(line)
        
    
    
    def per_base_quality(self, df_0):
        if(self.input_suffix == "fastq"):
            #df_0 = self.df_orig["quality"]
            df_0["quality"] = df_0["quality"].apply(lambda x: self.string_to_ascii_array(x))
            df_0 = df_0["quality"].str.split(",", expand=True).rename(columns = lambda x: "bp_" + str(x))
            df_0 = df_0.apply(pd.to_numeric)
            df_0.loc["avg"] = df_0.select_dtypes(pd.np.number).sum() / df_0.shape[0] #doing on select_dtypes means we only consider numerics.  
            
            stats_df = pd.DataFrame(df_0.loc["avg"]).transpose()
            stats_df.loc["SD"] = df_0.select_dtypes(pd.np.number).std()
            stats_df.loc["Median"] = df_0.select_dtypes(pd.np.number).median()
            stats_df.loc["Max"] = df_0.select_dtypes(pd.np.number).max()
            stats_df.loc["Min"] = df_0.select_dtypes(pd.np.number).min()
            stats_df.loc["Q1"] = df_0.select_dtypes(pd.np.number).quantile(0.25)
            stats_df.loc["Q3"] = df_0.select_dtypes(pd.np.number).quantile(0.75)
            #new_name = self.output_prefix + "_per_base_quality_report.csv"
            stats_df.to_csv(new_name, mode="w+", header=False, index=False)
        else:
            print("can't run per-base quality on this file. it's not as FASTQ (.fastq)")
            
    def per_sequence_quality(self, df_0, location, graph_title):
        if(self.input_suffix == "fastq"):
            #number of reads vs mean quality of base pairs in the read
            #df_0 = pd.DataFrame(self.df_orig["quality"])
            df_0["len"] = df_0["quality"].apply(lambda x: len(x))
            df_0["quality"] = df_0["quality"].apply(lambda x: self.avg_ascii(x))
            bins = range(0, 100, 10)
            hist = df_0.hist(column="quality", bins = bins)#, by = "len") This is actually calling matplotlib
            new_title = "Reads vs Average read quality of:" + graph_title
            for line in hist.flatten():
                line.set_title(new_title)
                line.set_xlabel("Average Read Quality")
                line.set_ylabel("Number of reads in file")
            new_name =  location + "_per_seq_quality_report.csv"
            print("report location:", new_name)
            print("hist location:", location + "_hist.jpg")
            plt.savefig(location + "_hist.jpg")
            df_0.to_csv(new_name, mode = "w+", header=False, index=False)
            
            
        else:
            print("can't run per-read quality on this file. it's not as FASTQ (.fastq)")
    
    def contig_stats(self, df_0):
        if(self.input_suffix == "fasta"):
            #get N50 from the contig
            #N50 is data pulled from contig lengths assembled end-to-end
            #df_0 = pd.DataFrame(self.df_orig["seq"])
            df_0["len"] = df_0["seq"].apply(lambda x: len(x))
            df_0.sort_values(by=["len"], ascending = False, inplace = True)
            total_contig_length = df_0["len"].sum()
            print(df_0)
            print("total contig length:", total_contig_length)
            halfway_point = total_contig_length / 2
            print("halfway:", halfway_point)
            cur_sum = total_contig_length
            n50_val = 0
            l50_val = 0
            count = 0
            for index, item in df_0.iterrows():
                prev_cur_sum = cur_sum
                cur_sum -= item["len"]
                l50_val += 1
                print(l50_val, prev_cur_sum, " -> " , cur_sum , "(", item["len"], ")")
                if(cur_sum < halfway_point):
                    print("n50 found:", item["len"])
                    n50_val = item["len"]
                    break
            return [n50_val, l50_val]
        else:
            print("not a fasta.  not running n50")
            return 0
                
if __name__ == "__main__":
    plt.ioff()
    
    op_mode = sys.argv[1]
    if(op_mode == "paired"):
    
        input_1 = sys.argv[2]
        input_2 = sys.argv[3]
        
        qc_1 = sys.argv[4]
        qc_2 = sys.argv[5]
        qc_singletons = sys.argv[6]
        
        file_destination = sys.argv[7]
        
        input_1_stats_obj = read_quality_metrics(input_1)
        input_2_stats_obj = read_quality_metrics(input_2)
        qc_1_stats_obj = read_quality_metrics(qc_1)
        qc_2_stats_obj = read_quality_metrics(qc_2)
        qc_singletons_stats_obj = read_quality_metrics(qc_singletons)
        
        input_df = pd.concat([input_1_stats_obj.df_orig, input_2_stats_obj.df_orig])
        input_df.index = range(input_df.shape[0])
        input_1_stats_obj.per_sequence_quality(input_df, file_destination + "/input", "raw sequences")
     
        
        qc_df = pd.concat([qc_1_stats_obj.df_orig, qc_2_stats_obj.df_orig, qc_singletons_stats_obj.df_orig])
        qc_df.index = range(qc_df.shape[0])
        qc_1_stats_obj.per_sequence_quality(qc_df, file_destination + "/quality_filter", "post quality-filter")

    elif(op_mode == "single"):
        input_1 = sys.argv[2]
        qc_1 = sys.argv[3]
        file_destination = sys.argv[4]
        input_1_stats_obj = read_quality_metrics(input_1)
        qc_1_stats_obj = read_quality_metrics(qc_1)
        
        input_1_stats_obj.per_sequence_quality(input_1_stats_obj.df_orig, file_destination + "/input", "raw sequences")
        qc_1_stats_obj.per_sequence_quality(qc_1_stats_obj.df_orig, file_destination + "/quality_filter", "post quality-filter")

   