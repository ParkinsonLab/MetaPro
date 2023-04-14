#this code is supposed to replace samtools fastq, because we need it to filter out data based on a criteria
#This code takes in a samfile from a BWA run, and bisects the data for the purposes of filtering the data.
#The use is to remove reads that align to some sort of contamination (host data, vectors)

import os
import sys
import pandas as pd
from datetime import datetime as dt
import re


def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df

def import_samfile(samfile, filter_stringency):
    no_hit_list = []
    hit_list = []
    all_queries = []
    len_chars= ["M","I","S","=","X"]
    read_count = 0
    with open(samfile, "r") as sam_in:
        for line in sam_in:
            if(line.startswith("@") or len(line) <= 1):
                continue
                
            else:
                read_count += 1
                line_parts = line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and store:
                query = line_parts[0]                        #  queryID= contig/readID,
                db_match = line_parts[2]                     #  geneID, and a
                flag = bin(int(line_parts[1]))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                                                            #  where each bit is a flag for a specific descriptor.
                all_queries.append(query)
                if (flag[8]=="1"): 
                    no_hit_list.append(query)
                else:
                    hit_list.append(query)
                    
                    
    hit_set = set(hit_list)
    no_hit_set = set(no_hit_list)
    
    common_set = set(hit_list).intersection(no_hit_list)
    final_no_hit_list = []
    final_hit_list = []
    
    #there will be reads where the one half of the pair is rejected, and the other half is accepted.  We deal with the straddlers here.
    if(filter_stringency == "high"):
        #take only things that are ensured to be in the no-hit-list.  no overlaps
        final_no_hit_list = [x for x in no_hit_set if x not in common_set] #done this way so that ordering doesn't matter.  issues with set subtraction
        final_hit_list = list(hit_set)
    elif(filter_stringency == "low"):
        final_hit_list = [x for x in hit_set if x not in common_set]
        final_no_hit_list = list(no_hit_set)
    else: #bypass for singletons
        final_hit_list = list(hit_set)
        final_no_hit_list = list(no_hit_set)
    
    all_queries = list(set(all_queries))
    print("query count:", len(all_queries))    
    print(dt.today(), "SAMFILE reads processed:", read_count)
    return final_hit_list, final_no_hit_list
  


if __name__ == "__main__":
    data_style = sys.argv[1]
    filter_stringency = sys.argv[2]
    samfile = sys.argv[3]
        
    if(data_style == "paired"):
        pair_1_raw = sys.argv[4]
        pair_2_raw = sys.argv[5]
        
        pair_1_pass = sys.argv[6]
        pair_2_pass = sys.argv[7]
        
        pair_1_reject = sys.argv[8]
        pair_2_reject = sys.argv[9]
        
        print(dt.today(), "importing data")
        hit_list, no_hit_list = import_samfile(samfile, filter_stringency)
        
        pair_1_df = import_fastq(pair_1_raw)
        pair_2_df = import_fastq(pair_2_raw)
        

        print(dt.today(), "dividing data")
        pair_1_filter_pass = pair_1_df[pair_1_df["ID"].isin(no_hit_list)]
        pair_2_filter_pass = pair_2_df[pair_2_df["ID"].isin(no_hit_list)]
        
        pair_1_filter_reject = pair_1_df[~pair_1_df["ID"].isin(no_hit_list)]
        pair_2_filter_reject = pair_2_df[~pair_2_df["ID"].isin(no_hit_list)]
        
        
        pair_1_raw_list = pair_1_df["ID"].tolist()
        pair_1_pass_list = pair_1_filter_pass["ID"].tolist()
        pair_1_reject_list = pair_1_filter_reject["ID"].tolist()
        
        pair_2_raw_list = pair_2_df["ID"].tolist()
        pair_2_pass_list = pair_2_filter_pass["ID"].tolist()
        pair_2_reject_list = pair_2_filter_reject["ID"].tolist()
        
        p1_total_reads = len(pair_1_raw_list)
        p1_pass_reads = len(pair_1_pass_list)
        p1_reject_reads = len(pair_1_reject_list)
        
        p2_total_reads = len(pair_2_raw_list)
        p2_pass_reads = len(pair_2_pass_list)
        p2_reject_reads = len(pair_2_reject_list)
        
        print("P1 | raw:", p1_total_reads, "host:", p1_reject_reads, "no-host:", p1_pass_reads, "diff:", p1_total_reads - p1_pass_reads - p1_reject_reads)
        print("P2 | raw:", p2_total_reads, "host:", p2_reject_reads, "no-host:", p2_pass_reads, "diff:", p2_total_reads - p2_pass_reads - p2_reject_reads)
            
        print(dt.today(), "paired: exporting reads")
        
        pair_1_filter_pass["ID"] = "@" + pair_1_filter_pass["ID"]
        pair_2_filter_pass["ID"] = "@" + pair_2_filter_pass["ID"]
        
        pair_1_filter_reject["ID"] = "@" + pair_1_filter_reject["ID"]
        pair_2_filter_reject["ID"] = "@" + pair_2_filter_reject["ID"]
        
        pair_1_filter_pass.to_csv(pair_1_pass, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_filter_pass.to_csv(pair_2_pass, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        pair_1_filter_reject.to_csv(pair_1_reject, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_filter_reject.to_csv(pair_2_reject, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        print(dt.today(), "done BWA filter")
                
    elif(data_style == "single"):
        singletons_raw = sys.argv[4]
        singletons_pass = sys.argv[5]
        singletons_reject = sys.argv[6]
        
        print(dt.today(), "importing singleton data")
        hit_list, no_hit_list = import_samfile(samfile, "single")
        singletons_df = import_fastq(singletons_raw)
        print(dt.today(), "dividing singletons")
        singletons_filter_pass = singletons_df[singletons_df["ID"].isin(no_hit_list)]
        singletons_filter_reject = singletons_df[singletons_df["ID"].isin(hit_list)]
        
        singletons_filter_pass["ID"] = "@" + singletons_filter_pass["ID"]
        singletons_filter_reject["ID"] = "@" + singletons_filter_reject["ID"]
        
        print(dt.today(), "exporting singletons")
        singletons_filter_pass.to_csv(singletons_pass, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        singletons_filter_reject.to_csv(singletons_reject, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        print(dt.today(), "bwa filter singletons done!")
                