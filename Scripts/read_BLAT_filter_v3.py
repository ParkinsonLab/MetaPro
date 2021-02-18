# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 23:05:20 2019

@author: Billy
"""

#!/usr/bin/env python

import sys
import pandas as pd
from datetime import datetime as dt


def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df


    
def import_blat(BLAT_tab_file):
    contaminated_seqs = []
    query_seq = ""

    with open(BLAT_tab_file, "r") as tabfile:
        for line in tabfile:
            if len(line) < 2:
                continue
            line_parts = line.split("\t")
            if query_seq == line_parts[0]:
                continue
            else:
                query_seq = line_parts[0]
                contaminated_seqs.append(query_seq)
    return contaminated_seqs

if __name__ == "__main__":
    data_style              = sys.argv[1]
    operating_mode          = sys.argv[2]
    
    if(data_style == "paired"):
        pair_1_raw_file     = sys.argv[3]
        pair_2_raw_file     = sys.argv[4]
        
        pair_1_blat         = sys.argv[5]
        pair_2_blat         = sys.argv[6]
        
        pair_1_accepted     = sys.argv[7]
        pair_2_accepted     = sys.argv[8]
        
        pair_1_rejected     = sys.argv[9]
        pair_2_rejected     = sys.argv[10]
        
        pair_1_blat_list = import_blat(pair_1_blat)
        pair_2_blat_list = import_blat(pair_2_blat)
        
        pair_1_raw_df = import_fastq(pair_1_raw_file)
        pair_2_raw_df = import_fastq(pair_2_raw_file)
        
        #grab the list of reads based off of our filter rules: logic OR if its to be loose (because we take the inverse), or logic AND if we're being conversative
        #BLAT looks for the stuff we want filtered.  in instances where 1 read is filtered, and the other isn't: low stringency lets the filtered read go.  high stringency arrests the unfiltered read
        common_list = []
        if(operating_mode == "low"):
            common_list = list(set(pair_1_blat_list).intersection(set(pair_2_blat_list)))
        elif(operating_mode == "high" or operating_mode == "HIGH"):
            common_list = list(set().union(pair_1_blat_list,pair_2_blat_list))
            
        print("pair 1 blat:", len(pair_1_blat_list), "pair 2 blat:", len(pair_2_blat_list), "common:", len(common_list))
            
        pair_1_accepted_df = pair_1_raw_df[~pair_1_raw_df["ID"].isin(common_list)]
        pair_2_accepted_df = pair_2_raw_df[~pair_2_raw_df["ID"].isin(common_list)]
        
        pair_1_accepted_id = sorted(pair_1_accepted_df["ID"].tolist())
        pair_2_accepted_id = sorted(pair_2_accepted_df["ID"].tolist())
        
        
        
        if(pair_1_accepted_id != pair_2_accepted_id):
            sys.exit("accepted DFs do not all have the same ID.  there's a problem with the data")
        else:
            print(dt.today(), "BLAT filtering OK.  pair 1 IDs match pair 2")
        
        pair_1_rejected_df = pair_1_raw_df[pair_1_raw_df["ID"].isin(common_list)]
        pair_2_rejected_df = pair_2_raw_df[pair_2_raw_df["ID"].isin(common_list)]
        
        # export it
        
        pair_1_accepted_df["ID"] = "@" + pair_1_accepted_df["ID"]
        pair_2_accepted_df["ID"] = "@" + pair_2_accepted_df["ID"]
        pair_1_rejected_df["ID"] = "@" + pair_1_rejected_df["ID"]
        pair_2_rejected_df["ID"] = "@" + pair_2_rejected_df["ID"]

        pair_1_accepted_df.to_csv(pair_1_accepted, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_accepted_df.to_csv(pair_2_accepted, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        pair_1_rejected_df.to_csv(pair_1_rejected, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_rejected_df.to_csv(pair_2_rejected, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        
    elif(data_style == "single"):
        singletons_raw_seqs = sys.argv[3]
        singletons_blat     = sys.argv[4]
        singletons_accepted = sys.argv[5] #data that passes through the filter
        singletons_rejected = sys.argv[6] #data that gets caught by the filter
        
        singletons_blat_list = import_blat(singletons_blat)
        singletons_raw_df = import_fastq(singletons_raw_seqs)
        
        singletons_accepted_df = singletons_raw_df[~singletons_raw_df["ID"].isin(singletons_blat_list)] 
        singletons_rejected_df = singletons_raw_df[singletons_raw_df["ID"].isin(singletons_blat_list)]
        
        singletons_accepted_df["ID"] = "@" + singletons_accepted_df["ID"]
        singletons_rejected_df["ID"] = "@" + singletons_rejected_df["ID"]
        
        singletons_accepted_df.to_csv(singletons_accepted, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        singletons_rejected_df.to_csv(singletons_rejected, sep = "\n", mode = "w", header = False, index = False, quoting = 3)