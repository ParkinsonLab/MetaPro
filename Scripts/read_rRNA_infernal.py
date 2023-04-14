<<<<<<< HEAD
import os
import os.path
import sys
import pandas as pd


#This module takes in the Output report of the infernal tool, and the fastq it scanned.
#The goal is to bisect the fastq into 2 piles:  entries with IDs that were located by Infernal (rRNA)
#-> and entries that weren't (mRNA), and then export them

def extract_rRNA_ID(inf_file):
    ID_list = set()
    inf_list = open(inf_file, mode='r')
    
    for item in inf_list:
        if(not item.startswith("#")):
            # needs the @ at the start, for a match with the FASTQ's IDs
            ID_list.add("@" + item.split()[0])
=======
#This script is meant to sort the paired reads into sections, based on its mismatch.
#there's an 
#This script replaces the orphanizer. 
#this replaces the read_rRNA_infernal 

import pandas as pd
from datetime import datetime as dt
import os
import sys

def import_infernal_rRNA(inf_file):
    ID_list = set()
    inf_list = open(inf_file, mode='r')
    line_count = 0
    for item in inf_list:
        line_count += 1
        if(not item.startswith("#")):
            # needs the @ at the start, for a match with the FASTQ's IDs
            fastq_list = item.split()
            #if(len(fastq_list) < 2):
            #print("count:", line_count, "fastq id split:", fastq_list)
            fastq_id = item.split()[2]
            
            #print("FASTQ ID located:", fastq_id)
            ID_list.add("@" + fastq_id)
>>>>>>> db_shrink
        elif(len(item) == 2):
            break
    #pandas can't deal with sets, but we only need unique elements        
    return list(ID_list)
<<<<<<< HEAD

def filter_rRNA(rRNA_ID_list, fastq_sequence, mRNA_loc, rRNA_loc, file_name):
    #import the fastq as a DF
<<<<<<< HEAD
    fastq_df = pd.read_csv(fastq_sequence, sep="\n", skip_blank_lines=False, header=None, names = [None])
=======
    fastq_df = pd.read_csv(fastq_sequence, sep="\n", skip_blank_lines=False, header=None, names = [None], quoting=3)
>>>>>>> 4d5286c... committing final-ish code.
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    #doing it inline so we don't create another DF
    #mRNA segment
    mRNA_export = os.path.join(mRNA_loc, file_name + "_mRNA.fastq")
    rRNA_export = os.path.join(rRNA_loc, file_name + "_rRNA.fastq")
<<<<<<< HEAD
    fastq_df[~fastq_df["ID"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False)
    fastq_df[fastq_df["ID"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False)
=======
    fastq_df[~fastq_df["ID"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False, quoting = 3)
    fastq_df[fastq_df["ID"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False, quoting = 3)
>>>>>>> 4d5286c... committing final-ish code.
    
    
if __name__ == "__main__":
    
    inf_file = sys.argv[1] #infernal
    fastq_sequence = sys.argv[2]
    mRNA_location = sys.argv[3]
    rRNA_location = sys.argv[4]
    segment_root_name = sys.argv[5]
    ID_list = extract_rRNA_ID(inf_file)
    filter_rRNA(ID_list, fastq_sequence, mRNA_location, rRNA_location, segment_root_name)
    
=======
    
def import_fastq(read_file):
    raw_df = pd.read_csv(read_file, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting = 3)
    fastq_df = pd.DataFrame(raw_df.values.reshape(int(len(raw_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    return fastq_df

def import_barrnap(barrnap_file):
    #If there's a better way to do this, I'd like to see it.
    ID_list = set()
    barrnap_list = open(barrnap_file, mode='r')
    
    for item in barrnap_list:
        if(not item.startswith("#")):
            ID_list.add("@" + item.split("\t")[0])
    #pandas can't deal with sets, but we only need unique elements        
    return list(ID_list)


if __name__ == "__main__":
    
    operating_mode = sys.argv[1]        #option: liberal constraints OR, or conversative constraints AND
    data_style = sys.argv[2]            #option: paired, or single
    
    
    data_flag = "single"
    if(data_style == "PAIRED" or data_style == "paired"):
        data_flag = "paired"
        
    if(data_flag == "paired"):
        
        pair_1_inf_file = sys.argv[3]       #the post-filter infernal_file
        pair_2_inf_file = sys.argv[4]       #IN
        
        pair_1_barrnap_file = sys.argv[5]   #the post-filter barrnap file
        pair_2_barrnap_file = sys.argv[6]   #IN
    
        pair_1_raw_file = sys.argv[7]       #the input data, to populate
        pair_2_raw_file = sys.argv[8]       #IN
    
        pair_1_accepted_file = sys.argv[9]  #the read files that get accepted
        pair_2_accepted_file = sys.argv[10]  #OUT
    
        pair_1_rejected_file = sys.argv[11]  #the read files that get rejected 
        pair_2_rejected_file = sys.argv[12]  #OUT
    
    else:
        single_inf_file = sys.argv[3]
        single_barrnap_file = sys.argv[4]
        single_raw_file = sys.argv[5]
        single_accepted_file = sys.argv[6]
        single_rejected_file = sys.argv[7]
    #-------------------------------------------------------
    
    if(data_flag == "paired"):
        
        pair_1_id_list = import_infernal_rRNA(pair_1_inf_file)
        pair_2_id_list = import_infernal_rRNA(pair_2_inf_file)
        
        pair_1_raw_df = import_fastq(pair_1_raw_file)
        pair_2_raw_df = import_fastq(pair_2_raw_file)
        
        pair_1_barrnap_list = import_barrnap(pair_1_barrnap_file)
        pair_2_barrnap_list = import_barrnap(pair_2_barrnap_file)
        
        pair_1_id_list += pair_1_barrnap_list
        pair_2_id_list += pair_2_barrnap_list
        
        #grab the common IDs.  They'll be the ones we want. The ones that meet the filter
        if(operating_mode == "low"):     #liberal:  only things that are included on both sides is rRNA
            common_id_list = list(set(pair_1_id_list).intersection(pair_2_id_list)) 
        elif(operating_mode == "high"):  #conservative: anything in this unioned list is rRNA
            common_id_list = list(set().union(pair_1_id_list, pair_2_id_list))
        
        #if the read is not inside the rRNA list, we take it.
        pair_1_accepted_df = pair_1_raw_df[~pair_1_raw_df["ID"].isin(common_id_list)]
        pair_2_accepted_df = pair_2_raw_df[~pair_2_raw_df["ID"].isin(common_id_list)]
        
        pair_1_accepted_size = pair_1_accepted_df.shape[0]
        pair_1_raw_size = pair_1_raw_df.shape[0]
        print(pair_1_accepted_size, "vs:", len(common_id_list))
        
        pair_1_rejected_df = pair_1_raw_df[pair_1_raw_df["ID"].isin(common_id_list)]
        pair_2_rejected_df = pair_2_raw_df[pair_2_raw_df["ID"].isin(common_id_list)]
        
        
        pair_1_accepted_df.to_csv(pair_1_accepted_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_accepted_df.to_csv(pair_2_accepted_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        pair_1_rejected_df.to_csv(pair_1_rejected_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        pair_2_rejected_df.to_csv(pair_2_rejected_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
    else:
        single_id_list = import_infernal_rRNA(single_inf_file)
        single_barrnap_list = import_barrnap(single_barrnap_file)
        single_raw_df = import_fastq(single_raw_file)
        
        single_id_list += single_barrnap_list
        
        single_accepted_df = single_raw_df[~single_raw_df["ID"].isin(single_id_list)]
        single_rejected_df = single_raw_df[single_raw_df["ID"].isin(single_id_list)]
        single_accepted_df.to_csv(single_accepted_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        single_rejected_df.to_csv(single_rejected_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        
>>>>>>> db_shrink
