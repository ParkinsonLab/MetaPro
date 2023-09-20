#This script is meant to sort the paired reads into sections, based on its mismatch.
#there's an 
#This script replaces the orphanizer. 
#this replaces the read_rRNA_infernal 

#Sept 13, 2023: exports the junk as ID lists. but it won't bother exporting infernal-only mRNA.
#Also: we don't really care about what tool removes what.  If we needed it, that could be reverse-engineered.  



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
        elif(len(item) == 2):
            break
    #pandas can't deal with sets, but we only need unique elements        
    return list(ID_list)
    
def import_fastq(read_file):
    raw_df = pd.read_csv(read_file, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting = 3)
    fastq_df = pd.DataFrame(raw_df.values.reshape(int(len(raw_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    return fastq_df





if __name__ == "__main__":
    
    operating_mode = sys.argv[1]        #option: liberal constraints OR, or conversative constraints AND
    data_style = sys.argv[2]            #option: paired, or single
    
    
    data_flag = "single"
    if(data_style == "PAIRED" or data_style == "paired"):
        data_flag = "paired"
        
    if(data_flag == "paired"):
        
        pair_1_inf_file = sys.argv[3]       #the post-filter infernal_file
        pair_2_inf_file = sys.argv[4]       #IN
        
        pair_1_raw_file = sys.argv[5]       #the input data, to populate
        pair_2_raw_file = sys.argv[6]       #IN
    
        pair_1_accepted_file = sys.argv[7]  #the read files that get accepted
        pair_2_accepted_file = sys.argv[8]  #OUT
    
        pair_1_rejected_file = sys.argv[9]  #the read files that get rejected 
        pair_2_rejected_file = sys.argv[10]  #OUT
    
    else:
        single_inf_file = sys.argv[3]
        single_raw_file = sys.argv[4]
        single_accepted_file = sys.argv[5]
        single_rejected_file = sys.argv[6]
    #-------------------------------------------------------
    
    if(data_flag == "paired"):
        
        pair_1_id_list = import_infernal_rRNA(pair_1_inf_file)
        pair_2_id_list = import_infernal_rRNA(pair_2_inf_file)
        
        pair_1_raw_df = import_fastq(pair_1_raw_file)
        pair_2_raw_df = import_fastq(pair_2_raw_file)

        
        #grab the common IDs between 1 and 2.  They'll be the ones we want. The ones that meet the filter
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
        single_raw_df = import_fastq(single_raw_file)
        
        
        single_accepted_df = single_raw_df[~single_raw_df["ID"].isin(single_id_list)]
        single_rejected_df = single_raw_df[single_raw_df["ID"].isin(single_id_list)]
        single_accepted_df.to_csv(single_accepted_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        single_rejected_df.to_csv(single_rejected_file, sep = "\n", mode = "w", header = False, index = False, quoting = 3)
        
        