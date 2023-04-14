#This script is meant to sort the paired reads into sections, based on its mismatch.
#Its logic is based on pre-and-post data

import pandas as pd
from datetime import datetime as dt
import os
import sys

def import_IDs_only(read_file):
    raw_df = pd.read_csv(read_file, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting = 3)
    fastq_df = pd.DataFrame(raw_df.values.reshape(int(len(raw_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    id_list = fastq_df["ID"].tolist()
    
    return id_list
    
def import_fastq(read_file):
    raw_df = pd.read_csv(read_file, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting = 3)
    fastq_df = pd.DataFrame(raw_df.values.reshape(int(len(raw_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    return fastq_df


if __name__ == "__main__":

    pair_1_post_file = sys.argv[1]      #the post-filter reads.  Reads purged of whatever it was supposed to purge
    pair_2_post_file = sys.argv[2]      #IN
    
    pair_1_raw_file = sys.argv[3]       #the input data, to populate
    pair_2_raw_file = sys.argv[4]       #IN
    
    pair_1_accepted_file = sys.argv[5]  #the read files that get accepted
    pair_2_accepted_file = sys.argv[6]  #OUT
    
    pair_1_rejected_file = sys.argv[7]  #the read files that get rejected 
    pair_2_rejected_file = sys.argv[8]  #OUT
    
    operating_mode  = sys.argv[9]       #option: either OR (loose), or AND (tight, conservative)
    
    #-------------------------------------------------------
    pair_1_id_list = import_IDs_only(pair_1_in_file)
    pair_2_id_list = import_IDs_only(pair_2_in_file)
    
    pair_1_raw_df = import_fastq(pair_1_raw_file)
    pair_2_raw_df = import_fastq(pair_2_raw_file)
    
    common_id_list = list(set(pair_1_id_list) & set(pair_2_id_list))
    
    pair_1_accepted_df = pair_1_raw_df[pair_1_raw_df["ID"].isin(common_id_list)]
    pair_2_accepted_df = pair_2_raw_df[pair_2_raw_df["ID"].isin(common_id_list)]
    
    pair_1_rejected_df = pair_1_raw_df[~pair_1_raw_df["ID"].isin(common_id_list)]
    pair_2_rejected_df = pair_2_raw_df[~pair_2_raw_df["ID"].isin(common_id_list)]
    
    
    
    
    
    