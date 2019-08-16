import os
import os.path
import sys
import pandas as pd


#This module takes in the Output report of the barrnap tool, and the fastq it scanned.
#The goal is to bisect the fastq into 2 piles:  entries with IDs that were located by Infernal (rRNA)
#-> and entries that weren't (mRNA), and then export them

def extract_rRNA_ID(barrnap_file):
    #If there's a better way to do this, I'd like to see it.
    ID_list = set()
    barrnap_list = open(barrnap_file, mode='r')
    
    for item in barrnap_list:
        if(not item.startswith("#")):
            ID_list.add("@" + item.split("\t")[0])
    #pandas can't deal with sets, but we only need unique elements        
    return list(ID_list)

def filter_rRNA(rRNA_ID_list, fastq_sequence, mRNA_loc, rRNA_loc, file_name):
    #import the fastq as a DF
    fastq_df = pd.read_csv(fastq_sequence, sep="\n", skip_blank_lines=False, quoting=3, header=None, names = [None])
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    #doing it inline so we don't create another DF
    #mRNA segment
    mRNA_export = os.path.join(mRNA_loc, file_name + "_mRNA.fastq")
    rRNA_export = os.path.join(rRNA_loc, file_name + "_rRNA.fastq")
    fastq_df[~fastq_df["ID"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False, quoting = 3)
    fastq_df[fastq_df["ID"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False, quoting = 3)
    
    
if __name__ == "__main__":
    
    barrnap_file = sys.argv[1]
    fastq_sequence = sys.argv[2]
    mRNA_location = sys.argv[3]
    rRNA_location = sys.argv[4]
    segment_name = sys.argv[5]
    ID_list = extract_rRNA_ID(barrnap_file)
    filter_rRNA(ID_list, fastq_sequence, mRNA_location, rRNA_location, segment_name)
    