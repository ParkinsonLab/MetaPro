import os
import os.path
import sys
import pandas as pd
import time

#This module takes in the Output report of the barrnap tool, and the fastq it scanned.
#The goal is to bisect the fastq into 2 piles:  entries with IDs that were located by Barrnap (rRNA)
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

def extract_rRNA_ID_fasta(barrnap_file):
    #If there's a better way to do this, I'd like to see it.
    ID_list = set()
    barrnap_list = open(barrnap_file, mode='r')
    
    for item in barrnap_list:
        if(not item.startswith("#")):
            ID_list.add(">" + item.split("\t")[0])
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
    fastq_df[fastq_df["ID"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False, quoting = 3)
    fastq_df[~fastq_df["ID"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False, quoting = 3)
    #writing mRNA last, so we can use it to check the progress of the Barrnap step
    
def filter_rRNA_fasta(rRNA_ID_list, fasta_sequence, mRNA_loc, rRNA_loc, file_name):
    fasta_df = pd.read_csv(fasta_sequence, error_bad_lines=False, header=None, sep="\n")  # import the fasta
    fasta_df.columns = ["row"]
    #There's apparently a possibility for NaNs to be introduced in the raw fasta.  We have to strip it before we process (from DIAMOND proteins.faa)
    fasta_df.dropna(inplace=True)
    new_df = pd.DataFrame(fasta_df.loc[fasta_df.row.str.contains('>')])  # grab all the IDs
    new_df.columns = ["names"]
    new_data_df = fasta_df.loc[~fasta_df.row.str.contains('>')]  # grab the data
    new_data_df.columns = ["data"]
    fasta_df = new_df.join(new_data_df, how='outer')  # join them into a 2-col DF
    fasta_df["names"] = fasta_df.fillna(method='ffill')  # fill in the blank spaces in the name section
    fasta_df.dropna(inplace=True)  # remove all rows with no sequences
    fasta_df.index = fasta_df.groupby('names').cumcount()  # index it for transform
    temp_columns = fasta_df.index  # save the index names for later
    fasta_df = fasta_df.pivot(values='data', columns='names')  # pivot
    fasta_df = fasta_df.T  # transpose
    fasta_df["sequence"] = fasta_df[fasta_df.columns[:]].apply(lambda x: "".join(x.dropna()), axis=1)  # consolidate all cols into a single sequence
    fasta_df.drop(temp_columns, axis=1, inplace=True)
    
    #doing it inline so we don't create another DF
    #mRNA segment
    mRNA_export = os.path.join(mRNA_loc, file_name + "_mRNA.fasta")
    rRNA_export = os.path.join(rRNA_loc, file_name + "_rRNA.fasta")
    fasta_df[fasta_df["names"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False, quoting = 3)
    fasta_df[~fasta_df["names"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False, quoting = 3)

if __name__ == "__main__":
    
    barrnap_file = sys.argv[1]
    input_sequence = sys.argv[2]
    mRNA_location = sys.argv[3]
    rRNA_location = sys.argv[4]
    segment_name = sys.argv[5]
    
    fasta_sequence = ""
    fastq_sequence = ""
    if(input_sequence.endswith(".fastq")):
        fastq_sequence = input_sequence
        ID_list = extract_rRNA_ID(barrnap_file)
        filter_rRNA(ID_list, fastq_sequence, mRNA_location, rRNA_location, segment_name)
    
    else:
        fasta_sequence = input_sequence
        ID_list = extract_rRNA_ID_fasta(barrnap_file)
        filter_rRNA_fasta()