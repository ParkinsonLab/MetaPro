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

def import_fasta(file_name_in):
    fasta_df = pd.read_csv(file_name_in, error_bad_lines=False, header=None, sep="\n")  # import the fasta
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
    return fasta_df

if __name__ == "__main__":
    input_file = sys.argv[1]            #IN: host-free seqs
    BLAT_tab_file = sys.argv[2]         #IN: BLAT report
    output_file = sys.argv[3]           #OUT: final host-free seqs
    contaminant_output_file = sys.argv[4]#OUT: host-data that was stripped
    
    #import the sequence
    if(input_file.endswith(".fastq")):
        input_seqs = import_fastq(input_file)#SeqIO.to_dict(SeqIO.parse(input_file, "fastq"))
    elif(input_file.endswith(".fasta")):
        input_seqs = import_fasta(input_file)
    #with open(BLAT_tab_file, "r") as tabfile:
    contaminated_seqs = []
    query_seq = ""
    #we don't know enough about BLATOUT files.  we assume we just want the IDs, but we aren't entirely sure
    #what the rules are here.  So we can leave this inefficient loop in for now. (feb 06, 2019)
    #pick out the results from the BLATOUT report
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
        
    
    #The BLAT report spits out what was aligned as human.  
    #so we take the NOT of the intersection to get host-free
    #split the input sequence by host/host-free and export
    host_free_df = input_seqs[~input_seqs.ID.isin(contaminated_seqs)].to_csv(output_file, sep = '\n', mode = "w+", header = False, index = False, quoting = 3)
    host_only_df = input_seqs[input_seqs.ID.isin(contaminated_seqs)].to_csv(contaminant_output_file, sep = '\n', mode = "w+", header = False, index = False, quoting = 3)
   