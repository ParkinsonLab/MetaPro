#Jan 29, 2019:
#-------------------------------------
#This module takes in all final-result fastqs and fastas, and gives us a report on the change in reads throughout the process.
#It serves as a pipeline estimator.
#It'll generate a report at the end of each stage.
#The table is .tsv

import pandas as pd
import os
import sys

def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
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
    

    pre_operation_seq_file = sys.argv[1]    #IN: file before the procedure
    post_operation_seq_file = sys.argv[2]   #IN: file after the procedure
    summary_table_file = sys.argv[3]        #OUT: report table to be generated
    if not(summary_table_file.endswith(".tsv")):
        summary_table_file += ".tsv"
    
    #import the read files
    pre_df = None
    if(pre_operation_seq_file.endswith(".fa") or pre_operation_seq_file.endswith(".fasta") or pre_operation_seq_file.endswith(".faa")):
        pre_df = import_fasta(pre_operation_seq_file)
        
    elif(pre_operation_seq_file.endswith("fq") or pre_operation_seq_file.endswith("fastq")):
        pre_df = import_fastq(pre_operation_seq_file)
        
    post_df = None
    if(post_operation_seq_file.endswith(".fa") or post_operation_seq_file.endswith(".fasta") or post_operation_seq_file.endswith(".faa")):
        post_df = import_fasta(post_operation_seq_file)
        
    elif(post_operation_seq_file.endswith("fq") or post_operation_seq_file.endswith("fastq")):
        post_df = import_fastq(post_operation_seq_file)
        
    #do math (get counts, get diff)    
    pre_count = pre_df.shape[0]
    post_count = post_df.shape[0]
    diff = pre_count - post_count
    diff_percent = round(((diff / pre_count) * 100), 2)
    
    #export it!
    data_table = open(summary_table_file, "w")
    data_table.write("number of reads before operation\t" + "number of reads after operation\n")
    data_table.write(str(pre_count) + "\t" + str(post_count) + "\n")
    data_table.write("difference\t" + str(pre_count - post_count) + "\n")
    data_table.write("diff percent\t" + str(diff_percent))
    data_table.close()
    
    
    
        

