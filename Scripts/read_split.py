#!/usr/bin/env python
# Oct 16, 2017
# -------------------------------------------------
# this module looks to be responsible for splitting up the fasta or fastq.

# Nov 21, 2017
# -------------------------------------- 
# upon further inspection, this thing really does just split up the fastq, and fasta.  but with unnecessary logic
# According to BJ:  splitting fasta isn't necessary anymore.  It was a leftover piece of code.  
# This thing will only split Fastq until we see a need for something else

import math as m
import os
import os.path
import sys
import pandas as pd
from datetime import datetime as dt

def split_fastq(file_name_in, file_name_out, chunks):#split_count = 4):
    print(dt.today(), "FASTQ file name in:", file_name_in)
    #FASTQ has 4 lines per entry.
    file_base_name = os.path.splitext(file_name_in)[0]
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    #At this point, we've already got the number of reads.
    #chunks = m.ceil(len(fastq_df) / split_count) #how many sequences each split file will have
    #print("total df length:", len(fastq_df))
    print("chunk size:", chunks)
    if(chunks < 1):
        print(dt.today(), "chunks is set to a default of 10000.  originally:", chunks)
        chunks = 10000
    #if(chunks < 1):
    #    print("split count too large. not enough info to split")
    #    chunks = 1
    #    print("new split count:", split_count)
        
    #for i in range(0, split_count):
    index_count = 0
    while(True):
        print("working on segment :", index_count +1, "of fastq splitter")
        #fancy naming
        new_file_name = file_name_out + "_" + str(index_count) + ".fastq"
        if(split_count == 1):
            new_file_name = file_name_out + ".fastq"
        #split file by selective selection, and writing
        start_index = int(index_count * chunks)
        end_index = int(((index_count+1) * chunks))
        #if(chunks == 1):
        #    end_index += 1 #override on splits that only have 1
        index_count += 1
        if not(fastq_df.iloc[start_index:end_index, :].empty):
            fastq_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize = chunks, mode = "w+", index=False, sep='\n', header=False, quoting = 3)
        else:
            print("empty frame detected.  no sense in running the rest of the fastq splitter")
            break

def split_fasta(file_name_in, file_name_out, chunks):#split_count = 4):
    #modded to take in fixed chunks
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
    # At this point, we've already got the number of reads.
    #chunks = m.ceil(len(fasta_df) / split_count)  # how many sequences each split file will have
    print("total df length:", len(fasta_df))
    print("chunk size:", chunks)
    if (chunks < 1):
        print("split count too large. not enough info to split")
        chunks = 1
        print("new split count:", split_count)


    #for i in range(0, split_count):
    index_count = 0
    while(True):
        print("working on segment :", index_count + 1, "of FASTA splitter" )
        # fancy naming
        new_file_name = file_name_out + "_" + str(index_count) + ".fasta"
        
        # split file by selective selection, and writing
        start_index = int(index_count * chunks)
        end_index = int(((index_count + 1) * chunks))
        # if(chunks == 1):
        #    end_index += 1 #override on splits that only have 1
        index_count += 1
        if not (fasta_df.iloc[start_index:end_index, :].empty):
            fasta_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize=chunks, mode="w+", sep='\n', header=False)
        else:
            print("empty frame detected.  no sense in running the rest of FASTA splitter")
            break
    
if __name__ == "__main__":
    if len(sys.argv) == 4:
        input_file = sys.argv[1]
        output_name = sys.argv[2]
        split_count = int(sys.argv[3])
        if(split_count < 2):
            print("not splitting file, just moving it")
            split_count = 1
            
        input_extension = os.path.splitext(input_file)[1]
        if ((input_extension == ".fastq") or (input_extension ==".fq")):
            split_fastq(input_file, output_name, split_count)
        elif input_extension == ".fasta" or input_extension == ".ffn" or input_extension == ".fna" or input_extension == ".faa" or input_extension == ".fa":
            split_fasta(input_file, output_name, split_count)
        else:
            print("This is not a fastq or fasta.  Not splitting the file")
    else:
        print("only: ", len(sys.argv), "number of args.  not running")
        print("wrong number of args to file splitter.  not running")
    sys.exit()