#!/usr/bin/env python
#This is supposed to figure out what's leftover as pair 1, pair 2, contigs, and orphans
#It was originally named Map read contigs
# code seems to:
# 1) sift through the SAM file, and divides it into 2 sets:  stuff that's been mapped, and stuff that's not mapped
# for a significant performance gain, this code can be redesigned to run only 1 type of fastq.  
# There's nothing really stopping it, besides the formation of the contig-components map, which can be separated
import os.path
import sys
import pandas as pd

if __name__ == "__main__":
    read_path           = sys.argv[1]   #fastq
    sam_path            = sys.argv[2]   #sam
    output_path         = sys.argv[3]   #just a location
    
    print("Reads:", read_path)
    print("SAM:", sam_path)
    print("Output Location:", output_path)
    
    #imports all fastq files into df
    read_df = pd.read_csv(read_path, header=None, names=[None], sep ='\n', skip_blank_lines = False)
    read_df = pd.DataFrame(read_df.values.reshape(int(len(read_df) / 4), 4))
    
    read_df.columns = ["ID", "sequence", "junk", "quality"]
    
    
    #import the SAM file
    read_sam_df = pd.read_csv(sam_path, error_bad_lines=False, header=None, sep="\t")
    read_sam_df.iloc[:, 1] = read_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])

    #select mapped and unmapped slices from reads
    mapped_read_sam_df = "@" + read_sam_df.loc[read_sam_df.iloc[:, 1] == "0"].iloc[:, 0]          # grab the first column of all the rows in the original sam df with flag 0.
                                                                                                        # -> first column is IDs
    mapped_read_sam_df = mapped_read_sam_df.drop_duplicates()                                       # remove all the duplicates (this may happen)
    mapped_read_sam_df.columns = ["ID"]                                                               # label the single column to be "ID"
   
    # 1 is unmapped, 0 is mapped
    unmapped_read_sam_df = "@" + read_sam_df.loc[read_sam_df.iloc[:, 1] == "1"].iloc[:, 0]        # do the same for unmapped -> flag 1
    unmapped_read_sam_df = unmapped_read_sam_df.drop_duplicates()
    unmapped_read_sam_df.columns = ["ID"]
    
    unmapped_read_sam_df = unmapped_read_sam_df[~unmapped_read_sam_df.isin(mapped_read_sam_df)] # then, take only the keys found uniquely in unmapped
    
    #------------------------------
    #write it
    read_name = os.path.basename(read_path)
    # export only the rows that are unmapped
    read_df[read_df.ID.isin(unmapped_read_sam_df)].to_csv(os.path.join(output_path, read_name), sep='\n', mode ="w+", header=False, index=False)
    # then export only the rows that are mapped
    read_df[read_df.ID.isin(mapped_read_sam_df)].to_csv(os.path.join(output_path, os.path.splitext(read_name)[0] + "_mapped" + os.path.splitext(read_name)[1]), sep='\n', mode ="w+", header=False, index=False)
