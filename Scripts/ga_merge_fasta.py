import os
import sys
import time
from datetime import datetime as dt
import pandas as pd

#this code will merge all FASTAs with the same source 
#it's meant to condense all of the leftover FASTA reads from running them in different segments of chocophlan via BWA
#notes: overlaps in IDs are gauranteed to be the same data.  So we don't even need to check if a read ID has been used before

#dec 09, 2021: No, this isn't necessarily correct. We should be taking the intersection of reads. 
#slowdown is from BLAT essentially doing all the work twice.


def import_fastas(dir_list, segment_name):
    #collect all the reads into 1 dict.  Store the IDs separately.
    read_dict = dict()
    reads_list = list()
    sample_count = 0
    
    for item in dir_list:
        if(item.startswith(segment_name)):
            
            sample_read_list = list()
            ID = ""
            seq = ""
            file_to_open = os.path.join(in_dir, item)
            print("working on:", file_to_open)
            with open(file_to_open, "r") as in_fasta:
                
                for line in in_fasta:
                    if(line.startswith(">")):
                        if(ID == ""):
                            ID = line.strip("\n")
                        else:
                            #new line. store the seq
                            read_dict[ID] = seq
                            sample_read_list.append(ID)
                            ID = line.strip("\n")
                            seq = ""
                    else:
                        seq += line

            #store the last entry
            #sample_dict[ID] = seq
            read_dict[ID] = seq
            sample_read_list.append(ID)
            print("reads found:", len(sample_read_list))
            reads_list.append(sample_read_list)
            print(dt.today(), "main collection:", len(reads_list))
            
            sample_count += 1

    return read_dict, reads_list


if __name__ == "__main__":
    print(dt.today(), "Starting GA fasta merge")
    in_dir = sys.argv[1] #eg: GA_BWA/
    segment_name = sys.argv[2]  #eg: pair_1_00, contigs_1, singletons_9111
    out_dir = sys.argv[3]
    out_path = os.path.join(out_dir, segment_name + ".fasta")
    dir_list = os.listdir(in_dir)
    
    print(dt.today(), "starting import")
    read_dict, reads_list = import_fastas(dir_list, segment_name)
    print(dt.today(), "getting read intersection")
    for item in reads_list:
        print("number of reads:", len(item)) 
    intersected_reads = set.intersection(*map(set, reads_list)) #get the intersection of all reads
    print("number of intersected reads:", len(intersected_reads))


    print(dt.today(), "starting export")
    #export it                    
    with open(out_path, "w") as out_fasta:
        for item in intersected_reads:
            out_line = item + "\n"
            out_fasta.write(out_line)
            out_line = read_dict[item]
            out_fasta.write(out_line)
    
    print(dt.today(), "done!")