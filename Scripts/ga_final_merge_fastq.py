#use this to merge leftover fastq
#take in fastas, resolve, and get fastq
#produces 2 or 4 files




import os
import sys
import time
from datetime import datetime as dt
import pandas as pd


def import_raw_fastq(raw_dir, context):
    fastq_file = os.path.join(raw_dir, context + ".fastq")
    fastq_df = pd.read_csv(fastq_file, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df 

def import_fasta(dir, extension, context):
    read_set = set()
    dir_list = os.listdir(dir)
    for file_name in dir_list:
        if(file_name.endswith(extension)):
            if(file_name.startswith(context)):
                
                full_path = os.path.join(dir, file_name)
                print("full path:", full_path)
                with open(full_path, "r") as fasta_in:
                    id = ""
                    seq = ""
                    for raw_line in fasta_in:
                        line = raw_line.strip("\n")
                        if(line.startswith(">")):
                            id = line.strip(">")
                            read_set.add(id)
    return read_set


def merge_fasta(dir, extension, context):
    fasta_dict = dict()
    dir_list = os.listdir(dir)
    for file_name in dir_list:
        if(file_name.endswith(extension)):
            if(file_name.startswith(context)):
                
                full_path = os.path.join(dir, file_name)
                print("full path:", full_path)
                with open(full_path, "r") as fasta_in:
                    id = ""
                    seq = ""
                    for raw_line in fasta_in:
                        line = raw_line.strip("\n")
                        if(line.startswith(">")):
                            if(id != ""):
                                fasta_dict[id] = seq
                            seq = ""
                            id = line#.strip(">")
                        else:
                            seq += line
    return fasta_dict
             

if __name__ == "__main__":
    #designed to merge the leftover reads
    raw_dir = sys.argv[1]   #assemble-contigs
    dmd_dir = sys.argv[2]   #GA_dmd
    op_mode = sys.argv[3]
    export_dir = sys.argv[4]
    dir_dict = dict()
    #dir_dict["bwa"] = bwa_dir
    #dir_dict["blat"] = blat_dir
    dir_dict["dmd"] = dmd_dir
    fasta_dict = dict() #stores all the dicts of all contexts
    fastq_dict = dict()
    context_list = ["singletons"]
    if(op_mode == "pair" or op_mode == "paired"):
        context_list = ["singletons", "pair_1", "pair_2"]
        
    for context in context_list:
        
        fasta_dict[context] = import_fasta(dmd_dir, "fasta", context) 
            
    for context in context_list:
        fastq_dict[context] = import_raw_fastq(raw_dir, context)
    
    #start singletons
    s_reads = fasta_dict["singletons"]
    s_df = fastq_dict["singletons"]
    
    
    #resolve orphans
    #collect common reads. move to their own thing
    if(op_mode == "paired" or op_mode == "pair"):
        print(dt.today(), "handling paired")
        p1_reads = fasta_dict["pair_1"]
        p2_reads = fasta_dict["pair_2"]
        
        p1_only_reads = p1_reads - p2_reads
        p2_only_reads = p2_reads - p1_reads
        common_reads = p1_reads.intersection(p2_reads)
        
        p1_df = fastq_dict["pair_1"]
        p2_df = fastq_dict["pair_2"]
        p1_common_df = p1_df[p1_df["ID"].isin(common_reads)]
        p2_common_df = p2_df[p2_df["ID"].isin(common_reads)]
        
        p1_only_df = p1_df[~p1_df["ID"].isin(common_reads)]
        p2_only_df = p2_df[~p2_df["ID"].isin(common_reads)]
        
        p1_export_path = os.path.join(export_dir, "leftover_pair_1.fastq")
        p2_export_path = os.path.join(export_dir, "leftover_pair_2.fastq")
        p1_common_df["ID"] = "@" + p1_common_df["ID"]
        p2_common_df["ID"] = "@" + p2_common_df["ID"]
        
        p1_common_df.to_csv(p1_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
        p2_common_df.to_csv(p2_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
        
        s_reads.union(p1_only_reads)
        s_reads.union(p2_only_reads)
        
        s_df = pd.concat([s_df, p1_only_df, p2_only_df])
        
        
    
    #export s fastq
    print(dt.today(), "exporting singletons")
    s_export_df = s_df[s_df["ID"].isin(s_reads)]
    s_export_df["ID"] = "@" + s_export_df["ID"]
    s_export_path = os.path.join(export_dir, "leftover_singletons.fastq")
    s_export_df.to_csv(s_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
    
    
    #handle contigs
    contig_dict = merge_fasta(dmd_dir, "fasta", "contig") 
    
        
    contig_export_path = os.path.join(export_dir, "leftover_contigs.fasta")
   
    print(dt.today(), "exporting to:", contig_export_path)
    
    
    with open(contig_export_path, "w") as out_file:
        for item in contig_dict.keys():
            out_line = item + "\n"
            out_file.write(out_line)
            out_line = contig_dict[item] + "\n"
            out_file.write(out_line)
            