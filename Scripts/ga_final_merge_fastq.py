#use this to merge leftover fastq
#take in fastas, resolve, and get fastq
#produces 2 or 4 files

#mar 11, 2025: remake this code with multithreading.
#g

import os
import sys
import time
from datetime import datetime as dt
import pandas as pd
import multiprocessing  as mp
import threading as th
from math import ceil, log2
import queue
from concurrent.futures import ThreadPoolExecutor

class import_reads:
    def __init__ (self):
        self.set_dict = dict()
        self.lock = th.RLock()
        self.result_queue = queue.Queue()

    def merge_sets(self, set_1, set_2):
        return set_1.union(set_2)

    def import_single_fasta(self, fasta_path, set_list):
        read_set = set()

        with open(fasta_path, "r") as fasta_in:
            id = ""
            seq = ""
            for raw_line in fasta_in:
                line = raw_line.strip("\n")
                if(line.startswith(">")):
                    id = line.strip(">")
                    read_set.add(id)

        with self.lock:
            set_list.append(read_set)


    def import_sections(self, dir, extension, context):
        read_set = set()
        dir_list = os.listdir(dir)
        file_list = list()
        set_list = list()
        thread_list = list()
        
        
        for file_name in dir_list:
            if(file_name.endswith(extension)):
                if(file_name.startswith(context)):
                    full_path = os.path.join(dir, file_name)
                    file_list.append(full_path)
                    import_thread = th.Thread(target = self.import_single_fasta, args = (full_path, set_list))
                    thread_list.append(import_thread)
                    import_thread.start()

        for thread in thread_list:
            thread.join()
        


        merge_levels = ceil(log2(len(set_list)))
        
        #we've got all the sets generated.  now merge
        for level in range(merge_levels):
            thread_list = list()
            next_sets = list()
            results = [None] * (len(set_list) // 2 + len(set_list) % 2)
            for i in range(0, len(set_list) - 1, 2):
                thread_index = i // 2
                thread = th.Thread(
                    target = lambda index, set_1, set_2: results.__setitem__(index, self.merge_sets(set_1, set_2)),
                    args = (thread_index, set_list[i], set_list[i+1])
                )
                thread_list.append(thread)

                thread.start()
                #print(dt.today(), context, "thread[", str(i), "] launched")

            if(len(set_list) % 2 == 1):
                results[-1] = set_list[-1]

            for thread in thread_list:
                thread.join()
            set_list = results

        return set_list[0]

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
             

        

def import_raw_fastq(raw_dir, context, thread_lock, fastq_dict):
    fastq_file = os.path.join(raw_dir, context + ".fastq")
    fastq_df = pd.read_csv(fastq_file, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    #
    # 
    with thread_lock:
        fastq_dict[context] = fastq_df
    # return fastq_df 





if __name__ == "__main__":
    #designed to merge the leftover reads
    raw_dir = os.path.abspath(sys.argv[1])   #assemble-contigs
    dmd_dir = os.path.abspath(sys.argv[2])   #GA_dmd
    op_mode = sys.argv[3]
    export_dir = os.path.abspath(sys.argv[4])
    dir_dict = dict()
    #dir_dict["bwa"] = bwa_dir
    #dir_dict["blat"] = blat_dir
    dir_dict["dmd"] = dmd_dir
    ID_dict = dict() #stores all the dicts of all contexts
    fastq_dict = dict()
    import_obj = import_reads()
    thread_list = list()
    thread_lock = th.RLock()
    
    context_list = ["singletons"]
    if(op_mode == "pair" or op_mode == "paired"):
        context_list.append("pair_1")
        context_list.append("pair_2")
        
    for context in context_list:
        print(dt.today(), "start ID fetch:", context)
        ID_dict[context] = import_obj.import_sections(dmd_dir, "fasta", context) 
        print(dt.today(), "start fastq import:", context)
        fastq_thread = th.Thread(target = import_raw_fastq, args = (raw_dir, context, thread_lock, fastq_dict))
        fastq_thread.start()
        thread_list.append(fastq_thread)
        #print(dt.today(), "done fastq import:", context)
   
    for item in thread_list:
        item.join()
        
    print(dt.today(), "imports done")
    #start singletons
    s_reads = ID_dict["singletons"]
    s_df = fastq_dict["singletons"]
    
    
    #resolve orphans
    #collect common reads. move to their own thing
    if(op_mode == "paired" or op_mode == "pair"):
        print(dt.today(), "handling paired")
        p1_reads = ID_dict["pair_1"]
        p2_reads = ID_dict["pair_2"]
        
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
        print(dt.today(), "exporting paired")
        p1_common_df.to_csv(p1_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
        p2_common_df.to_csv(p2_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
        print(dt.today(), "done exporting paired")
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
            
    print(dt.today(), "done")