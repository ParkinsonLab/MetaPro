#merges all genes.
#converts genes to proteins
#merges all proteins
import os
import sys
import time
from datetime import datetime as dt
import multiprocessing as mp
import psutil as psu
import math
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio.Seq import Seq


def convert_seq_wrapper(dict_in, dict_out, list_of_keys):
    for item in list_of_keys:
        seq = dict_in[item]
        dict_out[item] = convert_seq(seq)


def convert_seq(seq):
    #tried to make our own, but there's too many considerations. just use BioSeq's.  We've at least eliminated all the indexing and descriptions
    my_seq = Seq(seq)
    final_string = my_seq.translate(stop_symbol = "")
    return str(final_string)

def merge_fasta(dir, extension, context, dir_name, master_gene_dict):
    dir_list = os.listdir(dir)
    #print(dir_list)
    for file_name in dir_list:
    
        if(file_name.endswith(extension)):
            if(file_name.startswith(context)):
                #print(dt.today(), "working on:", dir_name, file_name)
                file_path = os.path.join(dir, file_name)
                
                with open(file_path, "r") as fasta_in:
                    id = ""
                    seq = ""
                    for line in fasta_in:
                        cleaned_line = line.strip("\n")
                        if(line.startswith(">")):
                            if(seq != ""):
                                
                                master_gene_dict[id] = seq
                                
                                #if(extension == ".faa"):
                                    #print(id, ":", seq)
                                    #time.sleep(2)
                                seq = ""
                            id = cleaned_line#.split(" ")[0]
                            #print("new id:", id)
                            #time.sleep(1)
                        else:
                            seq += cleaned_line
                    master_gene_dict[id] = seq
    #Scrub the dict:
    if '' in master_gene_dict.keys():
        print("empty string found")
        del master_gene_dict['']
                    
                    
            
                    
                    
                    
    
def export_fasta(fasta_dict, export_path):
    
    with open(export_path, "w") as out_file:
        
        for item in fasta_dict.keys():
            #print("writing:", [item])
            out_file.write(item + "\n")
            #print("writing:", fasta_dict[item])
            out_file.write(fasta_dict[item] + "\n")
            #time.sleep(1)

def mem_checker(threshold):
    #threshold is a percentage
    mem = psu.virtual_memory()
    available_mem = mem.available
    total_mem = mem.total
    
    available_pct = 100 * available_mem / total_mem
    
    if(float(available_pct) <= float(threshold)):
        return False
    else:
        return True
            

if __name__ == "__main__":
    bwa_dir = sys.argv[1]
    blat_dir = sys.argv[2]
    dmd_dir = sys.argv[3]
    export_dir = sys.argv[4]
    dir_dict = dict()
    dir_dict["bwa"] = bwa_dir
    dir_dict["blat"] = blat_dir
    dir_dict["dmd"] = dmd_dir
    
    context_list = ["pair", "singleton", "contigs"]
    gene_list = ["bwa", "blat"]
    
    master_gene_dict = dict()
    print(dt.today(), "starting")
    for context in context_list:
        for item in gene_list:
            print(dt.today(), "working on:", context, item)
            merge_fasta(dir_dict[item], ".fna", context, item, master_gene_dict)
            
    #sys.exit("stop")
    print(dt.today(), "number of genes:", len(master_gene_dict.keys()))
    key_list = list(master_gene_dict.keys())
    #print(">gi|421913313|ref|NZ_CANR01000173.1|:66-995|573|g__Klebsiella.s__Klebsiella_pneumoniae|UniRef90_K4S416|UniRef50_K4S416")
    #print(master_gene_dict[">gi|421913313|ref|NZ_CANR01000173.1|:66-995|573|g__Klebsiella.s__Klebsiella_pneumoniae|UniRef90_K4S416|UniRef50_K4S416"])
    
    master_prot_dict = dict()
    manager = mp.Manager()
    mgr_prot_dict = manager.dict()
    
    for context in context_list:
        merge_fasta(dir_dict["dmd"], ".faa", context, "dmd", master_prot_dict)
    print(dt.today(), "DMD proteins", len(master_prot_dict.keys()))
    
    work_divider = math.ceil(mp.cpu_count() * 0.9)
    index_offset = math.ceil(len(key_list)/work_divider)
    print("key list:", len(key_list))
    print("offset used:", index_offset)
    process_list = list()
    for i in range(0, work_divider):
        job_submitted = False
        while(job_submitted == False):
            if(mem_checker(5)):
                if(i == work_divider-1):
                    #print(dt.today(), "final launch:", i*index_offset, "->", (i+1)*index_offset, "but really, :")
                    process = mp.Process(target = convert_seq_wrapper, args = (master_gene_dict, mgr_prot_dict, key_list[i*index_offset:]))
                    
                    process.start()
                    job_submitted = True
                else:
                    #print(dt.today(), "launching:", i*index_offset, "->", (i+1)*index_offset)
                    process = mp.Process(target = convert_seq_wrapper, args = (master_gene_dict, mgr_prot_dict, key_list[i*index_offset:(i+1)*index_offset]))
                    process.start()
                    job_submitted = True
            else:
                time.sleep(5)
                print(dt.today(), "delaying conversion job:", i*index_offset, "->", (i+1)*index_offset, end = "\r")
        process_list.append(process)
        
        
    print(dt.today(), "all processes launched")
    for item in process_list:
        item.join()
    print(dt.today(), "all converters finished")
    prot_dict = dict(mgr_prot_dict)
    master_prot_dict.update(prot_dict)
    
    print(dt.today(), "final prot keys:", len(master_prot_dict.keys()))
    #for item in master_gene_dict.keys():
    #    master_prot_dict[item] = convert_seq(master_gene_dict[item])
    #print(dt.today(), "number of prots:", len(master_prot_dict.keys()))
    
    export_path = os.path.join(export_dir, "all_proteins.faa")
    export_fasta(master_prot_dict, export_path)
    print(dt.today(), "exported")
    
    
    