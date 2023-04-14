import os
import sys
import time 
from datetime import datetime as dt
import multiprocessing as mp
import threading as th


#merges gene maps. removes/resolves dupes


def merge_map(file_dir, master_dict, context):
    #will auto-resolve paired maps too if left on pair_
    #master_dict = dict()
    map_count = 0
    gene_set = set()
    map_list = os.listdir(file_dir)
    for file_name in map_list:
        if(file_name.endswith("_gene_map.tsv")):
            if(file_name.startswith(context)):
                full_path = os.path.join(file_dir, file_name)
                #print(dt.today(), "opening:", full_path)
                map_count += 1
                line_count = 0
                with open(full_path, "r") as map_in:
                    for line in map_in:
                        inner_dict = dict()
    
                        line_count += 1
                        line_split = line.strip("\n").split("\t")
                        gene_name = line_split[0]
                        gene_length = line_split[1]
                        
                        read_list = line_split[3:]
                        for read in read_list:
                            read_split = read.split("<")
                            read_id = read_split[0]
                            score = read_split[1].split(">")[1]
                            inner_dict["score"] = score
                            inner_dict["gene"] = gene_name
                            #gene_set.add(inner_dict["gene"])
                            inner_dict["gene_length"] = gene_length
                            #print("Score:", score)
                            #print("gene name", gene_name)
                            #print("length:", gene_length)
                            try:
                                read_entry = master_dict[read_id]
                                if(read_entry["score"] < score):
                                    master_dict[read_id] = inner_dict
                                    #print(dt.today(), context, "overwriting read:", read_id)
                                    #time.sleep(10)
                                
                            except KeyError:
                                master_dict[read_id] = inner_dict
                                #print(dt.today(), context, "inserted new read:", read_id)
                                #print(master_dict[read_id]["gene"])
                    #print("line count:", line_count)
                                
                print("number of maps scanned:", map_count, end = "\r")
    print("number of maps scanned:", map_count)
    
    #print("number of unique genes:", len(gene_set))
    return master_dict
    
def transform_map(master_dict):
    #reorients the data into gene-based index
    new_map_dict = dict()
    
    read_count = 0
    gene_name_set = set()
    for read in master_dict.keys():
        #print("read:", read)
        read_count += 1
        read_data = master_dict[read]
        gene_name = read_data["gene"]
        gene_name_set.add(gene_name)
        gene_length = read_data["gene_length"]
        
        try:
            old_entry = new_map_dict[gene_name]
            old_entry["reads"].add(read)
            new_map_dict[gene_name] = old_entry
            #print("updating:", gene_name, new_map_dict[gene_name]["gene_length"], len(new_map_dict[gene_name]["reads"]))
            #time.sleep(1)
        except KeyError:
            inner_dict = dict()
            inner_dict["reads"] = set([read])
            inner_dict["gene_length"] = gene_length
            new_map_dict[gene_name] = inner_dict
            #print("new gene:", gene_name, new_map_dict[gene_name]["gene_length"], len(new_map_dict[gene_name]["reads"]))
    
    #for item in new_map_dict:
    #    print(item, new_map_dict[item])
    
    return new_map_dict
    
def export_map(master_dict, export_path):
    with open(export_path, "w") as out_file:
        for item in master_dict.keys():
            #print("export gene:", item, master_dict[item]["gene_length"], len(master_dict[item]["reads"]))
            inner_data = master_dict[item]
            out_line = item + "\t" + str(inner_data["gene_length"]) + "\t" + str(len(inner_data["reads"]))
            for read in inner_data["reads"]:
                out_line += "\t" + read
                
            out_line += "\n"
            out_file.write(out_line)
            
            
    
                                    

if __name__ == "__main__":
    
    bwa_dir = sys.argv[1]
    blat_dir = sys.argv[2]
    dmd_dir = sys.argv[3]
    export_dir = sys.argv[4]
    
    dir_dict = dict()
    dir_dict["bwa"] = bwa_dir
    dir_dict["blat"] = blat_dir
    dir_dict["dmd"] = dmd_dir
    print(dt.today(), "start map merge")
    manager = mp.Manager()
    mgr_read_dict = manager.dict()
    thread_list = list()
    
    command_args = ["pair_", "singletons", "contigs"]
    """
    for dir in dir_dict.keys():
        for item in command_args:
            thread = th.Thread(target = merge_map, args = (dir_dict[dir], mgr_read_dict, item))
            thread.start()
            thread_list.append(thread)
    
    for item in thread_list:
        item.join()
    read_dict = dict(mgr_read_dict)
    """
    
    read_dict = dict()
    for dir in dir_dict.keys():
        for item in command_args:
            print(dir)
            merge_map(dir_dict[dir], read_dict, item)
    
    
    
    print("number of reads:", len(read_dict.keys()))
    
    
    
    print(dt.today(), "start transform")
    
    gene_dict = transform_map(read_dict)
    print("number of genes:", len(gene_dict.keys()))
    
    #for item in gene_dict:
    #    print(item, gene_dict[item])
    #sys.exit("stop")
    
    print(dt.today(), "start export")
    export_path = os.path.join(export_dir, "gene_map.tsv")
    export_map(gene_dict, export_path)
    print(dt.today(), "end")