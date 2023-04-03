#parse the wevote results. get the taxa. get the class. pull the libs from our chocophlan
#march 10, 2023:  install a threshold for reads.  There is noise in taxa.
#only include taxa for which there is a 1%< occurence

import os
import sys
import time 
from datetime import datetime as dt

def import_nodes(nodes_file):
    nodes_dict = dict()
    with open(nodes_file, "r") as nodes_in:
        for line in nodes_in:
            cleaned_line = line.strip("\n").split("\t")
            taxid = cleaned_line[0]
            rank = cleaned_line[4]
            nodes_dict[taxid] = rank
            
    return nodes_dict

def import_wevote(wevote_file, existence_percent):
    #filter the taxa. only include taxa that reps 1% or more
    read_counter = 0
    taxa_tally_dict = dict()
    unique_taxa = set()
    tally_dict = dict()
    
    #read the file
    with open(wevote_file, "r") as wevote_in:
        for line in wevote_in:
            read_counter += 1
            cleaned_line = line.strip("\n").split("\t")
            final_taxa = cleaned_line[-1]
            #unique_taxa.add(final_taxa)
            if(final_taxa in taxa_tally_dict):
                taxa_tally_dict[final_taxa] += 1
            else:
                taxa_tally_dict[final_taxa] = 1
                

    #figure out percentages
    for taxa in taxa_tally_dict:
        rep_val = taxa_tally_dict[taxa] * 100/ read_counter
        if(rep_val >= existence_percent):
            unique_taxa.add(taxa)
            tally_dict[taxa] = rep_val
    return unique_taxa, tally_dict
    
def import_taxa_class_map(taxa_class_map_path):
    taxa_class_dict = dict()
    
    with open(taxa_class_map_path, "r") as taxa_class_in:
        for line in taxa_class_in:
            cleaned_line = line.strip("\n").split("\t")
            taxa = cleaned_line[0]
            class_level = cleaned_line[1].split("|")
            class_taxa = class_level[1]
            taxa_class_dict[taxa] = class_taxa
    return taxa_class_dict

def export_lines(lib_file_path, nodes_dict, class_item, yes_count, no_count):
    if("1236" in lib_file_path):
        print(lib_file_path)
    exist_flag = "no"
    if(os.path.exists(lib_file_path)):
        exist_flag = "yes"
        yes_count += 1
    else:
        no_count += 1
        
    out_file.write(exist_flag +"|"+nodes_dict[class_item] + "|" +  class_item + ".fasta" + "|" + lib_file_path + "\n")



if __name__ == "__main__":
    wevote_file_path = sys.argv[1]
    taxa_class_map_path = sys.argv[2]
    nodes_file = sys.argv[3]
    export_lib_file = sys.argv[4]
    reject_lib_file = sys.argv[5]
    lib_root_path = sys.argv[6]
    exist_percent = float(sys.argv[7])
    
    unique_taxa, tally_dict = import_wevote(wevote_file_path, exist_percent)
    taxa_class_dict = import_taxa_class_map(taxa_class_map_path)
    nodes_dict = import_nodes(nodes_file)
    
    for item in tally_dict:
        print(item, tally_dict[item])
    
    class_set = set()
    reject_set = set()
    for item in unique_taxa:
        try:
            class_set.add(taxa_class_dict[item])
        except KeyError:
            reject_set.add(item)
            
    yes_count = 0
    no_count = 0
    no_find_count = 0
    with open(export_lib_file, "w") as out_file:
        for class_taxa_item in sorted([int(i) for i in class_set]):
            try:
                class_item = str(class_taxa_item)
                lib_file_path = os.path.join(lib_root_path, class_item + ".fasta")
                if(class_item == "1236" or class_item == "1760" or class_item == "28211"):
                    #bypasser for our specific split libs for these taxa
                    for i in range(0, 3):
                        lib_file_path = os.path.join(lib_root_path, class_item + "_" + str(i) + ".fasta")
                        export_lines(lib_file_path, nodes_dict, class_item, yes_count, no_count)
                        
                elif(class_item == "91061"):
                    #bypasser for our specific split libs for these taxa
                    for i in range(0, 2):
                        lib_file_path = os.path.join(lib_root_path, class_item + "_" + str(i) + ".fasta")
                        export_lines(lib_file_path, nodes_dict, class_item, yes_count, no_count)
                else:
                    export_lines(lib_file_path, nodes_dict, class_item, yes_count, no_count)
                        
                
                
            except KeyError:
                no_find_count += 1
                out_file.write("can't find" + "|" + class_item + ".fasta" + "\n")
    
    with open(reject_lib_file, "w") as out_file:
        
        for reject_taxa_item in sorted([int(i) for i in reject_set]):
            reject_item = str(reject_taxa_item)
            try:
                out_file.write(nodes_dict[reject_item] + "|" + reject_item + "\n")
            except KeyError:
                out_file.write("can't find" + "|" + reject_item + "\n")
                
    print("no:", no_count, "| yes:", yes_count, "| no-find:", no_find_count)
    