#Makes the new taxa table summary.
#takes in taxa_classifications from taxonomic_annotations/final_results
#and the contig read count from BWA-ing raw data onto contig segments.

#tallies all the new data.
import os
import sys
from datetime import datetime as dt
import multiprocessing as mp
import pandas as pd

def import_names(names_file):
    #makes 2 dicts: names to node, node to name
    #names = {}
    names_1 = dict()
    names_2 = dict()
    with open(names_file, "r") as infile:
        for line in infile:
            cols = line.split("\t|\t")
            taxid = cols[0]
            name = cols[1]
            if "scientific name" in cols[3]:
                names_1[taxid] = name
                names_2[name] = taxid
        else:
            names_1[0] = "unclassified"
            names_2["unclassified"]  = 0
    return names_1, names_2

def import_nodes(nodes_file):
    nodes = dict()
    with open(nodes_file, "r") as infile:
        for line in infile:
            cols = line.split("\t|\t")
            taxid = cols[0]
            parent = cols[1]
            rank = cols[2]
            nodes[taxid] = (parent, rank)
        else:
            nodes["0"] = ("0", "unclassified")
    return nodes
            
            
def import_names_classification(taxa_class_file):
    #we're using the constrain classification file.  no need to import names/nodes, or rebuild taxa trees
    read_taxa_dict = dict()
    with open(taxa_class_file, "r") as taxa_class_report:
        for line in taxa_class_report:
            cleaned_line = line.strip("\n")
            line_list = cleaned_line.split("\t")
            if(line_list[0] == "U"):
                continue
            else:
                read = line_list[1]
                taxa_tree = line_list[2]
                if(read in read_taxa_dict):
                    print(dt.today(), "multiple taxa per read.  this is bad")
                    sys.exit("BAD")
                else:
                    read_taxa_dict[read] = taxa_tree
    return read_taxa_dict
            
def import_contig_read_count(contig_read_count_file):
    contig_read_count_dict = dict()
    with open(contig_read_count_file, "r") as contig_read_count_report:
        for line in contig_read_count_report:
            cleaned_line = line.strip("\n")
            line_list = cleaned_line.split("\t")
            contig_name = line_list[0]
            read_count = float(line_list[1])  
            contig_read_count_dict[contig_name] = read_count
    return contig_read_count_dict

def import_read_in_contig_v2(read_in_contig_file):
    #takes in the new, better contig_read file, which is a lookup yes/no for each read
    read_contig_lookup_dict = dict()
    with open(read_in_contig_file, "r") as read_in_contig_report:
        for line in read_in_contig_report:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            read = line_split[0]
            is_contig = line_split[1]
            read_contig_lookup_dict[read] = is_contig
    return read_contig_lookup_dict
    
def convert_read_count(read, contig_read_count_dict, read_contig_lookup_dict):
    #convert the read into a count.  If it's in a contig, take the contig-counted number
    read_count = 1
    if(read.startswith("gene")):
        if(read in contig_read_count_dict):
            read_count = contig_read_count_dict[read]
        else:
            read_count = 0
    else:
        is_contig = read_contig_lookup_dict[read]
        if(is_contig == "yes"):
            read_count = 0
        
    return read_count
    
def import_life_list(life_list_file):
    life_dict = dict()
    with open(life_list_file, "r") as life_report:
        for line in life_report:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            category = line_split[0]
            names = line_split[1:]
            for name in names:
                life_dict[name] = category
    return life_dict      
    
def convert_taxa(x, life_dict):
    taxa_part = x.split(";")
    taxa_part.reverse()
    
    category = "unclassified"
    
    for item in taxa_part:
        if(item == "root"):
            break
        taxa_name = item.split("_")[1]
        if(taxa_name in life_dict):
            category = life_dict[taxa_name]
            break
    
    
        
    return category

if __name__ == "__main__":
    #names_file = sys.argv[1]
    #nodes_file = sys.argv[2]
    names_class_file = sys.argv[1] #constrain classification
    contig_read_count_file = sys.argv[2] #the contig read count file
    read_in_contig_file = sys.argv[3] # the better contig file
    #life_list_file = sys.argv[4]
    #export_report = sys.argv[5]
    export_report = sys.argv[4]

    print(dt.today(), "importing data")
    #life_dict = import_life_list(life_list_file)
    read_taxa_dict = import_names_classification(names_class_file)
    contig_read_count_dict = import_contig_read_count(contig_read_count_file)
    read_in_contig_dict = import_read_in_contig_v2(read_in_contig_file)
    
    print(dt.today(), "final assembly")
    
    
    final_df = pd.DataFrame.from_dict(read_taxa_dict, orient = "index")
    final_df.columns = ["taxa"]
    final_df["read"] = final_df.index
    final_df.reset_index(inplace = True, drop = True)
    #transform the reads into read counts, then groupby, sum, and provide the total read count per taxa
    #final_df["taxa"] = final_df["taxa"].apply(lambda x: convert_taxa(x, life_dict))
    final_df["read"] = final_df["read"].apply(lambda x: convert_read_count(x, contig_read_count_dict, read_in_contig_dict))
    final_df = final_df.groupby(["taxa"]).sum()
    
    final_df.to_csv(export_report, sep = "\t")
    
    #print(final_df)
    print(dt.today(), "taxa report: done")
    
    