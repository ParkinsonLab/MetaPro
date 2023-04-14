#This filters our genes that didn't quite make it after the conversion of contig-segments to reads made them disappear.
#it'll remove the genes that didn't make it.

#use the final_gene_map

import sys
import pandas as pd
import time

def import_gene_map(gene_map_file):
    gene_list = []
    with open(gene_map_file, "r") as gene_map_in:
        skip_header = True
        for line in gene_map_in:
            if(skip_header):
                skip_header = False
                continue
            else:    
                cleaned_line = line.strip("\n")
                line_split = cleaned_line.split("\t")
                gene_name = line_split[0]
                gene_list.append(gene_name)
    return gene_list
     
if __name__ == "__main__":
    EC_file = sys.argv[1]       #protein_ECs.all
    gene_map_file = sys.argv[2] #final_gene_map.tsv
    export_ec_file = sys.argv[3]
    
    gene_list = import_gene_map(gene_map_file)
    
    #sys.exit("parser")
    
    EC_df = pd.read_csv(EC_file, sep = "\t", header = None, names = ["gene", "EC_count", "EC"])
    #gene_map_df = pd.read_csv(gene_map_file, sep = "\t", usecols = ["gene", "length", "count"])
    
    
    #print(gene_map_df)
    #print(gene_map_df["gene"])
    print(EC_df["gene"])
    EC_df = EC_df[EC_df["gene"].isin(gene_list)]
    #print(EC_df)
    
    EC_df.to_csv(export_ec_file, sep = "\t", header = None, index = None)