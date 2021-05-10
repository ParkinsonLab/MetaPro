#!/usr/bin/env python
#this code merges all of the WEVOTE input files into 1 neat table, so WEVOTE can use items
#Why a V2?  because we need it to accept more than the prescribed number of input files originally laid out.
#We also wanted to make it as explicit as possible, so it's easy to understand.

import sys
import os
from datetime import datetime as dt

def import_classifications(input_file, contig2read_map):
    #list of dictionaries
    Input_classification = {}
    with open(input_file, "r") as Input:
        for line in Input:
            columns = line.split("\t")
            if columns[0] == "score":
                continue
            Seq_ID = columns[1]
            Tax_ID = columns[2].strip("\n")
            if Seq_ID in contig2read_map:
                for Read_ID in contig2read_map[Seq_ID]:
                    if Read_ID not in Input_classification:
                        Input_classification[Read_ID] = Tax_ID
            else:
                Input_classification[Seq_ID] = Tax_ID
        #Input_classifications.append(Input_classification)
        return Input_classification


def combine_results(Input_classifications, ):
    #combine everything.
    Combined_Classification = {}
    for Index, Classification in enumerate(Input_classifications):
        Combined_set = set(Combined_Classification)
        Classification_set = set(Classification)   
        for Read in Combined_set.intersection(Classification_set):
            Combined_Classification[Read].append(Classification[Read])
            
        for Read in Combined_set.symmetric_difference(Classification_set):
            if Read in Combined_Classification:
                Combined_Classification[Read].append("0")
            else:
                if Index > 0:
                    Combined_Classification[Read] = ["0"]
                    for x in range(Index - 1):
                        Combined_Classification[Read].append("0")
                    Combined_Classification[Read].append(Classification[Read])
                else:
                    Combined_Classification[Read] = [Classification[Read]]

    return Combined_Classification


def import_contig2read_map(file_name):
    #import the contig2read map
    contig2read_map = {}
    with open(file_name, "r") as mapping:
        for line in mapping:
            if len(line) > 5:
                entry = line.strip("\n").split("\t")
                contig2read_map[entry[0]] = entry[2:]

    return contig2read_map

def export_wevote_result(Combined_Classification, output_file):
    with open(output_file, "w") as out:
        for Read in Combined_Classification:
            out.write(Read + "," + ",".join(Combined_Classification[Read]) + "\n")

if __name__ == "__main__":
    #why did we input the same ga_taxon file twice? so we know it's being used 
    #twice, as a hack to give it a heavier weight to WEVOTE, yet still allow this code to change
    #if we wanted to add a truly separate taxa result to the collection.
    
    contig2read_file = sys.argv[1]
    output_file = sys.argv[2]
    ga_taxon_results_0 = sys.argv[3] 
    ga_taxon_results_1 = sys.argv[4]
    ga_taxon_results_2 = sys.argv[5]
    kaiju_results = sys.argv[6]
    centrifuge_results = sys.argv[7]
    
    #import the contig map
    contig2read_map = dict()
    if(os.path.exists(contig2read_file)):
        contig2read_map = import_contig2read_map(contig2read_file)
    else:
        print(dt.today(), "no contigs to consider.  continuing")
    
    #import each classification result
    Input_classifications = []
    if(os.path.exists(ga_taxon_results_0)):
        Input_classifications.append(import_classifications(ga_taxon_results_0, contig2read_map))
    if(os.path.exists(ga_taxon_results_1)):
        Input_classifications.append(import_classifications(ga_taxon_results_1, contig2read_map))
    if(os.path.exists(ga_taxon_results_2)):
        Input_classifications.append(import_classifications(ga_taxon_results_2, contig2read_map))
    Input_classifications.append(import_classifications(kaiju_results, contig2read_map))
    Input_classifications.append(import_classifications(centrifuge_results, contig2read_map))

    #combine them all for Wevote
    Combined_Classification = combine_results(Input_classifications)
    
    #export it
    export_wevote_result(Combined_Classification, output_file)
    