#!/usr/bin/env python
#what's changed?  DETECT-2 no longer needs to be filtered.  

import os.path
import sys
from datetime import datetime as dt
import multiprocessing as mp

def create_swissprot_map(SWISS_PROT_MAP):
    mapping_dict = {}
    final_mapping_dict = dict()
    with open(SWISS_PROT_MAP, "r") as mapping:
        for line in mapping.readlines():
            line_as_list = line.split("\t")
            mapping_dict[line_as_list[0]] = set(line_as_list[2:])
    
    
    for ec in mapping_dict.keys():
        gene_list = mapping_dict[ec]
        
        for item in gene_list:
            if(item in final_mapping_dict):
                final_mapping_dict[item].append(ec)
            else:
                final_mapping_dict[item] = [ec]
    
    return final_mapping_dict
    
#not used anymore. We don't filter DETECT-2
def filter_and_export_detect_predictions(detect_file, detect_ECs):
    with open(detect_file, "r") as topred:
        with open(detect_ECs, "w") as cutoff:
            for line in topred.readlines():
                line_as_list = line.split("\t")
                if line_as_list[2] == "probability":
                    continue
                if float(line_as_list[2]) >= 0.2 and int(line_as_list[3]) > 5:
                    cutoff.write(line)
                    
                    
def filter_and_export_priam_predictions(priam_file, priam_ECs):
    with open(priam_file, "r") as ECs:
        with open(priam_ECs, "w") as processedECs:
            ID = None
            EC = None
            for line in ECs.readlines():
                if line.startswith(">"):
                    ID = line.split(" ")[0][1:].strip("\n")
                    continue
                elif line == "\n":
                    continue
                elif ID and not line.startswith("#"):
                    EC = line.split(" ")[0]
                    processedECs.write("\t".join([ID, EC]))
                    processedECs.write("\n")
                    ID = None
                    EC = None
                    continue 


def filter_and_export_diamond_predictions(diamond_file, diamond_ECs):
    with open(diamond_file, "r") as blastout:
        with open(diamond_ECs, "w") as ecout:
            for line in blastout.readlines():
                line_as_list = line.strip().split("\t")
                for EC in mapping_dict:
                    if line_as_list[1] in mapping_dict[EC]:
                        ecout.write("\t".join([line_as_list[0], EC + "\n"]))


def import_detect_ec(detect_fbeta_file, gene_ec_dict):
    #DETECT-2's fbeta file gets left as-is.
    #gene_ec_dict = dict()
    with open(detect_fbeta_file, "r") as detect_fbeta:
        for line in detect_fbeta:
            list_line = line.split("\t")
            if(list_line[2] == "probability"):
                continue
            else:
                key = list_line[0]
                EC_val = list_line[1]
                if(key in gene_ec_dict):
                    gene_ec_dict[key] += EC_val
                else:
                    gene_ec_dict[key] = [EC_val]
    #return gene_ec_dict
    
def import_priam_ec(priam_sequence_ec, gene_ec_dict):
    #gene_ec_dict = dict()
    query_name = "None"
    ec_list = []
    with open(priam_sequence_ec, "r") as priam_ec:
        for line in priam_ec:
        
            #print(line)
            if(line == "\n"):
                continue
            
            elif(line.startswith(">")):
                #new line
                if(len(ec_list) > 0):
                    if(query_name is "None"):
                        print(dt.today(), "This shouldn't happen.  full EC list.  no query")
                    else:
                        #print("entry inserted:", query_name, ec_list)
                        gene_ec_dict[query_name] = ec_list
                query_name = line.strip(">")
                query_name = query_name.strip("\n")
            else:
                if not(line.startswith("#")):
                    list_line = line.split("\t")
                    ec = list_line[0]
                    probability = float(list_line[1])
                    if(probability >= 0.5):
                        ec_list.append(ec)
                        #print("INSERTED!", query_name, line)
            if(query_name is "None"):
                print(dt.today(), "This shouldn't be happening.  a line was skipped")
    #return gene_ec_dict
    
def import_diamond_ec(diamond_proteins_blastout, swissprot_map_dict, gene_length_dict, gene_ec_dict):
    #gene_ec_dict = dict()
    with open(diamond_proteins_blastout, "r") as diamond_ec:
        for item in diamond_ec:
            #note that this DIAMOND call has special parameters (for some reason)
            list_line = item.split("\t")
            #print(list_line)
            query_name = list_line[0].strip("\n")
            query_name = query_name.strip(">")
            swissprot_name = list_line[1]
            alignment_length = int(list_line[2])
            query_length = 0
            if(query_name in gene_length_dict):
                query_length = gene_length_dict[query_name]
                
            bitscore = float(list_line[8])
            coverage_percentage = 100 * 3 * alignment_length / query_length
            identity_percent = float(list_line[11].strip("\n"))
            
            if(query_length >= 100):
                if(bitscore >= 60):
                    #take this hit
                    if(swissprot_name in swissprot_map_dict):
                        EC_val = swissprot_map_dict[swissprot_name]
                        if(query_name in gene_ec_dict):
                            gene_ec_dict[query_name] += EC_val
                        else:
                            gene_ec_dict[query_name] = EC_val
                else:
                    EC_val = ["0.0.0.0"]
            else:
                if((identity_percent >= 85) and (coverage_percent >= 65)):
                    #take this hit
                    if(swissprot_name in swissprot_map_dict):
                        EC_val = swissprot_map_dict[swissprot_name]
                        if(query_name in gene_ec_dict):
                            gene_ec_dict[query_name] += EC_val
                        else:
                            gene_ec_dict[query_name] = EC_val
                else:
                    EC_val = ["0.0.0.0"]
        #return gene_ec_dict
                    
            
def import_gene_map(gene_map_file):
    gene_length_dict = dict()
    with open(gene_map_file, "r") as gene_map:
        for line in gene_map:
            list_line = line.split("\t")
            gene_name = list_line[0]
            gene_length = int(list_line[1])
            gene_length_dict[gene_name] = gene_length
    return gene_length_dict          
            
            
        

if __name__ == "__main__":
    detect_file     = sys.argv[1]
    priam_file      = sys.argv[2]
    diamond_file    = sys.argv[3]
    SWISS_PROT      = sys.argv[4]
    SWISS_PROT_MAP  = sys.argv[5]
    gene_map_file   = sys.argv[6]
    Output_file     = sys.argv[7]
    
    # print(detect_file)
    # print(priam_file)
    # print(diamond_file)
    # print(SWISS_PROT)
    # print(SWISS_PROT_MAP)
    # print(gene_map_file)
    # print(Output_file)
    
    #Input_Name = os.path.splitext(os.path.basename(Input_File))[0]
    # detect_dir = os.path.dirname(detect_file)
    # priam_dir = os.path.dirname(priam_file)
    # diamond_dir = os.path.dirname(diamond_file)
    
    # detect_ECs = os.path.join(detect_dir, Input_Name + ".toppred.cutoff")
    # priam_ECs = os.path.join(priam_dir, Input_Name + ".ECs")
    # diamond_ECs = os.path.join(diamond_dir, Input_Name + ".ECs")
    
    
    gene_length_dict = import_gene_map(gene_map_file)
    
    ec_process_list = []
    #-----------------------------------------
    # import the data
    manager = mp.Manager()
    diamond_ec_manager_dict = manager.dict()
    swissprot_map_dict = create_swissprot_map(SWISS_PROT_MAP)
    
    for item in swissprot_map_dict:
        print(item, swissprot_map_dict[item])
    diamond_ec_process = mp.Process(
        target = import_diamond_ec, 
        args = (diamond_file, swissprot_map_dict, gene_length_dict, diamond_ec_manager_dict)
    )
    diamond_ec_process.start()
    ec_process_list.append(diamond_ec_process)
    
    priam_ec_manager_dict = manager.dict()
    priam_ec_process = mp.Process(
        target = import_priam_ec, 
        args = (priam_file, priam_ec_manager_dict)
    )
    priam_ec_process.start()
    ec_process_list.append(priam_ec_process)
    
    detect_ec_manager_dict = manager.dict()
    detect_ec_process = mp.Process(
        target = import_detect_ec,
        args = (detect_file, detect_ec_manager_dict)
    )
    detect_ec_process.start()
    ec_process_list.append(detect_ec_process)
    
    for item in ec_process_list:
        item.join()
    ec_process_list[:] = []
    
    
    diamond_ec_dict = dict(diamond_ec_manager_dict)
    priam_ec_dict = dict(priam_ec_manager_dict)
    detect_ec_dict = dict(detect_ec_manager_dict)
    
    #print(priam_ec_dict)
    
    #--------------------------------------------------
    # merge the results
    diamond_keys = set(diamond_ec_dict.keys())
    priam_keys = set(priam_ec_dict.keys())
    common_keys = diamond_keys & priam_keys
    
    #for item in diamond_keys:
    #    print(item)
        
    #print("========================")
    #for item in priam_keys:
    #    print(item)
    # for item in priam_keys:
        # print(item, priam_ec_dict[item])
        
    # print("======================================================")
    # for item in detect_ec_dict.keys():
        # print(item, detect_ec_dict[item])
    
    #take the intersection of the ECs found for each gene between priam and diamond
    common_dict = dict()
    for item in common_keys:
        priam_ec_set = set(priam_ec_dict[item])
        diamond_ec_set = set(diamond_ec_dict[item])
        
        common_intersection_set = diamond_ec_set.intersection(priam_ec_set)
        
        common_dict[item]= list(common_intersection_set)
        #print("priam:", len(priam_ec_set), "diamond:", len(diamond_ec_set), "intersection:", len(common_intersection_set))
    
    for key in detect_ec_dict.keys():
    
        
        if(key in common_dict):
            detect_ec_list = detect_ec_dict[key]
            #print("before:", len(common_dict[key]))
            common_dict[key] += detect_ec_list
            common_dict[key] = list(set(common_dict[key])) #gets rid of dupes, and turns it back into a list
            #print("after:", len(common_dict[key]))
            #print("=-=========-=-=-=-=-==========")
        else:
            common_dict[key] = detect_ec_dict[key]
            #print("new:", len(common_dict[key]))
    #common_dict = sorted(common_dict)      
    
    #----------------------------------------------
    #export the final ec list
    with open(Output_file, "w") as ec_out:
        for item in sorted(common_dict.keys()):
            
            ec_list = common_dict[item]
            ec_string = ""
            for ec in ec_list:
                actual_ec = ec.strip(" ")
                ec_string += actual_ec + "|"
            ec_string = ec_string[:-1]
            ec_out.write(item + "\t" + str(len(ec_list)) + "\t" + ec_string + "\n")
    #export the final ec list
        
    
    
    
    
    
