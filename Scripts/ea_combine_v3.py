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
                    #if(key == "BAH62832.1"):
                    #    print(dt.today(), key, "key already present.  appending", EC_val)
                    
                    new_list = gene_ec_dict[key] #Apparently manager dicts can't handle appends.  this is eye-opening....
                    new_list.append(EC_val)
                    gene_ec_dict[key] = new_list
                    #if(key == "BAH62832.1"):
                    #    print(key, gene_ec_dict[key])
                else:
                    #if(key == "BAH62832.1"):
                    #    print(dt.today(), "new gene.  adding")
                    gene_ec_dict[key] = [EC_val]
                    #if(key == "BAH62832.1"):
                    #    print(key, gene_ec_dict[key])
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
                        ec_list[:] = []
                        
                query_name = line.strip(">")
                query_name = query_name.strip("\n")
                query_list = query_name.split(" ")
                query_name = query_list[0]
                
                #print("query:", query_name)
            else:
                if not(line.startswith("#")):
                    #if(query_name == "BAH62832.1"):
                    #    print("line is valid:", line)
                    list_line = line.split("\t")
                    ec = list_line[0].strip("\n")
                    ec = ec.strip(" ")
                    probability = float(list_line[1])
                    if(probability >= 0.5):
                        ec_list.append(ec)
                        #if(query_name == "BAH62832.1"):
                        #    print("INSERTED!", query_name, line)
                    else:
                        #if(query_name == "BAH62832.1"):
                        print("line valid but not inserted")
                #else:
                    #if(query_name == "BAH62832.1"):
                    #print("this line is skipped:", line)
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
            coverage_percent = 100 * 3 * alignment_length / query_length
            identity_percent = float(list_line[11].strip("\n"))
            
            if(query_length >= 100):
                if(bitscore >= 60):
                    #take this hit
                    if(swissprot_name in swissprot_map_dict):
                        EC_val = swissprot_map_dict[swissprot_name]
                        if(query_name in gene_ec_dict):
                            old_entry = gene_ec_dict[query_name]
                            old_entry += EC_val
                            gene_ec_dict[query_name] = old_entry
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
            
            
def debug_export(ec_dict, file_name):
    #used specifically to debug.  exports the dictionary
    with open(file_name, "w") as debug_out:
        for key in ec_dict:
            line_to_write = str(key)
            ec_list = ec_dict[key]
            for item in ec_list:
                line_to_write += "\t" + str(item)
            line_to_write += "\n"
            debug_out.write(line_to_write)

if __name__ == "__main__":
    detect_file     = sys.argv[1]
    priam_file      = sys.argv[2]
    diamond_file    = sys.argv[3]
    SWISS_PROT      = sys.argv[4]
    SWISS_PROT_MAP  = sys.argv[5]
    gene_map_file   = sys.argv[6]
    Output_file     = sys.argv[7]
    #tester_file     = sys.argv[8]
    
    # priam_ec_dict = dict()
    # import_priam_ec(priam_file, priam_ec_dict)
    
    #for item in priam_ec_dict:
    #    print(item, priam_ec_dict[item])
    
    
    gene_length_dict = import_gene_map(gene_map_file)
    
    ec_process_list = []
    #-----------------------------------------
    # import the data
    manager = mp.Manager()
    diamond_ec_manager_dict = manager.dict()
    swissprot_map_dict = create_swissprot_map(SWISS_PROT_MAP)
    
    #for item in swissprot_map_dict:
    #   print(item, swissprot_map_dict[item])
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
    
    print(dt.today(), "waiting for import jobs to finish")
    for item in ec_process_list:
        item.join()
    ec_process_list[:] = []
    print(dt.today(), "import jobs finished")
    
    
    
    
    diamond_ec_dict = dict(diamond_ec_manager_dict) #key: gene.  value: list of ECs
    priam_ec_dict = dict(priam_ec_manager_dict)     #key: gene.  value: list of ECs
    detect_ec_dict = dict(detect_ec_manager_dict)   #key: gene.  value: list of ECs
    print(dt.today(), "finished converting dicts")
    
    
    #--------------------------------------------------
    # merge the results.  get the genes for which both tools came back with an EC
    diamond_keys = set(diamond_ec_dict.keys())
    priam_keys = set(priam_ec_dict.keys())
    common_keys = diamond_keys & priam_keys
    
    print(dt.today(), "finished getting common genes")
    #take the intersection of the ECs found for each gene between priam and diamond
    common_dict = dict()
    for item in common_keys:
        priam_ec_set = set(priam_ec_dict[item])
        diamond_ec_set = set(diamond_ec_dict[item])
        combined_ec = list(diamond_ec_set.intersection(priam_ec_set))
        common_dict[item] = combined_ec
    
    
    
    for key in detect_ec_dict.keys():
    #then take all of DETECT's results
        if(key in common_dict):
            detect_ec_list = detect_ec_dict[key]
            old_entry = common_dict[key]
            old_entry += detect_ec_list
            common_dict[key] = old_entry
            common_dict[key] = list(set(common_dict[key])) #gets rid of dupes, and turns it back into a list
        else:
            common_dict[key] = detect_ec_dict[key]
     
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
        
    
    
    
    
    
