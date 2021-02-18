#!/usr/bin/env python
#what's changed?  DETECT-2 no longer needs to be filtered.  
#we now use the ENZYME database to filter out bad multi-EC annotations.
#in v5: it now informs where the annotation came from
import os.path
import sys
from datetime import datetime as dt
import multiprocessing as mp
import time 
class ec:
    def __init__(self, ec, prob):
        self.ec = ec
        self.prob = prob
        

def import_freq_file(freq_file):
    #import this multi-EC protein frequency file from Nirvana
    freq_dict = dict()
    double_count = 0
    
    with open (freq_file, "r") as freq_report:
        skip_first_line = True
        for line in freq_report:
            if(skip_first_line):
                skip_first_line = False
                continue
            else:   
                cleaned_line = line.strip("\n")
                line_split = cleaned_line.split("\t")
                EC_split = line_split[0].split("_")
                EC_0 = EC_split[0]
                EC_1 = EC_split[1]
                if("n" in EC_0):
                    continue
                if("n" in EC_1):
                    continue
                #print("EC 0:", EC_0, "EC 1:", EC_1)
                protein_count = line_split[1]
                ec_set = set([EC_0])
                ec_set.add(EC_1)
                ec_set = sorted(ec_set)
                ec_key = ec_set[0] + "_" + ec_set[1]
                #print("ec set:", ec_set)
                #print("ec key:", ec_key)
                #time.sleep(1)
                if(ec_key in freq_dict):
                    double_count += 1
                    print(dt.today(), "we've got a double.  taking whichever's bigger")
                    old_count = freq_dict[ec_key]
                    if(old_count >= protein_count):
                        continue
                    else:
                        freq_dict[ec_key] = protein_count
                else:
                    freq_dict[ec_key] = protein_count
    print("double-count:", double_count)
    return freq_dict

def create_swissprot_map(SWISS_PROT_MAP):
    
    mapping_dict = {} #key: EC, val: swissprot proteins
    final_mapping_dict = dict() #key: swissprot, val: ECs
    with open(SWISS_PROT_MAP, "r") as mapping:
        for line in mapping.readlines():
            cleaned_line = line.strip("\n")
            line_as_list = cleaned_line.split("\t")
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


def import_detect_ec(detect_fbeta_file, gene_ec_dict, prob_dict):
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
                probability = list_line[2]
                if(key in gene_ec_dict):
                    
                    new_list = gene_ec_dict[key] #Apparently manager dicts can't handle appends.  this is eye-opening....
                    
                    if(EC_val in new_list):
                        print(dt.today(), "this is a problem.  the same EC is annotated twice")
                        sys.exit("killer")
                    
                    new_list.append(EC_val)
                    prob_key = key + "_" + EC_val
                    prob_dict[prob_key] = probability
                    
                    gene_ec_dict[key] = new_list
                    
                
                else:
                    gene_ec_dict[key] = [EC_val]
                    prob_key = key + "_" + EC_val
                    prob_dict[prob_key] = probability
                   
def import_priam_ec_v2(priam_sequence_ec, gene_ec_dict, prob_ec_dict):
    #just uses E-values.  This one's the low-quality one that uses a cutoff of E-5, regardless if PRIAM counted it
    query_name = "None"
    ec_list = []
    line_count = 0
    with open(priam_sequence_ec, "r") as priam_ec:
        for line in priam_ec:
            line_count += 1
            #print(line)
            if(line == "\n"):
                continue
            
            elif(line.startswith(">")):
                #new line
                if(len(ec_list) > 0):
                    if(query_name is "None"):
                        print(dt.today(), "This shouldn't happen.  full EC list.  no query")
                    else:
                        gene_ec_dict[query_name] = ec_list
                        #print("===================================")
                        ec_list[:] = []
                        
                query_name = line.strip(">")
                query_name = query_name.strip("\n")
                query_list = query_name.split(" ")
                query_name = query_list[0]
                #print(dt.today(), "working on", query_name)
                
            else:
                list_line = line.split("\t")
                ec = list_line[0].strip("\n")
                ec = ec.strip(" ")
                ec = ec.strip("#")
                e_value = float(list_line[2])
                #print(dt.today(), "PRIAM: looking at e-value:", e_value)
                if(e_value < 1e-5):
                    #print("evalue passed:", ec)
                    #time.sleep(0.5)
                    ec_list.append(ec)
                    prob_key = query_name + "_" + ec
                    if(prob_key in prob_ec_dict):
                        continue
                    else:
                        prob_ec_dict[prob_key] = e_value
                    
            if(query_name is "None"):
                print(dt.today(), "This shouldn't be happening.  a line was skipped")
    
def import_priam_ec(priam_sequence_ec, gene_ec_dict, prob_ec_dict):
    #gene_ec_dict = dict()
    query_name = "None"
    ec_list = []
    line_count = 0
    with open(priam_sequence_ec, "r") as priam_ec:
        for line in priam_ec:
            line_count += 1
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
                        
                        prob_key = query_name + "_" + ec
                        if(prob_key in prob_ec_dict):
                            continue
                        else:
                            prob_ec_dict[prob_key] = probability
                        
                    else:
                        print("line valid but not inserted")
                
            if(query_name is "None"):
                print(dt.today(), "This shouldn't be happening.  a line was skipped")
    #return gene_ec_dict

def import_diamond_ec_v2(diamond_proteins_blastout, swissprot_map_dict, gene_length_dict, gene_ec_lq_dict, gene_ec_hq_dict):
    #This one uses E-values
    with open(diamond_proteins_blastout, "r") as diamond_ec:
        for item in diamond_ec:
            #note that this DIAMOND call has special parameters (for some reason)
            list_line = item.split("\t")
            #print(list_line)
            query_name = list_line[0].strip("\n")
            query_name = query_name.strip(">")
            swissprot_name = list_line[1]
            

            e_value = float(list_line[7])
            #print(dt.today(), "DIAMOND e-value:", e_value)
            if(e_value <= 1e-5):
                #take this hit
                if(swissprot_name in swissprot_map_dict):
                    EC_val = swissprot_map_dict[swissprot_name]
                    if(query_name in gene_ec_lq_dict):
                        old_entry = gene_ec_lq_dict[query_name]
                        old_entry += EC_val
                        gene_ec_lq_dict[query_name] = old_entry
                    else:
                        gene_ec_lq_dict[query_name] = EC_val
                else:
                    print(dt.today(), "LQ import diamond V2:", swissprot_name, "not found in swissprot map.  but it's fine")
                    gene_ec_lq_dict[query_name] = ["0.0.0.0"]
            
            if(e_value <= 1e-10):
                if(swissprot_name in swissprot_map_dict):
                    EC_val = swissprot_map_dict[swissprot_name]
                    if(query_name in gene_ec_hq_dict):
                        old_entry = gene_ec_hq_dict[query_name]
                        old_entry += EC_val
                        gene_ec_hq_dict[query_name] = old_entry
                    else:
                        gene_ec_hq_dict[query_name] = EC_val
                else:
                    print(dt.today(), "HQ import diamond V2:", swissprot_name, "not found in swissprot map.  but it's fine")
                    gene_ec_hq_dict[query_name] = ["0.0.0.0"]
            else:   
                gene_ec_hq_dict[query_name] = ["0.0.0.0"]
                    
                    
                    
            
    
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
                        print(dt.today(), "import diamond V1:", swissprot_name, "not found in swissprot map. but it's fine")
                        EC_val = ["0.0.0.0"]
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


def extract_value(x):
    #used for sorting the ECs tied with prob
    split = x.split("_")
    return float(split[0]) 


def get_intersection(common_keys, diamond_ec_dict, priam_ec_dict):
    common_dict = dict()
    for item in common_keys:
        priam_ec_set = set(priam_ec_dict[item])
        diamond_ec_set = set(diamond_ec_dict[item])
        combined_ec = list(diamond_ec_set.intersection(priam_ec_set))
        if(len(combined_ec) > 0):
            common_dict[item] = combined_ec    
    return common_dict
    
def merge_everything(detect_ec_dict, common_dict, detect_prob_dict, priam_prob_dict, freq_dict):    
    #sys.exit("pre-death")
    #key is gene
    final_dict = dict()
    origin_dict = dict() #key: gene, val: origin of match
    for key in detect_ec_dict.keys():
    #then take all of DETECT's results
        detect_ec_list = sorted(detect_ec_dict[key])
        if(key in common_dict):
            #arrange everything by high prob of ECs, but separate (DETECT, PRIAM/diamond)
            #then take the top 2.  if they're found in the ENZYME freq file, it's valid. else just take the highest prob EC
            #if there's more than 1 detect EC, we'll take those.  disregard the others.
            #if the pair exists, take the pair.  if they don't take the top DETECT hit
            if(len(detect_ec_list) > 1):
                pair_list = []
                prob_list = []
                for ec in detect_ec_list:
                    prob_key = key + "_" + ec
                    prob = detect_prob_dict[prob_key]
                    prob_list.append(prob + "_" + ec)
                prob_list = sorted(prob_list, key = extract_value)
                prob_list.reverse()
                prob_ec_0 = prob_list[0].split("_")[1]
                prob_ec_1 = prob_list[1].split("_")[1]
                pair_list.append(prob_ec_0)
                pair_list.append(prob_ec_1)
                pair_list = sorted(pair_list)
                ec_pair_key = pair_list[0] + "_" + pair_list[1]
                if(ec_pair_key in freq_dict):
                    final_dict[key] = pair_list
                    origin_dict[key] = "detect_multi"
                else:
                    final_dict[key] = [prob_ec_0]
                    origin_dict[key] = "detect_multi"
                    
                
            elif(len(detect_ec_list) == 1):
                #only 1 detect EC.  use the detect EC, take the top priam/detect
                prob_list = []
                pair_list = []
                for ec in common_dict[key]:
                    prob_key = key + "_" + ec
                    prob = priam_prob_dict[prob_key]
                    prob_list.append(str(prob) + "_"+ ec)
                prob_list = sorted(prob_list, key = extract_value)
                prob_list.reverse()
                ec_0 = prob_list[0].split("_")[1]
                pair_list.append(ec_0)
                pair_list += detect_ec_list
                pair_list = sorted(pair_list)
                ec_pair_key = pair_list[0] + "_" + pair_list[1]
                if(ec_pair_key in freq_dict):
                    final_dict[key] = pair_list
                    origin_dict[key] = "detect_priam_multi"
                else:
                    final_dict[key] = detect_ec_list
                    origin_dict[key] = "detect_priam_multi"
            
        else:
            if(len(detect_ec_list) > 1):
                pair_list = []
                prob_list = []
                for ec in detect_ec_list:
                    prob_key = key + "_" + ec
                    prob = detect_prob_dict[prob_key]
                    prob_list.append(str(prob) + "_" + ec)
                prob_list = sorted(prob_list, key = extract_value)
                prob_list.reverse()
                prob_ec_0 = prob_list[0].split("_")[1]
                prob_ec_1 = prob_list[1].split("_")[1]
                pair_list.append(prob_ec_0)
                pair_list.append(prob_ec_1)
                pair_list = sorted(pair_list)
                ec_pair_key = pair_list[0] + "_" + pair_list[1]
                if(ec_pair_key in freq_dict):
                    final_dict[key] = pair_list
                    origin_dict[key] = "detect_multi"
                else:
                    print(dt.today(), "detect pair failed", ec_pair_key)
                    final_dict[key] = [prob_ec_0]
                    origin_dict[key] = "detect"
            elif(len(detect_ec_list) == 1):
                final_dict[key] = detect_ec_list
                origin_dict[key] = "detect"
   
    #dealing with stuff only in priam/diamond
    for key in common_dict.keys():
        if(key in final_dict):
            continue
        else:
            
            common_ec_list = common_dict[key]
            if(len(common_ec_list) > 1):
                pair_list = []
                prob_list = []
                for ec in common_ec_list:
                    prob_key = key + "_" + ec
                    prob = priam_prob_dict[prob_key]
                    prob_list.append(str(prob) + "_"+ ec)
                prob_list = sorted(prob_list, key = extract_value)
                prob_list.reverse()
                pair_list.append(prob_list[0].split("_")[1])
                pair_list.append(prob_list[1].split("_")[1])
                pair_list = sorted(pair_list)
                ec_pair_key = pair_list[0] + "_" + pair_list[1]
                
                if(ec_pair_key in freq_dict):
                    final_dict[key] = pair_list
                    origin_dict[key] = "priam_multi"
                else:
                    print(dt.today(), "pair failed", ec_pair_key)
                    final_dict[key] = [pair_list[0]]
                    origin_dict[key] = "priam"
            elif(len(common_ec_list) == 1):
                final_dict[key] = common_ec_list
                origin_dict[key] = "priam"
    return final_dict, origin_dict
    
def export_report(final_dict, Output_file):
    with open(Output_file, "w") as ec_out:
        for item in sorted(final_dict.keys()):
            skip_line = False
            ec_list = final_dict[item]
            ec_string = ""
            for ec in ec_list:
                actual_ec = ec.strip(" ")
                ec_string += actual_ec + "|"
            ec_string = ec_string[:-1]
            if(len(ec_list) == 0):
                continue
            if(not skip_line):     
                ec_out.write(item + "\t" + str(len(ec_list)) + "\t" + ec_string + "\n")  

def export_simple(origin_dict, output_file):
    with open(output_file + "_origin", "w") as origin_out:
        for item in origin_dict.keys():
            out_line = item + "\t" + origin_dict[item] + "\n"
            origin_out.write(out_line)
            
    
def thing_in_question(inspect_key, detect_ec_dict, priam_ec_lq_dict, priam_ec_hq_dict, diamond_ec_lq_dict, diamond_ec_hq_dict):
    #inspect_key = "gi|490134789|ref|NZ_KB822565.1|:c146428-142376|97138|g__Clostridium.s__Clostridium_sp_ASF356|UniRef90_unknown|UniRef50_K4KWD0"
    #inspect_key = "gi|498503065|ref|NZ_KB822571.1|:503009-504169|1235803|g__Parabacteroides.s__Parabacteroides_sp_ASF519|UniRef90_R6X568|UniRef50_E8QZN3"
    #inspect_key = "gi|483984714|ref|NZ_KB892660.1|:c111601-110108|33035|g__Blautia.s__Blautia_producta|UniRef90_unknown|UniRef50_R5TWD5"
    print("thing in question:")
    print(inspect_key)
    if(inspect_key in priam_ec_lq_dict):
        print("PRIAM LQ:", priam_ec_lq_dict[inspect_key])
    else:
        print("PRIAM HQ: none")
    if(inspect_key in priam_ec_hq_dict):
        print("PRIAM HQ:", priam_ec_hq_dict[inspect_key])
    else:
        print("PRIAM HQ: none")
        
    if(inspect_key in diamond_ec_hq_dict):
        print("DIAMOND HQ:", diamond_ec_hq_dict[inspect_key])
    else:
        print("DIAMOND HQ: none")        
    if(inspect_key in diamond_ec_lq_dict):
        print("DIAMOND LQ:", diamond_ec_lq_dict[inspect_key])
    else:
        print("DIAMOND LQ: none")
        
        
    if(inspect_key in detect_ec_dict):
        print("DETECT:", detect_ec_dict[inspect_key])
    else:
        print("DETECT: none")
    
    
    
if __name__ == "__main__":
    detect_file     = sys.argv[1]
    priam_file      = sys.argv[2]
    diamond_file    = sys.argv[3]
    SWISS_PROT_MAP  = sys.argv[4]
    gene_map_file   = sys.argv[5]
    freq_file       = sys.argv[6]
    hq_output       = sys.argv[7]
    lq_output       = sys.argv[8]
    
    
    
    gene_length_dict = import_gene_map(gene_map_file)
    freq_dict = import_freq_file(freq_file)
    print(dt.today(), "size of freq dict:", len(freq_dict))
    
    ec_process_list = []
    #-----------------------------------------
    # import the data
    manager = mp.Manager()
    diamond_ec_lq_mgr_dict = manager.dict()
    diamond_ec_hq_mgr_dict = manager.dict()
    swissprot_map_dict = create_swissprot_map(SWISS_PROT_MAP)
        
    diamond_ec_process = mp.Process(
        target = import_diamond_ec_v2, 
        args = (diamond_file, swissprot_map_dict, gene_length_dict, diamond_ec_lq_mgr_dict, diamond_ec_hq_mgr_dict)
    )
    diamond_ec_process.start()
    ec_process_list.append(diamond_ec_process)
    
    priam_ec_hq_mgr_dict = manager.dict()
    priam_prob_hq_mgr_dict = manager.dict()
    priam_ec_hq_process = mp.Process(
        target = import_priam_ec, 
        args = (priam_file, priam_ec_hq_mgr_dict, priam_prob_hq_mgr_dict)
    )
    priam_ec_hq_process.start()
    ec_process_list.append(priam_ec_hq_process)
    
    priam_ec_lq_mgr_dict = manager.dict()
    priam_prob_lq_mgr_dict = manager.dict()
    priam_ec_lq_process = mp.Process(
        target = import_priam_ec_v2, 
        args = (priam_file, priam_ec_lq_mgr_dict, priam_prob_lq_mgr_dict)
    )
    priam_ec_lq_process.start()
    ec_process_list.append(priam_ec_lq_process)
    
    detect_ec_mgr_dict = manager.dict()
    detect_prob_mgr_dict = manager.dict()
    detect_ec_process = mp.Process(
        target = import_detect_ec,
        args = (detect_file, detect_ec_mgr_dict, detect_prob_mgr_dict)
    )
    detect_ec_process.start()
    ec_process_list.append(detect_ec_process)
    
    print(dt.today(), "waiting for import jobs to finish")
    for item in ec_process_list:
        item.join()
    ec_process_list[:] = []
    print(dt.today(), "import jobs finished")
    
    
    
    
    diamond_ec_lq_dict = dict(diamond_ec_lq_mgr_dict) #key: gene.  value: list of ECs
    diamond_ec_hq_dict = dict(diamond_ec_hq_mgr_dict) #key: gene.  value: list of ECs
    priam_ec_lq_dict = dict(priam_ec_lq_mgr_dict)     #key: gene.  value: list of ECs
    priam_ec_hq_dict = dict(priam_ec_hq_mgr_dict)     #key: gene.  value: list of ECs
    detect_ec_dict = dict(detect_ec_mgr_dict)   #key: gene.  value: list of ECs
    detect_prob_dict = dict(detect_prob_mgr_dict) #key: gene_EC. val: prob
    priam_prob_lq_dict = dict(priam_prob_lq_mgr_dict)#key: gene_EC, val: prob
    priam_prob_hq_dict = dict(priam_prob_hq_mgr_dict)#key: gene_EC, val: prob
    print(dt.today(), "finished converting dicts")
    
    
    #testing constructs
    #inspect_key = "gi|555554272|ref|NZ_KI530591.1|:88165-88941|248039|g__Mucispirillum.s__Mucispirillum_schaedleri|UniRef90_V2Q7Z1|UniRef50_E4TJW6"
    #thing_in_question(inspect_key, detect_ec_dict, priam_ec_lq_dict, priam_ec_hq_dict, diamond_ec_lq_dict, diamond_ec_hq_dict)
    #inspect_key = "gi|555554272|ref|NZ_KI530591.1|:87704-88150|248039|g__Mucispirillum.s__Mucispirillum_schaedleri|UniRef90_V2Q981|UniRef50_B0TXG6"
    #print("============================================")
    #thing_in_question(inspect_key, detect_ec_dict, priam_ec_lq_dict, priam_ec_hq_dict, diamond_ec_lq_dict, diamond_ec_hq_dict)
    
    #sys.exit("charized")
    
    #--------------------------------------------------
    # merge the results.  get the genes for which both tools came back with an EC
    diamond_keys_lq = set(diamond_ec_lq_dict.keys())
    priam_keys_lq = set(priam_ec_lq_dict.keys())
    common_keys_lq = diamond_keys_lq & priam_keys_lq
    
    diamond_keys_hq = set(diamond_ec_hq_dict.keys())
    priam_keys_hq = set(priam_ec_hq_dict.keys())
    common_keys_hq = diamond_keys_hq & priam_keys_hq
    
    
    print(dt.today(), "finished getting common genes")
    #take the intersection of the ECs found for each gene between priam and diamond
    common_lq_dict = get_intersection(common_keys_lq, diamond_ec_lq_dict, priam_ec_lq_dict)
    common_hq_dict = get_intersection(common_keys_hq, diamond_ec_hq_dict, priam_ec_hq_dict)

    final_hq_dict, origin_hq_dict = merge_everything(detect_ec_dict, common_hq_dict, detect_prob_dict, priam_prob_hq_dict, freq_dict)   
    final_lq_dict, origin_lq_dict = merge_everything(detect_ec_dict, common_lq_dict, detect_prob_dict, priam_prob_lq_dict, freq_dict)            
    
    #----------------------------------------------
    #export the final ec list
    
    export_report(final_hq_dict, hq_output)
    export_report(final_lq_dict, lq_output)
    #export the final ec list
    export_simple(origin_hq_dict, hq_output)
    export_simple(origin_lq_dict, lq_output)
    
    
    
    
    
