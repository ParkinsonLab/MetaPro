#Oct 20, 2020:
#This script makes the contig map.
#-> it makes sure only 1 read-ID is represented.
#-> it also separates what's used by contig assembly and exports those reads all at once

import os
import sys
import pandas as pd
import re
from datetime import datetime as dt
import time

def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df


def get_match_score(cigar_segment):
    CIGAR = re.split("([MIDNSHPX=])", cigar_segment) # Split CIGAR string into list, placing
    CIGAR = CIGAR[:-1]                      #lop off the empty char artifact from the split
    position_count = 0                      #position counter, because the CIGAR string is split into alternating segments of <length><Label>, 
    length = 0
    matched = 0
    segment_length = 0
    for item in CIGAR:
        if((position_count %2) == 0):       #every even position (starting from 0) is going to be a length
            segment_length = int(item)
        elif((position_count %2) == 1):     #every odd position is going to be a label
            length += segment_length
            if(item == "M"):
                matched += segment_length
        position_count += 1
    if(length == 0):
        return 0
        
    match_score = 100 * (matched / length)
    
    return match_score

def add_new_read(status, contig_id, AS_score):
    inner_dict = dict()
    inner_dict["status"] = status
    inner_dict["contig"] = contig_id
    inner_dict["score"] = AS_score
    return inner_dict

def import_samfile(samfile):
    #only one read per contig
    read_status_dict = dict() #stores match state, and top hit if applicable
      
    contig_read_dict = dict()
    
    
    with open(samfile, "r") as sam_in:
        for line in sam_in:
            if line.startswith("@") or len(line) <= 1:
                continue
            else:
                line_split = line.split("\t")
                match_score = get_match_score(line_split[5]) #decode cigar segment
                read_id = line_split[0]
                contig_id = line_split[2]
                flag_data = line_split[1]
                
                AS_score = 0
                for item in line_split:
                    if(item.startswith("AS")):
                        AS_score = int(item.split(":")[-1])
                
                flag_bin = bin(int(flag_data))[2:].zfill(11)
                if(flag_bin[8] == "1"):
                    #not matched to anything
                    if(read_id in read_status_dict):
                        continue
                    else:
                        inner_dict = add_new_read("u", contig_id, AS_score)
                        read_status_dict[read_id] = inner_dict
                else:
                    if(read_id in read_status_dict):
                        #repeat read
                        if(read_status_dict[read_id]["status"] == "c"):
                            prev_contig_match = read_status_dict[read_id]["contig"]
                            if(prev_contig_match == contig_id):
                                #matches to same contig
                                continue
                            else:
                                prev_AS_score = read_status_dict[read_id]["score"]
                                if(AS_score > prev_AS_score):
                                    #new match is better.
                                    contig_read_dict[prev_contig_match].remove(read_id)
                                    if(len(contig_read_dict[prev_contig_match]) == 0):
                                        del contig_read_dict[prev_contig_match]
                                        
                                    if(contig_id in contig_read_dict):
                                        contig_read_dict[contig_id].append(read_id)
                                    else:
                                        contig_read_dict[contig_id] = [read_id]
                                    
                                    read_status_dict[read_id]["contig"] = contig_id
                                    read_status_dict[read_id]["score"] = AS_score
                                    #time.sleep(1)
                        else:
                            #repeat read previously unmatched (happens in paired)
                            inner_dict = add_new_read("c", contig_id, AS_score)
                            read_status_dict[read_id] = inner_dict
                            if(contig_id in contig_read_dict):
                                contig_read_dict[contig_id].append(read_id)
                            else:
                                contig_read_dict[contig_id] = [read_id]
                                   
                    else:
                        #fresh read
                        inner_dict = add_new_read("c", contig_id, AS_score)
                        read_status_dict[read_id] = inner_dict
                        
                        if(contig_id in contig_read_dict):
                            #old contig
                            contig_read_dict[contig_id].append(read_id)
                        else:
                            contig_read_dict[contig_id] = [read_id]
                            
    return contig_read_dict, read_status_dict

def merge_dict(paired_dict, singletons_dict, op_mode):
    
    final_dict = singletons_dict
    for item in paired_dict:
        if(item in final_dict):
            if(op_mode == "read_status"):
                print("this shouldn't happen")
                sys.exit("deather")
            else:
                final_dict[item] += paired_dict[item]
            #print("new:", final_dict[item])
        else:
            #new entry
            final_dict[item] = paired_dict[item]
    return final_dict

def export_contig_map(map_out_name, final_contig_read_dict):
    with open(map_out_name, "w") as map_out:
        for contig in final_contig_read_dict:
            out_line = contig + "\t" + str(len(final_contig_read_dict[contig]))
            for read in final_contig_read_dict[contig]:
                out_line += "\t" + read
            out_line += "\n"
            
            map_out.write(out_line)
            
def pull_unmapped_reads(read_id, read_status_dict):
    if(read_id in read_status_dict):
        read_status = read_status_dict[read_id]["status"]
        if(read_status == "c"):
            return False
        elif(read_status == "u"):
            return True
        else:
            return False
    else:
        print(dt.today(), "read ID not found in status dict. this shouldn't happen")
        sys.exit("what??")
        
def export_unmapped_reads(raw_read_location, export_location, sample_file, sample_read_status_dict):
    export_path = os.path.join(raw_read_location, sample_file + ".fastq")
    sample_df = import_fastq(export_path)
    unmapped_sample_df = sample_df[sample_df["ID"].apply(lambda x: pull_unmapped_reads(x, sample_read_status_dict))]
    unmapped_sample_df["ID"] = "@" + unmapped_sample_df["ID"]
    export_sample_path = os.path.join(export_location, sample_file + ".fastq")
    unmapped_sample_df.to_csv(export_sample_path, header = False, index = False, mode = "w", sep = "\n", quoting = 3)
       

if __name__ == "__main__":
    operating_mode = sys.argv[1]
    raw_read_location = sys.argv[2]
    export_location = sys.argv[3]
    singletons_sam = sys.argv[4]
    singletons_contig_read_dict, singletons_read_status_dict = import_samfile(singletons_sam)
    paired_contig_read_dict = dict()
    paired_read_status_dict = dict()
    
    if(operating_mode == "paired"):
        paired_sam = sys.argv[5]
        paired_contig_read_dict, paired_read_status_dict = import_samfile(paired_sam)
    
        singletons_reads = singletons_read_status_dict.keys()
        paired_reads = paired_read_status_dict.keys()
        common_keys = list(set(singletons_reads) & set(paired_reads))
        if(len(common_keys) > 0):
            print(dt.today(), "singleton keys in paired.  this shouldn't happen")
            sys.exit("death")
        
    final_contig_read_dict = merge_dict(paired_contig_read_dict, singletons_contig_read_dict, "contigs")
    contig_map_out = os.path.join(export_location, "contig_map.tsv")
    export_contig_map(contig_map_out, final_contig_read_dict)
    
    export_unmapped_reads(raw_read_location, export_location, "singletons", singletons_read_status_dict)
    
    if(operating_mode == "paired"):
        export_unmapped_reads(raw_read_location, export_location, "pair_1", paired_read_status_dict)
        export_unmapped_reads(raw_read_location, export_location, "pair_2", paired_read_status_dict)
            
            