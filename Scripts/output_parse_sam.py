import sys
import re
from datetime import datetime as dt
def import_sam_file(sam_file):
    unmapped_set = set()
    read_to_contig_dict = dict() #key: read, val: dict: count division, contigs
    
    db_match_dict = dict()
    len_chars= ["M","I","S","=","X"] 
    reads_with_multiple_contigs = 0
    with open(sam_file, "r") as sam_report:
        for line in sam_report:
            if(line.startswith("@")):
                continue
            else:
                cleaned_line = line.strip("\n")
                line_list = cleaned_line.split("\t")
                
                read_query = line_list[0]
                cigar_part = line_list[1]
                db_match = line_list[2]
                flag = bin(int(cigar_part))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                if(flag[8] == "1"):
                    unmapped_set.add(read_query)
                    continue
                else:
                    
                    CIGAR = re.split("([MIDNSHPX=])", line_list[5]) # Split CIGAR string into list, placing
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
                    
                    
                    if matched<length*0.9:                          # If alignment is <90% matched:
                        unmapped_set.add(read_query)                         # add it to the unmapped set and
                        continue                                    # skip to the next query.
                    else:
                        if(read_query in read_to_contig_dict):
                            reads_with_multiple_contigs += 1
                            inner_dict = read_to_contig_dict[read_query]
                            
                            old_set = inner_dict["contigs"]
                            old_set.add(db_match)
                            inner_dict["contigs"] = old_set
                            inner_dict["count"] = 1/len(old_set)
                            read_to_contig_dict[read_query] = inner_dict
                        else:
                            inner_dict = dict()
                            inner_dict["contigs"] = set([db_match])
                            inner_dict["count"] = 1
                            read_to_contig_dict[read_query] = inner_dict
                    
                        if(db_match in db_match_dict):
                            old_set = db_match_dict[db_match]
                            old_set.add(read_query)
                            db_match_dict[db_match] = old_set
                        else:
                            new_set = set([read_query])
                            db_match_dict[db_match] = new_set
                            
                                
    print("reads belonging to multiple_contigs:", reads_with_multiple_contigs)
    return db_match_dict, read_to_contig_dict
if __name__ == "__main__":
    sam_file = sys.argv[1]
    contig_read_count_export = sys.argv[2]
    reads_in_contig_export = sys.argv[3]
    contig_read_map_export = sys.argv[4]
    
    print(dt.today(), "started import")
    db_match_dict, read_to_contig_dict = import_sam_file(sam_file)
    print(dt.today(), "started final tally")
    final_dict = dict()
    
    for read in read_to_contig_dict:
        inner_dict = read_to_contig_dict[read]
        contig_set = inner_dict["contigs"]
        per_read_count = inner_dict["count"]
        for contig in contig_set:
            #if(contig == "gene_15931|GeneMark.hmm|96_nt|-|373|468>NODE_9971_length_468_cov_44.139241_g9455_i0"):
            #    print("contig found")

            if(contig in final_dict):
                final_dict[contig] += per_read_count
            else:
                final_dict[contig] = per_read_count
        
        
        
        
        
    print(dt.today(), "started export")
    #contig -> read 
    with open(contig_read_map_export, "w") as contig_map_out:
        for item in db_match_dict:
            true_read_count = 0
            if(item in final_dict):
                true_read_count = final_dict[item]
            else:
                print("This shouldn't happen.  contig in db match not found in final_dict")
                sys.exit("db_match error")
            
            out_line = item + "\t" + str(true_read_count)
            
            read_list = list(db_match_dict[item])
            for read in read_list:
                out_line += "\t" + read 
            out_line += "\n"    
            contig_map_out.write(out_line)
    
    #contig -> read_count
    with open(contig_read_count_export, "w") as out_file:
        for contig in final_dict:
            out_line = contig + "\t" + str(final_dict[contig]) + "\n"
            out_file.write(out_line)
            
    #list of reads in contigs    
    with open(reads_in_contig_export, "w") as read_to_contig_out:
        for read in read_to_contig_dict:
            read_to_contig_out.write(read + "\n")
    print(dt.today(), "finished")
    """
    with open(other_export_file, "w") as out_file:
        for contig in db_match_dict:
            out_line = contig
            reads = db_match_dict[contig]
            for read in reads:
                out_line += "\t" + read
            out_line += "\n"
            out_file.write(out_line)
    
    with open(export_file, "w") as out_file:
        for read in read_to_contig_dict:
            out_line = read
            contig_matches = read_to_contig_dict[read]
            for item in contig_matches:
                out_line += "\t" + item
            out_line += "\n"
            out_file.write(out_line)
    """
