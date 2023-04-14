#this creates a lookup table of all reads in the sample, and says if it's in a contig.

import sys
import os
from datetime import datetime as dt
import time
#this solely relies on sort working.  Have all the reads sorted, and walk through.

def import_reads_in_contig(contig_reads_file):
    reads_list = list()
    with open(contig_reads_file, "r") as contig_reads:
        for line in contig_reads:
            cleaned_line = line.strip("\n")
            reads_list.append(cleaned_line)
            
    reads_list = sorted(reads_list)
    return reads_list

def import_raw_fastq(fastq_file):
    read_id_list = []
    with open(fastq_file, "r") as fastq:
        line_count = 0
        for line in fastq:
            if(line_count % 4 == 0):
                read_id = line.strip("@")
                read_id = read_id.strip("\n")
                read_id = read_id.split(" ")[0]
                read_id_list.append(read_id)
            line_count += 1
    read_id_list = sorted(read_id_list)
    return read_id_list

if __name__ == "__main__":
    contig_reads_file = sys.argv[1]
    raw_fastq_file = sys.argv[2]
    export_report = sys.argv[3]
    print(dt.today(), "importing")
    reads_in_contig_list = import_reads_in_contig(contig_reads_file)
    reads_list = import_raw_fastq(raw_fastq_file)
    print(dt.today(), "processing and exporting")
    contig_index = 0
    final_dict = dict()
    with open(export_report, "w") as out_file:
        for item in reads_list:
            if(contig_index >= len(reads_in_contig_list)):
                out_line = item + "\t" + "no" + "\n"
            else:
                
                selected_contig_read = reads_in_contig_list[contig_index]
                out_line = ""
                if(item != selected_contig_read):
                    #time.sleep(0.01)
                    #print(dt.today(), item, "!=", selected_contig_read)
                    out_line = item + "\t" + "no" + "\n"
                    
                else:
                    
                    out_line = item + "\t" + "yes" + "\n"
                    #time.sleep(0.01)
                    #print(dt.today(),item, "==", selected_contig_read, "moving on")
                    #print("------------------------------------------------------------------")
                    #if(contig_index == len(reads_in_contig_list)):
                    #    contig_index = len(reads_in_contig_list)
                    #else:
                    contig_index += 1
            out_file.write(out_line)
    print(dt.today(), "done")
    
    
    
        
        
    
