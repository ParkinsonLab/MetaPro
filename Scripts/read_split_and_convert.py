#splits a fastq and converts it into fasta segments.
#used specifically for rRNA, to reduce the number of files

import os
import sys
import time

def import_fastq(fastq_file):
    fastq_dict = dict()
    line_count = 0
    seq = ""
    fastq_ID = ""
    with open(fastq_file, "r") as fastq_in:
        for line in fastq_in:
            if(line_count % 4 == 0):
                if(line_count != 0):
                    fastq_dict[fastq_ID] = seq
                    seq = ""

                fastq_ID = ">" + line.strip("\n").split(" ")[0].strip("@")
                

                line_count = 1
            else:
                if(line_count == 1):
                    seq = line.strip("\n")
                    line_count += 1
                else:
                    line_count += 1
        
        #final entry
        fastq_dict[fastq_ID] = seq
    return fastq_dict
            





if __name__ == "__main__":
    fastq_in_file = sys.argv[1]
    fasta_out_dir = sys.argv[2]
    read_split_count = int(sys.argv[3])


    fastq_dict = import_fastq(fastq_in_file)
    sorted_keys = list(sorted(fastq_dict.keys()))
    
    fastq_basename = os.path.basename(fastq_in_file).split(".")[0]
    file_split_count = 0
    for limit in range(0, len(sorted_keys), read_split_count):

        export_name = fastq_basename + "_" + str(file_split_count) + ".fasta"
        export_file_name = os.path.join(fasta_out_dir, export_name)
        file_split_count += 1
        with open(export_file_name, "w") as out_file:
            if(limit + read_split_count > len(sorted_keys)):
                #print("final batch")
                for i in range(limit, len(sorted_keys)):
                    out_file.write(sorted_keys[i] + "\n")
                    out_file.write(fastq_dict[sorted_keys[i]] + "\n")
            else:    
                #print("normal batch")
                for i in range(limit, limit + read_split_count):
                    out_key = sorted_keys[i]
                    out_seq = fastq_dict[sorted_keys[i]]
                    

                    out_file.write(out_key + "\n")
                    out_file.write(out_seq + "\n")
