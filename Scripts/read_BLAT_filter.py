<<<<<<< HEAD
#!/usr/bin/env python

import sys
from Bio import SeqIO

input_file = sys.argv[1]
input_seqs = SeqIO.to_dict(SeqIO.parse(input_file, "fastq"))
BLAT_tab_file = sys.argv[2]
output_file = sys.argv[3]
output_seqs = []
output_file_made = False
contaminat_output_file = sys.argv[4]
contaminat_output_seqs = []
contaminat_output_file_made = False

with open(BLAT_tab_file, "r") as tabfile:
    contaminated_seqs = []
    query_seq = ""
=======
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 23:05:20 2019

@author: Billy
"""

#!/usr/bin/env python

import sys
import pandas as pd
from datetime import datetime as dt


def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    return fastq_df



if __name__ == "__main__":
    input_file = sys.argv[1]            #IN: host-free seqs
    BLAT_tab_file = sys.argv[2]         #IN: BLAT report
    output_file = sys.argv[3]           #OUT: final host-free seqs
    contaminant_output_file = sys.argv[4]#OUT: host-data that was stripped
    
    #import the sequence
    if(input_file.endswith(".fastq")):
        input_seqs = import_fastq(input_file)#SeqIO.to_dict(SeqIO.parse(input_file, "fastq"))
    else:
        print("only FASTQs supported.  program ending")
    #with open(BLAT_tab_file, "r") as tabfile:
    contaminated_seqs = []
    query_seq = ""
    #we don't know enough about BLATOUT files.  we assume we just want the IDs, but we aren't entirely sure
    #what the rules are here.  So we can leave this inefficient loop in for now. (feb 06, 2019)
    #pick out the results from the BLATOUT report
    tabfile = open(BLAT_tab_file, "r")
>>>>>>> 4d5286c... committing final-ish code.
    for line in tabfile:
        if len(line) < 2:
            continue
        line_parts = line.split("\t")
        if query_seq == line_parts[0]:
            continue
        else:
            query_seq = line_parts[0]
<<<<<<< HEAD
            contaminated_seqs.append(query_seq)
            
    for seq in input_seqs:
        if seq not in contaminated_seqs:
            output_seqs.append(input_seqs[seq])
            if len(output_seqs) > 100000:
                if not output_file_made:
                    output_file_made = True
                    with open(output_file, "w") as outfile:
                        SeqIO.write(output_seqs, outfile, "fastq")
                else:
                    with open(output_file, "a") as outfile:
                        SeqIO.write(output_seqs, outfile, "fastq")
                output_seqs = []
        else:
            contaminat_output_seqs.append(input_seqs[seq])
            if len(contaminat_output_seqs) > 100000:
                if not contaminat_output_file_made:
                    contaminat_output_file_made = True
                    with open(contaminat_output_file, "w") as outfile:
                        SeqIO.write(contaminat_output_seqs, outfile, "fastq")
                else:
                    with open(contaminat_output_file, "a") as outfile:
                        SeqIO.write(contaminat_output_seqs, outfile, "fastq")
                contaminat_output_seqs = []
                
                
    if len(output_seqs) > 0:
        if not output_file_made:
            output_file_made = True
            with open(output_file, "w") as outfile:
                SeqIO.write(output_seqs, outfile, "fastq")
        else:
            with open(output_file, "a") as outfile:
                SeqIO.write(output_seqs, outfile, "fastq")
    if len(contaminat_output_seqs) > 0:
        if not contaminat_output_file_made:
            contaminat_output_file_made = True
            with open(contaminat_output_file, "w") as outfile:
                SeqIO.write(contaminat_output_seqs, outfile, "fastq")
        else:
            with open(contaminat_output_file, "a") as outfile:
                SeqIO.write(contaminat_output_seqs, outfile, "fastq")
=======
            if not(query_seq.startswith("@")):
                query_seq = "@" + str(query_seq)
            contaminated_seqs.append(query_seq)
    tabfile.close()
    
    #The BLAT report spits out what was aligned as human.  
    #so we take the NOT of the intersection to get host-free
    #split the input sequence by host/host-free and export
    host_free_df = input_seqs[~input_seqs.ID.isin(contaminated_seqs)]
    host_free_df.to_csv(output_file, sep = '\n', mode = "w+", header = False, index = False, quoting = 3)
    
    host_only_df = input_seqs[input_seqs.ID.isin(contaminated_seqs)]
    host_only_df.to_csv(contaminant_output_file, sep = '\n', mode = "w+", header = False, index = False, quoting = 3)
   
>>>>>>> 4d5286c... committing final-ish code.
