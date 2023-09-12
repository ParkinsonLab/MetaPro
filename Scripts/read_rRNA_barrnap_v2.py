import os
import os.path
import sys
import pandas as pd
import time

#This module takes in the Output report of the barrnap tool, and the fastq it scanned.
#The goal is to bisect the fastq into 2 piles:  entries with IDs that were located by Barrnap (rRNA)
#-> and entries that weren't (mRNA), and then export them

#Aug 14, 2023: Changed to work only with FASTA, and not FASTQ.

def extract_rRNA_ID(barrnap_file):
    ID_list = set()
    barrnap_list = open(barrnap_file, mode='r')
    
    for item in barrnap_list:
        if(not item.startswith("#")):
            ID = ">" + item.split("\t")[0]
            #ID_list.add("@" + item.split("\t")[0])
            ID_list.add(ID)
    return list(ID_list)

def import_fasta(fasta_file):
    fasta_dict = dict()
    with open(fasta_file, "r") as fasta_in:
        fasta_ID = ""
        fasta_seq = ""
        for line in fasta_in:
            if(line.startswith (">")):
                if(fasta_ID != ""):
                    fasta_dict[fasta_ID] = fasta_seq
                    fasta_seq = ""
                fasta_ID = line.strip("\n")
            else:
                fasta_seq += line

    return fasta_dict


if __name__ == "__main__":
    
    barrnap_file    = sys.argv[1]
    fasta_seq_file  = sys.argv[2]
    mRNA_file   = sys.argv[3]
    junk_file  = sys.argv[4]

    fasta_dict = import_fasta(fasta_seq_file)
    fasta_ID_list = extract_rRNA_ID(barrnap_file)

    mRNA_export_path = mRNA_file
    junk_export_path = junk_file

    with open(mRNA_export_path, "w") as mRNA_out:
        for item in fasta_dict.keys():
            if(item in fasta_ID_list):
                mRNA_out.write(item + '\n')
                mRNA_out.write(fasta_dict[item])

    with open(junk_export_path, "w") as junk_out:
        for item in fasta_dict.keys():
            if not (item in fasta_ID_list):
                junk_out.write(item + "\n")
                junk_out.write(fasta_dict[item])

