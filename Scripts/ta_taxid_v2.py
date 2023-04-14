#!/usr/bin/env python
#This program extracts TaxID from a gene map file, and the accension file, and dumps it to one.
import sys
import os
from datetime import datetime as dt
def import_gene_map(gene2read_map_in):
    gene2read_dict = {}

    with open(gene2read_map_in, "r") as gene2read:
        for line in gene2read:
            if len(line) > 5:
                entry = line.split("\t")
                gene2read_dict[entry[0]] = entry[3:]
    return gene2read_dict
    
def import_genes(gene2read_dict):
    #collection of genes, extracted.
    #change this, later.  doesn't need a 2nd loop
    genes = {}

    for gene in gene2read_dict:
        try:
            accession = gene.split("|")[3].split(".")[0]
        except:
            accession = gene.split(".")[0]
        if accession in genes:
            genes[accession].append(gene)
        else:
            genes[accession] = [gene]
    return genes

def import_accession(accession2taxid_map_in):
    accession2taxid_dict = {}

    with open(accession2taxid_map_in, "r") as accession2taxid:
        for line in accession2taxid:
            columns = line.split("\t")
            if columns[0] in genes:
                accession2taxid_dict[columns[0]] = columns[2].strip("\n")
    return accession2taxid_dict

def get_read_taxa(gene2read_dict, genes, accession2taxid_dict):
    read2taxid_dict = {}

    for accession in genes:
        for gene in genes[accession]:
            for read in gene2read_dict[gene]:
                read = read.strip("\n")
                if read not in read2taxid_dict and accession in accession2taxid_dict:
                    read2taxid_dict[read] = accession2taxid_dict[accession]
                    
    return read2taxid_dict

if __name__ == "__main__":

    gene2read_map_in = sys.argv[1]
    accession2taxid_map_in = sys.argv[2]
    read2taxid_map_out = sys.argv[3]
    
    accession_map_size = os.path.getsize(accession2taxid_map_in)
    gene_map_size = os.path.getsize(gene2read_map_in)
    print("gene map size:", gene_map_size)
    print("accession map size:", accession_map_size)
    if(gene_map_size < 1):
        print(dt.today(), "no gene map available. exiting")
        sys.exit()
    else:
        print("things are running")
        gene2read_dict = import_gene_map(gene2read_map_in)
        genes = import_genes(gene2read_dict)
        accession2taxid_dict = import_accession(accession2taxid_map_in)

        read2taxid_dict = get_read_taxa(gene2read_dict, genes, accession2taxid_dict)

        

        #Write the results out
        with open(read2taxid_map_out, "w") as read2taxid:
            for read in read2taxid_dict:
                read2taxid.write("C" + "\t" + read + "\t" + read2taxid_dict[read] + "\n")