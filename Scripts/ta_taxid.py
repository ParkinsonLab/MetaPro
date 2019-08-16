#!/usr/bin/env python
#This program extracts TaxID from a gene map file, and the accension file, and dumps it to one.
import sys

gene2read_map_in = sys.argv[1]
accession2taxid_map_in = sys.argv[2]
read2taxid_map_out = sys.argv[3]

gene2read_dict = {}

with open(gene2read_map_in, "r") as gene2read:
    for line in gene2read:
        if len(line) > 5:
            entry = line.split("\t")
            gene2read_dict[entry[0]] = entry[3:]
            
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

accession2taxid_dict = {}

with open(accession2taxid_map_in, "r") as accession2taxid:
    for line in accession2taxid:
        columns = line.split("\t")
        if columns[0] in genes:
            accession2taxid_dict[columns[0]] = columns[2].strip("\n")

read2taxid_dict = {}

for accession in genes:
    for gene in genes[accession]:
        for read in gene2read_dict[gene]:
            read = read.strip("\n")
            if read not in read2taxid_dict and accession in accession2taxid_dict:
                read2taxid_dict[read] = accession2taxid_dict[accession]

#Write the results out
with open(read2taxid_map_out, "w") as read2taxid:
    for read in read2taxid_dict:
        read2taxid.write("C" + "\t" + read + "\t" + read2taxid_dict[read] + "\n")