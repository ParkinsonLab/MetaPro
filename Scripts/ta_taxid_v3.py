#Sept 30, 2021:
#---------------------------------------
#This version is made to handle the changed formatting of Humann3's chocophlan database.
#Differences from humann2: The geneID has changed. The taxid is the first section in the gene ID. 
#followed by the whole taxa tree.
#goal here is to get every read with a taxid.  making mods to compensate for new design + keep old function

#!/usr/bin/env python

#!/usr/bin/env python
#This program extracts TaxID from a gene map file, and the accension file, and dumps it to one.
import sys
import os
from datetime import datetime as dt
import time

def import_gene_map(gene2read_map_in):
    gene2read_dict = {}
    check_db_flag = True
    with open(gene2read_map_in, "r") as gene2read:
        for line in gene2read:
            if len(line) > 5:
                entry = line.split("\t")
                gene2read_dict[entry[0]] = entry[3:]
                if(check_db_flag):
                    gene_id = entry[0]
                    gene_id_split = gene_id.split("|")
                    #humann3's chocophlan changed their format
                    if(len(gene_id_split) == 5):
                        check_db_flag = False
                        db_mode = "h3"
                    #humann2's format used to have 9 sections
                    elif(len(gene_id_split) == 9):
                        db_mode = "h2"
                        check_db_flag = False
                    #else it's a protein and it's inconclusive
                        
    
    
    return gene2read_dict, db_mode
    
def import_genes(gene2read_dict, op_mode):
    #collection of genes, extracted.
    #change this, later.  doesn't need a 2nd loop
    genes = {}
    accessioned_genes = {}
    
    if(op_mode == "h3"):
        #h3 mode, where the taxid is in the name, and there's no accession
        for gene in gene2read_dict:
            if("|k__" in gene):
                taxid_part = gene.split("|")[0]
                taxid = taxid_part.split("_")[0]
                if(taxid.isnumeric()):
                    if(taxid in genes):
                        genes[taxid].append(gene)
                        
                        #print(dt.today(), "adding existing gene to taxa_gene map[", taxid, "]",  gene)
                    else:
                        genes[taxid] = [gene]
                        #print(dt.today(), "adding new gene to taxa_gene map[", taxid, "]",  gene)
                else:
                    print(dt.today(), "this isn't supposed to happen. taxid not numeric. something wrong with format")
                    print(taxid, gene)
                    time.sleep(50)
                    sys.exit()
            else:            
                try:
                    accession = gene.split("||")[1]
                    print(dt.today(), "special protein found:", accession)
                    if(accession in accessioned_genes):
                        accessioned_genes[accession].append(gene)
                    else:
                        accessioned_genes[accession] = [gene]
                except:
                    accession = gene.split(".")[0]
                    #print(dt.today(), "regular protein found", accession)
                    if(accession in accessioned_genes):
                        accessioned_genes[accession].append(gene)
                    else:
                        accessioned_genes[accession] = [gene]
    else:
        #h2 mode
        for gene in gene2read_dict:
            try:
                accession = gene.split("|")[3].split(".")[0]
            except:
                accession = gene.split(".")[0]
            if accession in accessioned_genes:
                accessioned_genes[accession].append(gene)
            else:
                accessioned_genes[accession] = [gene]
    return genes, accessioned_genes

def import_accession(accession2taxid_map_in, accessioned_genes):
    accession2taxid_dict = {}

    with open(accession2taxid_map_in, "r") as accession2taxid:
        for line in accession2taxid:
            columns = line.split("\t")
            if columns[0] in accessioned_genes:
                accession2taxid_dict[columns[0]] = columns[2].strip("\n")
    return accession2taxid_dict

def get_read_taxa(gene2read_dict, accessioned_genes, taxa_genes, accession2taxid_dict):
    #there's 2 categories: genes assigned an accession only (from NCBI NR, and the updated chocophlan with a taxa in its geneID)
    read2taxid_dict = {}
    print(dt.today(), "converting accession-based entries")
    for accession in accessioned_genes:
        for gene in accessioned_genes[accession]:
            for read in gene2read_dict[gene]:
                read = read.strip("\n")
                if read not in read2taxid_dict and accession in accession2taxid_dict:
                    read2taxid_dict[read] = accession2taxid_dict[accession]
                    
    print(dt.today(), "converting taxa-based entries")
    for taxid in taxa_genes:
        for gene in taxa_genes[taxid]:
            for read in gene2read_dict[gene]:
                read = read.strip("\n")
                if(read not in read2taxid_dict):
                    read2taxid_dict[read] = taxid
    print(dt.today(), "done")
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
        print(dt.today(), "importing gene map")
        gene2read_dict, op_mode = import_gene_map(gene2read_map_in)
        print(dt.today(), "extracting taxa/accessions from gene map: ", op_mode)
        taxa_genes, accessioned_genes = import_genes(gene2read_dict, op_mode)
        print(dt.today(), "importing accession2taxid")
        accession2taxid_dict = import_accession(accession2taxid_map_in, accessioned_genes)
        
        print(dt.today(), "assigning reads to taxa from accession + gene ID")
        read2taxid_dict = get_read_taxa(gene2read_dict, accessioned_genes, taxa_genes, accession2taxid_dict)

        

        #Write the results out
        print(dt.today(), "exporting")
        with open(read2taxid_map_out, "w") as read2taxid:
            for read in read2taxid_dict:
                read2taxid.write("C" + "\t" + read + "\t" + read2taxid_dict[read] + "\n")
                
        print(dt.today(), "ALL DONE: ta_taxid_v3.py")