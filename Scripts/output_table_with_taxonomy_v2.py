#!/usr/bin/env python
#Mar 12, 2020:
#This exports the taxonomy table.  It now also takes into account the gene segments.

import sys
from datetime import datetime as dt

def translate_contig_segement_map(contig_segment_dict, contig_map_dict):
    for item in contig_segment_dict:
        contig_name = item.split("|")[1]
        contig_segment_percent = contig_segment_dict[item]
        contig_reads = contig_map_dict[contig_name]
        contig_segment_dict[item] = int(round(float(contig_reads) * contig_segment_percent))
    
def import_contig_map(contig_map_path):
    contig_map_dict = dict()
    with open(contig_map_path, "r") as contig_map:
        
        for line in contig_map:
            line_split = line.split("\t")
            contig_name = line_split[0]
            read_count = line_split[1]
            contig_map_dict[contig_name] = read_count
    return contig_map_dict

def import_gene_report(gene_report_path):
    #There's a better solution than to import the gene report + contig map, and that's to make a contig segment -> read count map.  But... baby steps, here.
    contig_segment_dict = dict()
    with open(gene_report_path, "r") as gene_report:
        contig_name = ""
        contig_segment_name_list = []
        start_of_loop = True
        contig_segment_sum = 0
        
        for line in gene_report:
            
            cleaned_line = line.strip("\n")
            cleaned_line = cleaned_line.split(" ")
            
            cleaned_line = [i for i in cleaned_line if i]
            if(len(cleaned_line) > 1):
                #print(cleaned_line)
                if(cleaned_line[0] == "FASTA"):
                    if(start_of_loop == False):
                        for item in contig_segment_name_list:
                            #print("contig segment length:", contig_segment_dict[item], "contig segment sum:", contig_segment_sum)
                            
                            contig_segment_dict[item] = contig_segment_dict[item] / contig_segment_sum
                        contig_segment_name_list[:] = []
                    contig_name = cleaned_line[3]
                    #print("contig:", contig_name)
                    start_of_loop = False
                    contig_segment_sum = 0
                if(cleaned_line[0] == "Model" or cleaned_line[0] == "Gene" or cleaned_line[0] == "#" or cleaned_line[0] == "Predicted"):
                    continue
                if(cleaned_line[0].isdigit()):
                    #print(cleaned_line)
                    contig_segment_name = "gene_" + cleaned_line[0] + "|" + contig_name
                    #print("contig segment name:", contig_segment_name)
                    contig_segment_sum += int(cleaned_line[4])
                    #print("contig segment sum:", contig_segment_sum)
                    contig_segment_name_list.append(contig_segment_name)
                    contig_segment_dict[contig_segment_name] = int(cleaned_line[4])
        for item in contig_segment_name_list:
        #    print("contig segment length:", contig_segment_dict[item], "contig segment sum:", contig_segment_sum)
            contig_segment_dict[item] = contig_segment_dict[item] / contig_segment_sum
        contig_segment_name_list[:] = []
        
    return contig_segment_dict

def import_classification_table(read2taxonomy):
    read2taxonomy_dict = {}
    with open(read2taxonomy, "r") as infile:
        for line in infile:
            cols = line.split("\t")
            read = cols[1]
            read2taxonomy_dict[read] = cols[2].strip("\n")

    return read2taxonomy_dict
    
def import_gene_map(gene2read):
    mapped_reads = 0
    gene2read_dict = {}
    with open(gene2read, "r") as infile:
        for line in infile:
            cols = line.split("\t")
            gene = cols[0]
            gene_len = cols[1]
            reads = []
            for read in cols[3:]:
                reads.append(read.strip("\n"))
            mapped_reads += len(reads)
            if gene in gene2read_dict:
                gene2read_dict[gene][1].extend(reads)
            else:
                gene2read_dict[gene] = (gene_len, reads)    
             
    return gene2read_dict
    
def import_gene_ec_map(gene2EC):    
    #this should export the EC info as 1 large string
    gene2EC_dict = {}
    with open(gene2EC, "r") as infile:
        for line in infile:
            cols = line.split("\t")
            gene = cols[0]
            EC_string = cols[2].strip("\n")
            gene2EC_dict[gene] = EC_string
    return gene2EC_dict
    
if __name__ == "__main__":
    gene2read               = sys.argv[1] #IN: the gene_map (genes to reads)
    gene2EC                 = sys.argv[2] #IN: proteins.EC all
    read2taxonomy           = sys.argv[3] #IN: constrain classification.tsv
    assembly_gene_report    = sys.argv[4] #IN: the MGM gene report
    contig_map              = sys.argv[5] #IN: contig map
    taxa_table              = sys.argv[6] #OUT: taxa table

    contig_map_dict = import_contig_map(contig_map)
    contig_segment_dict = import_gene_report(assembly_gene_report)
    translate_contig_segement_map(contig_segment_dict, contig_map_dict)

    read2taxonomy_dict = import_classification_table(read2taxonomy)

    gene2read_dict = import_gene_map(gene2read)
    
    gene2EC_dict = import_gene_ec_map(gene2EC)

    
    #for item in read2taxonomy_dict:
    #    print("keys:", item)
    #sys.exit("break")
            
    RPKM_dict = {}

    RPKM_taxonomy_dict = {}
    for gene in gene2read_dict:
        #RPKM_div = ((float(gene2read_dict[gene][0])/float(1000))*(mapped_reads/float(1000000)))
        gene_EC_val = ""
        
        
        if(gene in gene2EC_dict):
            gene_EC_val = gene2EC_dict[gene]
            if(gene_EC_val == ""):
                gene_EC_val = "0.0.0.0"
            print("gene_EC_val from list:", gene_EC_val)

        else:
            gene_EC_val = "0.0.0.0"
            
        taxon_count_dict = {} #key: taxon.  val: the read counts associated with the taxon (to accomodate for not knowing the reads inside the contig segments)
        for read in gene2read_dict[gene][1]:
            real_read = ""
            if(read.startswith("gene")):
                transformed_read_list = read.split("|")
                tail = transformed_read_list[-1].split(">")[-1]
                header = transformed_read_list[0]
                real_read = header + "|" + tail
            else:
                real_read = read
            
            if(read in read2taxonomy_dict):
                taxon = read2taxonomy_dict[read]
                if taxon in taxon_count_dict:
                    #taxon_count_dict[taxon].append(read)
                    if(real_read.startswith("gene")):
                        taxon_count_dict[taxon] += contig_segment_dict[real_read]
                    else:      
                        taxon_count_dict[taxon] += 1
                else:
                    if(real_read.startswith("gene")):
                        taxon_count_dict[taxon] = contig_segment_dict[real_read]
                    else:
                        taxon_count_dict[taxon] = 1
                    #taxon_count_dict[taxon].append(read)
                    #taxon_count_dict[taxon] = [read]
            else:
                print(dt.today(), read, "in gene map:", gene, "but no TA association")
        for taxon in taxon_count_dict:
            RPKM_taxonomy_dict[gene+"||"+taxon] = [gene, taxon, gene2read_dict[gene][0], taxon_count_dict[taxon], gene_EC_val]
            #RPKM_taxonomy_dict[gene+"||"+taxon].append(len(taxon_count_dict[taxon])/RPKM_div)

    with open(taxa_table, "w") as taxa_table_out:
        taxa_table_out.write("GeneID\tTaxonomy\tLength\tReads\tEC#\n")
        for entry in RPKM_taxonomy_dict:
            taxa_table_out.write("\t".join(str(x) for x in RPKM_taxonomy_dict[entry]) + "\n")

