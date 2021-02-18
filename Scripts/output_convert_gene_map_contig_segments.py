#converts the gene portions of the gene map back into constituent reads

import sys
import time



def import_contig_map(contig_map_file):
    contig_dict = dict()
    contig_count_dict = dict()
    
    with open(contig_map_file, "r") as contig_map_in:
        for line in contig_map_in:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            contig = line_split[0]
            read_count = line_split[1]
            reads = line_split[2:]
            contig_count_dict[contig] = read_count
            if(contig in contig_dict):
                print("contig already rep'd.  This shouldn't happen")
                sys.exit("duplicate contigs in contig map?????")
            else:
                contig_dict[contig] = reads
    return contig_dict, contig_count_dict

def import_gene_map(gene_map_file): 
    gene_dict = dict()
    gene_length_dict = dict()
    with open(gene_map_file, "r") as gene_map_in:
        for line in gene_map_in:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            gene = line_split[0]
            length = line_split[1]
            reads = line_split[3:]
            reads = [x for x in reads if x]
            #print("reads:", reads)
            
            if(gene in gene_dict):
                print("This shouldn't happen.  gene already scanned???")
                sys.exit("Bad gene map")
            else:
                gene_dict[gene] = reads
                gene_length_dict[gene] = length
    return gene_dict, gene_length_dict
            
            

if __name__ == "__main__":
    gene_map_file = sys.argv[1]
    contig_map_file = sys.argv[2]
    gene_map_export = sys.argv[3]
    
    contig_dict, contig_count_dict = import_contig_map(contig_map_file)
    gene_dict, gene_length_dict = import_gene_map(gene_map_file)
    
    final_count_dict = dict()
    final_map_dict = dict()
    
    for gene in gene_dict:
        cleaned_reads = gene_dict[gene]
        final_reads = []
        read_count = 0
        for read in cleaned_reads:
            
            if(read.startswith("gene")):
                if(read in contig_dict):
                    if(gene in final_count_dict):
                        #print("append contig count:", read, contig_count_dict[read])
                        read_count += float(contig_count_dict[read])
                        #time.sleep(1)
                    else:
                        #print("new contig count:", read, contig_count_dict[read])
                        read_count = float(contig_count_dict[read])
                        #time.sleep(1)
                        
                    contig_reads = contig_dict[read]
                    for contig_read in contig_reads:
                        final_reads.append(contig_read)
                else:
                    #print("contig didn't pass QC in contig map formation.  No reads assigned")
                    if(gene in final_count_dict):
                        read_count += 0
                    else:
                        read_count = 0
            else:
                final_reads.append(read)
                if(gene in final_map_dict):
                    read_count += 1
                else:
                    read_count = 1
            final_map_dict[gene] = final_reads
            #time.sleep(1)
            #print("===========================================")
            #print(gene, cleaned_reads, "looking at:", read)
            #print("read count:",read_count)
        final_count_dict[gene] = read_count   
        
                    
    with open(gene_map_export, "w") as out_file:
        header = "gene" + "\t" + "length" + "\t" + "count" + "\t" + "reads" + "\n"
        out_file.write(header)
        for item in final_map_dict:
            gene_length = str(gene_length_dict[item])
            gene_count = str(final_count_dict[item])
            if(float(gene_count) > 0): 
                out_line = item + "\t" + gene_length + "\t" + gene_count 
                reads = final_map_dict[item]
                for read in reads:
                    out_line += "\t" + read 
                out_line += "\n"
                out_file.write(out_line)
            
            
                
    