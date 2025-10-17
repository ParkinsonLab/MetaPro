#this code just makes the final stats reports
#not much logic here
#changed to now operate on a single batch.  expected to loop and append to file

import sys
import os
import pandas as pd
from datetime import datetime as dt
def fastq_count(item):
    lines = 0
    if(os.path.exists(item)):
        with open(item, "r") as infile:
            for line in infile:
                lines += 1
            
    if(lines % 4 != 0):
        print(dt.today(), "Don't use this number. Error in fastq counting in:", item, lines / 4)
        
    return lines / 4
    


def annotated_count(map):
    annotated_mRNA = 0
    read_list = []

    genes = 0
    with open(map, "r") as infile:
        for line in infile:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            reads = line_split[3:]
            read_count_entry = line_split[2]
            
            genes += 1
            if (not read_list):
                read_list = reads
            else:
                read_list += reads
    #there will be some reads that annotate to multiple genes.  This only counts unique reads.            
    annotated_mRNA = len(set(read_list))            
    
    return annotated_mRNA, genes

#format's changed.  this needs changing too
def ec_count(map):
    ecs = set()
    with open(map, "r") as infile:
        for line in infile:
            #ecs.add(line.split("\t")[2].strip())
            ec_line_list = line.split("\t")
            ec_portion = ec_line_list[2].strip("\n")
            ec_list = ec_portion.split("|")
            for ec in ec_list:
                ecs.add(ec)
    return len(ecs)

def check_paired_data(p1, p2, message):
    p1_count = fastq_count(p1)
    p2_count = fastq_count(p2)
    if(p1_count == p2_count):
        return True
    else:
        print(dt.today(), "bad data in:", message)
        sys.exit()


if __name__ == "__main__":
    raw_sequence = sys.argv[1]
    s_seq = sys.argv[2]
    p1_seq = sys.argv[3]
    p2_seq = sys.argv[4]    
    sample_header = sys.argv[5] #qf, hosts, vec, etc.etc.
    gene_map_location   = sys.argv[6]
    ec_location         = sys.argv[7]
    output_file         = sys.argv[8]
    operating_mode      = sys.argv[9]

    qc_s = ""
    qc_p1 = ""
    qc_p2 = ""

    qc_p1_unique = ""
    qc_p2_unique = ""
    qc_s_unique = ""
    host_p1 = ""
    host_p2 = ""
    host_2 = ""
    vectors_p1 = ""
    vectors_p2 = ""
    vectors_s = ""
    rRNA_p1 = ""
    rRNA_p2 = ""
    rRNA_s = ""
    mRNA_p1 = ""
    mRNA_p2 = ""
    mRNA_s = ""
    
    headings = []
    data = []

    if(sample_header == "raw"):
        headings.append("Total reads")
        raw_sequence_count = fastq_count(raw_sequence)
        data.append(str(int(raw_sequence_count)))
    elif(sample_header == "qf"):
        #qc_s = ""
        #if(operating_mode == "single"):
        #    qc_s = #os.path.join(quality_location, "singletons_hq.fastq")
        #else:
        #    qc_s       = os.path.join(quality_location, "singletons_with_duplicates.fastq")
        qc_s = s_seq
        qc_p1 = p1_seq#os.path.join(quality_location, "pair_1_match.fastq")
        qc_p2 = p2_seq#os.path.join(quality_location, "pair_2_match.fastq")

        check_paired_data(qc_p1, qc_p2, "quality")
    elif(sample_header == "qf_u"):
        #qc_p1_unique    = os.path.join(quality_location, "pair_1.fastq")
        #qc_p2_unique    = os.path.join(quality_location, "pair_2.fastq")
        #qc_s_unique     = os.path.join(quality_location, "singletons.fastq")
    
        qc_s_unique = s_seq
        qc_p1_unique = p1_seq
        qc_p2_unique = p2_seq

        headings.append("High quality reads")
        quality_sequence_count = fastq_count(qc_p1) + fastq_count(qc_s)
        data.append(str(int(quality_sequence_count)))

    elif(sample_header.contains("host")):
        #host_p1         = os.path.join(host_location, "pair_1_full_hosts.fastq")
        #host_p2         = os.path.join(host_location, "pair_2_full_hosts.fastq")
        #host_s          = os.path.join(host_location, "singletons_full_hosts.fastq")    

        host_s = s_seq
        host_p1 = p1_seq
        host_p2 = p2_seq

    elif(sample_header == "vec"):
        #vectors_p1      = os.path.join(vectors_location, "pair_1_full_vectors.fastq")
        #vectors_p2      = os.path.join(vectors_location, "pair_2_full_vectors.fastq")
        #vectors_s       = os.path.join(vectors_location, "singletons_full_vectors.fastq")
        vectors_s = s_seq
        vectors_p1 = p1_seq
        vectors_p2 = p2_seq



    elif(sample_header == "rRNA"):
        #rRNA_p1         = os.path.join(repop_location, "pair_1_rRNA.fastq")
        #rRNA_p2         = os.path.join(repop_location, "pair_2_rRNA.fastq")
        #rRNA_s          = os.path.join(repop_location, "singletons_rRNA.fastq")
        rRNA_s = s_seq
        rRNA_p1 = p1_seq
        rRNA_p2 = p2_seq


    elif(sample_header == "mRNA"):
        #mRNA_p1         = os.path.join(repop_location, "pair_1.fastq")
        #mRNA_p2         = os.path.join(repop_location, "pair_2.fastq")
        #mRNA_s          = os.path.join(repop_location, "singletons.fastq")
        mRNA_s = s_seq
        mRNA_p1 = p1_seq
        mRNA_p2 = p2_seq

        headings.append("Annotated mRNA reads")
        annotated_mRNA_count, genes_count = annotated_count(gene_map_location)
        data.append(str(int(annotated_mRNA_count)))
    elif(sample_header == "prot"):

    
    

        
        
    #gene_to_read_map = os.path.join(gene_map_location, "gene_map.tsv")
    lq_ec_map = os.path.join(ec_location, "lq_proteins.ECs_All")
    ec_map = os.path.join(ec_location, "proteins.ECs_All")
    
    
    
    
    check_paired_data(host_p1, host_p2, "host")
    check_paired_data(rRNA_p1, rRNA_p2, "rRNA+tRNA")
    check_paired_data(mRNA_p1, mRNA_p2, "putative_mRNA")
    check_paired_data(vectors_p1, vectors_p2, "vectors")
    

    

    

    

    headings.append("% high quality")
    quality_sequence_pct = quality_sequence_count / raw_sequence_count
    data.append("%.2f" % (quality_sequence_pct*100))
    
    headings.append("host reads found in sample")
    host_read_counts = fastq_count(host_p1) + fastq_count(host_s)
    data.append(str(int(host_read_counts)))
    
    headings.append("% host reads in sample")
    host_pct = host_read_counts / raw_sequence_count
    data.append("%.2f" % (host_pct * 100))
    
    headings.append("vector reads found in sample")
    vectors_read_counts = fastq_count(vectors_p1) + fastq_count(vectors_s)
    data.append(str(int(vectors_read_counts)))
    
    headings.append("% vector reads in sample")
    vectors_pct = vectors_read_counts / raw_sequence_count
    data.append("%.2f" % (vectors_pct * 100))

    headings.append("rRNA + tRNA reads")
    rRNA_sequence_count = fastq_count(rRNA_p1) + fastq_count(rRNA_s)
    data.append(str(int(rRNA_sequence_count)))

    headings.append("% rRNA + tRNA reads")
    rRNA_sequence_pct = rRNA_sequence_count / raw_sequence_count
    data.append("%.2f" % (rRNA_sequence_pct*100))

    headings.append("Putative mRNA reads")
    mRNA_sequence_count = fastq_count(mRNA_p1) + fastq_count(mRNA_s)
    data.append(str(int(mRNA_sequence_count)))

    headings.append("% putative mRNA reads")
    mRNA_sequence_pct = mRNA_sequence_count / raw_sequence_count
    data.append("%.2f" % (mRNA_sequence_pct*100))

    

    headings.append("% of putative mRNA reads annotated")
    annotated_mRNA_pct = annotated_mRNA_count / mRNA_sequence_count
    data.append("%.2f" % (annotated_mRNA_pct*100))

    headings.append("Unique transcripts")
    data.append(str(int(genes_count)))

    headings.append("High-Quality unique enzymes")
    unique_ec_count = ec_count(ec_map)
    data.append(str(int(unique_ec_count)))
    
    headings.append("Low-Quality unique enzymes")
    unique_ec_count = ec_count(lq_ec_map)
    data.append(str(int(unique_ec_count)))

    with open(output_file, "a") as outfile:
        outfile.write("\t".join(headings))
        outfile.write("\n")
        outfile.write("\t".join(data))
