#This code performs all of the merging needed for GA to be completed.
#note: the contig gene->read map does not need converting
#oct 24, 2020: The paired gene map needs reconciliation
import os
import sys
import multiprocessing as mp
from datetime import datetime as dt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import pandas as pd

def export_proteins(diamond_proteins_file, gene_trans_dict, final_proteins):
    with open(diamond_proteins_file, "a") as diamond_proteins:
        for item in gene_trans_dict:
            SeqIO.write(item, diamond_proteins, "fasta")
            
def scrub_duplicates(fasta_in, fasta_out):
    #scrubs a fasta for duplicate entries.  
    imported_seqs = dict()
    header = "None"
    seq_body = "None"
    skip_flag = False
    #inspect_key = ">gi|483984714|ref|NZ_KB892660.1|:c111601-110108|33035|g__Blautia.s__Blautia_producta|UniRef90_unknown|UniRef50_R5TWD5"
    with open(fasta_in, "r") as in_seq:
        for line in in_seq:
            cleaned_line = line.strip("\n")
            
            if(cleaned_line.startswith(">")):
                if(header == "None"):
                    header = cleaned_line
                if(header in imported_seqs):
                    #header already imported.  skip it
                    header = "None"
                    skip_flag = True
                else:
                    skip_flag = False
                    if(seq_body != "None"):
                        #gene data gathered
                        imported_seqs[header] = seq_body
                        #if(header == inspect_key):
                        #    print(dt.today(), "inspect key found:", header)#, seq_body)
                            
                        
                        seq_body = "None"
                        header = cleaned_line
                        
            else:
                if(skip_flag):
                    seq_body = "None"
                    
                else:
                    if(seq_body == "None"):
                        seq_body = cleaned_line
                    else:
                        seq_body += cleaned_line
    
    with open(fasta_out, "w") as out_seq:
        for item in imported_seqs:
            out_line = item + "\n" + imported_seqs[item] + "\n"
            out_seq.write(out_line)    
            
def convert_genes_to_proteins(mapped_gene_file):#, section, gene_trans_dict):
    # WRITE OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs.  It's for a downstream tool.     
    # (.faa; fasta-format):
    # convert previously mapped genes to the AA format.
    #apparently seqIO can't handle duplicates in its file... it's kinda bullshit.
    
    
    gene_seqs = SeqIO.index(mapped_gene_file,"fasta")           # key=geneID, value=SeqRecord
    #inspect_key = "gi|483984714|ref|NZ_KB892660.1|:c111601-110108|33035|g__Blautia.s__Blautia_producta|UniRef90_unknown|UniRef50_R5TWD5"
    
    #gene_trans= []
    gene_trans_dict = dict()
    for gene in gene_seqs:                                  # Take each BWA&BLAT-aligned genes

        translated_seqRecord = SeqRecord(seq= gene_seqs[gene].seq.translate(stop_symbol=""), id= gene_seqs[gene].id, description= gene_seqs[gene].description)
        gene_id = gene
        gene_trans_dict[gene_id] = translated_seqRecord.seq                                                        #  and translate its SeqRecord sequence to aa.
        if(len(gene_trans_dict[gene_id]) * 3 > len(gene_seqs[gene].seq)):
            print(dt.today(), gene_id, "protein longer than bp seq.  protein;", len(gene_trans_dict[gene_id]), "bp:", len(gene_seqs[gene].seq))
            sys.exit("this is real bad")

    return gene_trans_dict
        

def concatenate_gene_maps_v2(path, gene_map_dict, context):
    #merge all gene maps given a directory
    #the gene length stays constant.  The number of reads do not.
    path_contents = os.listdir(path)
    duplicates_found = False
    #gene_map_dict = dict()
    for content in path_contents:
        item = os.path.join(path, content)
        base_name = os.path.basename(item)
        #print("base name:", base_name)
        if(base_name.endswith("_gene_map.tsv")):
            #print("context:", context, "checking:", item)
            #time.sleep(1)
            if(base_name.startswith(context)):
                #print("context:", context, "looking at:", item)
                with open(item, "r") as split_gene_map:
                    for line in split_gene_map:
                        line_split = line.strip("\n").split("\t")
                        number_of_reads = int(line_split[2])
                        gene_length = int(line_split[1])
                        gene_name = line_split[0]
                        reads = line_split[3:]
                  
                        if(gene_name in gene_map_dict):
                            old_pack = gene_map_dict[gene_name]
                            old_reads = old_pack[2:]
                            old_gene_length = int(old_pack[0])
                            
                            new_reads = list(set(old_reads + reads))
                            if(len(new_reads) < (len(old_reads) + len(reads))):
                                duplicates_found = True
                                print("---------------------------------")
                                print(dt.today(), gene_name, "duplicates removed", len(new_reads), "<-", len(old_reads) + len(reads))
                                print(gene_name, "new:", new_reads)
                                print(gene_name, "old:", old_reads)
                                print(gene_name, "added:", reads)
                                #time.sleep(3)
                            #print(gene_name, "old reads:", len(old_reads), "current reads:", len(reads), "combined:", len(new_reads))
                            #time.sleep(2)
                            
                            new_number_of_reads = len(new_reads)
                            final_gene_length = gene_length
                            if(int(old_gene_length) != int(gene_length)):
                                print("gene lengths for same gene in table do not match. This shouldn't happen")
                                print("taking larger length:", old_gene_length, gene_length)
                                print("for:", gene_name)
                                if(int(old_gene_length) > int(gene_length)):
                                    final_gene_length = old_gene_length
                                else:
                                    final_gene_length = gene_length
                            
                            gene_map_dict[gene_name] = [final_gene_length, new_number_of_reads] + new_reads
                            #print("merged entry:", gene_name, gene_map_dict[gene_name])
                        else:
                            
                            gene_map_dict[gene_name] = [gene_length, number_of_reads] + reads
    if(duplicates_found):
        print(dt.today(),"duplicates found")
    else:
        print(dt.today(), path, context, "OK! no duplicates found")
        
    
def merge_fastas(path_0, path_1, section, header, extension, bypass_ID_check = True):
    #This walks through each pairing, and checks for duplicates.
    #it has to be unique
    IDs_used = list()
    skip_this_line = False
    path_contents = os.listdir(path_0)
    final_fasta_file = os.path.join(path_1, header + "_" + section + extension)
    with open(final_fasta_file, "w") as final_fasta:
        for content in path_contents:
            item = os.path.join(path_0, content)
            
            if(content.startswith(section)) and (content.endswith(extension)):
                print("fasta merge:", item)
                if(not bypass_ID_check):
                    print(item, "bypassing ID check for merging raw data")
                time_spent_checking = 0
                time_writing = 0
                IDs_skipped = False
                ID_skip_count = 0
                with open(item, "r") as sample_fasta:
                    for line in sample_fasta:
                        if(line.startswith(">")):
                            ID = line
                            start_time = time.time()
                            if(bypass_ID_check):
                                if(ID in IDs_used):
                                    skip_this_line = True
                                    IDs_skipped = True
                                    ID_skip_count += 1
                                    #print("ID skipped:", item)
                                else:
                                    IDs_used.append(ID)
                                    skip_this_line = False
                                end_time = time.time()
                                check_time = end_time - start_time
                                time_spent_checking += check_time
                                #print("check time:", end_time - start_time)
                        if not(skip_this_line):
                            start_write = time.time()
                            final_fasta.write(line)
                            end_write = time.time()
                            time_writing = end_write - start_write
                        #else:
                        #    print(dt.today(), ID, "already written.  skipping")
                print(item, "time spent checking for overlap:", time_spent_checking)
                print(item, "time writing:", time_writing)
                if(IDs_skipped):
                    print(item, "IDs were skipped:", ID_skip_count)
                print("==================================")
    print(dt.today(), "find it at:", final_fasta_file)
    return final_fasta_file
    
    
def export_gene_map(gene_map, export_path, header = None):
    final_gene_map_file = ""
    if(header == None):
    
        final_gene_map_file = os.path.join(export_path, "gene_map.tsv")
    else:
        final_gene_map_file = os.path.join(export_path, header + "_gene_map.tsv")
    with open(final_gene_map_file, "w") as gene_map_out:
        #out_line = "geneID" + "\t" + "gene length" + "\t" + "Reads"
        for gene_name_key in gene_map:
            gene_entry = gene_map[gene_name_key]
            gene_ID = gene_name_key
            gene_length = gene_entry[0]
            number_of_reads = gene_entry[1]
            
            out_line = gene_ID + "\t" + str(gene_length) + "\t" + str(number_of_reads)
            for item in gene_entry[2:]:
                out_line += "\t" + item
            out_line += "\n"
            gene_map_out.write(out_line)
            
def make_merge_fasta_process(process_store, path_0, path_1, header, tail):
    section = ["pair_1", "pair_2", "contigs", "singletons"]
    
    for item in section:
        merge_fastas(path_0, path_1, item, header, tail)
        #merge_process = mp.Process(target = merge_fastas, args = (path_0, path_1, item, header, tail))
        #merge_process.start()
        #process_store.append(merge_process)
        #final_fasta_file = os.path.join(path_1, header + "_" + item + tail)

def make_merge_leftover_fasta_process(process_store, path_0, path_1, header, tail):
    section = ["contigs", "singletons"]
    
    for item in section:
        merge_fastas(path_0, path_1, item, header, tail, False)
        #merge_process = mp.Process(target = merge_fastas, args = (path_0, path_1, item, header, tail))
        #merge_process.start()
        #process_store.append(merge_process)
        #final_fasta_file = os.path.join(path_1, header + "_" + item + tail)
        
        
def merge_all_proteins(path, gene_transcripts_dict, export_path):
    #this just exports.
    all_proteins_path = os.path.join(export_path, "all_proteins.faa")
    skip_this_line = False
    IDs_used = list()
    with open(all_proteins_path, "w") as proteins_out:
        for item in os.listdir(path):
            
            if(item.endswith(".faa")):
                full_path = os.path.join(os.path.abspath(path), item)
                ID = "None"
                with open(full_path, "r") as dmd_proteins:
                    for line in dmd_proteins:
                    
                        if(line.startswith(">")):
                            ID = line
                            if(ID in IDs_used):
                                #print(dt.today(), "protein already written. skipping")
                                skip_this_line = True
                            else:
                                IDs_used.append(ID)
                                skip_this_line = False
                                #print("IDs used:", IDs_used)
                                #time.sleep(0.5)
                            
                        if not skip_this_line:
                            proteins_out.write(line)
                        else:
                            print("not writing:", ID)
                            #print(line)
                            #time.sleep(0.5)
        for gene in gene_transcripts_dict:
            #print("EXPORTING gene transcript:", gene, gene_transcripts_dict[gene])
            out_line = ">" + gene + "\n" 
            out_line += str(gene_transcripts_dict[gene]) + "\n"
            
            proteins_out.write(out_line)
        #for gene in gene_transcripts_list:
        #    out_line = str(gene.id) + "\n" + str(gene.seq) + "\n"
        #    proteins_out.write(out_line)
        #SeqIO.write(gene_transcripts_list, proteins_out, "fasta")
        
def handle_final_proteins(final_path, export_path):
    #just a function dump so we can do this step in parallel with other things

    #scrub for duplicates (because seqIO index doesn't like duplicates)
    bwa_blat_proteins_file = merge_fastas(final_path, final_path, "all", "BWA_BLAT_proteins", ".fna")
    unique_fna_file = "unique_" + os.path.basename(bwa_blat_proteins_file)
    real_unique_fna_file = os.path.join(final_path, unique_fna_file)
    scrub_duplicates(bwa_blat_proteins_file, real_unique_fna_file)
    print(dt.today(), real_unique_fna_file, "should have no duplicates")
    #sys.exit("Duplicates removed")
    gene_transcripts_dict = convert_genes_to_proteins(real_unique_fna_file)
    #print("GENE TRANSCRIPT:", len(gene_transcripts_list))
    #for item in gene_transcripts_list:
    #    print(item)
    
    merge_all_proteins(final_path, gene_transcripts_dict, export_path)
    
                    

                
def merge_dicts(list_of_dicts):
    #gene dict is: key: gene_name, val: [gene length, number of reads, reads]
    final_dict = {}
    
    for item in list_of_dicts:
        #if final dict is empty
        if(not bool(final_dict)):
            print("final dict empty")
            final_dict = item
            continue
        else:
            print("final dict not empty")
            for gene in item:
                if(gene in final_dict):
                    gene_length = item[gene][0]
                    
                    reads = item[gene][2:]
                    old_reads = final_dict[gene][2:]
                    new_reads = list(set(reads + old_reads))
                    new_entry = [gene_length, len(new_reads)] + new_reads
                    
                    if(len(new_reads) < (len(old_reads) + len(reads))):
                        print(gene, "dupes removed | old:", len(old_reads), "added:", len(reads), "combined:", len(new_reads))
                        print("new:", new_reads)
                        print("old:", old_reads)
                        print("added:", reads)
                        print("=====================================")
                    #else:
                    #    print(gene, "no dupes found | old:", len(old_reads), "added:", len(reads), "combined:", len(new_reads))
                    final_dict[gene] = new_entry
                    #time.sleep(1)
                    
                else:
                    #add new gene to map
                    final_dict[gene] = item[gene]
    return final_dict
    
def reconcile_paired_gene_map(pair_1_gene_map, pair_2_gene_map):
    #take all reads that mapped between these 2 sets.
    #for the purpose of final export
    all_reads = []
    for item in pair_1_gene_map:
        reads = pair_1_gene_map[item][2:]
        all_reads += reads
        
    for item in pair_2_gene_map:
        reads = pair_2_gene_map[item][2:]
        all_reads += reads
        
    all_reads = list(set(all_reads))
    return all_reads
    
def reconcile_paired_gene_map_v2(pair_1_gene_map, pair_2_gene_map, message):
    #we need to rebuild the gene maps of paired data
    read_details_dict = dict()
    final_gene_map = dict()
    repeat_reads = 0
    disagreements = 0
    new_read = 0
    print(dt.today(), "working on:", message)
    for p1_gene in pair_1_gene_map:
        reads = pair_1_gene_map[p1_gene][2:]
        gene_length = pair_1_gene_map[p1_gene][0]
        for read in reads:
            
            if(read in read_details_dict):
                print(dt.today(), "a repeat read in pair1.  This shouldn't happen")
                sys.exit("death")
            else:
                score = 0
                real_read_name = "none"
                if("<bitscore>" in read):
                    real_read_name = read.split("<bitscore>")[0]
                    score = read.split("<bitscore>")[1]
                elif("<AS_score>" in read):
                    real_read_name = read.split("<AS_score>")[0]
                    score = read.split("<AS_score>")[1]
                
                inner_dict = dict()
                inner_dict["gene"] = p1_gene
                inner_dict["score"] = score
                inner_dict["gene_length"] = gene_length
                read_details_dict[real_read_name] = inner_dict
            
    for p2_gene in pair_2_gene_map:
        gene_length = pair_2_gene_map[p2_gene][0]
        reads = pair_2_gene_map[p2_gene][2:]
        for read in reads:
            score = 0
            real_read_name = "none"
            if("<bitscore>" in read):
                real_read_name = read.split("<bitscore>")[0]
                score = read.split("<bitscore>")[1]
            elif("<AS_score>" in read):
                real_read_name = read.split("<AS_score>")[0]
                score = read.split("<AS_score>")[1]
            
            if(real_read_name in read_details_dict):
                repeat_reads += 1
                #read already processed
                old_details = read_details_dict[real_read_name]
                old_score = old_details["score"]
                if(old_score < score):
                    #disagreement found
                    disagreements += 1
                    inner_dict = dict()
                    inner_dict["gene"] = p2_gene
                    inner_dict["score"] = score
                    inner_dict["gene_length"] = gene_length
                    read_details_dict[real_read_name] = inner_dict
                    
               
            else:
                #new p2-only read
                new_read += 1
                inner_dict = dict()
                inner_dict["gene"] = p2_gene
                inner_dict["score"] = score
                inner_dict["gene_length"] = gene_length
                read_details_dict[real_read_name] = inner_dict
        
    #make the new gene map
    all_paired_reads = []
    for read in read_details_dict:
        inner_dict = read_details_dict[read]
        gene_name = inner_dict["gene"]
        gene_length = inner_dict["gene_length"]
        all_paired_reads.append(read)
        if(gene_name in final_gene_map):
            old_reads = final_gene_map[gene_name][2:]
            old_gene_count = final_gene_map[gene_name][1]
            old_reads.append(read)
            #print(dt.today(), "old final gene map:", gene_name, final_gene_map[gene_name])
            final_gene_map[gene_name] = [gene_length, old_gene_count + 1] +  old_reads
            #print(dt.today(), "new final_gene_map:", gene_name, final_gene_map[gene_name])
            
            #
        else:
            final_gene_map[gene_name] = [gene_length, 1, read]
            #print("new gene entry:", gene_name, final_gene_map[gene_name])
            
            #time.sleep(1)
    print(message, "repeat:", repeat_reads)
    print(message, "disagreements:", disagreements)
    print(message, "new reads:", new_read)
    
    if(message == "cs"):
        for item in final_gene_map:
            if(type(final_gene_map[item]) is not list):
                print(message, item,  type(final_gene_map[item]))
                time.sleep(1)
            
    return final_gene_map, all_paired_reads
        
                

def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df    
    
    
def export_leftover_paired_reads(all_paired_reads, pair_1_df, pair_2_df, export_path):
    pair_1_leftover = pair_1_df[~pair_1_df["ID"].isin(all_paired_reads)]
    pair_1_leftover["ID"] = "@" + pair_1_leftover["ID"]
    
    
    pair_2_leftover = pair_2_df[~pair_2_df["ID"].isin(all_paired_reads)]
    pair_2_leftover["ID"] = "@" + pair_2_leftover["ID"]
    
    pair_1_export_path = os.path.join(export_path, "GA_leftover_pair_1.fastq")
    pair_2_export_path = os.path.join(export_path, "GA_leftover_pair_2.fastq")
    
    pair_1_leftover.to_csv(pair_1_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
    pair_2_leftover.to_csv(pair_2_export_path, mode = "w", sep= "\n", header = False, index = False, quoting = 3)
    
    
def make_second_merge_process(process_store, path_0, path_1, header, tail):
    section = ["BWA_annotated", "BLAT_annotated"]
    for item in section:
        merge_fastas(path_0, path_1, item, header, tail)
        #merge_process = mp.Process(target = merge_fastas, args = (path_0, path_1, item, header, tail))
        #merge_process.start()
        #process_store.append(merge_process)    
        
        
def final_gene_map_merge(gene_map_0, gene_map_1):
    final_gene_map = gene_map_0
    for gene in gene_map_1:
        if(gene in final_gene_map):
            new_reads = gene_map_1[gene][2:]
            old_reads = gene_map_0[gene][2:]
            gene_length = final_gene_map[gene][0]
            combined_reads = old_reads + new_reads
            final_gene_map[gene] = [gene_length, len(combined_reads)] + combined_reads
            
        else:
            final_gene_map[gene] = gene_map_1[gene]
            
    return final_gene_map
    
    
    
if __name__ == "__main__":
    assemble_path           = sys.argv[1]
    bwa_path                = sys.argv[2]
    blat_path               = sys.argv[3]
    diamond_path            = sys.argv[4]
    final_path              = sys.argv[5]
    export_path             = sys.argv[6]
    operating_mode          = sys.argv[7]
    
    if(operating_mode == "single"):
        manager = mp.Manager()
        mgr_bwa_singletons_gene_map = manager.dict()
        mgr_bwa_contig_gene_map     = manager.dict()

        mgr_blat_singletons_gene_map= manager.dict()
        mgr_blat_contig_gene_map    = manager.dict()
        
        mgr_dia_singletons_gene_map = manager.dict()
        mgr_dia_contig_gene_map     = manager.dict()
        
        process_store = []

        #merge the annotated genes (for protein translation), and leftover contig + singletons
        make_merge_leftover_fasta_process(process_store, diamond_path, export_path, "GA_leftover", ".fasta")  #leftover reads
        make_merge_fasta_process(process_store, blat_path, final_path, "BLAT_annotated", ".fna")    #BLAT genes
        make_merge_fasta_process(process_store, bwa_path, final_path, "BWA_annotated", ".fna")      #BWA genes
        make_merge_fasta_process(process_store, diamond_path, final_path, "dmd", ".faa")            #DIAMOND proteins
        
        #merge the gene maps, but group them by data context (pair1, pair2, singletons, contigs) for each tool (BWA, BLAT, DIAMOND)
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        #----------------------------------------------------------------------
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        #----------------------------------------------------------------------
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        
        
        print(dt.today(), "concatenating gene maps")
        for item in process_store:
            item.join()
        process_store[:] = []
        
        #-----------------------------------------------------------------------
        #convert to a standard dict
        bwa_singleton_gene_map  = dict(mgr_bwa_singletons_gene_map)
        bwa_contig_gene_map     = dict(mgr_bwa_contig_gene_map)
        
        blat_singleton_gene_map = dict(mgr_blat_singletons_gene_map)
        blat_contig_gene_map    = dict(mgr_blat_contig_gene_map)
        
        dia_singleton_gene_map  = dict(mgr_dia_singletons_gene_map)
        dia_contig_gene_map     = dict(mgr_dia_contig_gene_map)
        
        singleton_gene_map_list = [bwa_singleton_gene_map, blat_singleton_gene_map, dia_singleton_gene_map]
        contig_gene_map_list    = [bwa_contig_gene_map, blat_contig_gene_map, dia_contig_gene_map]
        
        #merge the gene maps by category
        singletons_gene_map = merge_dicts(singleton_gene_map_list)
        contig_gene_map     = merge_dicts(contig_gene_map_list)
        
        #merge all gene maps
        final_gene_map, other_reads = reconcile_paired_gene_map_v2(contig_gene_map, singletons_gene_map, "cs")    
        
        print(dt.today(), "done converting")
        
        #secondary combine on the genes (BWA and BLAT)
        make_second_merge_process(process_store, final_path, final_path, "all", ".fna")
        
        for item in process_store:
            item.join()
        process_store[:] = []
        
        #convert genes to proteins, and merge with diamond's export
        handle_final_proteins(final_path, export_path)
        
        export_gene_map(final_gene_map, export_path)
        
        print(dt.today(), "We're at the end")
        
        
    elif(operating_mode == "paired"):
        pair_1_raw_df = import_fastq(os.path.join(assemble_path, "pair_1.fastq"))
        pair_2_raw_df = import_fastq(os.path.join(assemble_path, "pair_2.fastq"))
        
        manager = mp.Manager()
        mgr_bwa_pair_1_gene_map     = manager.dict()
        mgr_bwa_pair_2_gene_map     = manager.dict()
        mgr_bwa_singletons_gene_map = manager.dict()
        mgr_bwa_contig_gene_map     = manager.dict()

        mgr_blat_pair_1_gene_map    = manager.dict()
        mgr_blat_pair_2_gene_map    = manager.dict()
        mgr_blat_singletons_gene_map= manager.dict()
        mgr_blat_contig_gene_map    = manager.dict()
        
        mgr_dia_pair_1_gene_map     = manager.dict()
        mgr_dia_pair_2_gene_map     = manager.dict()
        mgr_dia_singletons_gene_map = manager.dict()
        mgr_dia_contig_gene_map     = manager.dict()
        
        
        process_store = []
        
        
        #merge the annotated genes (for protein translation), and leftover contig + singletons
        make_merge_leftover_fasta_process(process_store, diamond_path, export_path, "GA_leftover", ".fasta")  #leftover reads
        make_merge_fasta_process(process_store, blat_path, final_path, "BLAT_annotated", ".fna")    #BLAT genes
        make_merge_fasta_process(process_store, bwa_path, final_path, "BWA_annotated", ".fna")      #BWA genes
        make_merge_fasta_process(process_store, diamond_path, final_path, "dmd", ".faa")            #DIAMOND proteins
        
        
        #merge the gene maps, but group them by data context (pair1, pair2, singletons, contigs) for each tool (BWA, BLAT, DIAMOND)
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_pair_1_gene_map, "pair_1"))
        process.start()
        process_store.append(process)
        
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_pair_2_gene_map, "pair_2"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (bwa_path, mgr_bwa_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        #----------------------------------------------------------------------
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_pair_1_gene_map, "pair_1"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_pair_2_gene_map, "pair_2"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (blat_path, mgr_blat_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        #----------------------------------------------------------------------
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_pair_1_gene_map, "pair_1"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_pair_2_gene_map, "pair_2"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_singletons_gene_map, "singletons"))
        process.start()
        process_store.append(process)
        
        process = mp.Process(target = concatenate_gene_maps_v2, args = (diamond_path, mgr_dia_contig_gene_map, "contig"))
        process.start()
        process_store.append(process)
        
        
        
        print(dt.today(), "concatenating gene maps")
        for item in process_store:
            item.join()
        process_store[:] = []
        
        #-----------------------------------------------------------------------
        #convert to a standard dict
        bwa_pair_1_gene_map     = dict(mgr_bwa_pair_1_gene_map)
        bwa_pair_2_gene_map     = dict(mgr_bwa_pair_2_gene_map)
        bwa_singleton_gene_map  = dict(mgr_bwa_singletons_gene_map)
        bwa_contig_gene_map     = dict(mgr_bwa_contig_gene_map)
        
        blat_pair_1_gene_map    = dict(mgr_blat_pair_1_gene_map)
        blat_pair_2_gene_map    = dict(mgr_blat_pair_2_gene_map)
        blat_singleton_gene_map = dict(mgr_blat_singletons_gene_map)
        blat_contig_gene_map    = dict(mgr_blat_contig_gene_map)
        
        dia_pair_1_gene_map     = dict(mgr_dia_pair_1_gene_map)
        dia_pair_2_gene_map     = dict(mgr_dia_pair_2_gene_map)
        dia_singleton_gene_map  = dict(mgr_dia_singletons_gene_map)
        dia_contig_gene_map     = dict(mgr_dia_contig_gene_map)
        
        
        pair_1_gene_map_list    = [bwa_pair_1_gene_map, blat_pair_1_gene_map, dia_pair_1_gene_map]
        pair_2_gene_map_list    = [bwa_pair_2_gene_map, blat_pair_2_gene_map, dia_pair_2_gene_map]
        singleton_gene_map_list = [bwa_singleton_gene_map, blat_singleton_gene_map, dia_singleton_gene_map]
        contig_gene_map_list    = [bwa_contig_gene_map, blat_contig_gene_map, dia_contig_gene_map]
        
        #merge the gene maps by category
            
        
        
        pair_1_gene_map     = merge_dicts(pair_1_gene_map_list)
        pair_2_gene_map     = merge_dicts(pair_2_gene_map_list)
        singletons_gene_map = merge_dicts(singleton_gene_map_list)
        contig_gene_map     = merge_dicts(contig_gene_map_list)
        
        #reconcile the paired reads (For export)
        paired_gene_map, all_paired_reads = reconcile_paired_gene_map_v2(pair_1_gene_map, pair_2_gene_map, "paired only")
        export_leftover_paired_reads(all_paired_reads, pair_1_raw_df, pair_2_raw_df, export_path)
        
        
        #merge all gene maps
        contigs_and_singletons_gene_map, other_reads = reconcile_paired_gene_map_v2(contig_gene_map, singletons_gene_map, "cs")    
        
        final_gene_map = final_gene_map_merge(contigs_and_singletons_gene_map, paired_gene_map)
        
        print(dt.today(), "done converting")
        
        
        #secondary combine on the genes (BWA and BLAT)
        make_second_merge_process(process_store, final_path, final_path, "all", ".fna")
        
        for item in process_store:
            item.join()
        process_store[:] = []
        
        
        #convert genes to proteins, and merge with diamond's export
        handle_final_proteins(final_path, export_path)
        
        export_gene_map(final_gene_map, export_path)
        
        print(dt.today(), "We're at the end")
        

            