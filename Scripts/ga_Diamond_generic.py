#!/usr/bin/env python
#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - added .strip("\n") to extracting from the gene2read_file
# - fixed sort function to return float
# - reversed the sort direction to get largest scores at top
# - pulled gene_map() out of the readtype loop
# - rewrote gene_map()
# - removed checks to see if contig/read was already assigned to same gene (shouldn't happen)
# - removed check to see if query=read was already mapped by DMD (shouldn't happen)
# - multiplied align_len*3 to convert aa->nt
# - changed align_len to float(align_len) (to make sure it's a number)
# - stores aligned proteins in new dict prot2read_map, instead of gene2read_map:
#   readability esp. during WRITE OUTPUT: BWA&BLAT&DMD-aligned;
#   also only check new DMD alignements against protein list for duplicates
# - paired-end DMD outputs are combined and analyzed at the same time
# - fixed WRITE OUTPUT: non-BWA&BLAT&DMD-aligned: changed "mapped_reads.add" to "break"

# NOTE:
# - Filenames for unmerged paired-end reads must be specified last.
# - Sometimes reads will be matched to multiple genes. This occurs between
#   BLAT contig alignement and DMD contig alignment.
# - DMDpp only takes the read<->gene match with the top score. For duplicate
#   top scores, it only takes the one that shows up first in the initial sort.
# - The paired-end DMD outputs are combined and analyzed at the same time, in order to
#   rank the matches between different ends of the same read, and not double count that read.
# - Sometimes a read is aligned via multiple contigs and annoteted only according
#   to the best matched contig only. Theoretically, this breaks up the other contigs.

#OTHER NOTES: Feb 27, 2019
# - this is a messy Piece-of-shit code that needs to be written without those dumb import loops.
# Aug 26, 2019:  it's much cleaner now, but there's a possibility for empty-files to be created, which interferes with the error-check of the pipe.


import os
import os.path
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime as dt
from shutil import copyfile
import pandas as pd
import itertools
import multiprocessing as mp


#####################################
# FUNCTIONS:
def check_file_safety(file_name):
    if(os.path.exists(file_name)):
        if(os.path.getsize(file_name) == 0):
            #copyfile(gene2read_file, new_gene2read_file)
            print(file_name, "is unsafe.  skipping")
            return False
            #sys.exit(DMD_tab_file_1 + " -> DMD tab file 1 is empty.  aborting")
        else:
            print(file_name, "exists and is safe")
            return True
    else:
        #sys.exit(DMD_tab_file_1 + " -> DMD tab file 1 is missing.  aborting")
        print(file_name, "is missing. skipping")
        return False



# read .dmdout file and acquire read length:
def get_dmd_hit_details(dmd_out_file, reads_in):
    seqrec = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:]) 
#   List of lists [dmdout info, read length]
    # get info from .blatout file:
    print (dt.today(), 'Reading ' + str(os.path.basename(dmd_out_file)) + '.')
    with open(dmd_out_file,"r") as dmd_out:
        hits= []                                        # List of lists containing .blatout fields.
        for line in dmd_out:                            # In the .dmdout file:
            if len(line)>=2:                            # If length of line >= 2,
                info_list= line.strip("\n").split("\t") #  make a list from .dmdout tab-delimited fields,
                read_id = info_list[0]                     #  get the queryID= contig/readID, -> 
                info_list.append(len(seqrec[read_id].seq))#  add query length (int) to end of list,
                hits.append(info_list)                  #  and append the list to the hits list.
    print(dt.today(), "finished read-aligning")
    return hits

# sort function: sort by score:
# (12th field of the .dmdout file)
def sortbyscore(line):
    return line[11]


def initial_dmdout_filter(hits):
    print(dt.today(), "working on initial dmdout filter")
    #sorts and sifts out stuff from the dmdout
    sorted_hits = sorted(hits, key=sortbyscore, reverse=True)     # sort alignment list by high score:
    del hits
        
    # DMD threshold:
    identity_cutoff= 85
    length_cutoff= 0.65
    score_cutoff= 60
    
    # tracking DMD-assigned & unassigned:
    #query2prot_map = defaultdict(set)        # Dict of DMD-aligned contig/readID<->protID(s).
    #queryIScontig = {}                       # Dict of DMD-aligned contig/readID<->contig? (True/False).
    unmapped = set()                         # Set of unmapped contig/readIDs.
    
    query_details_dict = dict() #dict of dict: level 1: query.  level 2: contig-bool, protein
    #Only going to be at most 1 protein per query. 
    #Also query is a catch-all term for read-ID (because read-ID is assumed to come from the raw data, but this one also includes contig-IDs.  We're keeping the term to avoid confusion)
    
    # loop through sorted DMD-aligned reads/contigs:
    for line in sorted_hits:
        inner_dict = dict() 
        
        # extract & store data:
        query = line[0]                      # query= contig/readID
        db_match = line[1]                   # proteinID
        seq_identity = float(line[2])        # sequence identity
        align_len = 3*int(line[3])           # alignment length (aa->nt)
        score = float(line[11])              # score
        seq_len = int(line[12])              # query sequence length
        
        # is this alignment the highest-score match for that query?
        
        #print("---------------------------------------------")
        if query in query_details_dict:
            continue                        # skip to next qurey as this alignment will have a lower score,
        
            
        # is query a contig?:
        contig_flag = False 
        if query in contig2read_map:        # If query is found in contig list,
            contig_flag = True                    #  then mark as contig,
        
        # does query alignment meet threshold?: # (identity, length, score):
        if seq_identity<=identity_cutoff or align_len<=seq_len*length_cutoff or score<=score_cutoff:
            unmapped.add(query)             # if threshold is failed, add to unmapped set and
            continue                        # skip to the next query.
   
        # store info for queries that remain:
        inner_dict["is_contig"] = contig_flag
        inner_dict["protein"] = db_match
        query_details_dict[query] = inner_dict
        
            
    print(dt.today(), "done initial filtering")        
    return unmapped, query_details_dict
    
# add DMD-aligned reads that meet threshold
# to the aligned gene/protID<->readID(s) dict:
def form_prot_map(hits, mapped_reads, contig2read_map, prot2read_map): 
    #contig2read map is external.  should contain upstream runs' contigs
    #prot2read map is external.  should contain prior calls' proteins
    
    unmapped, query_details_dict = initial_dmdout_filter(hits) #unmapped, query2prot_map, queryIScontig = initial_dmdout_filter(hits)
    # EXPAND HERE to deal with multiple high-score genes for a read.
    print(dt.today(), "finished initial sorting")  
    
    #checking
    if(len(prot2read_map) == 0):
        print(dt.today(), "prot2read map is empty")
    else:
        print(dt.today(), "prot2read map is not empty")
        
    if(len(contig2read_map) == 0):
        print(dt.today(), "contig2read map is empty")
    else:
        print(dt.today(), "contig2read map not empty")
        

      
    # DEBUG:
    contigread_inmapped = 0
    contigread_inmapped_inprot = 0
    contigread_inprot = 0
    read_inmapped = 0
    read_inmapped_inprot = 0
    read_inprot = 0
    
    query_keys = list(query_details_dict.keys())#query2prot_map.keys())
    print(dt.today(), "read ID 2 prot map keys:", len(query_keys))
 
    print(dt.today(), "sifting through read-ID -> protein map")
    # FINAL remaining DMD-aligned queries:
    #for query in query2prot_map:                        # contig_ID or read_ID.  DIAMOND spits out queries.
    for query in query_keys:
        db_match = query_details_dict[query]["protein"]#list(query2prot_map[query])[0]        # geneID (pull out of 1-element set)
        contig_flag = query_details_dict[query]["is_contig"]#queryIScontig[query]                    # contig?

        # RECORD alignment:
        if contig_flag:                                      
            #print(dt.today(), query, "is part of a contig")
            for read in contig2read_map[query]:    
                if read in mapped_reads:                # (Track how many contig reads have already been mapped by DMD---i.e., through other contigs---
                    contigread_inmapped+= 1             
                    #print(dt.today(), query, "is also mapped")
                    if read in prot2read_map[db_match]: #  and how many of those were already assigned to this particular prot---i.e., contigs aligned to same prot.)
                        contigread_inmapped_inprot+= 1  
                        
                elif read in prot2read_map[db_match]:   # (Check to see if any of the contig reads are already assigned to this particular prot, but not by DMD.)
                    contigread_inprot+= 1              
                    
                if read not in mapped_reads:            #  if not already assigned by DMD to a diff prot***, then and append their readIDs to aligned prot<->read dict, and mark them as assigned by DMD.
                    #print(dt.today(), "contig:", query, ":", read, "read part of contig, but not mapped yet.  adding contig-read")
                    prot2read_map[db_match].append(read) 
                    mapped_reads.add(read)                
                
        
        else:
            #print(dt.today(), query, "is not a contig")
            # DEBUG:
            if query in mapped_reads:                   # (Check to see if the read has already been mapped
                read_inmapped+= 1                       #  by DMD---it shouldn't be---and
                if query in prot2read_map[db_match]:    #  and to the same prot, for that matter.
                    read_inmapped_inprot+= 1
            elif query in prot2read_map[db_match]:      # (Check to see if the read is already assigned to this
                read_inprot+= 1                         #  particular prot, but not by DMD---it shouldn't be.)
            
            if query in mapped_reads:
                print("=====================================")
                print(dt.today(), "read has been mapped already.  THIS SHOULDN'T BE HAPPENING")
                print(query)
            else:    
                prot2read_map[db_match].append(query)       #  append its readID to aligned prot<->read dict,
                mapped_reads.add(query)                     #  and mark it as assigned by DMD.
            #print(dt.today(), "added to mapped reads")
        
        #print("========================================================")
        # *** This deals with reads that show up in multiple contigs.
        # Just use the read alignment from the contig that had the best alignment score.
        # This could result in "broken" contigs...

    # DEBUG (for this datatype):
    print ('no. contig reads already mapped by DMD= ' + str(contigread_inmapped))
    print ('no. contig reads mapped by DMD to same prot= ' + str(contigread_inmapped_inprot))
    print ('no. contig reads mapped by NOT DMD to same prot= ' + str(contigread_inprot) + ' (should be 0)')
    print ('no. reads already mapped by DMD= ' + str(read_inmapped) + ' (should be 0)')
    print ('no. reads already mapped by DMD to same prot= ' + str(read_inmapped_inprot) + ' (should be 0)')
    print ('no. reads already mapped by NOT DMD to same prot= ' + str(read_inprot) + ' (should be 0)')

    # Remove contigs/reads previously added to the unmapped set but later found to have a mapping:
    # This prevents re-annotation by a later program.
    # Such queries end up in the unmapped set when they DMD-aligned to multiple prots, where one
    # alignment is recorded, while the other alignments fail the "on-the-fly" alignment-threshold filter.
    print ('umapped no. (before double-checking mapped set)= ' + str(len(unmapped)))
    #for query in query2prot_map:                        # Take all contigs/reads to be mapped and
    for query in query_keys:
        try:                                            #  if they exist in the unmapped set, then
            unmapped.remove(query)                      #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))

    # return unmapped set:
    return unmapped

# WRITE OUTPUT: rewrite gene<->read map file to include DMD-aligned:
# [BWA&BLAT&DMD-aligned geneID, length, #reads, readIDs ...]
def write_proteins_genemap(prot2read_map, Prot_DB, new_gene2read_file, prot_file):
    reads_count= 0
    proteins= []
    unique_reads_set = set()
    print(dt.today(), "starting to append to new gene map")
    
    #with open(new_gene2read_file,"a") as out_map:               
    with open(new_gene2read_file,"a") as out_map:               

        # write proteins:
        for record in SeqIO.parse(Prot_DB,"fasta"):         # Loop through SeqRec of all prot in PROTdb:
                                                            #  (PROTdb is needed to get the aa sequence.)
            if record.id in prot2read_map:                  #  If PROTdb prot is one of the matched proteins,
                proteins.append(record)                     #  append the SeqRec to proteins list (for next file), and
                out_map.write(record.id + "\t" + str(len(record.seq)*3) + "\t" + str(len(prot2read_map[record.id]))) #multiplied by 3 because proteins come in amino acids (groups of 3)
                                                            #  write [aligned protID, length (in nt), #reads, ...],
                for read in prot2read_map[record.id]:
                    out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                    reads_count+= 1
                else:
                    out_map.write("\n")                     #  and a new line character.
                    
    with open(prot_file,"w") as out_prot:
        SeqIO.write(proteins, out_prot, "fasta")            # and aligned proteins aa seqs.
    


def check_prot2read_map(prot2read_map):
    read_count = 0
    unique_reads = set()
    for key in prot2read_map:
        read_list = prot2read_map[key]
        read_count += len(read_list)
        unique_reads.update(read_list)
    if(read_count != len(unique_reads)):
        print(dt.today(), "there's a problem with this prot map")
        return False
    else:
        print(dt.today(), "this prot map is OK")
        return True
def filter_consumed_reads(read_file, DMD_tab_file, output_file, mapped_reads, prev_mapping_count, contig2read_map, prot2read_map):    

    
    # process DIAMOND output: # readtype sets: contigs, merged:
    print(dt.today(), "starting seqIO on reads")
    read_seqs = SeqIO.index(read_file, os.path.splitext(read_file)[1][1:]) #import reads (index.   key: read_id -> val: seq)
    if(check_file_safety(DMD_tab_file)):
        print(dt.today(), DMD_tab_file, "is ok")
    else:
        print(dt.today(), DMD_tab_file, "is not ok.  killing program")
        sys.exit()
    print(dt.today(), "finished seqIO on reads")

    # read DMD output & get read/contig lengths:
    DMD_hits = get_dmd_hit_details(DMD_tab_file, read_seqs)     # Store info in DMD_hits (list of lists).

    # process DMD-aligned reads:
    
    unmapped_reads = form_prot_map(DMD_hits, mapped_reads, contig2read_map, prot2read_map)    # Store DMD-aligned contigs/reads in prot2read_map
    
    unmapped_len_before= len(unmapped_reads)            # DEBUG
    for read in read_seqs:                              # Take all non-BWA&BLAT-aligned contigs/reads (input to DMD)
        if read not in mapped_reads:                    #  that are still unmapped and
            unmapped_reads.add(read)                    #  add them to the unmapped_reads set (won't add duplicates).

    print ('no. additional contigs/reads completely unmapped by DMD= ' + str(len(unmapped_reads)-unmapped_len_before))



    # print no. aligned reads from current readtype set:
    print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file) + '\n')
    prev_mapping_count= len(mapped_reads)
    
    
def construct_contig2read_map(contig2read_file):
    print(dt.today(), "reading contig map")
    # make dict of contigID<->readsID(s):
    contig2read_map= {}                                 #Input: key->contig | val->reads
    with open(contig2read_file,"r") as mapping:
        for line in mapping:
            if len(line)>5:                             # line starts with 'NODE_'
                entry= line.strip("\n").split("\t")     # break tab-separated into list
                contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
    print(dt.today(), "finished reading contig map")
    return contig2read_map

    
def import_gene_map(gene2read_file):
    # make dict of BWA&BLAT-aligned geneID<->readID(s):
    mapped_reads = set()
    gene2read_map = dict()
    old_gene_map_reads_count = 0
    with open(gene2read_file,"r") as mapping:
        for line in mapping:
            if len(line)>5:                             # line at least 5 characeters?
                entry = line.strip("\n").split("\t")
                gene2read_map[entry[0]] = entry[3:]      # key=geneID, value=list of readID(s) (using [:] syntax ensures a list, even if one ele)
                mapped_reads.update(entry[3:])          # DEBUG
                old_gene_map_reads_count += len(entry[3:])
    return gene2read_map, mapped_reads, old_gene_map_reads_count
    
    
def write_unmapped_seqs(unmapped_reads, reads_in, reads_out):
    # WRITE unmerged1 OUTPUT: non-BWA&BLAT&DMD-aligned:
    if(len(unmapped_reads) == 0):
        print(dt.today(), "no unmapped reads found.  skipping")
    else:
        read_seqs = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:])
        unmapped_seqs= []
        for read in unmapped_reads:
            read_key = read
            #if(not read.startswith("@")):
            #    read_key = "@" + read
            if(read_key in read_seqs):
                unmapped_seqs.append(read_seqs[read_key])
            #else:
            #    print(read_key,"from",DMD_tab_file_1, DMD_tab_file_2,   "not found in:", read_file_1)
        with open(reads_out,"w") as outfile:
            SeqIO.write(unmapped_seqs, outfile, "fasta")
            
        
#####################################
if __name__ == "__main__":

    Prot_DB                 = sys.argv[1]   # INPUT: AA db used for DIAMOND alignement
    contig2read_file        = sys.argv[2]   # INPUT: [contigID, #reads, readIDs ...]
    new_gene2read_file      = sys.argv[3]   # OUTPUT: [BWA&BLAT&DMD-aligned gene/protID, length, #reads, readIDs ...]
    prot_file               = sys.argv[4]   # OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs (.faa; fasta-format)
    
    
    reads_in    = sys.argv[5]
    dmd_in      = sys.argv[6]
    reads_out   = sys.argv[7]
    
    input_safety = check_file_safety(reads_in) and check_file_safety(dmd_in)
    
    if(input_safety):
        #"global" vars
        contig2read_map = construct_contig2read_map(contig2read_file)   #Input: key->contig | val->reads
        gene2read_map = {}                                   #Input: key->geneID | val->reads
        mapped_reads = set()                                 # tracks all reads mapped.  
        prot2read_map = defaultdict(list)                    # dict of DMD-aligned protID<->readID(s) #  key=protID, value=list of readID(s)
        prev_mapping_count = 0
        
        
        #gene2read_map, mapped_reads, prev_mapping_count = import_gene_map(gene2read_file)
        process_store = []
        
        dmd_hits = get_dmd_hit_details(dmd_in, reads_in)
        unmapped_reads = form_prot_map(dmd_hits, mapped_reads, contig2read_map, prot2read_map)
        write_prot_map_process = mp.Process(target = write_proteins_genemap, args = (prot2read_map, Prot_DB, new_gene2read_file, prot_file))
        write_prot_map_process.start()
        process_store.append(write_prot_map_process)
        
        
        write_unmapped_reads_process = mp.Process(target = write_unmapped_seqs, args = (unmapped_reads, reads_in, reads_out))
        write_unmapped_reads_process.start()
        process_store.append(write_unmapped_reads_process)
        
        for item in process_store:
            item.join()
        process_store[:] = []
        
        
    else:
        if(os.path.exists(reads_in)):
            print("reads weren't annotated.  leftovers")
            copyfile(reads_in, reads_out)
        else:
            print("no reads remaining")
            open(reads_out, "a").close()
        
       
        
        