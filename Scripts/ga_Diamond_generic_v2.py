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
# - stores aligned proteins in new dict prot2read_map, instead of prot2read_map:
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
    
    hits = sorted(hits, key=sortbyscore, reverse = True)     # sort alignment list by high score:
    
    
    return hits

# sort function: sort by score:
#sort by bitscore
def sortbyscore(line):
    return line[11]



    
# add DMD-aligned reads that meet threshold
# to the aligned gene/protID<->readID(s) dict:
def form_prot_map(identity_cutoff, length_cutoff, score_cutoff, hits, contig2read_map): 
    # sort alignment list by high score:
    prot2read_map = dict()
    mapped = set()
    unmapped = set()                         # Set of unmapped contig/readIDs.
    query_details_dict = dict()

    # BLAT threshold:
    #identity_cutoff= 85
    #length_cutoff= 0.65
    #score_cutoff= 60

    repeat_reads = 0
    disagreements = 0
    one_sided_alignments = 0
    
    # loop through sorted BLAT-aligned reads/contigs:
    for line in hits:
        inner_details_dict = dict()
        # "ON-THE-FLY" filter:
        
        # extract & store data:
        query = line[0]                      # queryID= contig/readID
        db_match = line[1]                   # geneID
        seq_identity = float(line[2])        # sequence identity
        align_len = 3*int(line[3])             # alignment length
        score = float(line[11])              # score
        seq_len = int(line[12])              # query sequence length
        bitscore = float(line[11])
        
        # does query alignment meet threshold?:
        # (identity, length, score):
        if seq_identity<=identity_cutoff or align_len<=seq_len*length_cutoff or score<=score_cutoff:
            #unmapped.add(query)             # if threshold is failed, add to unmapped set and
            inner_details_dict = dict()
            inner_details_dict["match"] = False
            query_details_dict[query] = inner_details_dict
        else:
            
            # is query a contig?:
            if query in contig2read_map:        # If query is found in contig list,
                contig= True                    #  then mark as contig,
            else:                               #  if not,
                contig= False                   #  mark as not a contig.
                
                
            # is this alignment the highest-score match for that query?
            if query in query_details_dict: #query2gene_map:         # If alignment previously found for this contig/read,
                repeat_reads += 1
                inner_details_dict = query_details_dict[query]
                if(inner_details_dict["match"]):
                    old_bitscore = inner_details_dict["bitscore"]
                    disagreements += 1
                    if(old_bitscore < bitscore):
                        print("old:", old_bitscore, "new:", bitscore)
                        #update, but this shouldn't happen
                        inner_details_dict["gene"] = db_match
                        inner_details_dict["bitscore"] = bitscore
                        inner_details_dict["is_contig"] = contig
                else:
                    one_sided_alignments += 1
                    inner_details_dict["match"] = True
                    inner_details_dict["is_contig"] = contig
                    inner_details_dict["gene"] = db_match
                    inner_details_dict["bitscore"] = bitscore
                
                                            
                              # skip to the next query.
            else:
                #new hit
                inner_details_dict = dict()
                inner_details_dict["match"] = True
                inner_details_dict["is_contig"] = contig
                inner_details_dict["gene"] = db_match
                inner_details_dict["bitscore"] = bitscore
                query_details_dict[query] = inner_details_dict
            #queryIScontig[query] = contig        # Store contig (T/F) info.
            #query2gene_map[query].add(db_match) # Collect all aligned genes for contig/read;
                                                #  query2gene_map[query] is a set, although there should be
                                                #  no more than one gene alignement getting through filter.
            
    # EXPAND HERE to deal with multiple high-score genes for a read.
    print(dt.today(), "repeat reads:", repeat_reads)
    print(dt.today(), "disagreements:", disagreements)
    print(dt.today(), "one-sided alignments:", one_sided_alignments)
    
    # FINAL remaining BLAT-aligned queries:
    for query in query_details_dict:#query2gene_map:                        # contig/readID
        inner_dict = query_details_dict[query]
        if(inner_dict["match"]):
            mapped.add(query)
            db_match = inner_dict["gene"] 
            contig = inner_dict["is_contig"]
            bitscore = inner_dict["bitscore"]
            
            # RECORD alignment:
            if contig:                                      # If query is a contig, then
                for read in contig2read_map[query]:         #  take all reads making up that contig and
                    read_line = read + "<bitscore>" + str(bitscore)
                    if(db_match in prot2read_map):
                        prot2read_map[db_match].append(read_line)#  append their readIDs to aligned gene<->read dict,
                    else:
                        prot2read_map[db_match] = [read_line]
                        
            else:
                read_line = query + "<bitscore>" + str(bitscore)
                if(db_match in prot2read_map):
                    prot2read_map[db_match].append(read_line)
                else:
                    prot2read_map[db_match] = [read_line]

        
        else:
            unmapped.add(query)
    return unmapped, mapped, prot2read_map

# WRITE OUTPUT: rewrite gene<->read map file to include DMD-aligned:
# [BWA&BLAT&DMD-aligned geneID, length, #reads, readIDs ...]
def write_proteins_genemap(prot2read_map, Prot_DB, new_prot2read_file, prot_file):
    reads_count= 0
    proteins= []
    unique_reads_set = set()
    print(dt.today(), "starting to append to new gene map")
    
    #with open(new_gene2read_file,"a") as out_map:               
    with open(new_prot2read_file,"a") as out_map:               

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
            
        with open(reads_out,"w") as outfile:
            SeqIO.write(unmapped_seqs, outfile, "fasta")
            
        
#####################################
if __name__ == "__main__":
    identity_cut            = sys.argv[1]
    length_cut              = sys.argv[2]
    score_cut               = sys.argv[3]

    Prot_DB                 = sys.argv[4]   # INPUT: AA db used for DIAMOND alignement
    contig2read_file        = sys.argv[5]   # INPUT: [contigID, #reads, readIDs ...]
    new_gene2read_file      = sys.argv[6]   # OUTPUT: [BWA&BLAT&DMD-aligned gene/protID, length, #reads, readIDs ...]
    prot_file               = sys.argv[7]   # OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs (.faa; fasta-format)
    
    
    reads_in    = sys.argv[8]
    dmd_in      = sys.argv[9]
    reads_out   = sys.argv[10]
    
    identity_cutoff = int(identity_cut)
    length_cutoff = float(length_cut)
    score_cutoff = int(score_cut)
    
    input_safety = check_file_safety(reads_in) and check_file_safety(dmd_in)
    
    if(input_safety):
        #"global" vars
        contig2read_map = construct_contig2read_map(contig2read_file)   #Input: key->contig | val->reads
        
        
        dmd_hits = get_dmd_hit_details(dmd_in, reads_in)
        unmapped_reads, mapped_reads, prot2read_map = form_prot_map(identity_cutoff, length_cutoff, score_cutoff, dmd_hits, contig2read_map)
        write_proteins_genemap(prot2read_map, Prot_DB, new_gene2read_file, prot_file)
        
        
        write_unmapped_seqs(unmapped_reads, reads_in, reads_out)
        
        
    else:
        if(os.path.exists(reads_in)):
            print("reads weren't annotated.  leftovers")
            copyfile(reads_in, reads_out)
        else:
            print("no reads remaining")
            open(reads_out, "a").close()
        
       
        
        