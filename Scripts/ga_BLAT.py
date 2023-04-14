#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - added .strip("\n") to extracting from the gene2read_file
# - fixed sort function to return float
# - reversed the sort direction to get largest scores at top
# - pulled gene_map() out of the readtype loop
# - rewrote gene_map()
# - removed checks to see if contig/read was already assigned to same gene (shouldn't happen)
# - removed check to see if query=read was already mapped by BLAT (shouldn't happen)
# - changed align_len to float(align_len) (to make sure it's a number)
# - paired-end BLAT outputs are combined and analyzed at the same time

# NOTE:
# - Filenames for unmerged paired-end reads must be specified last.
# - Any read aligned by BWA will not be input into BLAT. Even BWA-aligned contigs
#   have unique reads (not found in other contigs).
# - BLATpp only takes the read<->gene match with the top score. For duplicate
#   top scores, it only takes the one that shows up first in the initial sort.
# - The paired-end BLAT outputs are combined and analyzed at the same time, in order to
#   rank the matches between different ends of the same read, and not double count that read.
# - Sometimes a read is aligned via multiple contigs and annoteted only according
#   to the best matched contig only. Theoretically, this breaks up the other contigs.

import os
import os.path
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO

DNA_DB= sys.argv[1]                 # INPUT: DNA db used for BLAT alignement
contig2read_file= sys.argv[2]       # INPUT: [contigID, #reads, readIDs ...]
gene2read_file= sys.argv[3]         # INPUT: [BWA-aligned geneID, length, #reads, readIDs ...]
                                    # ->OUTPUT: [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
gene_file= sys.argv[4]              # OUTPUT: BWA&BLAT-aligned geneIDs and aa seqs (.fna; fasta-format)
new_gene2read_file = sys.argv[5]    # OUTPUT: new gene2read_file instead of overwriting the old map
new_gene_file = sys.argv[6]         # OUTPUT new gene.fna 
# make dict of contigID<->readsID(s):
contig2read_map= {}
with open(contig2read_file,"r") as mapping:
    for line in mapping:
        if len(line)>5:                             # line starts with 'NODE_'
            entry= line.strip("\n").split("\t")     # break tab-separated into list
            contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)

# make dict of BWA-aligned geneID<->readID(s):
BWAreads= []                                        # DEBUG
gene2read_map= defaultdict(list)                    # Dict of BWA&BLAT-aligned geneID<->readID(s)
with open(gene2read_file,"r") as mapping:           #  initialized w BWA-alignments.
    for line in mapping:
        if len(line)>5:                             # line at least 5 characeters?
            entry= line.strip("\n").split("\t")
            gene2read_map[entry[0]]= entry[3:]      # key=geneID, value=list of readID(s)
                                                    # (using [:] syntax ensures a list, even if one ele)
            BWAreads.extend(entry[3:])              # DEBUG

# tracking BLAT-assigned:
mapped_reads= set()                                 # tracks BLAT-assigned reads
prev_mapping_count= 0

# DEBUG:
if len(set(BWAreads))==len(BWAreads):
    print ('BWA-aligned reads are all unique.\n')
else:
    print ('There are repeating BWA-aligned reads:')
    print ('no. unique reads= ' + str(len(set(BWAreads))))
    print ('no. total reads= ' + str(len(BWAreads)) + '\n')
BWAreads_count= Counter(BWAreads)                   # dict of read<->no. of contigs


#####################################
# FUNCTIONS:

# read .blatout file and acquire read length:
def read_aligned(tsv, seqrec):                          # List of lists [blatout info, read length]=
                                                        #  read_aligned(.blatout file, dict contig/readID<->SeqRecord)
    # get info from .blatout file:
    print ('Reading ' + str(os.path.basename(tsv)) + '.')
    with open(tsv,"r") as tabfile:
        hits= []                                        # List of lists containing .blatout fields.
        for line in tabfile:                            # In the .blatout file:
            if len(line)>=2:                            # If length of line >= 2,
                info_list= line.strip("\n").split("\t") #  make a list from .blatout tab-delimited fields,
                query= info_list[0]                     #  get the queryID= contig/readID,
                info_list.append(len(seqrec[query].seq))#  add query length (int) to end of list,
                hits.append(info_list)                  #  and append the list to the hits list.

    # return info:
    return hits

# sort by score:
# (12th field of the .blatout file)
def sortbyscore(line):
    return float(line[11])

# add BLAT-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:
def gene_map(hits):                         # fail-mapped contig/readIDs=
                                            #  gene_map(list of list of blatout fields)
    # sort alignment list by high score:
    sorted_hits= sorted(hits, key=sortbyscore, reverse=True)
    del hits

    # BLAT threshold:
    identity_cutoff= 85
    length_cutoff= 0.65
    score_cutoff= 60

    # tracking BLAT-assigned & unassigned:
    query2gene_map= defaultdict(set)        # Dict of BLAT-aligned contig/readID<->geneID(s).
    queryIScontig= {}                       # Dict of BLAT-aligned contig/readID<->contig? (True/False).
    unmapped= set()                         # Set of unmapped contig/readIDs.

    # loop through sorted BLAT-aligned reads/contigs:
    for line in sorted_hits:
    
        # "ON-THE-FLY" filter:
        
        # extract & store data:
        query= line[0]                      # queryID= contig/readID
        db_match= line[1]                   # geneID
        seq_identity= float(line[2])        # sequence identity
        align_len= int(line[3])             # alignment length
        score= float(line[11])              # score
        seq_len= int(line[12])              # query sequence length
        
        # is this alignment the highest-score match for that query?
        if query in query2gene_map:         # If alignment previously found for this contig/read,
            continue                        # skip to next qurey as this alignment will have a lower score,
                                            # and don't add to unmapped set.
        # is query a contig?:
        if query in contig2read_map:        # If query is found in contig list,
            contig= True                    #  then mark as contig,
        else:                               #  if not,
            contig= False                   #  mark as not a contig.
                                        
        # does query alignment meet threshold?:
        # (identity, length, score):
        if seq_identity<=identity_cutoff or align_len<=seq_len*length_cutoff or score<=score_cutoff:
            unmapped.add(query)             # if threshold is failed, add to unmapped set and
            continue                        # skip to the next query.
        
        # store info for queries that remain:
        queryIScontig[query]= contig        # Store contig (T/F) info.
        query2gene_map[query].add(db_match) # Collect all aligned genes for contig/read;
                                            #  query2gene_map[query] is a set, although there should be
                                            #  no more than one gene alignement getting through filter.
        
    # EXPAND HERE to deal with multiple high-score genes for a read.
    
    # DEBUG:
    contigread_inmapped= 0
    contigread_inmapped_ingene= 0
    contigread_ingene= 0
    read_inmapped= 0
    read_inmapped_ingene= 0
    read_ingene= 0
    
    # FINAL remaining BLAT-aligned queries:
    for query in query2gene_map:                        # contig/readID
    
        db_match= list(query2gene_map[query])[0]        # geneID (pull out of 1-element set)
        contig= queryIScontig[query]                    # contig?
        
        # RECORD alignment:
        if contig:                                      # If query is a contig, then
            for read in contig2read_map[query]:         #  take all reads making up that contig and
            
                # DEBUG:
                if read in mapped_reads:                # (Track how many contig reads have already been
                    contigread_inmapped+= 1             #  mapped by BLAT---i.e., through other contigs---
                    if read in gene2read_map[db_match]: #  and how many of those were already assigned to this
                        contigread_inmapped_ingene+= 1  #  particular gene---i.e., contigs aligned to same gene.)
                elif read in gene2read_map[db_match]:   # (Check to see if any of the contig reads are already
                    contigread_ingene+= 1               #  assigned to this particular gene, but not by BLAT.)
                
                if read not in mapped_reads:            #  if not already assigned by BLAT to a diff gene***, then
                    gene2read_map[db_match].append(read)#  append their readIDs to aligned gene<->read dict,
                    mapped_reads.add(read)              #  and mark them as assigned by BLAT.
        elif not contig:                                # If query is a read, then
        
            # DEBUG:
            if query in mapped_reads:                   # (Check to see if the read has already been mapped
                read_inmapped+= 1                       #  by BLAT---it shouldn't be---and
                if query in gene2read_map[db_match]:    #  and to the same gene, for that matter.
                    read_inmapped_ingene+= 1
            elif query in gene2read_map[db_match]:      # (Check to see if the read is already assigned to this
                read_ingene+= 1                         #  particular gene, but not by BLAT---it shouldn't be.)
            
            gene2read_map[db_match].append(query)       #  append its readID to aligned gene<->read dict,
            mapped_reads.add(query)                     #  and mark it as assigned by BLAT.

        # *** This deals with reads that show up in multiple contigs.
        # Just use the read alignment from the contig that had the best alignment score.
        # This could result in "broken" contigs...

    # DEBUG (for this datatype):
    print ('no. contig reads already mapped by BLAT= ' + str(contigread_inmapped))
    print ('no. contig reads mapped by BLAT to same gene= ' + str(contigread_inmapped_ingene))
    print ('no. contig reads mapped by NOT BLAT to same gene= ' + str(contigread_ingene) + ' (should be 0)')
    print ('no. reads already mapped by BLAT= ' + str(read_inmapped) + ' (should be 0)')
    print ('no. reads already mapped by BLAT to same gene= ' + str(read_inmapped_ingene) + ' (should be 0)')
    print ('no. reads already mapped by NOT BLAT to same gene= ' + str(read_ingene) + ' (should be 0)')

    # Remove contigs/reads previously added to the unmapped set but later found to have a mapping:
    # This prevents re-annotation by a later program.
    # Such queries end up in the unmapped set when they BLAT-aligned to multiple genes, where one
    # alignment is recorded, while the other alignments fail the "on-the-fly" alignment-threshold filter.
    print ('umapped no. (before double-checking mapped set)= ' + str(len(unmapped)))
    for query in query2gene_map:                        # Take all contigs/reads to be mapped and
        try:                                            #  if they exist in the unmapped set, then
            unmapped.remove(query)                      #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))

    # return unmapped set:
    return unmapped

#####################################

# check number of readtype sets (file inputs)
read_sets = int((len(sys.argv)-7)/3)
if (len(sys.argv)-7) % 3 != 0:
    sys.exit('Incorrect number of readtype sets.')

# process BLAT output:
# readtype sets: contigs, merged:
for x in range(read_sets):
    read_file= sys.argv[3*x+7]      # INPUT: non-BWA-aligned contig/readIDs and seqs (.fasta)
    read_seqs= SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
                                    # dict of non-BWA-aligned read SeqRecords: key=contig/readID
                                    #  (second argument specifies filetype, e.g., "fasta")
    BLAT_tab_file= sys.argv[3*x+8]  # INPUT: BLAT-aligned contig/readIDs (.blatout)
    output_file= sys.argv[3*x+9]    # OUTPUT: non-BWA&BLAT-aligned contig/readIDs and seqs (.fasta)

    # read BLAT output & get read/contig lengths:
    BLAT_hits= read_aligned(BLAT_tab_file, read_seqs)   # Store info in BLAT_hits (list of lists).

    # process BLAT-aligned reads:
    unmapped_reads= gene_map(BLAT_hits)                 # Add BLAT-aligned contigs/reads to gene2read_map
                                                        #  (aligned geneID<->readID(s) dict),
                                                        #  and return a set of failed mapping readIDs.
    # add reads never mapped by BLAT in the
    # first place to the unmapped set:
    unmapped_len_before= len(unmapped_reads)            # DEBUG
    for read in read_seqs:                              # Take all non-BWA-aligned contigs/reads (input to BLAT)
        if read not in mapped_reads:                    #  that are still unmapped and
            unmapped_reads.add(read)                    #  add them to the unmapped_reads set (won't add duplicates).

    print ('no. additional contigs/reads completely unmapped by BLAT= ' + str(len(unmapped_reads)-unmapped_len_before))

    # WRITE OUTPUT: non-BWA&BLAT-aligned contig/readIDs:
    # and seqs (.fasta):
    unmapped_seqs= []                                   # Initialize list of SeqRecords.
    for read in unmapped_reads:                         # Put corresponding SeqRecords for unmapped_reads
        unmapped_seqs.append(read_seqs[read])           #  into unmapped_seqs
    with open(output_file,"w") as outfile:
        SeqIO.write(unmapped_seqs, outfile, "fasta")    #  and write it to file.

    # print no. aligned reads from current readtype set:
    print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file) + '\n')
    prev_mapping_count= len(mapped_reads)
'''
# process BLAT output:
# readtype sets: unmerged1, unmerged2
if numsets==4:
    
    # unmerged1:
    x= 2
    read_file_1= sys.argv[3*x+7]
    read_seqs_1= SeqIO.index(read_file_1, os.path.splitext(read_file_1)[1][1:])
    BLAT_tab_file_1= sys.argv[3*x+8]
    output_file_1= sys.argv[3*x+9]
    BLAT_hits_1= read_aligned(BLAT_tab_file_1, read_seqs_1) # Store info in BLAT_hits_1 (list of lists).

    # unmerged2:
    x= 3
    read_file_2= sys.argv[3*x+7]
    read_seqs_2= SeqIO.index(read_file_2, os.path.splitext(read_file_2)[1][1:])
    BLAT_tab_file_2= sys.argv[3*x+8]
    output_file_2= sys.argv[3*x+9]
    BLAT_hits_2= read_aligned(BLAT_tab_file_2, read_seqs_2) # Store info in BLAT_hits_2 (list of lists).

    # process BLAT-aligned reads together:
    BLAT_hits= BLAT_hits_1+BLAT_hits_2                  # Concatenate BLAT info lists,
    unmapped_reads= gene_map(BLAT_hits)                 #  add BLAT-aligned contigs/reads to gene2read_map
                                                        #  and return a set of failed mapping readIDs.
    # add reads never mapped by BLAT in the
    # first place to the unmapped set:
    unmapped_len_before= len(unmapped_reads)
    for read in read_seqs_1:
        if read not in mapped_reads:
            unmapped_reads.add(read)
    print ('(1st paired end) no. additional contigs/reads completely unmapped by BLAT= ' + str(len(unmapped_reads)-unmapped_len_before))

    # add reads never mapped by BLAT in the
    # first place to the unmapped set:
    unmapped_len_before= len(unmapped_reads)
    for read in read_seqs_2:
        if read not in mapped_reads:
            unmapped_reads.add(read)
    print ('(2nd paired end) no. additional contigs/reads completely unmapped by BLAT= ' + str(len(unmapped_reads)-unmapped_len_before) + ' (should be 0)')

    # WRITE unmerged1 OUTPUT: non-BWA&BLAT-aligned:
    unmapped_seqs= []
    for read in unmapped_reads:
        unmapped_seqs.append(read_seqs_1[read])
    with open(output_file_1,"w") as outfile:
        SeqIO.write(unmapped_seqs, outfile, "fasta")

    # WRITE unmerged2 OUTPUT: non-BWA&BLAT-aligned:
    unmapped_seqs= []
    for read in unmapped_reads:
        unmapped_seqs.append(read_seqs_2[read])
    with open(output_file_2,"w") as outfile:
        SeqIO.write(unmapped_seqs, outfile, "fasta")

    # print no. aligned reads from current readtype set:
    print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file_1))
    print ('  and ' + os.path.basename(read_file_2) + '\n')
    prev_mapping_count= len(mapped_reads)
'''
# WRITE OUTPUT: rewrite gene<->read mapfile to include BLAT-aligned:
# [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
reads_count= 0
genes= []
with open(new_gene2read_file,"w") as out_map:               # Delete old gene2read_file and write a new one.
    for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                        #  (DNA db is needed to get the sequence.)
        if record.id in gene2read_map:                  #  If DNA db gene is one of the matched genes,
            genes.append(record)                        #  append the SeqRec to genes list (for next file), and
            out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                        #  write [aligned geneID, length, #reads, ...],
            for read in gene2read_map[record.id]:
                out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                reads_count+= 1
            else:
                out_map.write("\n")                     #  and a new line character.

# WRITE OUTPUT: BWA&BLAT-aligned geneIDs and seqs (.fna; fasta-format):
# (this wasn't done in BWA post-processing)
with open(gene_file,"w") as outfile:
    SeqIO.write(genes, outfile, "fasta")

# print BWA+BLAT stats:
print (str(reads_count) + ' reads were mapped with BWA and BLAT.')
print ('Reads mapped to ' + str(len(genes)) + ' genes.')