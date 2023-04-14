#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - fixed 'while line.startswith("@"): continue' infinte loop
# - changed gene_map() to deal with thresholds and
#   multiple alignments before the "Recording" part.

# NOTES:
# - Some contigs/reads may match to the DB multiple times.
# - Filenames for unmerged paired-end reads must be specified last.

# "non-BWA-aligned" = unmapped should consist of:
# (1) contigs/reads not aligned by BWA at all
# (2) BWA-aligned contigs containing reads involved in other contigs
# (3) BWA-aligned contigs/reads that don't meet match % threshold
# (4) BWA-aligned reads that are part of a contig (THIS SHOULDT'T HAPPEN), but not the contig
# (5) BWA-aligned contigs/reads that align to multiple genes (PAIRED READS CONSIDERED TOGETHER)

import os
import os.path
import re
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO

DNA_DB= sys.argv[1]             # INPUT: DNA db used for BWA alignement
contig2read_file= sys.argv[2]   # INPUT: [contigID, #reads, readIDs ...]
gene2read_file= sys.argv[3]     # OUTPUT: [BWA-aligned geneID, length, #reads, readIDs ...]

# make initial dict of contigID<->readsID(s):
contig2read_map= {}
contig_reads= []                                    # list of just reads
with open(contig2read_file,"r") as mapping:
    for line in mapping:
        if len(line)>5:                             # line starts with 'NODE_'
            entry= line.strip("\n").split("\t")     # break tab-separated into list
            contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
            contig_reads.extend(entry[2:])          # append all the reads

# make new dict only of contigs with unique reads:
# (hard to tell w BWA what contigs match better, so for reads associated
# w multiple matched contigs, avoid choosing btw contigs for now.)
contig_reads_count= Counter(contig_reads)           # dict of read<->no. of contigs
contig2read_map_uniq= {}
contig_unique_reads= []                             # DEBUG
for contig in contig2read_map:
    for read in contig2read_map[contig]:            # If contig has
        if contig_reads_count[read]>1:              #  a read assoc. w multiple contigs
            break                                   #  then throw the contig away,
    else:
        contig2read_map_uniq[contig]= contig2read_map[contig]
                                                    #  else, store it in the unique dict.
        contig_unique_reads.extend(contig2read_map[contig]) # DEBUG

# tracking BWA-assigned:
gene2read_map= defaultdict(list)                    # dict of BWA-aligned geneID<->readID(s)
mapped_reads= set()                                 # tracks BWA-assigned reads
mapped_list= []
prev_mapping_count= 0

#####################################
# FUNCTION:
# add BWA-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:

# additional filtering steps:
# (1) only use contigs that contain unique reads (not shared w other contigs)
# (2) make sure BWA-aligned (.sam file has aligned & non-aligned data)
# (3) only use contigs/reads where >=90% length of the seq was matched

# CIGAR string describes how read aligns with the ref. Consists of >=1 components.
# Each component comprises an operator and no. bases which the op applies to.
#
# Operators:
# D	Deletion; the nucleotide is present in the reference but not in the read.
# H	Hard Clipping; the clipped nucleotides are not present in the read.
# I	Insertion; the nucleotide is present in the read  but not in the reference.
# M	Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
# N	Skipped region; a region of nucleotides is not present in the read.
# P	Padding; padded area in the read and not in the reference.
# S	Soft Clipping;  the clipped nucleotides are present in the read.
# =	Read Match; the nucleotide is present in the reference.
# X	Read Mismatch; the nucleotide is present in the reference.

def gene_map(sam):                                      # Set of unmapped contig/readIDs=
                                                        #  gene_map(BWA .sam file)
    # tracking BWA-assigned & unassigned:
    query2gene_map= defaultdict(set)                    # Dict of BWA-aligned contig/readID<->geneID(s).
    queryIScontig= {}                                   # Dict of BWA-aligned contig/readID<->contig? (True/False).
    unmapped= set()                                     # Set of unmapped contig/readIDs.

    len_chars= ["M","I","S","=","X"]                    # These particular CIGAR operations cause the
                                                        #  alignment to step along the query sequence.
                                                        # Sum of lengths of these ops=length of seq.

    # first, process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
        
            # "ON-THE-FLY" filter:
        
            # extract & store data:
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  go to the next query (restart for).
            line_parts= line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and store:
            query= line_parts[0]                        #  queryID= contig/readID,
            db_match= line_parts[2]                     #  geneID, and a
            flag= bin(int(line_parts[1]))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                                                        #  where each bit is a flag for a specific descriptor.

            # is query BWA-aligned?:
            if flag[8]=="1":                            # If contig/read is NOT BWA-ALIGNED (9th digit=1),
                unmapped.add(query)                     # add it to the unmapped set and
                continue                                # skip to the next query.
            
            # is query a contig made of unique reads?:
            if query in contig2read_map:                # If query is a contig (searches through keys)
                if query in contig2read_map_uniq:       #  and it's made of contig-unique reads,
                    contig= True                        #  then mark as contig and move on.
                else:                                   # Otherwise, contig contains non-unique reads,
                    unmapped.add(query)                 #  therefore add it to the unmapped set and
                    continue                            #  skip to the next query.
            else:
                contig= False                           # If query isn't a contig, just move on.
            
            # does query alignment meet threshold?:
            length= 0
            matched= 0
            CIGAR= re.split("([MIDNSHPX=])", line_parts[5]) # Split CIGAR string into list, placing
                                                            #  all chars w/in [...] into own field
                                                            #  (e.g., 9S41M50S->['9','S','41','M','50','S','']).
            for index in range(len(CIGAR))[:-1]:            # Loop CIGAR elements (last element=''),
                if CIGAR[index+1] in len_chars:             # Use CIGAR operations that step along the query seq,
                    length+= int(CIGAR[index])              #  to determine length of query.
                if CIGAR[index+1]=="M":                     # Use CIGAR match operation to
                    matched+= int(CIGAR[index])             #  determine no. nuclotides matched.
            if matched<length*0.9:                          # If alignment is <90% matched:
                unmapped.add(query)                         # add it to the unmapped set and
                continue                                    # skip to the next query.

            # store info for queries that remain:
            query2gene_map[query].add(db_match)             # Collect all aligned genes for contig/read.
            queryIScontig[query]= contig                    # Store contig (T/F) info.

    print ('Reading ' + str(os.path.basename(sam)) + '.')
    print ('no. queries (that meet initial threshold)= ' + str(len(query2gene_map)))
    
    # remove read queries that are part of contigs:
    # THOUGH THIS SHOULDN'T HAPPEN (and it doesn't).
    for query in query2gene_map:                        # Check all queries
        if queryIScontig[query]:                        #  (that are reads) to see
            if query in contig_reads:                   #  if they belong to any contigs, and if so
                del query2gene_map[query]               #  delete the read from the query list.
                del queryIScontig[query]                # Don't add it to the unmapped set, as it shouldn't be
                                                        #  annotated as a stand-alone read.

    print ('no. (after removal of read queries in contigs)= ' + str(len(query2gene_map)))

    # remove queries aligning to multiple genes:
    for query, genes in list(query2gene_map.items()):   # Check all queries to see
        if len(genes)>1:                                #  if they are aligend to multiple genes, and if so
            unmapped.add(query)                         #  add it to the unmapped set, then
            del query2gene_map[query]                   #  delete it from the query list.
            del queryIScontig[query]

    print ('no. (after removal of queries aligning to multiple genes)= ' + str(len(query2gene_map)))

    # FINAL remaining BWA-aligned queries:
    for query in query2gene_map:                        # contig/readID
    
        db_match= list(query2gene_map[query])[0]        # geneID (pull out of 1-element set)
        contig= queryIScontig[query]                    # contig?
    
        # RECORD alignments:
        if contig:                                      # If query is a contig, then
            gene2read_map[db_match].extend(contig2read_map_uniq[query])
                                                        #  add the readIDs of all the reads making up
                                                        #  that contig to the aligned gene<->read dict,
            for read in contig2read_map_uniq[query]:    #  and mark all those reads
                mapped_reads.add(read)                  #  as assigned by BWA.
                mapped_list.append(read)                # DEBUG (not a set so will double w paired reads)
        else:                                           # If query is a read,
            gene2read_map[db_match].append(query)       #  append its readID to aligned gene<->read dict,
            mapped_reads.add(query)                     #  and mark it as assigned by BWA.
            mapped_list.append(query)                   # DEBUG (not a set so will double w paired reads)

    # Remove contigs/reads previously added to the unmapped set but later found to have a mapping:
    # This prevents re-annotation by a later program.
    # Such queries end up in the unmapped set when they BWA-aligned to multiple genes, where one
    # alignment is recorded, while the other alignments fail the "on-the-fly" filter; or
    # when one end of an unmerged paired read aligns and the other end doesnt.
    print ('umapped no. (before double-checking mapped set)= ' + str(len(unmapped)))
    for query in query2gene_map:                        # Take all contigs/reads to be mapped and
        try:                                            #  if they exist in the unmapped set, then
            unmapped.remove(query)                      #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))


    # DEBUG (check so far):
    if len(set(mapped_list))==len(mapped_list):
        print ('So far, BWA-aligned reads are all unique.')
    else:
        print ('Repeating BWA-aligned reads thus far:')
        print ('no. unique reads= ' + str(len(set(mapped_list))))
        print ('no. total reads= ' + str(len(mapped_list)))
        print ('no reads in the set= ' + str(len(mapped_reads)))

    # return unmapped set:
    return unmapped

#####################################

# check number of readtype sets (file inputs)
numsets= (len(sys.argv)-4)/3

if numsets not in [2,4]:
    sys.exit('Incorrect number of readtype sets.')

# process BWA output:
# readtype sets: contigs, merged, unmerged1, unmerged2
for x in range(int((len(sys.argv)-4)/3)):
    read_file= sys.argv[3*x+4]      # INPUT: all contig/readIDs and seqs (.fasta)
    read_seqs= SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
                                    # dict of all read SeqRecords: key=contig/readID
                                    #  (second argument specifies filetype, e.g., "fasta")
    BWA_sam_file= sys.argv[3*x+5]   # INPUT: BWA-aligned&unaligned contig/readIDs (.sam)
    output_file= sys.argv[3*x+6]    # OUTPUT: non-BWA-aligned contig/readIDs and seqs (.fasta)

    # extraction of "non-BWA-aligned" and "BWA-aligned":
    if x!=3:                                        # Only do once for unmerged paired end (same .sam file).
        unmapped_reads= gene_map(BWA_sam_file)      # Store BWA-aligned contigs/reads in gene2read_map
                                                    #  (aligned geneID<->readID(s) dict),
                                                    #  and return a set of unmapped readIDs
    # WRITE OUTPUT: non-BWA-aligned contig/readIDs:
    # and seqs (.fasta):
    unmapped_seqs= []                               # Inintialize list of SeqRecords.
    for read in unmapped_reads:                     # Put corresponding SeqRecords for unmapped_reads
        if(read in read_seqs):
            unmapped_seqs.append(read_seqs[read])       #  into unmapped_seqs
        else:
            print("ignoring:", read, "can't find in read_seqs")
    with open(output_file,"w") as out:
        SeqIO.write(unmapped_seqs, out, "fasta")    #  and write it to file.

    # print no. aligned reads from current readtype set:
    print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file))
    if x!=2: print ('')
    prev_mapping_count= len(mapped_reads)

# WRITE OUTPUT: write gene<->read mapfile of BWA-aligned:
# [BWA-aligned geneID, length, #reads, readIDs ...]
reads_count= 0
genes= []
with open(gene2read_file,"w") as out_map:
    for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                        #  (DNA db is needed to get the sequence.)
        if record.id in gene2read_map:                  #  If DNA db gene is one of the matched genes,
            genes.append(record)                        #  append the SeqRec to genes list (NOT REALLY USED), and
            out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                        #  write [aligned geneID, length, #reads, ...],
            for read in gene2read_map[record.id]:
                out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                reads_count+= 1
            else:
                out_map.write("\n")                     #  and a new line character.

# print BWA stats:
print (str(reads_count) + ' reads were mapped with BWA.')
print ('Reads mapped to ' + str(len(genes)) + ' genes.')