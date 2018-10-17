#!/usr/bin/env python
#Jan 09, 2018
#-------------------------------------------------------
# The base version of this code comes from Jordan Ang
# This code tries to figure out what's been lost due to contig assembly
# may 07,2018
# We've got a competing version called "contig duplicate remover"
# code's been changed to only produce the missing contig map
#---------------------------------------------------------
# Now with some commenting!
# CHANGES:
# - rewrote contig_map(), also added a check to remove mapped reads from the unmapped set.

# NOTES:
# - For unmerged pairs, where only one of a pair is mapped to a contig, the other
#   will be discarded, and the read considered to be part of the contig.
# - It is possible for a read to align to multiple contigs.

import os
import os.path
import sys
from collections import defaultdict

contig_map_file = sys.argv[1]            # OUTPUT: contig<->read map: [contigID, #reads, readIDs ...]

# tracking BWA-aligned:
contig2read_map= defaultdict(set)   # dict of BWA-aligned contigID<->readID(s)
mapped_reads= set()                 # tracks BWA-assigned reads
mapped_list= []

#####################################
# FUNCTION:
def make_contig_map(sam, unmapped, mapped_reads, contig2read_map):

    # process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
        
            # valid line?
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  skip to the next line.
        
            # extract data:
            line_parts= line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and extract:
            readID= line_parts[0]                       #  the readID,
            flag= bin(int(line_parts[1]))[2:].zfill(11) #  the flag--converting to 11-digit binary:
                                                        #  each bit is a flag for a specific descriptor.
            # is read BWA-aligned?:
            if flag[8]=="1":                            # If read is NOT BWA-ALIGNED (9th digit=1),
                unmapped.add(readID)                    #  add it to the unmapped set and
                continue                                #  skip to the next read.
        
            # store info in map:
            contigID= line_parts[2]                     # Extract the contigID, and
            contig2read_map[contigID].add(readID)       #  add the readID to the contigID dict entry
                                                        #  (it's a set, so it won't store duplicates, though it
                                                        #  is possible for a read to align to multiple contigs),
            mapped_reads.add(readID)                    #  and mark the read as mapped to a contig.

    # Remove reads previously added to the unmapped set but later found to align to a contig:
    # Occurs for unmerged pairs, where only one of a pair is mapped to a contig, the other will
    # need to be removed from the unmapped set.
    print ('\n' + str(os.path.basename(sam)))
    print ('umapped no. (before double-checking mapped set)= ' + str(len(unmapped)))
    for read in mapped_reads:                           # Take all mapped reads and
        try:                                            #  if the exist in the unmapped set, then
            unmapped.remove(read)                       #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))

#####################################

# extraction of contig-aligned reads and unaligned reads:
unmapped_paired_reads= set()
unmapped_unpaired_reads= set()

for sam_file in sys.argv[2:]:
    make_contig_map(sam_file, unmapped_paired_reads, mapped_reads, contig2read_map)

# WRITE OUTPUT: contig<->read map
# [contigID, #reads in contig, read IDs...]
with open(contig_map_file,"w") as outfile:
    for contig in contig2read_map:                      # For each contig,
        outfile.write(contig + "\t" + str(len(contig2read_map[contig])))
                                                        # write [contigID, #reads in contig, ...]
        for read in list(contig2read_map[contig]):
            outfile.write("\t" + read)                  # [readIDs ...],
        else:
            outfile.write("\n")                         # and a new line character.
