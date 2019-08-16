#!/usr/bin/env python
#Supposedly, this code is to merge something?
#Need to find out what wevote does
# set intersection, and symmetric_difference are the exact opposites, which only serve to cover the entire region
# What we really want is a union, with no duplicates.
import sys

contig2read_file = sys.argv[1]
output_file = sys.argv[2]

contig2read_map = {}

with open(contig2read_file, "r") as mapping:
    for line in mapping:
        if len(line) > 5:
            entry = line.strip("\n").split("\t")
            contig2read_map[entry[0]] = entry[2:]

Input_classifications = []

for x in range((len(sys.argv) - 3)):
    Input_file = sys.argv[x + 3]
    Input_classification = {}
    with open(Input_file, "r") as Input:
        for line in Input:
            columns = line.split("\t")
            if columns[0] == "score":
                continue
            Seq_ID = columns[1]
            Tax_ID = columns[2].strip("\n")
            if Seq_ID in contig2read_map:
                for Read_ID in contig2read_map[Seq_ID]:
                    if Read_ID not in Input_classification:
                        Input_classification[Read_ID] = Tax_ID
            else:
                Input_classification[Seq_ID] = Tax_ID
        Input_classifications.append(Input_classification)

Combined_Classification = {}

#This loop seems to want to grab only the unique items
#"Read" is a key.  This is a dict
#There will be many files coming in.  This is just a weird (and error-prone) way to do merging.
for Index, Classification in enumerate(Input_classifications):
    Combined_set = set(Combined_Classification)
    Classification_set = set(Classification)   
    for Read in Combined_set.intersection(Classification_set):
        Combined_Classification[Read].append(Classification[Read])
    for Read in Combined_set.symmetric_difference(Classification_set):
        if Read in Combined_Classification:
            Combined_Classification[Read].append("0")
        else:
            if Index > 0:
                Combined_Classification[Read] = ["0"]
                for x in range(Index - 1):
                    Combined_Classification[Read].append("0")
                Combined_Classification[Read].append(Classification[Read])
            else:
                Combined_Classification[Read] = [Classification[Read]]

with open(output_file, "w") as out:
    for Read in Combined_Classification:
        out.write(Read + "," + ",".join(Combined_Classification[Read]) + "\n")