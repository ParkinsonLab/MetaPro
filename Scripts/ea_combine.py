#!/usr/bin/env python

import os.path
import sys

Input_File = sys.argv[1]
Input_Name = os.path.splitext(os.path.basename(Input_File))[0]
detect_file = sys.argv[2]
detect_dir = os.path.dirname(detect_file)
priam_file = sys.argv[3]
priam_dir = os.path.dirname(priam_file)
diamond_file = sys.argv[4]
diamond_dir = os.path.dirname(diamond_file)
SWISS_PROT = sys.argv[5]
SWISS_PROT_MAP = sys.argv[6]
Output_Dir = sys.argv[7]

mapping_dict = {}
with open(SWISS_PROT_MAP, "r") as mapping:
    for line in mapping.readlines():
        line_as_list = line.split("\t")
        mapping_dict[line_as_list[0]] = set(line_as_list[2:])

detect_ECs = os.path.join(detect_dir, Input_Name + ".toppred.cutoff")
with open(detect_file, "r") as topred:
    with open(detect_ECs, "w") as cutoff:
        for line in topred.readlines():
            line_as_list = line.split("\t")
            if line_as_list[2] == "probability":
                continue
            if float(line_as_list[2]) >= 0.2 and int(line_as_list[3]) > 5:
                cutoff.write(line)

priam_ECs = os.path.join(priam_dir, Input_Name + ".ECs")
with open(priam_file, "r") as ECs:
    with open(priam_ECs, "w") as processedECs:
        ID = None
        EC = None
        for line in ECs.readlines():
            if line.startswith(">"):
                ID = line.split(" ")[0][1:].strip("\n")
                continue
            elif line == "\n":
                continue
            elif ID and not line.startswith("#"):
                EC = line.split(" ")[0]
                processedECs.write("\t".join([ID, EC]))
                processedECs.write("\n")
                ID = None
                EC = None
                continue

diamond_ECs = os.path.join(diamond_dir, Input_Name + ".ECs")
with open(diamond_file, "r") as blastout:
    with open(diamond_ECs, "w") as ecout:
        for line in blastout.readlines():
            line_as_list = line.strip().split("\t")
            for EC in mapping_dict:
                if line_as_list[1] in mapping_dict[EC]:
                    ecout.write("\t".join([line_as_list[0], EC + "\n"]))

with open(os.path.join(Output_Dir, Input_Name + ".ECs_PB"), "w") as PB_out:
    with open(priam_ECs, "r") as priam_ECs_in:
        priam_preds = priam_ECs_in.readlines()
    with open(diamond_ECs, "r") as diamond_ECs_in:
        diamond_preds = diamond_ECs_in.readlines()
    PB_preds = []
    
    #take the ones common in BLAST (diamond) and PRIAM.  
    for priam_ec in priam_preds:
        if priam_ec in diamond_preds:
            PB_preds.append(priam_ec)
    PB_out.writelines(PB_preds)
with open(detect_ECs, "r") as detect_ECs_in:
    detect_preds = []
    for line in detect_ECs_in.readlines():
        line_as_list = line.split("\t")
        line_as_list = "\t".join(line_as_list[:2]) + "\n"
        detect_preds.append(line_as_list)
All_preds = set()

#DETECT-2 is designed to take all of the ECs.  not just one.  We will no longer discard the other ECs.  We must take them all from DETECT-2

#take everything from PRIAM. only take the unique ones (there's multiple of the same).  There will also be a score cutoff (it's the first number).  build in a feature to change the EC cutoffs
#There will also be multiple ECs from swissprot.  Take them all.  Also another cutoff for the BLAST results. (we need a number.) also make it parameterizable
for pred in detect_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
for pred in PB_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
with open(os.path.join(Output_Dir, Input_Name + ".ECs_All"), "w") as ec_out:
    ec_out.writelines(sorted(All_preds))
