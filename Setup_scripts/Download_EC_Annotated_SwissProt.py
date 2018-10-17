#!/usr/bin/env python
# Sept 22, 2018
# This code downloads and creates the swiss_prot_EC_map.tsv needed for the MetaPro Pipeline
# instructions:  Copy this file to a desired location and just run.  Code will get dumped in the same directory as this file
# Linux OS only, because of the wget calls

import subprocess

subprocess.call(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"])
subprocess.call(["gzip", "-d", "uniprot_sprot.fasta.gz"])
subprocess.call(["wget", "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"])

enzyme_dict = {}
EC_dict = {}
ID_EC = ""
Name_EC = ""
ID_EC_SwissProt = []

with open("enzyme.dat", "r") as enzyme_annotations:
	for line in enzyme_annotations.readlines():
		if line.startswith("ID"):
			if len(ID_EC_SwissProt) > 0:
				enzyme_dict[ID_EC] = ID_EC_SwissProt
				EC_dict[ID_EC] = Name_EC[:-1]
				ID_EC_SwissProt = []
				Name_EC = ""
			ID_EC = line[5:].strip()
		elif line.startswith("DE") and Name_EC == "":
			Name_EC = line[5:].strip()
		elif line.startswith("DR"):
			line = line[5:].split(";")
			for entry in line:
				entry = entry.strip()
				entry = entry[:entry.find(",")]
				if len(entry) > 1:
					ID_EC_SwissProt.append(entry)

with open("SwissProt_EC_Mapping.tsv", "w") as enzyme_map:
	for EC in sorted(enzyme_dict.keys()):
		enzyme_map.write(EC + "\t")
		enzyme_map.write(EC_dict[EC] + "\t")
		enzyme_map.write("\t".join(enzyme_dict[EC]) + "\n")

Swiss_Prot_Fasta = {}
ID_2_ID_Map = {}
Fasta_ID = ""
Fasta_ID_Trimmed = ""

with open("uniprot_sprot.fasta", "r") as Swiss_In:
	for line in Swiss_In.readlines():
		if line.startswith(">"):
			Fasta_ID = line.strip()
			Fasta_ID_Trimmed = Fasta_ID[4:Fasta_ID.find("|", 5)]
			Swiss_Prot_Fasta[Fasta_ID] = ""
			ID_2_ID_Map[Fasta_ID_Trimmed] = Fasta_ID
		elif Fasta_ID != "":
			Swiss_Prot_Fasta[Fasta_ID] += line.strip()

Annotated_SwissProt_IDs = []

for ID in enzyme_dict:
	Annotated_SwissProt_IDs.extend(enzyme_dict[ID])

Annotated_SwissProt_Fasta = {}

for ID in Annotated_SwissProt_IDs:
	try:
		Annotated_SwissProt_Fasta[ID_2_ID_Map[ID]] = Swiss_Prot_Fasta[ID_2_ID_Map[ID]]
	except:
		pass

with open("uniprot_sprot_annotated.fasta", "w") as Swiss_Out:
	for desc in Annotated_SwissProt_Fasta:
		Swiss_Out.write(desc.replace("sp|", "").replace("|", " ") + "\n")
		Swiss_Out.write(Annotated_SwissProt_Fasta[desc] + "\n")