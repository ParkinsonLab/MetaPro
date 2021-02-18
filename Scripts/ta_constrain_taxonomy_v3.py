#!/usr/bin/env python
#this code is supposed to constrain the taxonomy reporting down to a desired level.

import sys
import copy



#Read in library info "Nodes" and "Names"
#These 2 files are the taxonomic keys of every microbe we're got.  

def import_nodes(nodes_file):
    nodes = {}
    with open(nodes_file, "r") as infile:
        for line in infile:
            cols = line.split("\t|\t")
            taxid = cols[0]
            parent = cols[1]
            rank = cols[2]
            nodes[taxid] = (parent, rank)
        else:
            nodes["0"] = ("0", "unclassified")
    return nodes
        
def import_names(names_file):
    names = {}
    with open(names_file, "r") as infile:
        for line in infile:
            cols = line.split("\t|\t")
            taxid = cols[0]
            name = cols[1]
            if "scientific name" in cols[3]:
                names[taxid] = name
        else:
            names["0"] = "unclassified"
    return names
    
    
if __name__ == "__main__":
    target_rank = sys.argv[1]
    classification_file = sys.argv[2]
    nodes_file = sys.argv[3]
    names_file = sys.argv[4]
    consensus_classification_file = sys.argv[5]
        
    Ranks = ["unclassified", "no rank", "superkingdom", "phylum", "subphylum", "class", "subclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily","tribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies"]

    names = import_names(names_file)
    nodes = import_nodes(nodes_file)

    read_classifications = {}
    classification_count = {}
    # Go through the file, and convert the info into english
    #june 20, 2019:  Why was there a letter C in the file in the first place?? anyway, it's gone.  the columns are re-aligned. 
    with open(classification_file, "r") as infile:
        for line in infile:
            if len(line) < 1:
                continue
            cols = line.split("\t")
            read = cols[0]
            taxid = cols[1].strip("\n")
            if taxid != "0" and taxid in nodes:
                taxid2 = copy.deepcopy(taxid)
                taxonomy = ""
                nlevels = 0
                #species_group = "FALSE"
                if names[taxid2] == "root":
                    taxonomy = "root"
                else:
                    i = 0
                    while names[taxid2] != "root":
                        #print(Ranks.index(nodes[taxid2][1]))
                        if nlevels > 0:
                            taxonomy = ";"+taxonomy
                        i+=1
                        if i > 20:
                            taxonomy="no root" + taxonomy
                            break
                        #if nodes[taxid2][1] == "species group":
                        #    species_group = "TRUE"
                        #    taxonomy = "sg_"+names[taxid2]+taxonomy
                        #    nlevels = nlevels + 1
                        if nodes[taxid2][1] == "species":
                            taxonomy = "s_"+names[taxid2]+taxonomy
                            nlevels = 1
                        elif nodes[taxid2][1] == "genus":
                            taxonomy = "g_"+names[taxid2]+taxonomy
                            nlevels = 2
                        elif nodes[taxid2][1] == "family":
                            if nlevels!=0 and nlevels!=2:
                                if nlevels==1:
                                    taxonomy=";g_" + taxonomy
                            taxonomy = "f_"+names[taxid2]+taxonomy
                            nlevels = 3
                        elif nodes[taxid2][1] == "order":
                            if nlevels!=0 and nlevels!=3:
                                if nlevels==1:
                                    taxonomy=";f_;g_" + taxonomy
                                elif nlevels==2:
                                    taxonomy=";f_" + taxonomy
                            taxonomy = "o_"+names[taxid2]+taxonomy
                            nlevels = 4
                        elif nodes[taxid2][1] == "class":
                            if nlevels!=0 and nlevels!=4:
                                if nlevels==1:
                                    taxonomy=";o_;f_;g_" + taxonomy
                                elif nlevels==2:
                                    taxonomy=";o_;f_" + taxonomy
                                elif nlevels==3:
                                    taxonomy=";o_" + taxonomy
                            taxonomy = "c_"+names[taxid2]+taxonomy
                            nlevels = 5
                        elif nodes[taxid2][1] == "phylum":
                            if nlevels!=0 and nlevels!=5:
                                if nlevels==1:
                                    taxonomy=";c_;o_;f_;g_" + taxonomy
                                elif nlevels==2:
                                    taxonomy=";c_;o_;f_" + taxonomy
                                elif nlevels==3:
                                    taxonomy=";c_;o_" + taxonomy
                                elif nlevels==4:
                                    taxonomy=";c_" + taxonomy
                            taxonomy = "p_"+names[taxid2]+taxonomy
                            nlevels = 6
                        elif nodes[taxid2][1] == "superkingdom":
                            if nlevels!=0 and nlevels!=6:
                                if nlevels==1:
                                    taxonomy=";p_;c_;o_;f_;g_" + taxonomy
                                elif nlevels==2:
                                    taxonomy=";p_;c_;o_;f_" + taxonomy
                                elif nlevels==3:
                                    taxonomy=";p_;c_;o_" + taxonomy
                                elif nlevels==4:
                                    taxonomy=";p_;c_" + taxonomy
                                elif nlevels==5:
                                    taxonomy=";p_" + taxonomy
                            taxonomy = "k_"+names[taxid2]+taxonomy
                            nlevels = 7
                        taxid2 = nodes[taxid2][0]
                    else:
                        if nlevels!=0 and nlevels!=7:
                            if nlevels==1:
                                taxonomy="k_;p_;c_;o_;f_;g_" + taxonomy
                            elif nlevels==2:
                                taxonomy="k_;p_;c_;o_;f_;" + taxonomy
                            elif nlevels==3:
                                taxonomy="k_;p_;c_;o_;" + taxonomy
                            elif nlevels==4:
                                taxonomy="k_;p_;c_;" + taxonomy
                            elif nlevels==5:
                                taxonomy="k_;p_;" + taxonomy
                            elif nlevels==6:
                                taxonomy="k_;" + taxonomy
                        taxonomy = "root;" + taxonomy
                    taxonomy=taxonomy.replace(";;;;",";")
                    taxonomy=taxonomy.replace(";;;",";")
                    taxonomy=taxonomy.replace(";;",";")
                read_classifications[read] = (taxid, taxonomy)
            else:
                read_classifications[read] = ("0", names["0"])

    with open(consensus_classification_file, "w") as outfile:
        for read in read_classifications:
            if read_classifications[read][0] == "0":
                outfile.write("U\t" + read + "\t" + read_classifications[read][1] + "\n")
            else:
                outfile.write("C\t" + read + "\t" + read_classifications[read][1] + "\n")
