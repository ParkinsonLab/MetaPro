#!/usr/bin/env python

import sys
import matplotlib
from matplotlib import cm

# Two methods to define taxa in order of increasing priority: Cutoff or Tax ID list
cutoff = sys.argv[1]                    #IN: Proportion of annotated reads -> 0.01 from the example command list.
ID_list = sys.argv[2]                   #IN: literally a list of taxid we can optionally use.  If nothing, it's empty and it's fine.
nodes = sys.argv[3]                     #IN: nodes.dmp
names = sys.argv[4]                     #IN: names.dmp
gene2read = sys.argv[5]                 #IN: genes -> reads (from GA)
read2taxid = sys.argv[6]                #IN: reads -> taxid (from TA)
gene2EC = sys.argv[7]                   #IN: genes -> EC mapping (from EC, EC.All)
show_unclassified_flag = sys.argv[8]    #IN: either true or false.  Pull from config.  default: true
#raw_count = sys.argv[9]                 #OUT
RPKM = sys.argv[9]                     #OUT
cytoscape = sys.argv[10]                #OUT


show_unclassified = True
if(show_unclassified_flag == "No"):
    show_unclassified = False # temp, should be a user modifiable setting

rank_name = []
if ID_list == "None":
    rank_taxid = []
    print("no taxid supplied.  using blank")
else:
    print("using taxid supplied")
    rank_taxid = [ID.strip() for ID in ID_list.split(",")]

# parse nodes.dmp
nodes_dict = {}
with open(nodes, "r") as infile:
    for line in infile:
        cols = line.split("\t|\t")
        taxid = cols[0]
        parent = cols[1]
        nodes_dict[taxid] = parent

# parse names.dmp
names_dict = {}
with open(names, "r") as infile:
    for line in infile:
        cols = line.split("\t|\t")
        taxid = cols[0]
        name = cols[1]
        type = cols[3].strip("\t|\n")
        if type == "scientific name":
            names_dict[taxid] = name

# parse taxonomic annotation file
read2taxid_dict = {}
taxid_count = {"1": 0}
tax_annotated_count = 0
with open(read2taxid, "r") as infile:
    for line in infile:
        cols = line.split("\t")
        read = cols[0]
        taxid = cols[1].strip("\n")
        if taxid == "0":
            taxid = "1"
        if taxid != "1":
            tax_annotated_count += 1
        read2taxid_dict[read] = taxid
        if taxid in taxid_count:
            taxid_count[taxid] += 1
        else:
            taxid_count[taxid] = 1
    cutoff_count = tax_annotated_count * float(cutoff)

#generate a lineage dictionaries limited to only taxa found in read2taxid_dict
parent_list = {}
child_list = {}
combined_taxid = set()
combined_taxid.update(taxid_count)
combined_taxid.update(rank_taxid)
for taxid in combined_taxid:
    parent_list[taxid] = []
    child_taxid = taxid
    while True:
        try:
            parent_taxid = nodes_dict[child_taxid]
        except:
            parent_taxid = "1"
        if parent_taxid in taxid_count:
            parent_list[taxid].append(parent_taxid)
        if parent_taxid == "0" or parent_taxid == "1":
            break
        child_taxid = parent_taxid
for taxid in taxid_count:
    child_list[taxid] = set()
    for child_taxid in parent_list:
        if taxid in parent_list[child_taxid]:
            child_list[taxid].add(child_taxid)

# adds child nodes to parent nodesif they do not meet cutoff or are not in ID list
# if given no list of taxids, given cutoff is used
if ID_list == "None":
    for taxid in taxid_count:
        processed_children = set()
        while True:
            for child_taxid in child_list[taxid]:
                if child_taxid == "1" or taxid_count[child_taxid] == 0 or taxid_count[child_taxid] > cutoff_count:
                    processed_children.add(child_taxid)
                else:
                    shared_children = child_list[taxid].intersection(child_list[child_taxid])
                    for shared_taxid in shared_children:
                        if shared_taxid not in processed_children:
                            if shared_taxid == "1" or taxid_count[shared_taxid] == 0 or taxid_count[shared_taxid] > cutoff_count:
                                processed_children.add(shared_taxid)
                            else:
                                break
                    else:
                        parent_taxid = parent_list[child_taxid][0]
                        if taxid_count[child_taxid] < cutoff_count and parent_taxid != child_taxid:
                            taxid_count[parent_taxid] += taxid_count[child_taxid]
                            taxid_count[child_taxid] = 0
                        processed_children.add(child_taxid)
                    continue
            if len(processed_children) == len(child_list[taxid]):
                break
        if taxid_count[taxid] > 0:
            rank_taxid.append(taxid)
# if given a list of taxids
else:
    for taxid in taxid_count:
        if taxid not in rank_taxid:
            for x in range(len(parent_list[taxid])):
                parent_taxid = parent_list[taxid][x]
                if parent_taxid in rank_taxid:
                    taxid_count[parent_taxid] += taxid_count[taxid]
                    taxid_count[taxid] = 0
                    break

# modify read to taxid dict with updated taxids
for read in read2taxid_dict:
    taxid = read2taxid_dict[read]
    if taxid_count[read2taxid_dict[read]] == 0:
        for x in range(len(parent_list[taxid])):
            parent_taxid = parent_list[taxid][x]
            if taxid_count[parent_taxid] != 0:
                read2taxid_dict[read] = parent_taxid
                break

# sorting rank taxids on the basis of parentage
parent_sorted_rank_taxid = sorted(rank_taxid, key=lambda child: len(parent_list[child]))
try:
    del parent_sorted_rank_taxid[parent_sorted_rank_taxid.index("1")]
except:
    pass
sorting_list = []
while len(sorting_list) < len(parent_sorted_rank_taxid):
    reverse_list = sorting_list.copy()
    reverse_list.reverse()
    try:
        last = sorting_list[-1]
    except:
        sorting_list.append(parent_sorted_rank_taxid[0])
        last = sorting_list[-1]
    for taxid in parent_sorted_rank_taxid:
        parent_taxid = parent_list[taxid][0]
        if parent_taxid == last and taxid not in sorting_list:
            sorting_list.append(taxid)
            break
    else:
        for x in range(len(reverse_list)):
            for taxid in parent_sorted_rank_taxid:
                parent_taxid = parent_list[taxid][0]
                if parent_taxid == reverse_list[x] and taxid not in sorting_list:
                    sorting_list.append(taxid)
                    break
            else:
                continue
            break
        else:
            for taxid in parent_sorted_rank_taxid:
                if taxid not in sorting_list:
                    sorting_list.append(taxid)
                    break

# genereating rank names
rank_taxid = sorting_list
for taxid in rank_taxid:
    rank_name.append(names_dict[taxid])

if show_unclassified:
    rank_taxid.append("1")
    rank_name.append("Unclassified")

# parse gene annotations
mapped_reads = 0
gene2read_dict = {}
skip_header = True
with open(gene2read, "r") as infile:
    for line in infile:
        if(skip_header):
            skip_header = False
            continue
        cols = line.split("\t")
        gene = cols[0]
        gene_len = cols[1]
        reads = []
        for read in cols[3:]:
            reads.append(read.strip("\n"))
        mapped_reads += len(reads)
        if gene in gene2read_dict:
            gene2read_dict[gene][1].extend(reads)
        else:
            gene2read_dict[gene] = (gene_len, reads)

# parse EC annotations and give back an EC -> genes dict
EC2genes_dict = {}
with open(gene2EC, "r") as infile:
    for line in infile:
        if(len(line) > 1):
            cols = line.split("\t")
            gene = cols[0]
            EC = cols[2].strip("\n")
            if(len(EC) > 1):
                if EC in EC2genes_dict:
                    EC2genes_dict[EC].append(gene)
                else:
                    EC2genes_dict[EC] = [gene]

# Read count and RPKM Tables
#raw_count_dict = {}
RPKM_dict = {}
for gene in gene2read_dict:
    #raw_count_dict[gene] = [gene2read_dict[gene][0], len(gene2read_dict[gene][1])]
    RPKM_div = ((float(gene2read_dict[gene][0])/float(1000))*(mapped_reads/float(1000000)))
    RPKM_dict[gene] = [gene2read_dict[gene][0], len(gene2read_dict[gene][1])]
    EC_string = ""
    for EC in EC2genes_dict:
        if gene in EC2genes_dict[EC]:
            EC_string += EC + "|"
            #raw_count_dict[gene].append(EC)
            #RPKM_dict[gene].append(EC)
            #break
       
    #else:
        #raw_count_dict[gene].append("0.0.0.0")
        #RPKM_dict[gene].append("0.0.0.0")
    if(EC_string == ""):
        EC_string = "0.0.0.0|"
    EC_string = EC_string[:-1]
    RPKM_dict[gene].append(EC_string)    
    RPKM_dict[gene].append(len(gene2read_dict[gene][1])/RPKM_div)
    for taxa in rank_taxid:
        read_count = 0
        for read in gene2read_dict[gene][1]:
            try:
                if read2taxid_dict[read] == taxa:
                    read_count += 1
            except:
                pass
        else:
            #raw_count_dict[gene].append(read_count)
            RPKM_dict[gene].append(read_count / RPKM_div)

#with open(raw_count, "w") as raw_count_out:
#    raw_count_out.write("GeneID\tLength\tReads\tEC#\t" + "\t".join(str(x) for x in rank_name) + "\n")
#    for entry in raw_count_dict:
#        raw_count_out.write(entry + "\t" + "\t".join(str(x) for x in raw_count_dict[entry]) + "\n")
    #raw_count_out.write(",".join(str(x) for x in rank_taxid))
    #raw_count_out.write(",".join(str(x) for x in combined_taxid))



with open(RPKM, "w") as RPKM_out:
    RPKM_out.write("GeneID\tLength\tReads\tEC#\tRPKM\t" + "\t".join(str(x) for x in rank_name) + "\n")
    for entry in RPKM_dict:
        RPKM_out.write(entry + "\t" + "\t".join(str(x) for x in RPKM_dict[entry]) + "\n")
        #print("RPKM key:", entry)
    #raw_count_out.write(",".join(str(x) for x in rank_taxid))
    #RPKM_out.write(",".join(str(x) for x in rank_taxid))
    
# Cytoscape table
rank_colour = []
cs=cm.get_cmap("nipy_spectral", len(rank_taxid))
for i in range(cs.N):
    rgb = cs(i)[:3] # will return rgba, we take only first 3 so we get rgb
    rank_colour.append(matplotlib.colors.rgb2hex(rgb))



#1) ok, this export can't be done in a dictionary as it is.
#2) this table already sums up like-ECs, which is cool.
#3) nope, nvm.  we're already in the clear.   Why? because it collects all like-ECs already.   We don't need it separated by which gene it came from.


Cytoscape_dict = {}
for EC in EC2genes_dict:
    #print("EC:", EC)
    for entry in RPKM_dict:
        #if RPKM_dict[entry][2] == EC:
        ec_list = RPKM_dict[entry][2].split("|")
        for item in ec_list:
            if(item == EC):
                try:
                    for index, RPKM_val in enumerate(Cytoscape_dict[EC]):
                        Cytoscape_dict[EC][index] += RPKM_dict[entry][3 + index]
                except:
                    Cytoscape_dict[EC] = RPKM_dict[entry][3:]
    try:
        Cytoscape_dict[EC].append("piechart: attributelist=\"" + ",".join(str(x) for x in rank_name) + "\" colorlist=\"" + ",".join(str(x) for x in rank_colour) + ",#000000" + "\" showlabels=false\"")
    except:
        pass

#old single-EC format
# Cytoscape_dict = {}
# for EC in EC2genes_dict:
    # for entry in RPKM_dict:
        # if RPKM_dict[entry][2] == EC:
            # try:
                # for index, RPKM_val in enumerate(Cytoscape_dict[EC]):
                    # Cytoscape_dict[EC][index] += RPKM_dict[entry][3 + index]
            # except:
                # Cytoscape_dict[EC] = RPKM_dict[entry][3:]
    # try:
        # Cytoscape_dict[EC].append("piechart: attributelist=\"" + ",".join(str(x) for x in rank_name) + "\" colorlist=\"" + ",".join(str(x) for x in rank_colour) + ",#000000" + "\" showlabels=false\"")
    # except:
        # pass
# print("=============================================")        
# for item in Cytoscape_dict:
    # print(item, len(item))
        
        
with open(cytoscape, "w") as Cytoscape_out:
    Cytoscape_out.write("EC#\tRPKM\t" + "\t".join(str(x) for x in rank_name) + "\tOther\tPiechart\n")
    for entry in Cytoscape_dict:
        Cytoscape_out.write(entry + "\t" + "\t".join(str(x) for x in Cytoscape_dict[entry]) + "\n")
