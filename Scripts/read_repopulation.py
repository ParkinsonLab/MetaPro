#!/usr/bin/env python
#This file is just supposed to repopulate the mRNA with the right number of duplicates
#but in order to do so, it must find out what's mRNA.  
#There's better ways to do this, but we're told that it's such a minor thing, we should leave it be.
#It also sorta ruins the modularity aspect, since we now have to include more than the original 3 files, but.... yeah...
#We need a clever solution for this.
#Also, this needs to check for file existence.  If there's missing files in the arg, don't run
#also this thing is totally broken
import sys
import pandas as pd

def repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename):
    ref_file = pd.read_csv(ref_filename, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting=3)
    ref_df = pd.DataFrame(ref_file.values.reshape(int(len(ref_file)/4), 4))
    ref_df.columns = ["ID", "seq", "junk", "quality"]
    ref_df["ID"] = ref_df["ID"].apply(lambda x: x.split(" ")[0])

    mRNA_file = pd.read_csv(mRNA_filename, header=None, names=[None], sep = '\n', skip_blank_lines = False, quoting=3)
    mRNA_df = pd.DataFrame(mRNA_file.values.reshape(int(len(mRNA_file)/4), 4))
    mRNA_df.columns = ["ID", "seq", "junk", "quality"]
    mRNA_df["ID"] = mRNA_df["ID"].apply(lambda x: x.split(" ")[0])
    cluster_file = cluster_filename
    cluster_map = {}
    full_mRNA_file = output_filename

    reduplicated_ids = set()
    reduplicated_seqs = []

    #There's some things that the duplicate remover deems as a duplicate, but comes as a different ID.  
    #scanning the cluster like this means we keep it all.  We must scan the cluster file
    with open(cluster_file, "r") as clustr_read:
        rep = ""
        seq_id = ""
        for line in clustr_read:
            if line.startswith(">"):
                continue
            elif line.startswith("0"):
                rep = line[line.find(">") + 1:line.find("...")]
                seq_id = rep
                cluster_map[rep] = [seq_id]
            elif len(line) > 5:
                seq_id = line[line.find(">") + 1:line.find("...")]
                cluster_map[rep].append(seq_id)

    #puts the cluster and mRNA IDs together in a single list
    for sequence in mRNA_df["ID"]:
        sequence = sequence[1:] #removing @ from fastq
        if sequence in cluster_map:
            if len(cluster_map[sequence]) > 1:
                for seq_id in cluster_map[sequence]:
                    reduplicated_ids.add("@" + seq_id)
            else:
                reduplicated_ids.add("@" + sequence)
        else:
            reduplicated_ids.add("@" + sequence)

    #exports the full mRNA by fetching from ref
    ref_df[ref_df.ID.isin(sorted(reduplicated_ids))].to_csv(full_mRNA_file, sep = '\n', mode = "w+", header = False, index = False, quoting = 3)


if __name__ == "__main__":
    ref_filename = sys.argv[1]      #in: the file that the duplicates will come from
    mRNA_filename = sys.argv[2]     #in: the file that will contain the reads needing to be duplicated
    cluster_filename = sys.argv[3]  #in: the cluster file showing what was in-fact duplicated
    output_filename = sys.argv[4]   #out: the final output
    repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename)
