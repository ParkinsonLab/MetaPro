#this will reverse-engineer the host reads (unique) from the pipeline.
#"Where do I run this?" 
#-> on a node with singularity loaded.  use singularity as a virtual env, or if on niagara, load up anaconda3.  It needs pandas.
#"What do I supply the arguments?"
#-> the singletons from the host_read_filter/final_results (this is post-host)
#-> the singletons from before the host reads were filtered out (this is pre-host.  it'll be singletons.fastq, from quality_filter/final_results.  the smaller of the 2 singletons fastq files)
#-> some arbitrary output file with the extension ".fastq"

#The next step (after you run this code) will be to call the "read_repopulation.py" code
#you will supply <read_repopulation.py> with (and in the following order):

#-> the file with duplicate reads (singletons_with_duplicates.fastq, found in quality_filter/final_results)
#-> the host reads you extracted (the output of this code)
#-> a .clstr file (the one you want to get the duplicate reads from)
#-> the name of the output (with the extension ending in .fastq)

#Dec 14, 2020:  this also works for vectors, and anything that needs reversing.

import os
import sys
import pandas as pd

def import_fastq(fastq_file):
    ref_file = pd.read_csv(fastq_file, header = None, names = [None], sep = '\n', skip_blank_lines = False, quoting=3)
    ref_df = pd.DataFrame(ref_file.values.reshape(int(len(ref_file)/4), 4))
    ref_df.columns = ["ID", "seq", "junk", "quality"]
    ref_df["ID"] = ref_df["ID"].apply(lambda x: x.split(" ")[0])
    
    return ref_df


if __name__ == "__main__":
    post_host_file  = sys.argv[1]   #host_read_filter_final
    pre_host_file   = sys.argv[2]   #quality_filter
    output_file     = sys.argv[3]   #output
    
    post_host_df = import_fastq(post_host_file)
    pre_host_df = import_fastq(pre_host_file)
    print("=======================================")
    print(post_host_df)
    print("=======================================")
    print(pre_host_df)
    #grab the stuff that doesn't overlap.
    host_df = pre_host_df[~pre_host_df["ID"].isin(post_host_df["ID"])]
    print("=======================================")
    print(host_df)
    
    host_df.to_csv(output_file, sep = "\n", index = False, header = False, quoting = 3)