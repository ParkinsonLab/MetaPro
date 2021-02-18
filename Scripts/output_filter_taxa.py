#removes classifications where the entry is a contig segment that annotated to 0 reads, due to a lossy contig-segment map

import sys
import pandas as pd

def import_contig_segment(contig_segment_map):
    contig_dict = dict()
    with open(contig_segment_map, "r") as contig_in:
        for line in contig_in:
            cleaned_line = line.strip("\n")
            line_split = cleaned_line.split("\t")
            contig_name = line_split[0]
            read_count = line_split[1]
            contig_dict[contig_name] = (contig_name, read_count)
    contig_df = pd.DataFrame.from_dict(contig_dict, orient = "index", columns = ["contig", "count"])
    contig_df = contig_df.reset_index(drop = True)
    #print(contig_df)
    return contig_df
if __name__ == "__main__":
    constrain_file = sys.argv[1]
    contig_segment_map = sys.argv[2]
    export_constrain_file = sys.argv[3]
    
    constrain_df = pd.read_csv(constrain_file, sep = "\t", index_col = None, header = None, names = ["classified", "read", "taxa"])
    contig_segment_df = import_contig_segment(contig_segment_map)
    
    other_df = constrain_df[~constrain_df["read"].str.startswith("gene")]
    constrain_section = constrain_df[constrain_df["read"].str.startswith("gene")]
    #constrain_section.to_csv("stuff_that_starts_with_gene.tsv", sep = "\t")
    constrain_section = constrain_section[constrain_section["read"].isin(contig_segment_df["contig"])]
    #constrain_section.to_csv("stuff_that_is_only_inside_the_contig_segments.tsv", sep = "\t")
    
    
    final_df = pd.concat([other_df, constrain_section])
    final_df.to_csv(export_constrain_file, sep = "\t", header = None, index = None)