import os
import sys
import pandas as pd
#remove all contig segments not in the map

def import_contig_map_headers(contig_map):
    header_list = []
    #just grab the headers.  we don't need the reads
    with open(contig_map, "r") as map_file:
        for line in map_file:
            line_split = line.split("\t")
            header = line_split[0]
            if(not header_list):
                header_list = [">" + header]
            else:
                header_list.append(">" + header)
            
    #for item in header_list:
    #    print(item)
    return header_list
            
            
def import_fasta_plain(file_name_in):
    fasta_dict = {}
    ID = "0"
    with open(file_name_in) as fasta_file:
        for raw_line in fasta_file:
            line = raw_line.strip("\n")
            if(line.startswith(">")):
                ID = line
            else:
                #print("===============================================")
                #print("ID:", ID, "line:", line)
                if(ID not in fasta_dict):
                    fasta_dict[ID] = line
                else:
                    fasta_dict[ID] += line
    return fasta_dict

if __name__ == "__main__":
    contig_map_file = sys.argv[1]
    contigs_fasta = sys.argv[2]
    export_path = sys.argv[3]
    
    header_list = import_contig_map_headers(contig_map_file)
    contigs_dict = import_fasta_plain(contigs_fasta)
    
    contigs_df = pd.DataFrame.from_dict(contigs_dict, orient = "index", columns = ["read"])
    contigs_df["ID"] = contigs_df.index
    contigs_df.reset_index(drop = True, inplace = True)
    
    final_contigs_df = contigs_df[contigs_df["ID"].isin(header_list)]
    cols = ["ID", "read"]
    final_contigs_df = final_contigs_df[cols]
    print(final_contigs_df)
    final_contigs_df.to_csv(export_path, sep = "\n", quoting = 3, header = False, index = False, mode = "w")
    print("old:", contigs_df.shape)
    
    print("new", final_contigs_df.shape)
    
    
    