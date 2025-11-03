import os
import sys
from datetime import datetime as dt
from collections import defaultdict


def merge_map(file_dir, master_dict, context):
    """Merge gene maps, keeping highest scoring read assignments"""
    if not os.path.exists(file_dir):
        print(dt.today(), "Directory does not exist:", file_dir)
        return master_dict
        
    map_count = 0
    map_list = os.listdir(file_dir)
    
    for file_name in map_list:
        if file_name.endswith((".tsv")) and file_name.startswith(context):
            full_path = os.path.join(file_dir, file_name)
            map_count += 1
            
            with open(full_path, "r") as map_in:
                for line in map_in:
                    line = line.strip()
                    if not line:
                        continue
                        
                    line_split = line.split("\t")
                    if len(line_split) < 4:
                        continue
                        
                    gene_name = line_split[0]
                    gene_length = line_split[1]
                    read_list = line_split[3:]
                    
                    for read in read_list:
                        if not read or "<" not in read:
                            continue
                            
                        read_parts = read.split("<", 1)
                        read_id = read_parts[0]
                        
                        try:
                            score = int(read_parts[1].split(">")[1])
                        except (IndexError, ValueError):
                            score = 0
                        
                        # Keep higher scoring assignment
                        if read_id in master_dict:
                            if master_dict[read_id]["score"] < score:
                                master_dict[read_id] = {
                                    "score": score,
                                    "gene": gene_name,
                                    "gene_length": gene_length
                                }
                        else:
                            master_dict[read_id] = {
                                "score": score,
                                "gene": gene_name,
                                "gene_length": gene_length
                            }
            
            print("Maps scanned:", map_count, end="\r")
    
    print("Maps scanned:", map_count)
    return master_dict


def transform_map(master_dict, contig_dict):
    """Transform read-based map to gene-based map"""
    new_map_dict = defaultdict(lambda: {"reads": set(), "gene_length": ""})
    
    for read_id, read_data in master_dict.items():
        gene_name = read_data["gene"]
        gene_length = read_data["gene_length"]
        if(read_id in contig_dict):
            for contig_read in contig_dict[read_id]:
                new_map_dict[gene_name]["reads"].add(contig_read)
        else:
            new_map_dict[gene_name]["reads"].add(read_id)

        new_map_dict[gene_name]["gene_length"] = gene_length
    
    return dict(new_map_dict)


def export_map(master_dict, export_path):
    """Export gene map to TSV file"""
    with open(export_path, "w") as out_file:
        for gene_name, inner_data in master_dict.items():
            read_count = len(inner_data["reads"])
            gene_length = inner_data["gene_length"]
            
            out_line = gene_name + "\t" + gene_length + "\t" + str(read_count)
            
            for read_id in inner_data["reads"]:
                out_line += "\t" + read_id
            
            out_file.write(out_line + "\n")

def import_contig_map(contig_map):
    contig_dict = dict()
    with open(contig_map, "r") as contigs_in:
        for raw_line in contigs_in:
            line = raw_line.strip("\n")
            line_split = line.split("\t")
            contig = line_split[0]
            reads = set(line_split[2:])
            contig_dict[contig] = reads
    return contig_dict




if __name__ == "__main__":
    bt2_dir = sys.argv[1]
    dmd_dir = sys.argv[2]
    contig_map_in = sys.argv[3]
    final_map_out = sys.argv[4]

    file_list = os.listdir(bt2_dir)
    print("bt2 files:", file_list)
    print("dmd_files:", os.listdir(dmd_dir))
    
    dir_dict = {
        "bt2": bt2_dir,
        "dmd": dmd_dir
    }
    
    print(dt.today(), "start map merge")
    contig_dict = import_contig_map(contig_map_in)
    
    command_args = ["p", "s", "c"]
    read_dict = {}
    
    for dir_name in dir_dict:
        for context in command_args:
            print(dir_name)
            merge_map(dir_dict[dir_name], read_dict, context)
    
    print("number of reads:", len(read_dict))
    
    print(dt.today(), "start transform")
    gene_dict = transform_map(read_dict, contig_dict)
    print("number of genes:", len(gene_dict))
    
    print(dt.today(), "start export")
    
    export_path = final_map_out
    export_map(gene_dict, export_path)
    print(dt.today(), "end")