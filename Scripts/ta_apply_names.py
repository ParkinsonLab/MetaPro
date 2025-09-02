#!/usr/bin/env python
"""
Taxonomic classification script that converts taxid-based classifications
to human-readable taxonomic hierarchies constrained to a desired level.
"""

import sys
from datetime import datetime as dt
import copy
import time

def import_nodes(nodes_file):
    """
    Import taxonomic node information from NCBI nodes.dmp file.
    
    Returns:
        dict: Maps taxid -> (parent_taxid, rank)
    """
    nodes = {}
    with open(nodes_file, "r") as infile:
        for line in infile:
            cols = line.strip().split("\t|\t")
            if len(cols) >= 3:
                taxid = cols[0]
                parent = cols[1]
                rank = cols[2]
                nodes[taxid] = (parent, rank)
                #print("taxid:", taxid, "parent:", parent, "rank:", rank)
                #time.sleep(0.5)
    
    # Add unclassified entry
    nodes["0"] = ("0", "unclassified")
    return nodes


def import_names(names_file):
    """
    Import taxonomic names from NCBI names.dmp file.
    
    Returns:
        dict: Maps taxid -> scientific_name
    """
    names = {}
    with open(names_file, "r") as infile:
        for line in infile:
            cols = line.strip().split("\t|\t")
            if len(cols) >= 4:
                taxid = cols[0]
                name = cols[1]
                name_type = cols[3]
                if "scientific name" in name_type:
                    names[taxid] = name
    
    # Add unclassified entry
    names["0"] = "unclassified"
    return names


def get_rank_hierarchy():
    """
    Define the standard taxonomic rank hierarchy with prefixes.
    
    Returns:
        dict: Maps rank -> (prefix, level_number)
    """
    return {
        "superkingdom": ("k_", 1),
        "phylum": ("p_", 2),
        "class": ("c_", 3),
        "order": ("o_", 4),
        "family": ("f_", 5),
        "genus": ("g_", 6),
        "species": ("s_", 7)
    }


def build_taxonomy_path(taxid, nodes, names, max_iterations=50):
    """
    Build complete taxonomic path from taxid to root.
    
    Args:
        taxid: Starting taxonomic ID
        nodes: Node information dict
        names: Name information dict
        max_iterations: Maximum iterations to prevent infinite loops
        
    Returns:
        list: List of (taxid, name, rank) tuples from root to species
    """
    if taxid == "0" or taxid not in nodes:
        return []
    
    path = []
    current_taxid = taxid
    iterations = 0
    
    # Traverse up the taxonomic tree
    while current_taxid != "0" and names.get(current_taxid) != "root":
        if iterations > max_iterations:
            # Prevent infinite loops
            break
            
        if current_taxid in names and current_taxid in nodes:
            name = names[current_taxid]
            rank = nodes[current_taxid][1]
            path.append((current_taxid, name, rank))
            current_taxid = nodes[current_taxid][0]
        else:
            break
            
        iterations += 1
    
    # Reverse to get root-to-species order
    return list(reversed(path))


def format_taxonomy_string(taxonomy_path, rank_hierarchy):
    """
    Format taxonomy path into standardized string with proper hierarchy.
    
    Args:
        taxonomy_path: List of (taxid, name, rank) tuples
        rank_hierarchy: Dict mapping ranks to prefixes and levels
        
    Returns:
        str: Formatted taxonomy string
    """
    if not taxonomy_path:
        return "unclassified"
    
    # Create dict of ranks present in path
    path_ranks = {}
    for taxid, name, rank in taxonomy_path:
        if rank in rank_hierarchy:
            path_ranks[rank] = name
    
    # Build taxonomy string in hierarchical order
    taxonomy_parts = ["root"]
    
    for rank in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
        if rank in rank_hierarchy:
            prefix, level = rank_hierarchy[rank]
            if rank in path_ranks:
                taxonomy_parts.append(prefix + path_ranks[rank])
            else:
                # Add empty placeholder to maintain hierarchy
                taxonomy_parts.append(prefix)
    
    # Join and clean up the taxonomy string
    taxonomy = ";".join(taxonomy_parts)
    
    # Remove trailing empty ranks
    while taxonomy.endswith(";_"):
        taxonomy = taxonomy[:-2]
    
    return taxonomy


def process_classifications(classification_file, nodes, names):
    """
    Process the classification file and build taxonomies for each read.
    
    Args:
        classification_file: Path to classification input file
        nodes: Node information dict
        names: Name information dict
        
    Returns:
        dict: Maps read_id -> (taxid, taxonomy_string)
    """
    rank_hierarchy = get_rank_hierarchy()
    read_classifications = {}
    no_hit_taxid_count = 0
    
    with open(classification_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            cols = line.split("\t")
            if len(cols) < 2:
                continue
                
            read_id = cols[1]
            taxid = cols[2].strip()
            #print(dt.today(), "Working on read:", read_id, "taxid:", taxid)
            #time.sleep(1)
            
            if taxid == "0" or taxid not in nodes:
                if(taxid not in nodes):
                    print(dt.today(), taxid, "not in nodes")
                    no_hit_taxid_count += 1
                read_classifications[read_id] = ("0", "unclassified")
            else:
                taxonomy_path = build_taxonomy_path(taxid, nodes, names)
                taxonomy_string = format_taxonomy_string(taxonomy_path, rank_hierarchy)
                read_classifications[read_id] = (taxid, taxonomy_string)
    print("no-hit taxids:", no_hit_taxid_count)
    
    return read_classifications


def write_consensus_file(read_classifications, output_file):
    """
    Write the consensus classification file.
    
    Args:
        read_classifications: Dict of read_id -> (taxid, taxonomy)
        output_file: Path to output file
    """
    with open(output_file, "w") as outfile:
        for read_id, (taxid, taxonomy) in read_classifications.items():
            if taxid == "0":
                status = "U"  # Unclassified
            else:
                status = "C"  # Classified
            
            outfile.write(status + "\t" + read_id + "\t" + taxonomy + "\n")


def main():
    """Main function to orchestrate the taxonomic classification process."""
    if len(sys.argv) != 5:
        print("Usage: {} <classification_file> <nodes_file> <names_file> <output_file>".format(sys.argv[0]))
        sys.exit(1)
    
    classification_file = sys.argv[1]
    nodes_file = sys.argv[3]
    names_file = sys.argv[2]
    consensus_classification_file = sys.argv[4]
    
    print("Loading taxonomic database...")
    names = import_names(names_file)
    nodes = import_nodes(nodes_file)
    
    print("Processing classifications...")
    read_classifications = process_classifications(classification_file, nodes, names)
    
    print("Writing results...")
    write_consensus_file(read_classifications, consensus_classification_file)
    
    print("Classification complete. Processed {} reads.".format(len(read_classifications)))


if __name__ == "__main__":
    main()