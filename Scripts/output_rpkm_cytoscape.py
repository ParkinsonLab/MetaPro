#!/usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt

# --- Optimized I/O Functions ---

def import_names(names_file):
    """Parse names.dmp using Pandas to select and filter efficiently."""
    names_df = pd.read_csv(
        names_file, 
        sep=r'\t\|\t', 
        header=None, 
        engine='python',
        names=['taxid', 'name', 'unique_name', 'type'],
        dtype={'taxid': str}
    )
    # The last field on each line still carries a trailing "\t|" (or "\t|\n") from
    # the taxdump line terminator, since there's no closing "\t|\t" after the final
    # real field on the line -- strip that off before comparing.
    names_df['type'] = names_df['type'].str.replace(r'\t?\|\s*$', '', regex=True).str.strip()
    names_dict = names_df[names_df['type'] == "scientific name"].set_index('taxid')['name'].to_dict()
    return names_dict

def import_nodes(nodes_file):
    """Parse nodes.dmp using Pandas."""
    nodes_df = pd.read_csv(
        nodes_file, 
        sep=r'\t\|\t', 
        header=None, 
        engine='python',
        usecols=[0, 1],
        names=['taxid', 'parent_taxid'],
        dtype={'taxid': str, 'parent_taxid': str}
    )
    nodes_dict = nodes_df.set_index('taxid')['parent_taxid'].to_dict()
    return nodes_dict

def import_taxa(taxa_file, cutoff_arg):
    """Read the read2taxid file directly into a DataFrame for initial count/aggregation."""
    taxa_df = pd.read_csv(
        taxa_file, 
        sep="\t", 
        header=None, 
        usecols=[1, 2],
        names=['read', 'taxid'],
        dtype={'read': 'str', 'taxid': 'str'}
    )
    
    taxa_df['taxid'] = taxa_df['taxid'].replace('0', '1')
    
    read2taxid_dict = taxa_df.set_index('read')['taxid'].to_dict()
    
    taxid_count_series = taxa_df['taxid'].value_counts()
    taxid_count = taxid_count_series.to_dict()
    
    tax_annotated_count = taxid_count_series[taxid_count_series.index != '1'].sum()
    cutoff_count = tax_annotated_count * float(cutoff_arg)

    return read2taxid_dict, taxid_count, cutoff_count

def get_lineage_path(taxid, nodes_dict):
    """Generates the full lineage path from a taxid up to the root ('1')."""
    path = []
    current_taxid = taxid
    while current_taxid != '1' and current_taxid != '0':
        parent = nodes_dict.get(current_taxid)
        if parent is None or parent == current_taxid:
            break
        path.append(parent)
        current_taxid = parent
    return path

def aggregate_taxa_counts(taxid_count, nodes_dict, rank_taxid, cutoff_count):
    """Aggregates counts to parent nodes based on cutoff or target list."""
    
    count_df = pd.Series(taxid_count, name='count').to_frame()
    count_df.index.name = 'taxid'
    
    target_taxids = set(rank_taxid) if rank_taxid else set(count_df[count_df['count'] >= cutoff_count].index)
    target_taxids.add('1')
    
    taxids_to_aggregate = set(count_df.index) - target_taxids

    new_counts = count_df['count'].copy()
    
    for taxid in taxids_to_aggregate:
        if new_counts[taxid] == 0:
            continue
            
        path = get_lineage_path(taxid, nodes_dict)
        
        parent_found = None
        for parent_taxid in path:
            if parent_taxid in target_taxids:
                parent_found = parent_taxid
                break
        
        if parent_found:
            new_counts[parent_found] = new_counts.get(parent_found, 0) + new_counts[taxid]
            new_counts[taxid] = 0
            
    final_taxid_count = new_counts[new_counts > 0].to_dict()
    final_rank_taxid = sorted(list(set(taxid for taxid, count in final_taxid_count.items() if count > 0 and taxid in target_taxids)))

    return final_taxid_count, final_rank_taxid

def build_gene_taxa_tables(gene2read_dict, all_gene_reads, read2taxid_dict, EC_string_map, mapped_reads):
    """
    Shared groundwork for both RPKM and TPM: raw per-gene/per-taxon read counts
    (taxa_counts_pivot) and per-gene metadata -- length, total reads, EC string,
    and the RPKM divisor (gene_metadata_df). Computed once and handed to both
    calculate_rpkm() and calculate_tpm() so neither repeats this work.
    """
    gene_read_df = pd.DataFrame(all_gene_reads, columns=['GeneID', 'ReadID'])
    read_taxa_map_series = pd.Series(read2taxid_dict, name='TaxID')
    read_taxa_map_series.index.name = 'ReadID'
    merged_df = gene_read_df.merge(read_taxa_map_series.reset_index(), on='ReadID', how='left')

    taxa_counts_df = merged_df.groupby(['GeneID', 'TaxID']).size().reset_index(name='TaxaReadCount')
    taxa_counts_pivot = taxa_counts_df.pivot(index='GeneID', columns='TaxID', values='TaxaReadCount').fillna(0)

    gene_metadata_df = pd.DataFrame.from_dict({
        gene: {
            'Length': float(gene_len), 
            'Reads': len(reads),
            'EC#': EC_string_map.get(gene, "0.0.0.0"), 
            'RPKM_div': (float(gene_len) / 1000) * (mapped_reads / 1000000)
        } for gene, (gene_len, reads) in gene2read_dict.items()
    }, orient='index')
    gene_metadata_df.index.name = 'GeneID'

    return taxa_counts_pivot, gene_metadata_df

def calculate_rpkm(taxa_counts_pivot, gene_metadata_df, rank_taxid, names_dict):
    """
    Compute the RPKM table: overall per-gene RPKM (Reads / (Length_kb * mapped_reads_millions))
    plus a per-taxon RPKM breakdown. Returns the final, renamed output DataFrame
    (not yet written to disk).
    """
    RPKM_df = taxa_counts_pivot.merge(gene_metadata_df, on='GeneID', how='inner').reset_index()
    RPKM_df['RPKM'] = RPKM_df['Reads'] / RPKM_df['RPKM_div']

    taxa_cols_to_calc = [t for t in taxa_counts_pivot.columns if t in RPKM_df.columns]
    RPKM_df[taxa_cols_to_calc] = RPKM_df[taxa_cols_to_calc].div(RPKM_df['RPKM_div'], axis=0)

    rename_map = {t: names_dict.get(t, t) for t in rank_taxid}
    if "1" in rank_taxid:
        rename_map['1'] = 'Unclassified'

    final_taxa_cols = [t for t in rank_taxid if t in RPKM_df.columns]
    RPKM_out_df = RPKM_df[['GeneID', 'Length', 'Reads', 'EC#', 'RPKM'] + final_taxa_cols]
    RPKM_out_df.rename(columns=rename_map, inplace=True)

    return RPKM_out_df

def calculate_tpm(taxa_counts_pivot, gene_metadata_df, rank_taxid, names_dict):
    """
    Compute the TPM table: overall per-gene TPM plus a per-taxon TPM breakdown.
    TPM differs from RPKM in when the normalization happens:
      RPKM = reads / (length_kb * total_mapped_reads_millions)  -- one divisor, same for every gene's own reads
      TPM  = RPK   / (sum_of_RPK_over_all_genes / 1e6)          -- length-normalize first (RPK), then
                                                                    scale so all genes' TPM values sum to 1e6
    Built from taxa_counts_pivot/gene_metadata_df directly (both untouched by calculate_rpkm(),
    since .merge() returns a new object), so this cannot affect the RPKM table in any way.
    Returns the final, renamed output DataFrame (not yet written to disk).
    """
    TPM_df = taxa_counts_pivot.merge(gene_metadata_df, on='GeneID', how='inner').reset_index()

    taxa_cols_for_tpm = [t for t in taxa_counts_pivot.columns if t in TPM_df.columns]

    # RPK: reads per kilobase, length-normalized only (not yet scaled to sum to 1e6)
    TPM_df['RPK'] = TPM_df['Reads'] / (TPM_df['Length'] / 1000)
    TPM_df[taxa_cols_for_tpm] = TPM_df[taxa_cols_for_tpm].div(TPM_df['Length'] / 1000, axis=0)

    # Single scale factor shared by every gene/taxon column: sum of all genes' RPK, divided by 1e6
    tpm_scale_factor = TPM_df['RPK'].sum() / 1e6

    TPM_df['TPM'] = TPM_df['RPK'] / tpm_scale_factor
    TPM_df[taxa_cols_for_tpm] = TPM_df[taxa_cols_for_tpm].div(tpm_scale_factor)

    rename_map_tpm = {t: names_dict.get(t, t) for t in rank_taxid}
    if "1" in rank_taxid:
        rename_map_tpm['1'] = 'Unclassified'

    final_taxa_cols_tpm = [t for t in rank_taxid if t in TPM_df.columns]
    TPM_out_df = TPM_df[['GeneID', 'Length', 'Reads', 'EC#', 'TPM'] + final_taxa_cols_tpm]
    TPM_out_df.rename(columns=rename_map_tpm, inplace=True)

    return TPM_out_df

# --- Main Execution Block ---

if __name__ == "__main__":
    
    if len(sys.argv) < 11:
        print(f"Usage: python {sys.argv[0]} <cutoff> <ID_list> ... <cytoscape>")
        sys.exit(1)

    cutoff = sys.argv[1]
    ID_list = sys.argv[2]
    nodes = sys.argv[3]
    names = sys.argv[4]
    gene2read = sys.argv[5]
    read2taxid = sys.argv[6]
    gene2EC = sys.argv[7]
    show_unclassified_flag = sys.argv[8]
    RPKM_out_file = sys.argv[9]
    cytoscape_out_file = sys.argv[10]

    # Derived path, not a new CLI argument: keeps this script's call signature unchanged.
    _rpkm_root, _rpkm_ext = os.path.splitext(RPKM_out_file)
    TPM_out_file = _rpkm_root + ".TPM" + _rpkm_ext

    # --- Step 1: Import Data ---
    names_dict = import_names(names)
    nodes_dict = import_nodes(nodes)
    
    show_unclassified = (show_unclassified_flag != "No")
    rank_taxid_initial = [ID.strip() for ID in ID_list.split(",")] if ID_list != "None" else []
        
    read2taxid_dict, taxid_count, cutoff_count = import_taxa(read2taxid, cutoff)
    
    # --- Step 2: Lineage Simplification ---
    final_taxid_count, final_rank_taxid = aggregate_taxa_counts(
        taxid_count, 
        nodes_dict, 
        rank_taxid_initial, 
        cutoff_count
    )
    
    for read, taxid in read2taxid_dict.items():
        if final_taxid_count.get(taxid, 0) == 0:
            path = get_lineage_path(taxid, nodes_dict)
            for parent_taxid in path:
                if final_taxid_count.get(parent_taxid, 0) > 0:
                    read2taxid_dict[read] = parent_taxid
                    break
    
    rank_taxid = final_rank_taxid

    # --- Step 3: Efficient Sorting and Naming ---
    taxid_depth = {}
    for taxid in rank_taxid:
        taxid_depth[taxid] = len(get_lineage_path(taxid, nodes_dict))
    
    rank_taxid.sort(key=lambda t: (taxid_depth.get(t, 0), t))
    
    rank_name = [names_dict.get(taxid, "Unknown") for taxid in rank_taxid]

    if show_unclassified and "1" not in rank_taxid:
        rank_taxid.append("1")
        rank_name.append("Unclassified")
        
    # --- Step 4: Parse Gene and EC Annotations (Correct & Optimized) ---

    gene2read_dict = {}
    all_gene_reads = [] 
    mapped_reads = 0
    
    skip_header = True
    with open(gene2read, "r") as infile:
        for line in infile:
            if skip_header:
                skip_header = False
                continue
            cols = line.split("\t")
            gene = cols[0]
            gene_len = cols[1]
            reads = [read.strip("\n") for read in cols[3:]]
            
            mapped_reads += len(reads)
            if gene in gene2read_dict:
                gene2read_dict[gene][1].extend(reads)
            else:
                gene2read_dict[gene] = (gene_len, reads)
            
            all_gene_reads.extend([(gene, read) for read in reads])

    EC2genes_dict = {}
    with open(gene2EC, "r") as infile:
        for line in infile:
            line = line.strip()
            if line:
                cols = line.split("\t")
                gene = cols[0]
                EC = cols[1].strip()
                if len(EC) > 1:
                    EC2genes_dict.setdefault(EC, []).append(gene)

    EC_string_map = {}
    for gene in gene2read_dict: 
        EC_string = "|".join(EC for EC, genes in EC2genes_dict.items() if gene in genes)
        EC_string_map[gene] = EC_string if EC_string else "0.0.0.0"

    # --- Step 5/6: RPKM Calculation and Output ---

    taxa_counts_pivot, gene_metadata_df = build_gene_taxa_tables(
        gene2read_dict, all_gene_reads, read2taxid_dict, EC_string_map, mapped_reads
    )

    RPKM_out_df = calculate_rpkm(taxa_counts_pivot, gene_metadata_df, rank_taxid, names_dict)
    RPKM_out_df.to_csv(RPKM_out_file, sep="\t", index=False, float_format='%.6f')

    # --- Step 6b: TPM Calculation and Output ---

    TPM_out_df = calculate_tpm(taxa_counts_pivot, gene_metadata_df, rank_taxid, names_dict)
    TPM_out_df.to_csv(TPM_out_file, sep="\t", index=False, float_format='%.6f')


    # --- Step 7: Vectorized Cytoscape Table Generation (FINAL FIX applied here) ---
    
    rank_colour = []
    cs = plt.get_cmap("nipy_spectral", len(rank_taxid))
    for i in range(cs.N):
        rgb = cs(i)[:3]
        rank_colour.append(matplotlib.colors.rgb2hex(rgb))

    ec_gene_map = []
    for EC, genes in EC2genes_dict.items():
        for gene in genes:
            ec_gene_map.append({'EC#_single': EC, 'GeneID': gene})
    ec_gene_df = pd.DataFrame(ec_gene_map)

    # Use the renamed columns from RPKM_out_df
    # Note: RPKM_out_df contains 'RPKM' and the final renamed taxa columns
    cytoscape_cols_for_merge = [c for c in RPKM_out_df.columns if c not in ['GeneID', 'Length', 'Reads', 'EC#']]
    
    cytoscape_temp_df = ec_gene_df.merge(RPKM_out_df[['GeneID'] + cytoscape_cols_for_merge], on='GeneID', how='inner')

    # Summing all the RPKM columns based on what's left after merge
    cols_to_sum = [c for c in cytoscape_temp_df.columns if c not in ['EC#_single', 'GeneID']]
    Cytoscape_final_df = cytoscape_temp_df.groupby('EC#_single')[cols_to_sum].sum().reset_index()
    
    Cytoscape_final_df.rename(columns={'EC#_single': 'EC#'}, inplace=True)

    # 1. Add 'Other' and 'Piechart' columns
    Cytoscape_final_df['Other'] = 0.0
    
    rank_name_str = ",".join(str(x) for x in rank_name)
    rank_colour_str = ",".join(str(x) for x in rank_colour) + ",#000000"
    piechart_base = f"piechart: attributelist=\"{rank_name_str}\" colorlist=\"{rank_colour_str}\" showlabels=false\""
    
    Cytoscape_final_df['Piechart'] = piechart_base
    
    # 2. FINAL CRITICAL FIX: Build the column list based on columns PRESENT in the DataFrame
    
    # The list of RPKM columns the script *expects* in the output order:
    expected_rpkm_order = ['RPKM'] + rank_name 
    
    # The list of columns actually present in the DataFrame:
    present_cols = list(Cytoscape_final_df.columns)
    
    # Select only the columns that are PRESENT, in the expected order
    final_ordered_rpkm_cols = [c for c in expected_rpkm_order if c in present_cols]
    
    # The complete column list for the final write
    final_cytoscape_cols = ['EC#'] + final_ordered_rpkm_cols + ['Other', 'Piechart']
    
    # 3. Write the Cytoscape table
    Cytoscape_final_df[final_cytoscape_cols].to_csv(cytoscape_out_file, sep="\t", index=False, header=True, float_format='%.6f')