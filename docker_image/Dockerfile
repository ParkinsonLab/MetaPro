FROM parkinsonlab/mpro_stash:v1


#-------------------------------------------
# add the pipeline code -> prototype form.  need to make it pull from a single location in github for all time (maybe something like "latest" or something)
WORKDIR /pipeline
RUN mkdir /pipeline/Scripts && mkdir /pipeline/custom_databases && mkdir /pipeline/setup_scripts

ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/MetaPro.py /pipeline
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/MetaPro_paths.py /pipeline
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/MetaPro_commands.py /pipeline
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Config.ini /pipeline
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/Detect_2.2.8.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/Detect_2.2.9.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/assembly_deduplicate.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/assembly_make_contig_map.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/assembly_flush_bad_contigs.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ea_combine_v5.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ga_BLAT_generic_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ga_BWA_generic_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ga_Diamond_generic_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ga_Final_merge_v4.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_are_you_in_a_contig.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_contig_stats.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_convert_gene_map_contig_segments.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_data_change_metrics.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_EC_metrics.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_filter_ECs.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_filter_taxa.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_read_counts_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_read_quality_metrics.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_taxa_groupby.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_table_v3.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_parse_sam.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_reformat_rpkm_table.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/output_get_host_reads.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_BLAT_filter_v3.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_orphan.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_remove_tag.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_rRNA_barrnap.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_rRNA_infernal.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_repopulation.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_sam.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/bwa_read_sorter.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_sort.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/read_split.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_combine_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_constrain.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_constrain_taxonomy_v2.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_taxid.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_wevote_parser.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/remove_gaps_in_fasta.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/ta_contig_name_convert.py /pipeline/Scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Scripts/alter_taxa_for_krona.py /pipeline/Scripts


#--------------------------------------------
# get reference files
# 
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/custom_databases/EC_pathway.txt /pipeline/custom_databases
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/custom_databases/pathway_to_superpathway.csv /pipeline/custom_databases
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/custom_databases/FREQ_EC_pairs_3_mai_2020.txt /pipeline/custom_databases
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/custom_databases/SwissProt_EC_Mapping.tsv /pipeline/custom_databases


#--------------------------------------------
# get setup scripts
ADD https://raw.githubusercontent.com/billyc59/parkinson_int_meta_pipe/develop/Setup_scripts/Download_EC_Annotated_SwissProt.py /pipeline/setup_scripts


# get mgm license
ADD https://raw.githubusercontent.com/billytaj/parkinson_int_meta_pipe/develop/mgm_license/gm_key /pipeline

RUN chmod -R 555 /pipeline



WORKDIR /
CMD ["bash"]