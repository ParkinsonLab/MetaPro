
import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


#--------------------------------------------------------------------------------------        
class dir_path_obj:
#-------------------------------------------
# object for all paths



    def __init__ (self, config_path, args_pack):
        print("CHECKING CONFIG")
        if config_path:
            self.config = ConfigParser() #change this to ex
            self.config.read(config_path)
            print("USING CONFIG", config_path)
        else:
            print("no config found, defaulting")
            self.config = None

        self.raw_forward = args_pack["forward"]
        self.raw_reverse = args_pack["reverse"]
        self.raw_single = args_pack["single"]
        self.out_dir = args_pack["out_dir"]


        
        
    
        #--------------------------------------------------------------------------------------------
        # Labels.  
        # why? to change them during integration + new feature testing

        quality_filter_label_default                    = "quality_filter"
        host_filter_label_default                       = "host_filter"
        vector_filter_label_default                     = "vector_filter"
        rRNA_filter_label_default                       = "rRNA_filter"
        rRNA_filter_split_label_default                 = "rRNA_filter_split"
        rRNA_filter_convert_label_default               = "rRNA_filter_convert"
        rRNA_filter_barrnap_label_default               = "rRNA_filter_barrnap"
        rRNA_filter_barrnap_merge_label_default         = "rRNA_filter_barrnap_merge"
        rRNA_filter_barrnap_pp_label_default            = "rRNA_filter_barrnap_pp"
        rRNA_filter_infernal_label_default              = "rRNA_filter_infernal"
        rRNA_filter_infernal_prep_label_default         = "rRNA_filter_infernal_prep"
        rRNA_filter_splitter_label_default              = "rRNA_filter_splitter"
        rRNA_filter_post_label_default                  = "rRNA_filter_post"
        repop_label_default                             = "duplicate_repopulation"
        assemble_contigs_label_default                  = "assemble_contigs"
        destroy_contigs_label_default                   = "destroy_contigs"
        GA_pre_scan_label_default                       = "GA_pre_scan"
        GA_split_label_default                          = "GA_split"
        GA_BWA_label_default                            = "GA_BWA"
        GA_BWA_pp_label_default                         = "GA_BWA_pp"
        GA_BWA_merge_label_default                      = "GA_BWA_merge"
        GA_BLAT_label_default                           = "GA_BLAT"
        GA_BLAT_cleanup_label_default                   = "GA_BLAT_cleanup"
        GA_BLAT_cat_label_default                       = "GA_BLAT_cat"
        GA_BLAT_pp_label_default                        = "GA_BLAT_pp"
        GA_BLAT_merge_label_default                     = "GA_BLAT_merge"
        GA_DIAMOND_label_default                        = "GA_DMD"
        GA_DIAMOND_pp_label_default                     = "GA_DMD_pp"
        GA_final_merge_label_default                    = "GA_final_merge"
        taxon_annotation_label_default                  = "taxonomic_annotation"
        ec_annotation_label_default                     = "enzyme_annotation"
        ec_annotation_detect_label_default              = "enzyme_annotation_detect"
        ec_annotation_priam_label_default               = "enzyme_annotation_priam"
        ec_annotation_priam_split_label_default         = "enzyme_annotation_priam_split"
        ec_annotation_priam_cat_label_default           = "enzyme_annotation_priam_cat"
        ec_annotation_DIAMOND_label_default             = "enzyme_annotation_DMD"
        ec_annotation_pp_label_default                  = "enzyme_annotation_pp"
        output_label_default                            = "outputs"
        output_copy_gene_map_label_default              = "output_copy_gene_map"
        output_clean_EC_label_default                   = "output_clean_ec"
        output_copy_taxa_label_default                  = "output_copy_taxa"
        output_network_gen_label_default                = "output_network_generation"
        output_unique_hosts_singletons_label_default    = "output_unique_hosts_singletons"
        output_unique_hosts_pair_1_label_default        = "output_unique_hosts_pair_1"
        output_unique_hosts_pair_2_label_default        = "output_unique_hosts_pair_2"
        output_unique_vectors_singletons_label_default  = "output_unique_vectors_singletons"
        output_unique_vectors_pair_1_label_default      = "output_unique_vectors_pair_1"
        output_unique_vectors_pair_2_label_default      = "output_unique_vectors_pair_2"
        output_combine_hosts_label_default              = "output_combine_hosts"
        output_per_read_scores_label_default            = "output_per_read_scores"
        output_contig_stats_label_default               = "output_contig_stats"
        output_ec_heatmap_label_default                 = "output_ec_heatmap"
        output_taxa_groupby_label_default               = "output_taxa_groupby"
        output_read_count_label_default                 = "output_read_count"
        

        self.qc_label                   = self.value_assignment(self.config, "Labels", "quality_filter",                     quality_filter_label_default)
        self.host_label                      = self.value_assignment(self.config, "Labels", "host_filter",                        host_filter_label_default)
        self.vector_label                    = self.value_assignment(self.config, "Labels", "vector_filter",                      vector_filter_label_default)
        self.rRNA_filter_label                      = self.value_assignment(self.config, "Labels", "rRNA_filter",                        rRNA_filter_label_default)
        self.repop_label                            = self.value_assignment(self.config, "Labels", "repop",                              repop_label_default)
        self.assemble_contigs_label                 = self.value_assignment(self.config, "Labels", "assemble_contigs",                   assemble_contigs_label_default)
        self.destroy_contigs_label                  = self.value_assignment(self.config, "Labels", "destroy_contigs",                    destroy_contigs_label_default)
        self.GA_pre_scan_label                      = self.value_assignment(self.config, "Labels", "GA_pre_scan",                        GA_pre_scan_label_default)
        self.GA_split_label                         = self.value_assignment(self.config, "Labels", "GA_split",                           GA_split_label_default)
        self.GA_BWA_label                           = self.value_assignment(self.config, "Labels", "GA_BWA",                             GA_BWA_label_default)
        self.GA_BLAT_label                          = self.value_assignment(self.config, "Labels", "GA_BLAT",                            GA_BLAT_label_default)
        self.GA_DIAMOND_label                       = self.value_assignment(self.config, "Labels", "GA_DIAMOND",                         GA_DIAMOND_label_default)
        self.GA_final_merge_label                   = self.value_assignment(self.config, "Labels", "GA_final_merge",                     GA_final_merge_label_default)
        self.ta_label                               = self.value_assignment(self.config, "Labels", "ta",                                 taxon_annotation_label_default)
        self.ec_label                               = self.value_assignment(self.config, "Labels", "ec",                                 ec_annotation_label_default)
        self.output_label                           = self.value_assignment(self.config, "Labels", "outputs",                            output_label_default)
        
        self.qc_top_dir = os.path.join(self.out_dir, self.qc_label)
        self.qc_data_dir = os.path.join(self.qc_top_dir, "data")
        self.qc_export_dir = os.path.join(self.qc_top_dir, "final_results")
        self.qc_sort_dir = os.path.join(self.qc_data_dir, "0_sort")
        self.qc_ar_dir = os.path.join(self.qc_data_dir, "1_AR")
        self.qc_tag_dir = os.path.join(self.qc_data_dir, "2_tags")
        self.qc_merge_dir = os.path.join(self.qc_data_dir, "3_merge")
        self.qc_rem_dir = os.path.join(self.qc_data_dir, "4_qc_filter")
        self.qc_orphan_dir = os.path.join(self.qc_data_dir, "5_orphans")
        self.qc_dupes_dir = os.path.join(self.qc_data_dir, "6_dupes")

        self.host_top_dir = os.path.join(self.out_dir, self.host_label)
        self.host_data_dir = os.path.join(self.host_top_dir, "data")
        self.host_export_dir = os.path.join(self.host_top_dir, "final_results")
        self.host_bwa_dir = os.path.join(self.host_data_dir, "0_bwa_host")
        self.host_blat_dir = os.path.join(self.host_data_dir, "1_blat_host")

        self.vec_top_dir = os.path.join(self.out_dir, self.vector_label)
        self.vec_data_dir = os.path.join(self.vec_top_dir, "data")
        self.vec_export_dir = os.path.join(self.vec_top_dir, "final_results")
        self.vec_bwa_dir = os.path.join(self.vec_data_dir, "0_bwa_vec")
        self.vec_blat_dir = os.path.join(self.vec_data_dir, "1_blat_vec")

        self.rRNA_top_dir = os.path.join(self.out_dir, self.rRNA_filter_label)
        self.rRNA_data_dir = os.path.join(self.rRNA_top_dir, "data")
        self.rRNA_export_dir = os.path.join(self.rRNA_top_dir, "final_results")

        self.GA_BWA_top_dir = os.path.join(self.out_dir, self.GA_BWA_label)
        self.GA_BWA_data_dir = os.path.join(self.GA_BWA_top_dir, "data")
        self.GA_BWA_export_dir = os.path.join(self.GA_BWA_top_dir, "final_results")

        self.GA_BLAT_top_dir = os.path.join(self.out_dir, self.GA_BLAT_label)
        self.GA_BLAT_data_dir = os.path.join(self.GA_BLAT_top_dir, "data")
        self.GA_BLAT_export_dir = os.path.join(self.GA_BLAT_top_dir, "final_results")
