import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation

class mpro_dir:

    

    def make_dirs(self, path):
        if(not os.path.exists(path)):
            os.makedirs(path)
        else:
            if(not os.path.isdir(path)):
                os.makedirs(path)

    def make_dirs_from_list(self, key):
        for item in self.dir_dict[key]:
            print(dt.today(), "building dir:", item)
            self.make_dirs(self.dir_dict[item])
    
    def value_assignment(self, filetype, config, config_section, var_name, default):
        value = ""
        #print("CONFIG:", config)
        if config:
            if(config_section in config):
                if(var_name in config[config_section]):
                    value = config[config_section][var_name]
                    value = value.strip("\n")
                    if('"' in value):
                        value = value.strip('"')
                        #print("quotes cleaned")
                        #print(var_name, "found! using:", value)
                        #time.sleep(1)
                    print(var_name, "found! using:", value)
                    if(value == "0"):
                        print(var_name, "zero-setting detected: using default")
                else:
                    print(var_name, "no inner section found. using default", default)
                    value = default
            else:
                print("section:", config_section, "not found. using default for[", var_name, "]:", default)
                value = default
        else:
            print("no config, using default:", default)
            value = default


        if((filetype == "str") or (filetype == "path")):
            value = str(value)
        elif(filetype == "list"):
            
            value = value.strip(" ")
            list_str = sorted(value.split(","))
            value_list = list()
            for item in list_str:
                value_list.append(item)
                return value_list
        
        elif(filetype == "int"):
            value = int(value)

        elif((filetype == "float") or (filetype == "flt")):
            value = float(value)

        return value
        
    
    def get_label_dict(self):
        return self.label_dict

    def get_dir_dict(self):
        return self.dir_dict

    def __init__ (self, config_path, config_dict):
        
        self.config_dict = config_dict
        self.out_dir = self.config_dict["out_dir"]
        self.config_path = config_path
        if(not os.path.isabs(self.out_dir)):
            self.out_dir = os.path.abspath(self.out_dir)
            
        self.label_dict = dict()
        self.dir_dict = dict()
        self.dir_dict["main"] = self.out_dir


        if os.path.exists(self.config_path):
            config = ConfigParser() #change this to ex
            config.read(self.config_path)
            print("USING CONFIG", self.config_path)
        else:
            print("no config found, defaulting")
            config = None

        #--------------------------------------------------------------------------------------------
        # Labels.  
        # why? to change them during integration + new feature testing

        qf_label_default                    = "qf"
        host_filter_label_default                       = "host_filter"
        vector_filter_label_default                     = "vector_filter"
        rRNA_filter_label_default                       = "rRNA_filter"
        rRNA_filter_split_label_default                 = "rRNA_filter_split"
        rRNA_filter_convert_label_default               = "rRNA_filter_convert"
        rRNA_filter_bnap_label_default               = "rRNA_filter_bnap"
        rRNA_filter_bnap_merge_label_default         = "rRNA_filter_bnap_merge"
        rRNA_filter_bnap_pp_label_default            = "rRNA_filter_bnap_pp"
        rRNA_filter_infernal_label_default              = "rRNA_filter_infernal"
        rRNA_filter_infernal_prep_label_default         = "rRNA_filter_infernal_prep"
        rRNA_filter_splitter_label_default              = "rRNA_filter_splitter"
        rRNA_filter_post_label_default                  = "rRNA_filter_post"
        repop_label_default                             = "duplicate_repopulation"
        contigs_label_default                           = "assemble_contigs"
        GA_pre_scan_label_default                       = "GA_pre_scan"
        GA_split_label_default                          = "GA_split"
        GA_BT2_label_default                            = "GA_BT2"
        GA_BT2_pp_label_default                         = "GA_BT2_pp"
        GA_BT2_merge_label_default                      = "GA_BT2_merge"
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
        

        self.label_dict["qf"]                   = self.value_assignment("str", config, "Labels", "qf",                     qf_label_default)
        self.label_dict["host"]                      = self.value_assignment("list", config, "Labels", "host_filter",                        host_filter_label_default)
        self.label_dict["vec"]                    = self.value_assignment("str", config, "Labels", "vector_filter",                      vector_filter_label_default)
        self.label_dict["rRNA"]                      = self.value_assignment("str", config, "Labels", "rRNA_filter",                        rRNA_filter_label_default)
        self.label_dict["rRNA_split"]                = self.value_assignment("str", config, "Labels", "rRNA_filter_split",                  rRNA_filter_split_label_default)   
        self.label_dict["rRNA_convert"]              = self.value_assignment("str", config, "Labels", "rRNA_filter_convert",                rRNA_filter_convert_label_default)
        self.label_dict["rRNA_bnap"]              = self.value_assignment("str", config, "Labels", "rRNA_filter_bnap",                rRNA_filter_bnap_label_default)
        self.label_dict["rRNA_bnap_merge"]        = self.value_assignment("str", config, "Labels", "rRNA_filter_bnap_merge",          rRNA_filter_bnap_merge_label_default)
        self.label_dict["rRNA_bnap_pp"]           = self.value_assignment("str", config, "Labels", "rRNA_filter_bnap_pp",             rRNA_filter_bnap_pp_label_default)
        self.label_dict["rRNA_infernal"]             = self.value_assignment("str", config, "Labels", "rRNA_filter_infernal",               rRNA_filter_infernal_label_default)
        self.label_dict["rRNA_infernal_prep"]        = self.value_assignment("str", config, "Labels", "rRNA_filter_infernal_prep",          rRNA_filter_infernal_prep_label_default)
        self.label_dict["rRNA_splitter"]             = self.value_assignment("str", config, "Labels", "rRNA_filter_splitter",               rRNA_filter_splitter_label_default)
        self.label_dict["rRNA_post"]                 = self.value_assignment("str", config, "Labels", "rRNA_filter_post",                   rRNA_filter_post_label_default)
        self.label_dict["repop"]                            = self.value_assignment("str", config, "Labels", "repop",                              repop_label_default)
        self.label_dict["contigs"]                 = self.value_assignment("str", config, "Labels", "contigs",                                     contigs_label_default)
        self.label_dict["GA_pre_scan"]                      = self.value_assignment("str", config, "Labels", "GA_pre_scan",                        GA_pre_scan_label_default)
        self.label_dict["GA_split"]                         = self.value_assignment("str", config, "Labels", "GA_split",                           GA_split_label_default)
        self.label_dict["GA_BT2"]                           = self.value_assignment("str", config, "Labels", "GA_BT2",                             GA_BT2_label_default)
        self.label_dict["GA_BT2_pp"]                        = self.value_assignment("str", config, "Labels", "GA_BT2_pp",                          GA_BT2_pp_label_default)
        self.label_dict["GA_BT2_merge"]                     = self.value_assignment("str", config, "Labels", "GA_BT2_merge",                       GA_BT2_merge_label_default)
        self.label_dict["GA_BLAT"]                          = self.value_assignment("str", config, "Labels", "GA_BLAT",                            GA_BLAT_label_default)
        self.label_dict["GA_BLAT_cleanup"]                  = self.value_assignment("str", config, "Labels", "GA_BLAT_cleanup",                    GA_BLAT_cleanup_label_default)
        self.label_dict["GA_BLAT_cat"]                      = self.value_assignment("str", config, "Labels", "GA_BLAT_cat",                        GA_BLAT_cat_label_default)
        self.label_dict["GA_BLAT_pp"]                       = self.value_assignment("str", config, "Labels", "GA_BLAT_pp",                         GA_BLAT_pp_label_default)
        self.label_dict["GA_BLAT_merge"]                    = self.value_assignment("str", config, "Labels", "GA_BLAT_merge",                      GA_BLAT_merge_label_default)
        self.label_dict["GA_DMD"]                       = self.value_assignment("str", config, "Labels", "GA_DIAMOND",                         GA_DIAMOND_label_default)
        self.label_dict["GA_DMD_pp"]                    = self.value_assignment("str", config, "Labels", "GA_DIAMOND_pp",                      GA_DIAMOND_pp_label_default)
        self.label_dict["GA_final_merge"]                   = self.value_assignment("str", config, "Labels", "GA_final_merge",                     GA_final_merge_label_default)
        self.label_dict["TA"]                               = self.value_assignment("str", config, "Labels", "TA",                                 taxon_annotation_label_default)
        self.label_dict["EC"]                               = self.value_assignment("str", config, "Labels", "EC",                                 ec_annotation_label_default)
        self.label_dict["EC_detect"]                        = self.value_assignment("str", config, "Labels", "EC_detect",                          ec_annotation_detect_label_default)
        self.label_dict["EC_priam"]                         = self.value_assignment("str", config, "Labels", "EC_priam",                           ec_annotation_priam_label_default)
        self.label_dict["EC_priam_split"]                   = self.value_assignment("str", config, "Labels", "EC_priam_split",                     ec_annotation_priam_split_label_default)
        self.label_dict["EC_priam_cat"]                     = self.value_assignment("str", config, "Labels", "EC_priam_cat",                       ec_annotation_priam_cat_label_default)
        self.label_dict["EC_DIAMOND"]                       = self.value_assignment("str", config, "Labels", "EC_DIAMOND",                         ec_annotation_DIAMOND_label_default)
        self.label_dict["EC_pp"]                            = self.value_assignment("str", config, "Labels", "EC_pp",                              ec_annotation_pp_label_default)
        self.label_dict["out"]                           = self.value_assignment("str", config, "Labels", "outputs",                            output_label_default)
        self.label_dict["out_copy_gene_map"]             = self.value_assignment("str", config, "Labels", "output_copy_gene_map",               output_copy_gene_map_label_default)
        self.label_dict["out_clean_ec"]                  = self.value_assignment("str", config, "Labels", "output_clean_ec",                    output_clean_EC_label_default)
        self.label_dict["out_copy_taxa"]                 = self.value_assignment("str", config, "Labels", "output_copy_taxa",                   output_copy_taxa_label_default)
        self.label_dict["out_network_generation"]        = self.value_assignment("str", config, "Labels", "output_network_generation",          output_network_gen_label_default)
        self.label_dict["out_unique_hosts_singletons"]   = self.value_assignment("str", config, "Labels", "output_unique_hosts_singletons",     output_unique_hosts_singletons_label_default)
        self.label_dict["out_unique_hosts_pair_1"]       = self.value_assignment("str", config, "Labels", "output_unique_hosts_pair_1",         output_unique_hosts_pair_1_label_default)
        self.label_dict["out_unique_hosts_pair_2"]       = self.value_assignment("str", config, "Labels", "output_unique_hosts_pair_2",         output_unique_hosts_pair_2_label_default)
        self.label_dict["out_unique_vectors_singletons"] = self.value_assignment("str", config, "Labels", "output_unique_vectors_singletons",   output_unique_vectors_singletons_label_default)
        self.label_dict["out_unique_vectors_pair_1"]     = self.value_assignment("str", config, "Labels", "output_unique_vectors_pair_1",       output_unique_vectors_pair_1_label_default)
        self.label_dict["out_unique_vectors_pair_2"]     = self.value_assignment("str", config, "Labels", "output_unique_vectors_pair_2",       output_unique_vectors_pair_2_label_default)
        self.label_dict["out_combine_hosts"]             = self.value_assignment("str", config, "Labels", "output_combine_hosts",               output_combine_hosts_label_default)
        self.label_dict["out_per_read_scores"]           = self.value_assignment("str", config, "Labels", "output_per_read_scores",             output_per_read_scores_label_default)
        self.label_dict["out_contig_stats"]              = self.value_assignment("str", config, "Labels", "output_contig_stats",                output_contig_stats_label_default)
        self.label_dict["out_ec_heatmap"]                = self.value_assignment("str", config, "Labels", "output_ec_heatmap",                  output_ec_heatmap_label_default)
        self.label_dict["out_taxa_groupby"]              = self.value_assignment("str", config, "Labels", "output_taxa_groupby",                output_taxa_groupby_label_default)
        self.label_dict["out_read_count"]                = self.value_assignment("str", config, "Labels", "output_read_count",                  output_read_count_label_default)
        
        
        self.dir_dict["qf"] = os.path.join(self.out_dir, self.label_dict["qf"])
        self.dir_dict["qf_data"] = os.path.join(self.dir_dict["qf"], "data")
        self.dir_dict["qf_sort"] = os.path.join(self.dir_dict["qf_data"], "0_ID_sort")
        self.dir_dict["qf_adapt"] = os.path.join(self.dir_dict["qf_data"], "1_adapters")
        self.dir_dict["qf_tags"] = os.path.join(self.dir_dict["qf_data"], "2_tags")
        self.dir_dict["qf_merge"] = os.path.join(self.dir_dict["qf_data"], "3_merge")
        self.dir_dict["qf_hq"] = os.path.join(self.dir_dict["qf_data"], "4_hq")
        self.dir_dict["qf_orphan"] = os.path.join(self.dir_dict["qf_data"], "5_orphan")
        self.dir_dict["qf_dup"] = os.path.join(self.dir_dict["qf_data"], "6_dup")
        self.dir_dict["qf_export"] = os.path.join(self.dir_dict["qf"], "export")

        self.dir_dict["qf_list"] = ["qf", "qf_data", "qf_sort", "qf_adapt", "qf_tags", "qf_merge", "qf_hq", "qf_orphan", "qf_dup", "qf_export"]
        
        self.dir_dict["main_host"] = os.path.join(self.out_dir, "host")
        self.dir_dict["host_jobs"] = os.path.join(self.dir_dict["main_host"], "jobs")
        for item in self.label_dict["host"]:
            self.dir_dict[item + "_host"] = os.path.join(self.dir_dict["main_host"], item)
            self.dir_dict[item + "_host_data"] = os.path.join(self.dir_dict[item + "_host"], "data")
            self.dir_dict[item +"_host_scan"] = os.path.join(self.dir_dict[item + "_host_data"], "0_BT2_scan")
            self.dir_dict[item + "_host_export"] = os.path.join(self.dir_dict[item + "_host"], "export")

            self.dir_dict[item + "_dir_list"] = ["main_host", "host_jobs", item + "_host", item + "_host_data", item + "_host_scan", item + "_host_export"]

        self.dir_dict["vec"] = os.path.join(self.out_dir, self.label_dict["vec"])
        self.dir_dict["vec_data"] = os.path.join(self.dir_dict["vec"], "data")
        self.dir_dict["vec_scan"] = os.path.join(self.dir_dict["vec_data"], "0_BT2")
        self.dir_dict["vec_export"] = os.path.join(self.dir_dict["vec"], "export")

        self.dir_dict["vec_list"] = ["vec", "vec_data", "vec_scan", "vec_export"]

        
        self.dir_dict["rRNA"] = os.path.join(self.out_dir, self.label_dict["rRNA"])
        self.dir_dict["rRNA_data"] = os.path.join(self.dir_dict["rRNA"], "data")
        self.dir_dict["rRNA_convert"] = os.path.join(self.dir_dict["rRNA_data"], self.label_dict["rRNA_convert"])
        self.dir_dict["rRNA_jobs"] = os.path.join(self.dir_dict["rRNA"], "jobs")
        self.dir_dict["rRNA_mkrs"] = os.path.join(self.dir_dict["rRNA"], "mkrs")
        self.dir_dict["rRNA_split"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_split")
        self.dir_dict["rRNA_bnap"] = os.path.join(self.dir_dict["rRNA_data"], "bnap")
        self.dir_dict["rRNA_inf"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_inf")
        
        self.dir_dict["rRNA_export"] = os.path.join(self.dir_dict["rRNA"], "export")
        self.dir_dict["rRNA_mRNA"] = os.path.join(self.dir_dict["rRNA_export"], "mRNA")
        self.dir_dict["rRNA_other"] = os.path.join(self.dir_dict["rRNA_export"], "other")
        self.dir_dict["rRNA_list"] = ["rRNA", "rRNA_data", "rRNA_convert", "rRNA_jobs", "rRNA_mkrs", "rRNA_split", "rRNA_bnap", "rRNA_inf", "rRNA_export", "rRNA_mRNA", "rRNA_other"]

        self.dir_dict["repop"] = os.path.join(self.out_dir, self.label_dict["repop"])
        self.dir_dict["repop_data"] = os.path.join(self.dir_dict["repop"], "data")
        self.dir_dict["repop_export"] = os.path.join(self.dir_dict["repop"], "export")
        self.dir_dict["repop_list"] = ["repop", "repop_data", "repop_export"]

        self.dir_dict["contigs"] = os.path.join(self.out_dir, self.label_dict["contigs"])
        self.dir_dict["contigs_data"] = os.path.join(self.dir_dict["contigs"], "data")
        self.dir_dict["contigs_export"] = os.path.join(self.dir_dict["contigs"], "export")
        self.dir_dict["contigs_spades"] = os.path.join(self.dir_dict["contigs_data"], "0_spades")
        self.dir_dict["contigs_mgm"] = os.path.join(self.dir_dict["contigs_data"], "1_mgm")
        self.dir_dict["contigs_BT2"] = os.path.join(self.dir_dict["contigs_data"], "2_BT2")
        self.dir_dict["spades_transcripts"] = os.path.join(self.dir_dict["contigs_spades"], "transcripts.fasta")
        self.dir_dict["spades_done"] = os.path.join(self.dir_dict["contigs_spades"], "stage_7_terminate")
        
        self.dir_dict["contigs_list"] = ["contigs", "contigs_data", "contigs_spades", "contigs_mgm", "contigs_BT2", "contigs_export"]

        self.dir_dict["GA_ps"] = os.path.join(self.out_dir, self.label_dict["GA_pre_scan"])
        self.dir_dict["GA_ps_data"] = os.path.join(self.dir_dict["GA_ps"], "data")
        self.dir_dict["GA_ps_export"] = os.path.join(self.dir_dict["GA_ps"], "export")
        self.dir_dict["GA_ps_k2"] = os.path.join(self.dir_dict["GA_ps_data"], "k2")
        #self.dir_dict["GA_ps_wevote"] = os.path.join(self.dir_dict["GA_ps_data"], "wevote")
        self.dir_dict["GA_ps_list"] = ["GA_ps", "GA_ps_data", "GA_ps_export", "GA_ps_k2"]

        self.dir_dict["GA_split"] = os.path.join(self.out_dir, self.label_dict["GA_split"])

        self.dir_dict["GA_BT2"] = os.path.join(self.out_dir, self.label_dict["GA_BT2"])
        self.dir_dict["GA_BT2_data"] = os.path.join(self.dir_dict["GA_BT2"], "data")
        self.dir_dict["GA_BT2_jobs"] = os.path.join(self.dir_dict["GA_BT2"], "jobs")
        self.dir_dict["GA_BT2_mkrs"] = os.path.join(self.dir_dict["GA_BT2"], "markers")
        self.dir_dict["GA_BT2_run"] = os.path.join(self.dir_dict["GA_BT2_data"], "0_BT2")
        self.dir_dict["GA_BT2_pp"] = os.path.join(self.dir_dict["GA_BT2_data"], "1_pp")
        self.dir_dict["GA_BT2_export"] = os.path.join(self.dir_dict["GA_BT2"], "export")
        self.dir_dict["GA_BT2_list"] = ["GA_BT2", "GA_BT2_jobs", "GA_BT2_data", "GA_BT2_run", "GA_BT2_mkrs", "GA_BT2_pp", "GA_BT2_export"]

        self.dir_dict["GA_DMD"] = os.path.join(self.out_dir, self.label_dict["GA_DMD"])        
        self.dir_dict["GA_DMD_data"] = os.path.join(self.dir_dict["GA_DMD"], "data")
        
        self.dir_dict["GA_DMD_export"] = os.path.join(self.dir_dict["GA_DMD"], "export")
        self.dir_dict["GA_DMD_jobs"] = os.path.join(self.dir_dict["GA_DMD"], "jobs")
        self.dir_dict["GA_DMD_mkrs"] = os.path.join(self.dir_dict["GA_DMD"], "markers")
        self.dir_dict["GA_DMD_run"] = os.path.join(self.dir_dict["GA_DMD_data"], "0_dmd")
        self.dir_dict["GA_DMD_pp"] = os.path.join(self.dir_dict["GA_DMD_data"], "1_pp")
        self.dir_dict["GA_DMD_temp"] = os.path.join(self.dir_dict["GA_DMD_run"], "temp")
        self.dir_dict["GA_DMD_list"] = ["GA_DMD", "GA_DMD_data", "GA_DMD_jobs", "GA_DMD_mkrs", "GA_DMD_export", "GA_DMD_run", "GA_DMD_pp", "GA_DMD_temp"]

        self.dir_dict["GA_FM"] = os.path.join(self.out_dir, self.label_dict["GA_final_merge"])
        self.dir_dict["GA_FM_data"] = os.path.join(self.dir_dict["GA_FM"], "data")
        self.dir_dict["GA_FM_export"] = os.path.join(self.dir_dict["GA_FM"], "export")
        self.dir_dict["GA_FM_jobs"] = os.path.join(self.dir_dict["GA_FM"], "jobs")
        self.dir_dict["GA_FM_list"] = ["GA_FM", "GA_FM_data", "GA_FM_export", "GA_FM_jobs"]

        self.dir_dict["TA"] = os.path.join(self.out_dir, self.label_dict["TA"])
        self.dir_dict["TA_data"] = os.path.join(self.dir_dict["TA"], "data")
        self.dir_dict["TA_jobs"] = os.path.join(self.dir_dict["TA"], "jobs")
        self.dir_dict["TA_export"] = os.path.join(self.dir_dict["TA"], "export")
        self.dir_dict["TA_pull"] = os.path.join(self.dir_dict["TA_data"], "0_ga_extract")
        self.dir_dict["TA_kraken2"] = os.path.join(self.dir_dict["TA_data"], "1_kraken2")
        self.dir_dict["TA_wevote"] = os.path.join(self.dir_dict["TA_data"], "2_wevote")
        self.dir_dict["TA_list"] = ["TA", "TA_data", "TA_jobs", "TA_export", "TA_pull", "TA_kraken2", "TA_wevote"]

        self.dir_dict["EC"] = os.path.join(self.out_dir, self.label_dict["EC"])
        self.dir_dict["EC_data"] = os.path.join(self.dir_dict["EC"], "data")
        self.dir_dict["EC_jobs"] = os.path.join(self.dir_dict["EC"], "jobs")
        self.dir_dict["EC_export"] = os.path.join(self.dir_dict["EC"], "export")
        self.dir_dict["EC_detect"] = os.path.join(self.dir_dict["EC_data"], "0_detect")
        self.dir_dict["EC_priam"] = os.path.join(self.dir_dict["EC_data"], "1_priam")
        self.dir_dict["EC_DMD"] = os.path.join(self.dir_dict["EC_data"], "2_DMD")
        self.dir_dict["EC_list"] = ["EC", "EC_data", "EC_jobs", "EC_export", "EC_detect", "EC_priam", "EC_DMD"]

        self.dir_dict["out"] = os.path.join(self.out_dir, self.label_dict["out"])
        self.dir_dict["out_export"] = os.path.join(self.dir_dict["out"], "export")
        self.dir_dict["out_data"] = os.path.join(self.dir_dict["out"], "data")
        self.dir_dict["out_jobs"] = os.path.join(self.dir_dict["out"], "jobs")
        self.dir_dict["out_ng"] = os.path.join(self.dir_dict["out_data"], "metabolic_network")
        self.dir_dict["out_unique_hosts"] = os.path.join(self.dir_dict["out_data"], "unique_hosts")
        self.dir_dict["out_unique_vec"] = os.path.join(self.dir_dict["out_data"], "unique_vectors")
        self.dir_dict["out_heatmap"] = os.path.join(self.dir_dict["out_data"], "heatmap")
        
        self.dir_dict["out_list"] = ["out", "out_export", "out_data", "out_jobs", "out_ng", "out_unique_hosts", "out_unique_vec", "out_heatmap"]

