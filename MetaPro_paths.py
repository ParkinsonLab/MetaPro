#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

#April 12, 2021
#-------------------------------------------------------------------------
#MetaPro_paths.py
#This code manages all of the locations of files needed to run MetaPro. 
#It also controls the import of the configuration file.  
#Jan 27, 2025: moved around a bunch of things to make the paths more consistent/coherent.  New dir structure class!
import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


class time_obj:
    def compare_time(self, label):
        if(label == "qc"):
            return self.qc_end - self.qc_start
        elif(label == "host"):
            return self.host_end - self.host_start
        elif(label == "vec"):
            return self.vec_end - self.vec_start
        elif(label == "rRNA"):
            return self.rRNA_end - self.rRNA_start
        elif(label == "contigs"):


    def measure_time(self, label, start_or_end):
        if(start_or_end == "start"):
            if(label == "qc"):
                self.qc_start = time.time()
            elif(label == "host"):
                self.host_start = time.time()
            elif(label == "vec"):
                self.vec_start = time.time()
            elif(label == "rRNA"):
                self.rRNA_start = time.time()
            elif(label == "contigs"):
                self.contigs_start = time.time()
            elif(label == "repop"):
                self.repop_start = time.time()
            elif(label == "GA_BWA"):
                self.GA_BWA_start = time.time()
            elif(label == "GA_BLAT"):
                self.GA_BLAT_start = time.time()
            elif(label == "GA_DMD"):
                self.GA_DMD_start = time.time()
            elif(label == "TA"):
                self.TA_start = time.time()
            elif(label == "EC"):
                self.EC_start = time.time()
            elif(label == "out"):
                self.out_start = time.time()

        elif(start_or_end == "end"):
            if(label == "qc"):
                self.qc_end = time.time()
            elif(label == "host"):
                self.host_end = time.time()
            elif(label == "vec"):
                self.vec_end = time.time()
            elif(label == "rRNA"):
                self.rRNA_end = time.time()
            elif(label == "contigs"):
                self.contigs_end = time.time()
            elif(label == "repop"):
                self.repop_end = time.time()
            elif(label == "GA_BWA"):
                self.GA_BWA_end = time.time()
            elif(label == "GA_BLAT"):
                self.GA_BLAT_end = time.time()
            elif(label == "GA_DMD"):
                self.GA_DMD_end = time.time()
            elif(label == "TA"):
                self.TA_end = time.time()
            elif(label == "EC"):
                self.EC_end = time.time()
            elif(label == "out"):
                self.out_end = time.time()
            
            

    def __init__(self):
        self.time_now = dt.today()
        #timing vars
        self.start_time                     = time.time()
        self.end_time                       = 0
        self.qc_start                  = 0
        self.qc_end                    = 0
        
        self.host_start                     = 0
        self.host_end                       = 0

        
        self.vec_start                   = 0
        self.vec_end                     = 0

        self.rRNA_start              = 0  
        self.rRNA_end                = 0
    
        
        self.repop_start                    = 0
        self.repop_end                      = 0

        
        self.contigs_start         = 0
        self.contigs_end           = 0

        
       
        
        self.GA_BWA_start                   = 0
        self.GA_BWA_end                     = 0

        
        self.GA_BLAT_start                  = 0
        self.GA_BLAT_end                    = 0

        
        self.GA_DMD_start               = 0
        self.GA_DMD_end                 = 0

        
        self.TA_start                       = 0
        self.TA_end                         = 0
        
        self.EC_start                       = 0
        self.EC_end                         = 0

        self.out_start                  = 0
        self.out_end                    = 9

        



class dir_obj:

    def make_dirs(self, path):
        if(not os.path.exists(path)):
            os.makedirs(path)

    def value_assignment(self, config, config_section, var_name, default):
        value = ""
        #print("CONFIG:", config)
        if config:
            if(config_section in config):
                if(var_name in config[config_section]):
                    value = config[config_section][var_name]
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

        return value
    
    def get_label_dict(self):
        return self.label_dict

    def get_dir_dict(self):
        return self.dir_dict

    def __init__ (self, output_folder, config_path):

        self.out_dir = output_folder

        print("dir obj CHECKING CONFIG")
        self.config_path = config_path
        if(not os.path.isabs(self.config_path)):
            self.config_path = os.path.join(os.path.dirname(__file__), self.config_path)
        #print("full path:", self.config_path)

        if os.path.exists(self.config_path):
            config = ConfigParser() #change this to ex
            config.read(self.config_path)
            print("USING CONFIG", self.config_path)
        else:
            print("no config found, defaulting")
            config = None

        self.label_dict = dict()
        self.dir_dict = dict()
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
        contigs_label_default                           = "assemble_contigs"
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
        

        self.label_dict["qf"]                   = self.value_assignment(config, "Labels", "quality_filter",                     quality_filter_label_default)
        self.label_dict["host"]                      = self.value_assignment(config, "Labels", "host_filter",                        host_filter_label_default)
        self.label_dict["vec"]                    = self.value_assignment(config, "Labels", "vector_filter",                      vector_filter_label_default)
        self.label_dict["rRNA"]                      = self.value_assignment(config, "Labels", "rRNA_filter",                        rRNA_filter_label_default)
        self.label_dict["rRNA_split"]                = self.value_assignment(config, "Labels", "rRNA_filter_split",                  rRNA_filter_split_label_default)   
        self.label_dict["rRNA_convert"]              = self.value_assignment(config, "Labels", "rRNA_filter_convert",                rRNA_filter_convert_label_default)
        self.label_dict["rRNA_barrnap"]              = self.value_assignment(config, "Labels", "rRNA_filter_barrnap",                rRNA_filter_barrnap_label_default)
        self.label_dict["rRNA_barrnap_merge"]        = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_merge",          rRNA_filter_barrnap_merge_label_default)
        self.label_dict["rRNA_barrnap_pp"]           = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_pp",             rRNA_filter_barrnap_pp_label_default)
        self.label_dict["rRNA_infernal"]             = self.value_assignment(config, "Labels", "rRNA_filter_infernal",               rRNA_filter_infernal_label_default)
        self.label_dict["rRNA_infernal_prep"]        = self.value_assignment(config, "Labels", "rRNA_filter_infernal_prep",          rRNA_filter_infernal_prep_label_default)
        self.label_dict["rRNA_splitter"]             = self.value_assignment(config, "Labels", "rRNA_filter_splitter",               rRNA_filter_splitter_label_default)
        self.label_dict["rRNA_post"]                 = self.value_assignment(config, "Labels", "rRNA_filter_post",                   rRNA_filter_post_label_default)
        self.label_dict["repop"]                            = self.value_assignment(config, "Labels", "repop",                              repop_label_default)
        self.label_dict["contigs"]                 = self.value_assignment(config, "Labels", "contigs",                                     contigs_label_default)
        self.label_dict["destroy_contigs"]                  = self.value_assignment(config, "Labels", "destroy_contigs",                    destroy_contigs_label_default)
        self.label_dict["GA_pre_scan"]                      = self.value_assignment(config, "Labels", "GA_pre_scan",                        GA_pre_scan_label_default)
        self.label_dict["GA_split"]                         = self.value_assignment(config, "Labels", "GA_split",                           GA_split_label_default)
        self.label_dict["GA_BWA"]                           = self.value_assignment(config, "Labels", "GA_BWA",                             GA_BWA_label_default)
        self.label_dict["GA_BWA_pp"]                        = self.value_assignment(config, "Labels", "GA_BWA_pp",                          GA_BWA_pp_label_default)
        self.label_dict["GA_BWA_merge"]                     = self.value_assignment(config, "Labels", "GA_BWA_merge",                       GA_BWA_merge_label_default)
        self.label_dict["GA_BLAT"]                          = self.value_assignment(config, "Labels", "GA_BLAT",                            GA_BLAT_label_default)
        self.label_dict["GA_BLAT_cleanup"]                  = self.value_assignment(config, "Labels", "GA_BLAT_cleanup",                    GA_BLAT_cleanup_label_default)
        self.label_dict["GA_BLAT_cat"]                      = self.value_assignment(config, "Labels", "GA_BLAT_cat",                        GA_BLAT_cat_label_default)
        self.label_dict["GA_BLAT_pp"]                       = self.value_assignment(config, "Labels", "GA_BLAT_pp",                         GA_BLAT_pp_label_default)
        self.label_dict["GA_BLAT_merge"]                    = self.value_assignment(config, "Labels", "GA_BLAT_merge",                      GA_BLAT_merge_label_default)
        self.label_dict["GA_DIAMOND"]                       = self.value_assignment(config, "Labels", "GA_DIAMOND",                         GA_DIAMOND_label_default)
        self.label_dict["GA_DIAMOND_pp"]                    = self.value_assignment(config, "Labels", "GA_DIAMOND_pp",                      GA_DIAMOND_pp_label_default)
        self.label_dict["GA_final_merge"]                   = self.value_assignment(config, "Labels", "GA_final_merge",                     GA_final_merge_label_default)
        self.label_dict["TA"]                               = self.value_assignment(config, "Labels", "TA",                                 taxon_annotation_label_default)
        self.label_dict["EC"]                               = self.value_assignment(config, "Labels", "EC",                                 ec_annotation_label_default)
        self.label_dict["EC_detect"]                        = self.value_assignment(config, "Labels", "EC_detect",                          ec_annotation_detect_label_default)
        self.label_dict["EC_priam"]                         = self.value_assignment(config, "Labels", "EC_priam",                           ec_annotation_priam_label_default)
        self.label_dict["EC_priam_split"]                   = self.value_assignment(config, "Labels", "EC_priam_split",                     ec_annotation_priam_split_label_default)
        self.label_dict["EC_priam_cat"]                     = self.value_assignment(config, "Labels", "EC_priam_cat",                       ec_annotation_priam_cat_label_default)
        self.label_dict["EC_DIAMOND"]                       = self.value_assignment(config, "Labels", "EC_DIAMOND",                         ec_annotation_DIAMOND_label_default)
        self.label_dict["EC_pp"]                            = self.value_assignment(config, "Labels", "EC_pp",                              ec_annotation_pp_label_default)
        self.label_dict["out"]                           = self.value_assignment(config, "Labels", "outputs",                            output_label_default)
        self.label_dict["out_copy_gene_map"]             = self.value_assignment(config, "Labels", "output_copy_gene_map",               output_copy_gene_map_label_default)
        self.label_dict["out_clean_ec"]                  = self.value_assignment(config, "Labels", "output_clean_ec",                    output_clean_EC_label_default)
        self.label_dict["out_copy_taxa"]                 = self.value_assignment(config, "Labels", "output_copy_taxa",                   output_copy_taxa_label_default)
        self.label_dict["out_network_generation"]        = self.value_assignment(config, "Labels", "output_network_generation",          output_network_gen_label_default)
        self.label_dict["out_unique_hosts_singletons"]   = self.value_assignment(config, "Labels", "output_unique_hosts_singletons",     output_unique_hosts_singletons_label_default)
        self.label_dict["out_unique_hosts_pair_1"]       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_1",         output_unique_hosts_pair_1_label_default)
        self.label_dict["out_unique_hosts_pair_2"]       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_2",         output_unique_hosts_pair_2_label_default)
        self.label_dict["out_unique_vectors_singletons"] = self.value_assignment(config, "Labels", "output_unique_vectors_singletons",   output_unique_vectors_singletons_label_default)
        self.label_dict["out_unique_vectors_pair_1"]     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_1",       output_unique_vectors_pair_1_label_default)
        self.label_dict["out_unique_vectors_pair_2"]     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_2",       output_unique_vectors_pair_2_label_default)
        self.label_dict["out_combine_hosts"]             = self.value_assignment(config, "Labels", "output_combine_hosts",               output_combine_hosts_label_default)
        self.label_dict["out_per_read_scores"]           = self.value_assignment(config, "Labels", "output_per_read_scores",             output_per_read_scores_label_default)
        self.label_dict["out_contig_stats"]              = self.value_assignment(config, "Labels", "output_contig_stats",                output_contig_stats_label_default)
        self.label_dict["out_ec_heatmap"]                = self.value_assignment(config, "Labels", "output_ec_heatmap",                  output_ec_heatmap_label_default)
        self.label_dict["out_taxa_groupby"]              = self.value_assignment(config, "Labels", "output_taxa_groupby",                output_taxa_groupby_label_default)
        self.label_dict["out_read_count"]                = self.value_assignment(config, "Labels", "output_read_count",                  output_read_count_label_default)
        
        
        self.dir_dict["qf"] = os.path.join(self.out_dir, self.label_dict["quality_filter"])
        self.dir_dict["qf_data"] = os.path.join(self.dir_dict["gf"], "data")
        self.dir_dict["qf_sort"] = os.path.join(self.dir_dict["qf_data"], "0_ID_sort")
        self.dir_dict["qf_adapt"] = os.path.join(self.dir_dict["qf_data"], "1_adapters")
        self.dir_dict["qf_tags"] = os.path.join(self.dir_dict["qf_data"], "2_tags")
        self.dir_dict["qf_merge"] = os.path.join(self.dir_dict["qf_data"], "3_merge")
        self.dir_dict["qf_filter"] = os.path.join(self.dir_dict["qf_data"], "4_filter")
        self.dir_dict["qf_orphan"] = os.path.join(self.dir_dict["qf_data"], "5_orphan")
        self.dir_dict["qf_dup"] = os.path.join(self.dir_dict["qf_data"], "6_dup")
        self.dir_dict["qf_export"] = os.path.join(self.dir_dict["qf"], "export")

        self.dir_dict["qf_list"] = ["qf", "qf_data", "qf_sort", "qf_adapt", "qf_tags", "qf_merge", "qf_filter", "qf_orphan", "qf_dup", "qf_export"]

        self.dir_dict["host"] = os.path.join(self.out_dir, self.label_dict["host_filter"])
        self.dir_dict["host_data"] = os.path.join(self.dir_dict["host"], "data")
        self.dir_dict["host_scan"] = os.path.join(self.dir_dict["host_data"], "0_bt2_scan")
        self.dir_dict["host_export"] = os.path.join(self.dir_dict["host"], "export")

        self.dir_dict["host_list"] = ["host", "host_data", "host_scan", "host_export"]

        self.dir_dict["vec"] = os.path.join(self.out_dir, self.label_dict["vec"])
        self.dir_dict["vec_data"] = os.path.join(self.dir_dict["vec"], "data")
        self.dir_dict["vec_scan"] = os.path.join(self.dir_dict["vec_data"], "0_bt2")
        self.dir_dict["vec_export"] = os.path.join(self.dir_dict["vec"], "export")

        self.dir_dict["vec_list"] = ["vec", "vec_data", "vec_scan", "vec_export"]

        #we're keeping the split at barrnap.  Barrnap is O(n)
        self.dir_dict["rRNA"] = os.path.join(self.out_dir, self.label_dict["rRNA"])
        self.dir_dict["rRNA_data"] = os.path.join(self.dir_dict["rRNA"], "data")
        self.dir_dict["rRNA_jobs"] = os.path.join(self.dir_dict["rRNA"], "jobs")
        self.dir_dict["rRNA_split"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_split")
        self.dir_dict["rRNA_barrnap"] = os.path.join(self.dir_dict["rRNA_data"], "barrnap")
        self.dir_dict["rRNA_inf"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_inf")
        self.dir_dict["rRNA_s_fasta"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_s_fasta")
        self.dir_dict["rRNA_p1_fasta"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_p1_fasta")
        self.dir_dict["rRNA_p2_fasta"] = os.path.join(self.dir_dict["rRNA_data"], "rRNA_p2_fasta")
        self.dir_dict["rRNA_export"] = os.path.join(self.dir_dict["rRNA"], "export")
        self.dir_dict["rRNA_mRNA"] = os.path.join(self.dir_dict["rRNA_export"], "mRNA")
        self.dir_dict["rRNA_other"] = os.path.join(self.dir_dict["rRNA_export"], "other")
        self.dir_dict["rRNA_list"] = ["rRNA", "rRNA_data", "rRNA_jobs", "rRNA_split", "rRNA_barrnap", "rRNA_inf", "rRNA_s_fasta", "rRNA_p1_fasta", "rRNA_p2_fasta", "rRNA_export", "rRNA_mRNA", "rRNA_other"]

        self.dir_dict["repop"] = os.path.join(self.out_dir, self.label_dict["repop"])
        self.dir_dict["repop_data"] = os.path.join(self.dir_dict["repop"], "data")
        self.dir_dict["repop_export"] = os.path.join(self.dir_dict["repop"], "export")
        self.dir_dict["repop_list"] = ["repop", "repop_data", "repop_export"]

        self.dir_dict["contigs"] = os.path.join(self.out_dir, self.label_dict["contigs"])
        self.dir_dict["contigs_data"] = os.path.join(self.dir_dict["contigs"], "data")
        self.dir_dict["contigs_export"] = os.path.join(self.dir_dict["contigs"], "export")
        self.dir_dict["contigs_spades"] = os.path.join(self.dir_dict["data"], "0_spades")
        self.dir_dict["contigs_mgm"] = os.path.join(self.dir_dict["data"], "1_mgm")
        self.dir_dict["contigs_list"] = ["contigs", "contigs_data", "contigs_spades", "contigs_mgm", "contigs_export"]

        self.dir_dict["GA_BWA"] = os.path.join(self.out_dir, self.label_dict["GA_BWA"])
        self.dir_dict["GA_BWA_data"] = os.path.join(self.dir_dict["GA_BWA"], "data")
        self.dir_dict["GA_BWA_jobs"] = os.path.join(self.dir_dict["GA_BWA"], "jobs")
        self.dir_dict["GA_BWA_split"] = os.path.join(self.dir_dict["GA_BWA_data"], "0_split")
        self.dir_dict["GA_BWA_run"] = os.path.join(self.dir_dict["GA_BWA_data"], "1_BWA")
        self.dir_dict["GA_BWA_pp"] = os.path.join(self.dir_dict["GA_BWA_data"], "2_pp")
        self.dir_dict["GA_BWA_export"] = os.path.join(self.dir_dict["GA_BWA"], "export")
        self.dir_dict["GA_BWA_list"] = ["GA_BWA", "GA_BWA_jobs", "GA_BWA_data", "GA_BWA_split", "GA_BWA_run", "GA_BWA_pp", "GA_BWA_export"]

        self.dir_dict["GA_DMD"] = os.path.join(self.out_dir, self.label_dict["GA_DMD"])        
        self.dif_dict["GA_DMD_data"] = os.path.join(self.dir_dict["GA_DMD"], "data")
        self.dir_dict["GA_DMD_export"] = os.path.join(self.dir_dict["GA_DMD"], "export")
        self.dir_dict["GA_DMD_jobs"] = os.path.join(self.dir_dict["GA_DMD"], "jobs")
        self.dir_dict["GA_DMD_run"] = os.path.join(self.dir_dict["GA_DMD_data"], "0_dmd")
        self.dir_dict["GA_DMD_pp"] = os.path.join(self.dir_dict["GA_DMD_data"], "1_pp")
        self.dir_dict["GA_DMD_temp"] = os.path.join(self.dir_dict["GA_DMD_run"], "temp")
        self.dir_dict["GA_DMD_list"] = ["GA_DMD", "GA_DMD_data", "GA_DMD_jobs", "GA_DMD_export", "GA_DMD_run", "GA_DMD_pp", "GA_DMD_temp"]

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
        self.dir_dict["out_ng"] = os.path.join(self.dir_dict["data"], "metabolic_network")
        self.dir_dict["out_unique_hosts"] = os.path.join(self.dir_dict["data"], "unique_hosts")
        self.dir_dict["out_unique_vec"] = os.path.join(self.dir_dict["data"], "unique_vectors")
        self.dir_dict["out_heatmap"] = os.path.join(self.dir_dict["data"], "heatmap")
        
        self.dir_dict["out_list"] = ["out", "out_export", "out_data", "out_jobs", "out_ng", "out_unique_hosts", "out_unique_vec", "out_heatmap"]

        


class tool_path_obj:    
    def value_assignment(self, filetype, config, config_section, var_name, default):
        value = ""
        #print("CONFIG:", config)
        if config:
            if(config_section in config):
                if(var_name in config[config_section]):
                    value = config[config_section][var_name]
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
        
        elif(filetype == "int"):
            value = int(value)

        elif((filetype == "float") or (filetype == "flt")):
            value = float(filetype)

        return value
        
    def check_file_valid(self, file_0):
        #checks that a file is there, and larger than 0kb
        ok_flag = False
        if (os.path.exists(file_0)):
            if(os.path.getsize(file_0) > 0):
                ok_flag = True
        return ok_flag
        
    def check_dmd_valid(self):
        #if there's a .dmnd file
        #if it's sufficiently big
        dir_name = os.path.dirname(self.Prot_DB)
        basename = os.path.basename(self.Prot_DB)
        print("dir name:", dir_name)
        file_name = basename.split(".")[0]
        print("file name:", file_name)#, "extension:", extension)
        
        
        print(self.Prot_DB)
        if not(os.path.exists(self.Prot_DB)):
            sys.exit("file does not exists")
        else:
            print(dt.today(), self.Prot_DB, "exists")
        dmd_index_path = os.path.join(dir_name, file_name + ".dmnd")
        db_size = os.path.getsize(self.Prot_DB)
        
        if(os.path.exists(dmd_index_path)):
            if(os.path.getsize(dmd_index_path) >= db_size * 0.9):
                print(dt.today(), "DMD index is ok")

        
    #Mar 212, 2023: changed to be an external function because we now make custom DBs after a GA pre-scan
    def check_bwa_valid(self, DNA_DB):
        print(dt.today(), "BWA DB check:", DNA_DB)
        if(os.path.exists(DNA_DB)):
            if(os.path.isdir(DNA_DB)):
                
                file_list = os.listdir(DNA_DB)
                
                ok_flag = False
                file_count = 0
                for item in file_list:
                    
                    if(item.endswith(".fasta")):
                        #print("CHECKING:", item)
                        ext_0 = ".amb"
                        ext_1 = ".ann"
                        ext_2 = ".bwt"
                        ext_3 = ".pac"
                        ext_4 = ".sa"
                        
                        file_0 = os.path.join(DNA_DB, item + ext_0)
                        file_1 = os.path.join(DNA_DB, item + ext_1)
                        file_2 = os.path.join(DNA_DB, item + ext_2)
                        file_3 = os.path.join(DNA_DB, item + ext_3)
                        file_4 = os.path.join(DNA_DB, item + ext_4)
                        
                        ok_flag = self.check_file_valid(file_0)
                        ok_flag = self.check_file_valid(file_1)
                        ok_flag = self.check_file_valid(file_2)
                        ok_flag = self.check_file_valid(file_3)
                        ok_flag = self.check_file_valid(file_4)
                        
                        #print("OK FLAG:", ok_flag)
                        file_count += 1
                if(file_count == 0):
                    print(dt.today(), "Error: no fasta files found. BWA only accepts .fasta extensions")
                    sys.exit("empty BWA database")    
                if(not ok_flag):
                    sys.exit("BWA database has not been fully indexed. try reindexing the database")
                else:
                    print(dt.today(), "BWA database OK")
                    self.GA_DB_mode = "multi"
            else:
                if(self.DNA_DB.endswith(".fasta")):
                    ext_0 = ".amb"
                    ext_1 = ".ann"
                    ext_2 = ".bwt"
                    ext_3 = ".pac"
                    ext_4 = ".sa"
                    
                    file_0 = os.path.join(DNA_DB + ext_0)
                    file_1 = os.path.join(DNA_DB + ext_1)
                    file_2 = os.path.join(DNA_DB + ext_2)
                    file_3 = os.path.join(DNA_DB + ext_3)
                    file_4 = os.path.join(DNA_DB + ext_4)
                    
                    ok_flag = self.check_file_valid(file_0)
                    ok_flag = self.check_file_valid(file_1)
                    ok_flag = self.check_file_valid(file_2)
                    ok_flag = self.check_file_valid(file_3)
                    ok_flag = self.check_file_valid(file_4)

                    if(not ok_flag):
                        sys.exit("BWA Database FILE has not been fully indexed. Try reindexing the database")
                    else:
                        print(dt.today(), "BWA Database File OK")
                        self.GA_DB_mode = "single"
                else:
                    sys.exit("Error: no fasta file found. BWA only accepts .fasta extensions")
        else:
            exit_string = str(dt.today()) + " " + "Error: path does not exist. "+ DNA_DB
            sys.exit(exit_string)
        
            
        #if it's a fastq or fasta
        #if it's been indexed (all files present)

    def check_blat_valid(self, DNA_DB):
        #check that there's at least 1 fasta in the dict
        #but truth-be-told, this doesn't do anything, since BWA will use the same DB
        print(dt.today(), "BLAT DB check:", DNA_DB)
        if(os.path.isdir(DNA_DB)):
            file_list = os.listdir(DNA_DB)
            ok_flag = False
            file_count = 0
            for item in file_list:
                if(item.endswith(".fasta")):
                    ok_flag = self.check_file_valid(os.path.join(DNA_DB, item))
                    file_count += 1
            if(file_count == 0):
                print(dt.today(), "Error: no fasta file found.  BLAT accepts .fasta extensions only")
                sys.exit()
            if(ok_flag):
                print(dt.today(), "BLAT database OK")
                self.GA_DB_mode = "multi"
            else:
                sys.exit("Error with BLAT db. there's an empty fasta file")
        else:
            if(self.DNA_DB.endswith(".fasta")):
                ok_flag = self.check_file_valid(DNA_DB)
                if(ok_flag):
                    print(dt.today(), "BLAT Database file OK")
                    self.GA_DB_mode = "single"
                else:
                    sys.exit("Error with BLAT DB file. it's empty")
            else:
                print(dt.today(), "Error: no fasta file found.  BLAT accepts .fasta extensions only")
                sys.exit()

    def get_config_dict(self):
        return self.config_dict

    def __init__ (self, config_path):
        print("CHECKING CONFIG")
        self.config_path = config_path
        if(not os.path.isabs(self.config_path)):
            self.config_path = os.path.join(os.path.dirname(__file__), self.config_path)
        #print("full path:", self.config_path)

        if os.path.exists(self.config_path):
            config = ConfigParser() #change this to ex
            config.read(self.config_path)
            print("USING CONFIG", self.config_path)
        else:
            print("no config found, defaulting")
            config = None

        self.config_dict = dict()

        script_path             = "/pipeline/Scripts"
        tool_path               = "/pipeline_tools/"
        database_path           = self.value_assignment(config, "Databases", "database_path", "/project/j/jparkin/Lab_Databases")
        
        custom_database_path    = "/pipeline/custom_databases/"

        output_folder_default = "metapro_" + dt.today().strftime("%m%d%Y_%H%M%S")
        output_path_default = os.path.join(os.getcwd(), output_folder_default)
        self.output_path        = self.value_assignment(config, "Settings", "output_dir", output_path_default)


        
        #--------------------------------------------------
        # miscellaneous values
        BWA_mem_default = 50
        BLAT_mem_default = 10 #100MB
        DIAMOND_mem_default = 50 #60GB
        DETECT_mem_default = 50
        Infernal_mem_default = 50
        Barrnap_mem_default = 50
        repop_mem_default = 50
        TA_mem_threshold_default = 50
        BWA_pp_mem_default = 50
        BLAT_pp_mem_default = 50
        DIAMOND_pp_mem_default = 50
        GA_final_merge_mem_default = 5 
        EC_mem_threshold_default = 5
        
        cpu_default = int(math.ceil(os.cpu_count() * 0.75))
        rRNA_chunksize_default = 50000
        EC_chunksize_default = 50000
        GA_chunksize_default = 25000
        
        BWA_mem_footprint_default = 5
        BLAT_mem_footprint_default = 5
        DMD_mem_footprint_default = 10
        
        BWA_job_limit_default               = cpu_default
        BLAT_job_limit_default              = cpu_default
        DIAMOND_job_limit_default           = cpu_default
        DETECT_job_limit_default            = cpu_default
        Infernal_job_limit_default          = 1000
        Barrnap_job_limit_default           = 1000
        BWA_pp_job_limit_default            = cpu_default
        BLAT_pp_job_limit_default           = cpu_default
        DIAMOND_pp_job_limit_default        = cpu_default
        GA_final_merge_job_limit_default    = cpu_default
        repop_job_limit_default             = 1
        TA_job_limit_default                = cpu_default
        EC_job_limit_default                = cpu_default

        Barrnap_job_delay_default           = 5
        Infernal_job_delay_default          = 5
        BWA_job_delay_default               = 5
        BLAT_job_delay_default              = 5
        DIAMOND_job_delay_default           = 5
        DETECT_job_delay_default            = 5
        BWA_pp_job_delay_default            = 5
        BLAT_pp_job_delay_default           = 5
        DIAMOND_pp_job_delay_default        = 5
        GA_final_merge_job_delay_default    = 5
        repop_job_delay_default             = 10
        TA_job_delay_default                = 5
        EC_job_delay_default                = 1

        keep_all_default = "yes"
        keep_quality_default = "no"
        keep_host_default = "no"
        keep_repop_default = "no"
        keep_vector_default = "no"
        keep_assemble_contigs_default = "no"
        keep_rRNA_default = "no"
        keep_GA_BWA_default = "no"
        keep_GA_BLAT_default = "no"
        keep_GA_DIAMOND_default = "no"
        keep_GA_final_default = "no"
        keep_TA_default = "no"
        keep_EC_default = "no"
        keep_outputs_default = "no"
        
        filter_stringency_default = "high"

        
        BWA_cigar_default = 90
        BLAT_identity_default = 85
        BLAT_length_default = 0.65
        BLAT_score_default = 60
        DIAMOND_identity_default = 85
        DIAMOND_length_default = 0.65
        DIAMOND_score_default = 60
        
        #identity_cutoff= 85
        #length_cutoff= 0.65
        #score_cutoff= 60
        
        self.config_dict["target_rank"] = self.value_assignment("str", config, "Settings", "target_rank", "genus")
        self.config_dict["adapterremoval_minlength"]   = self.value_assignment("int", config, "Settings", "AdapterRemoval_minlength", 30)
        self.config_dict["show_unclassified"]          = self.value_assignment("str", config, "Settings", "Show_unclassified", "No")
        self.config_dict["bypass_log_name"]            = self.value_assignment("str", config, "Settings", "bypass_log_name", "bypass_log.txt")
        self.config_dict["debug_stop_flag"]           = self.value_assignment("str", config, "Settings", "debug_stop_flag", "none")
        self.config_dict["num_threads"]                = self.value_assignment("int", config, "Settings", "num_threads", os.cpu_count())
        if(self.config_dict["num_threads"] == 0):
            self.config_dict["num_threads"] = 1
        self.config_dict["taxa_exist_cutoff"]          = self.value_assignment("float", config, "Settings", "taxa_existence_cutoff", 0.1)
        self.config_dict["DNA_DB_mode"]                = self.value_assignment("str", config, "Settings", "DNA_DB_mode", "chocophlan") #used to indicate custom DB, or our grouped library
        #other setting is "custom"
        
        
        
        self.config_dict["RPKM_cutoff"]                = self.value_assignment("float", config, "Settings", "RPKM_cutoff", 0.01)
        self.config_dict["BWA_cigar_cutoff"]           = self.value_assignment("int", config, "Settings", "BWA_cigar_cutoff", BWA_cigar_default)
        self.config_dict["BLAT_identity_cutoff"]       = self.value_assignment("int", config, "Settings", "BLAT_identity_cutoff", BLAT_identity_default)
        self.config_dict["BLAT_length_cutoff"]         = self.value_assignment("int", config, "Settings", "BLAT_length_cutoff", BLAT_length_default)
        self.config_dict["BLAT_score_cutoff"]          = self.value_assignment("int", config, "Settings", "BLAT_score_cutoff", BLAT_score_default)            
        self.config_dict["DIAMOND_identity_cutoff"]    = self.value_assignment("int", config, "Settings", "DIAMOND_identity_cutoff", DIAMOND_identity_default)
        self.config_dict["DIAMOND_length_cutoff"]      = self.value_assignment("int", config, "Settings", "DIAMOND_length_cutoff", DIAMOND_length_default)
        self.config_dict["DIAMOND_score_cutoff"]       = self.value_assignment("int", config, "Settings", "DIAMOND_score_cutoff", DIAMOND_score_default)
        #-----------------------------------------------------------------------------------------------   

        
        self.config_dict["BWA_mem_footprint"]      = self.value_assignment("int", config, "Settings", "BWA_mem_footprint", BWA_mem_footprint_default)
        self.config_dict["BLAT_mem_footprint"]     = self.value_assignment("int", config, "Settings", "BLAT_mem_footprint", BLAT_mem_footprint_default)
        self.config_dict["DMD_mem_footprint"]      = self.value_assignment("int", config, "Settings", "DMD_mem_footprint", DMD_mem_footprint_default)
        

        #-------------------------------------------------------------------------------------------------
        self.config_dict["BWA_mem_threshold"]              = self.value_assignment("int", config, "Settings", "BWA_mem_threshold", BWA_mem_default)
        self.config_dict["BLAT_mem_threshold"]             = self.value_assignment("int", config, "Settings", "BLAT_mem_threshold", BLAT_mem_default)
        self.config_dict["DIAMOND_mem_threshold"]          = self.value_assignment("int", config, "Settings", "DIAMOND_mem_threshold", DIAMOND_mem_default)
        self.config_dict["DETECT_mem_threshold"]           = self.value_assignment("int", config, "Settings", "DETECT_mem_threshold", DETECT_mem_default)
        self.config_dict["Infernal_mem_threshold"]         = self.value_assignment("int", config, "Settings", "Infernal_mem_threshold", Infernal_mem_default)
        self.config_dict["Barrnap_mem_threshold"]          = self.value_assignment("int", config, "Settings", "Barrnap_mem_threshold", Barrnap_mem_default)
        self.config_dict["BWA_pp_mem_threshold"]           = self.value_assignment("int", config, "Settings", "BWA_pp_mem_threshold", BWA_pp_mem_default)
        self.config_dict["BLAT_pp_mem_threshold"]          = self.value_assignment("int", config, "Settings", "BLAT_pp_mem_threshold", BLAT_pp_mem_default)
        self.config_dict["DIAMOND_pp_mem_threshold"]       = self.value_assignment("int", config, "Settings", "DIAMOND_pp_mem_threshold", DIAMOND_pp_mem_default)
        self.config_dict["GA_final_merge_mem_threshold"]   = self.value_assignment("int", config, "Settings", "GA_final_merge_mem_threshold", GA_final_merge_mem_default)
        self.config_dict["TA_mem_threshold"]               = self.value_assignment("int", config, "Settings", "TA_mem_threshold", TA_mem_threshold_default)
        self.config_dict["repop_mem_threshold"]            = self.value_assignment("int", config, "Settings", "repop_mem_threshold", repop_mem_default)
        self.config_dict["EC_mem_threshold"]               = self.value_assignment("int", config, "Settings", "EC_mem_threshold", EC_mem_threshold_default)
            
        #-----------------------------------------------------------------------------------------------    
        
        self.config_dict["BWA_job_limit"]              = self.value_assignment("int", config, "Settings", "BWA_job_limit", BWA_job_limit_default)
        self.config_dict["BLAT_job_limit"]             = self.value_assignment("int", config, "Settings", "BLAT_job_limit", BLAT_job_limit_default)
        self.config_dict["DIAMOND_job_limit"]          = self.value_assignment("int", config, "Settings", "DIAMOND_job_limit", DIAMOND_job_limit_default)
        self.config_dict["DETECT_job_limit"]           = self.value_assignment("int", config, "Settings", "DETECT_job_limit", DETECT_job_limit_default)
        self.config_dict["Infernal_job_limit"]         = self.value_assignment("int", config, "Settings", "Infernal_job_limit", Infernal_job_limit_default)
        self.config_dict["Barrnap_job_limit"]          = self.value_assignment("int", config, "Settings", "Barrnap_job_limit", Barrnap_job_limit_default)
        self.config_dict["BWA_pp_job_limit"]           = self.value_assignment("int", config, "Settings", "BWA_pp_job_limit", BWA_pp_job_limit_default)
        self.config_dict["BLAT_pp_job_limit"]          = self.value_assignment("int", config, "Settings", "BLAT_pp_job_limit", BLAT_pp_job_limit_default)
        self.config_dict["DIAMOND_pp_job_limit"]       = self.value_assignment("int", config, "Settings", "DIAMOND_pp_job_limit", DIAMOND_pp_job_limit_default)
        self.config_dict["GA_final_merge_job_limit"]   = self.value_assignment("int", config, "Settings", "GA_final_merge_job_limit", GA_final_merge_job_limit_default)
        self.config_dict["TA_job_limit"]               = self.value_assignment("int", config, "Settings", "TA_job_limit", TA_job_limit_default)
        self.config_dict["repop_job_limit"]            = self.value_assignment("int", config, "Settings", "repop_job_limit", repop_job_limit_default)
        self.config_dict["EC_job_limit"]               = self.value_assignment("int", config, "Settings", "EC_job_limit", EC_job_limit_default)
        
        #------------------------------------------------------------------------
        
        self.config_dict["Infernal_job_delay"]         = self.value_assignment("int", config, "Settings", "Infernal_job_delay", Infernal_job_delay_default)
        self.config_dict["Barrnap_job_delay"]          = self.value_assignment("int", config, "Settings", "Barrnap_job_delay", Barrnap_job_delay_default)
        self.config_dict["BWA_job_delay"]              = self.value_assignment("int", config, "Settings", "BWA_job_delay", BWA_job_delay_default)
        self.config_dict["BLAT_job_delay"]             = self.value_assignment("int", config, "Settings", "BLAT_job_delay", BLAT_job_delay_default)
        self.config_dict["DIAMOND_job_delay"]          = self.value_assignment("int", config, "Settings", "DIAMOND_job_delay", DIAMOND_job_delay_default)
        self.config_dict["DETECT_job_delay"]           = self.value_assignment("int", config, "Settings", "DETECT_job_delay", DETECT_job_delay_default)
        self.config_dict["BWA_pp_job_delay"]           = self.value_assignment("int", config, "Settings", "BWA_pp_job_delay", BWA_pp_job_delay_default)
        self.config_dict["BLAT_pp_job_delay"]          = self.value_assignment("int", config, "Settings", "BLAT_pp_job_delay", BLAT_pp_job_delay_default)
        self.config_dict["DIAMOND_pp_job_delay"]       = self.value_assignment("int", config, "Settings", "DIAMOND_pp_job_delay", DIAMOND_pp_job_delay_default)
        self.config_dict["GA_final_merge_job_delay"]   = self.value_assignment("int", config, "Settings", "GA_final_merge_job_delay", GA_final_merge_job_delay_default)
        self.config_dict["TA_job_delay"]               = self.value_assignment("int", config, "Settings", "TA_job_delay", TA_job_delay_default)
        self.config_dict["repop_job_delay"]            = self.value_assignment("int", config, "Settings", "repop_job_delay", repop_job_delay_default)
        self.config_dict["EC_job_delay"]               = self.value_assignment("int", config, "Settings", "EC_job_delay", EC_job_delay_default)

        #------------------------------------------------------------------------------------------------
        self.config_dict["keep_all"]                   = self.value_assignment("str", config, "Settings", "keep_all", keep_all_default)
        self.config_dict["keep_quality"]               = self.value_assignment("str", config, "Settings", "keep_quality", keep_quality_default)
        self.config_dict["keep_host"]                  = self.value_assignment("str", config, "Settings", "keep_host", keep_host_default)
        self.config_dict["keep_vector"]                = self.value_assignment("str", config, "Settings", "keep_vector", keep_vector_default)
        self.config_dict["keep_rRNA"]                  = self.value_assignment("str", config, "Settings", "keep_rRNA", keep_rRNA_default)
        self.config_dict["keep_repop"]                 = self.value_assignment("str", config, "Settings", "keep_repop", keep_repop_default)
        self.config_dict["keep_assemble_contigs"]      = self.value_assignment("str", config, "Settings", "keep_assemble_contigs", keep_assemble_contigs_default)
        self.config_dict["keep_GA_BWA"]                = self.value_assignment("str", config, "Settings", "keep_GA_BWA", keep_GA_BWA_default)
        self.config_dict["keep_GA_BLAT"]               = self.value_assignment("str", config, "Settings", "keep_GA_BLAT", keep_GA_BLAT_default)
        self.config_dict["keep_GA_DIAMOND"]            = self.value_assignment("str", config, "Settings", "keep_GA_DIAMOND", keep_GA_DIAMOND_default)
        self.config_dict["keep_GA_final"]              = self.value_assignment("str", config, "Settings", "keep_GA_final", keep_GA_final_default)
        self.config_dict["keep_TA"]                    = self.value_assignment("str", config, "Settings", "keep_TA", keep_TA_default)
        self.config_dict["keep_EC"]                    = self.value_assignment("str", config, "Settings", "keep_EC", keep_EC_default)
        self.config_dict["keep_outputs"]               = self.value_assignment("str", config, "Settings", "keep_outputs", keep_outputs_default)
        

        
        
        #--------------------------------------------------------------------------------------
        
        self.config_dict["filter_stringency"]          = self.value_assignment("str", config, "Settings", "filter_stringency", filter_stringency_default)
        self.config_dict["GA_chunksize"]               = self.value_assignment("int", config, "Settings", "GA_chunk_size", GA_chunksize_default)
        self.config_dict["EC_chunksize"]               = self.value_assignment("int", config, "Settings", "EC_chunk_size", EC_chunksize_default)
        self.config_dict["rRNA_chunksize"]             = self.value_assignment("int", config, "Settings", "rRNA_chunk_size", rRNA_chunksize_default)
        
     
        
    
        
        
        #----------------------------------------------------------
        # Reference Databases
        # Note: default host is Mouse CDS
        
        #if config:
        self.config_dict["UniVec_Core"]        = self.value_assignment("path", config, "Databases", "UniVec_Core", os.path.join(database_path, "univec_core/UniVec_Core.fasta")) 
        self.config_dict["Adapter"]            = self.value_assignment("path", config, "Databases", "Adapter", os.path.join(database_path, "Trimmomatic_adapters/TruSeq3-PE-2.fa"))
        self.config_dict["Host"]               = self.value_assignment("path", config, "Databases", "Host",  os.path.join(database_path, "Mouse_cds/Mouse_cds.fasta"))
        self.config_dict["Rfam"]               = self.value_assignment("path", config, "Databases", "Rfam", os.path.join(database_path, "Rfam/Rfam.cm"))
        self.config_dict["DNA_DB"]             = self.value_assignment("path", config, "Databases", "DNA_DB", os.path.join(database_path, "ChocoPhlAn/ChocoPhlAn.fasta"))
        self.config_dict["source_taxa_DB"]     = self.value_assignment("path", config, "Databases", "source_taxa_db", os.path.join(database_path, "family_llbs"))
        self.config_dict["Prot_DB"]            = self.value_assignment("path", config, "Databases", "Prot_DB", os.path.join(database_path, "nr/nr"))
        self.config_dict["Prot_DB_reads"]      = self.value_assignment("path", config, "Databases", "Prot_DB_reads", os.path.join(database_path, "nr/nr"))
        self.config_dict["accession2taxid"]    = self.value_assignment("path", config, "Databases", "accession2taxid", os.path.join(database_path, "accession2taxid/accession2taxid"))
        self.config_dict["nodes"]              = self.value_assignment("path", config, "Databases", "nodes", os.path.join(database_path, "WEVOTE_db", "nodes.dmp"))
        self.config_dict["names"]              = self.value_assignment("path", config, "Databases", "names", os.path.join(database_path, "WEVOTE_db", "names.dmp"))
        self.config_dict["Kaiju_db"]           = self.value_assignment("path", config, "Databases", "Kaiju_db", os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi"))
        self.config_dict["Centrifuge_db"]      = self.value_assignment("path", config, "Databases", "Centrifuge_db", os.path.join(database_path, "centrifuge_db/nt"))
        self.config_dict["SWISS_PROT"]         = self.value_assignment("path", config, "Databases", "SWISS_PROT", os.path.join(database_path, "swiss_prot_db/swiss_prot_db"))
        self.config_dict["SWISS_PROT_map"]     = self.value_assignment("path", config, "Databases", "SWISS_PROT_map", os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv"))
        self.config_dict["PriamDB"]            = self.value_assignment("path", config, "Databases", "PriamDB", os.path.join(database_path, "PRIAM_db/"))
        self.config_dict["DetectDB"]           = self.value_assignment("path", config, "Databases", "DetectDB", os.path.join(database_path, "DETECTv2"))
        self.config_dict["WEVOTEDB"]           = self.value_assignment("path", config, "Databases", "WEVOTEDB", os.path.join(database_path, "WEVOTE_db/"))
        self.config_dict["EC_pathway"]         = self.value_assignment("path", config, "Databases", "EC_pathway", os.path.join(database_path, "EC_pathway.txt"))
        self.config_dict["path_to_superpath"]  = self.value_assignment("path", config, "Databases", "path_to_superpath", os.path.join(custom_database_path, "pathway_to_superpathway.csv"))
        self.config_dict["mgm_model"]          = self.value_assignment("path", config, "Databases", "MetaGeneMark_model", os.path.join(tool_path, "mgm/MetaGeneMark_v1.mod"))
        self.config_dict["enzyme_db"]          = self.value_assignment("path", config, "Databases", "enzyme_db", os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt"))
        self.config_dict["taxid_tree"]         = self.value_assignment("path", config, "Databases", "taxid_tree", os.path.join(custom_database_path, "taxid_trees", "family_tree.tsv"))
        self.config_dict["kraken2_db"]         = self.value_assignment("path", config, "Databases", "kraken2_db", os.path.join(custom_database_path, "kraken2_db"))

        #-------------------------------------------------------
        # test DBs
        self.GA_DB_mode = "multi" #by default for the new DB changes.
        self.check_dmd_valid()
        if(self.DNA_DB_mode == "custom"):
            self.check_bwa_valid(self.DNA_DB)
            self.check_blat_valid(self.DNA_DB)
        

        #----------------------------------------------------------
        # external tools
        
        #if config:
        self.config_dict["Python"]         = self.value_assignment("path", config, "Tools", "Python", "python3")
        self.config_dict["Java"]           = self.value_assignment("path", config, "Tools", "Java", "java -jar")
        self.config_dict["cdhit_dup"]      = self.value_assignment("path", config, "Tools", "cdhit_dup",  os.path.join(tool_path, "cdhit_dup/cd-hit-dup"))
        self.config_dict["AdapterRemoval"] = self.value_assignment("path", config, "Tools", "AdapterRemoval", os.path.join(tool_path, "adapterremoval/AdapterRemoval"))
        self.config_dict["vsearch"]        = self.value_assignment("path", config, "Tools", "vsearch", os.path.join(tool_path, "vsearch/vsearch"))
        self.config_dict["BWA"]            = self.value_assignment("path", config, "Tools", "BWA", os.path.join(tool_path, "BWA/bwa"))
        self.config_dict["SAMTOOLS"]       = self.value_assignment("path", config, "Tools", "SAMTOOLS", os.path.join(tool_path, "samtools/samtools"))
        self.config_dict["BLAT"]           = self.value_assignment("path", config, "Tools", "BLAT", os.path.join(tool_path, "PBLAT/pblat"))
        self.config_dict["DIAMOND"]        = self.value_assignment("path", config, "Tools", "DIAMOND", os.path.join(tool_path, "DIAMOND/diamond"))
        self.config_dict["Blastp"]         = self.value_assignment("path", config, "Tools", "Blastp", os.path.join(tool_path, "BLAST_p/blastp"))
        self.config_dict["Needle"]         = self.value_assignment("path", config, "Tools", "Needle", os.path.join(tool_path, "EMBOSS-6.6.0/emboss/stretcher"))
        self.config_dict["Makeblastdb"]    = self.value_assignment("path", config, "Tools", "Makeblastdb", os.path.join(tool_path, "BLAST_p/makeblastdb"))
        self.config_dict["Barrnap"]        = self.value_assignment("path", config, "Tools", "Barrnap", os.path.join(tool_path, "Barrnap/bin/barrnap"))
        self.config_dict["Infernal"]       = self.value_assignment("path", config, "Tools", "Infernal", os.path.join(tool_path, "infernal/cmsearch"))
        self.config_dict["Kaiju"]          = self.value_assignment("path", config, "Tools", "Kaiju", os.path.join(tool_path, "kaiju/kaiju"))
        self.config_dict["Centrifuge"]     = self.value_assignment("path", config, "Tools", "Centrifuge", os.path.join(tool_path, "centrifuge/centrifuge"))
        self.config_dict["Priam"]          = self.value_assignment("path", config, "Tools", "Priam", os.path.join(tool_path, "PRIAM_search/PRIAM_search.jar"))
        self.config_dict["Detect"]         = self.value_assignment("path", config, "Tools", "Detect", os.path.join(script_path, "Detect_2.2.10.py"))
        self.config_dict["BLAST_dir"]      = self.value_assignment("path", config, "Tools", "BLAST_dir", os.path.join(tool_path, "BLAST_p"))
        self.config_dict["WEVOTE"]         = self.value_assignment("path", config, "Tools", "WEVOTE", os.path.join(tool_path, "WEVOTE/WEVOTE"))
        self.config_dict["Spades"]         = self.value_assignment("path", config, "Tools", "Spades", os.path.join(tool_path, "SPAdes/bin/spades.py"))
        self.config_dict["MetaGeneMark"]   = self.value_assignment("path", config, "Tools", "MetaGeneMark", os.path.join(tool_path, "mgm/gmhmmp"))
        self.config_dict["kraken2"]        = self.value_assignment("path", config, "Tools", "kraken2", os.path.join(tool_path, "kraken2/kraken2"))
            
        #--------------------------------------------
        # Python scripts
        #if config:
        self.config_dict["sam_trimmer"]                = self.value_assignment("path", config, "code", "sam_trimmer", os.path.join(script_path, "read_sam.py"))
        self.config_dict["sort_reads"]                 = self.value_assignment("path", config, "code", "sort_reads", os.path.join(script_path, "read_sort.py"))
        self.config_dict["duplicate_repopulate"]       = self.value_assignment("path", config, "code", "duplicate_repopulation", os.path.join(script_path, "read_repopulation.py"))
        self.config_dict["orphaned_read_filter"]       = self.value_assignment("path", config, "code", "orphaned_read_filter", os.path.join(script_path, "read_orphan.py"))
        self.config_dict["remove_tag"]                 = self.value_assignment("path", config, "code", "remove_tag", os.path.join(script_path, "read_remove_tag.py"))
        self.config_dict["BLAT_Contaminant_Filter"]    = self.value_assignment("path", config, "code", "blat_contaminant_filter", os.path.join(script_path, "read_BLAT_filter_v3.py"))
        self.config_dict["File_splitter"]              = self.value_assignment("path", config, "code", "file_splitter", os.path.join(script_path, "read_split.py"))
        self.config_dict["barrnap_post"]               = self.value_assignment("path", config, "code", "barrnap_post", os.path.join(script_path, "read_rRNA_barrnap.py"))
        self.config_dict["rRNA_filter"]                = self.value_assignment("path", config, "code", "rRNA_filter", os.path.join(script_path, "read_rRNA_infernal.py"))
        self.config_dict["Map_contig"]                 = self.value_assignment("path", config, "code", "map_contig", os.path.join(script_path, "assembly_make_contig_map.py"))
        self.config_dict["flush_bad_contigs"]          = self.value_assignment("path", config, "code", "flush_bad_contigs", os.path.join(script_path, "assembly_flush_bad_contigs.py"))
        self.config_dict["contig_duplicate_remover"]   = self.value_assignment("path", config, "code", "contig_duplicate_remover", os.path.join(script_path, "assembly_deduplicate.py"))
        self.config_dict["Map_reads_gene_BWA"]         = self.value_assignment("path", config, "code", "ga_bwa_pp", os.path.join(script_path, "ga_BWA_generic_v2.py"))
        self.config_dict["Map_reads_gene_BLAT"]        = self.value_assignment("path", config, "code", "ga_blat_pp", os.path.join(script_path, "ga_BLAT_generic_v3.py"))
        self.config_dict["Map_reads_prot_DMND"]        = self.value_assignment("path", config, "code", "ga_dmd_pp", os.path.join(script_path, "ga_Diamond_generic_v2.py"))
        self.config_dict["GA_final_merge"]             = self.value_assignment("path", config, "code", "ga_final_merge", os.path.join(script_path, "ga_Final_merge_v4.py"))
        self.config_dict["GA_merge_fasta"]             = self.value_assignment("path", config, "code", "ga_merge_fasta", os.path.join(script_path, "ga_merge_fasta.py"))
        self.config_dict["GA_final_merge_fasta"]       = self.value_assignment("path", config, "code", "ga_final_merge_fasta", os.path.join(script_path, "ga_final_merge_fastq.py"))
        self.config_dict["GA_final_merge_proteins"]    = self.value_assignment("path", config, "code", "ga_final_merge_proteins", os.path.join(script_path, "ga_final_merge_proteins.py"))
        self.config_dict["GA_final_merge_maps"]        = self.value_assignment("path", config, "code", "ga_final_merge_maps", os.path.join(script_path, "ga_final_merge_map.py"))
        self.config_dict["EC_Annotation_Post"]         = self.value_assignment("path", config, "code", "ec_combine", os.path.join(script_path, "ea_combine_v5.py"))
        self.config_dict["Annotated_taxid"]            = self.value_assignment("path", config, "code", "ta_taxid", os.path.join(script_path, "ta_taxid_v3.py"))
        self.config_dict["Constrain_classification"]   = self.value_assignment("path", config, "code", "ta_constrain", os.path.join(script_path, "ta_constrain_taxonomy_v2.py"))
        self.config_dict["Classification_combine"]     = self.value_assignment("path", config, "code", "ta_combine", os.path.join(script_path, "ta_combine_v3.py"))
        self.config_dict["Wevote_parser"]              = self.value_assignment("path", config, "code", "ta_wevote", os.path.join(script_path, "ta_wevote_parser.py"))
        self.config_dict["taxa_table"]                 = self.value_assignment("path", config, "code", "output_taxa", os.path.join(script_path, "output_taxa_groupby.py"))
        self.config_dict["RPKM"]                       = self.value_assignment("path", config, "code", "output_rpkm", os.path.join(script_path, "output_table_v3.py"))
        self.config_dict["format_RPKM"]                = self.value_assignment("path", config, "code", "output_reformat", os.path.join(script_path, "output_reformat_rpkm_table.py"))
        self.config_dict["read_count"]                 = self.value_assignment("path", config, "code", "output_read_count", os.path.join(script_path, "output_read_counts_v2.py"))
        self.config_dict["read_quality_metrics"]       = self.value_assignment("path", config, "code", "output_qual", os.path.join(script_path, "output_read_quality_metrics.py"))
        self.config_dict["contig_stats"]               = self.value_assignment("path", config, "code", "output_contig_stats", os.path.join(script_path, "output_contig_stats.py"))
        self.config_dict["ec_heatmap"]                 = self.value_assignment("path", config, "code", "output_heatmap", os.path.join(script_path, "output_EC_metrics.py"))
        self.config_dict["data_change_metrics"]        = self.value_assignment("path", config, "code", "output_data_change", os.path.join(script_path, "output_data_change_metrics.py"))
        self.config_dict["get_unique_host_reads"]      = self.value_assignment("path", config, "code", "output_unique_hosts", os.path.join(script_path, "output_get_host_reads.py"))
        self.config_dict["remove_gaps_in_fasta"]       = self.value_assignment("path", config, "code", "remove_fasta_gaps", os.path.join(script_path, "remove_gaps_in_fasta.py"))
        self.config_dict["parse_sam"]                  = self.value_assignment("path", config, "code", "output_sam", os.path.join(script_path, "output_parse_sam.py"))
        self.config_dict["are_you_in_a_contig"]        = self.value_assignment("path", config, "code", "output_contig_check", os.path.join(script_path, "output_are_you_in_a_contig.py"))
        self.config_dict["convert_contig_segments"]    = self.value_assignment("path", config, "code", "output_convert_gene_map", os.path.join(script_path, "output_convert_gene_map_contig_segments.py"))
        self.config_dict["output_filter_taxa"]         = self.value_assignment("path", config, "code", "output_filter_taxa", os.path.join(script_path, "output_filter_taxa.py"))
        self.config_dict["output_filter_ECs"]          = self.value_assignment("path", config, "code", "output_filter_ec", os.path.join(script_path, "output_filter_ECs.py"))
        self.config_dict["bwa_read_sorter"]            = self.value_assignment("path", config, "code", "bwa_read_sorter", os.path.join(script_path, "bwa_read_sorter.py"))
        self.config_dict["ta_contig_name_convert"]     = self.value_assignment("path", config, "code", "ta_name_convert", os.path.join(script_path, "ta_contig_name_convert.py"))
        self.config_dict["GA_pre_scan_get_lib"]        = self.value_assignment("path", config, "code", "ga_pre_scan_get_lib", os.path.join(script_path, "ga_pre_scan_get_libs.py"))
        self.config_dict["GA_pre_scan_assemble_lib"]   = self.value_assignment("path", config, "code", "ga_pre_scan_assemble_lib", os.path.join(script_path, "ga_pre_scan_assemble_libs.py"))
        
