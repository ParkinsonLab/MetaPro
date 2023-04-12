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


import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation

class tool_path_obj:
    
    
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
                else:
                    print(var_name, "no inner section found. using default", default)
                    value = default
            else:
                print(config_section, "no section found, using default:", default)
                value = default
        else:
            print("no config, using default:", default)
            value = default

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

    def __init__ (self, config_path):
        print("CHECKING CONFIG")
        if config_path:
            config = ConfigParser() #change this to ex
            config.read(config_path)
            print("USING CONFIG", config_path)
        else:
            print("no config found, defaulting")
            config = None

        script_path             = "/pipeline/Scripts"
        tool_path               = "/pipeline_tools/"
        database_path           = "/project/j/jparkin/Lab_Databases/"
        custom_database_path    = "/pipeline/custom_databases/"
        

        
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
        
        self.target_rank                = self.value_assignment(config, "Settings", "target_rank", "genus")
        self.adapterremoval_minlength   = self.value_assignment(config, "Settings", "AdapterRemoval_minlength", 30)
        self.show_unclassified          = self.value_assignment(config, "Settings", "Show_unclassified", "No")
        self.bypass_log_name            = self.value_assignment(config, "Settings", "bypass_log_name", "bypass_log.txt")
        self.debug_stop_flag            = self.value_assignment(config, "Settings", "debug_stop_flag", "none")
        self.num_threads                = self.value_assignment(config, "Settings", "num_threads", os.cpu_count())
        if(self.num_threads == 0):
            self.num_threads = 1
        self.taxa_exist_cutoff          = self.value_assignment(config, "Settings", "taxa_existence_cutoff", 0.1)
        self.DNA_DB_mode                = self.value_assignment(config, "Settings", "DNA_DB_mode", "chocophlan") #used to indicate custom DB, or our grouped library
        #other setting is "custom"
        
        
        
        self.RPKM_cutoff                = self.value_assignment(config, "Settings", "RPKM_cutoff", 0.01)
        self.BWA_cigar_cutoff           = self.value_assignment(config, "Settings", "BWA_cigar_cutoff", BWA_cigar_default)
        self.BLAT_identity_cutoff       = self.value_assignment(config, "Settings", "BLAT_identity_cutoff", BLAT_identity_default)
        self.BLAT_length_cutoff         = self.value_assignment(config, "Settings", "BLAT_length_cutoff", BLAT_length_default)
        self.BLAT_score_cutoff          = self.value_assignment(config, "Settings", "BLAT_score_cutoff", BLAT_score_default)            
        self.DIAMOND_identity_cutoff    = self.value_assignment(config, "Settings", "DIAMOND_identity_cutoff", DIAMOND_identity_default)
        self.DIAMOND_length_cutoff      = self.value_assignment(config, "Settings", "DIAMOND_length_cutoff", DIAMOND_length_default)
        self.DIAMOND_score_cutoff       = self.value_assignment(config, "Settings", "DIAMOND_score_cutoff", DIAMOND_score_default)
        #-----------------------------------------------------------------------------------------------    
        
        self.BWA_mem_threshold              = self.value_assignment(config, "Settings", "BWA_mem_threshold", BWA_mem_default)
        self.BLAT_mem_threshold             = self.value_assignment(config, "Settings", "BLAT_mem_threshold", BLAT_mem_default)
        self.DIAMOND_mem_threshold          = self.value_assignment(config, "Settings", "DIAMOND_mem_threshold", DIAMOND_mem_default)
        self.DETECT_mem_threshold           = self.value_assignment(config, "Settings", "DETECT_mem_threshold", DETECT_mem_default)
        self.Infernal_mem_threshold         = self.value_assignment(config, "Settings", "Infernal_mem_threshold", Infernal_mem_default)
        self.Barrnap_mem_threshold          = self.value_assignment(config, "Settings", "Barrnap_mem_threshold", Barrnap_mem_default)
        self.BWA_pp_mem_threshold           = self.value_assignment(config, "Settings", "BWA_pp_mem_threshold", BWA_pp_mem_default)
        self.BLAT_pp_mem_threshold          = self.value_assignment(config, "Settings", "BLAT_pp_mem_threshold", BLAT_pp_mem_default)
        self.DIAMOND_pp_mem_threshold       = self.value_assignment(config, "Settings", "DIAMOND_pp_mem_threshold", DIAMOND_pp_mem_default)
        self.GA_final_merge_mem_threshold   = self.value_assignment(config, "Settings", "GA_final_merge_mem_threshold", GA_final_merge_mem_default)
        self.TA_mem_threshold               = self.value_assignment(config, "Settings", "TA_mem_threshold", TA_mem_threshold_default)
        self.repop_mem_threshold            = self.value_assignment(config, "Settings", "repop_mem_threshold", repop_mem_default)
        self.EC_mem_threshold               = self.value_assignment(config, "Settings", "EC_mem_threshold", EC_mem_threshold_default)
            
        #-----------------------------------------------------------------------------------------------    
        
        self.BWA_job_limit              = self.value_assignment(config, "Settings", "BWA_job_limit", BWA_job_limit_default)
        self.BLAT_job_limit             = self.value_assignment(config, "Settings", "BLAT_job_limit", BLAT_job_limit_default)
        self.DIAMOND_job_limit          = self.value_assignment(config, "Settings", "DIAMOND_job_limit", DIAMOND_job_limit_default)
        self.DETECT_job_limit           = self.value_assignment(config, "Settings", "DETECT_job_limit", DETECT_job_limit_default)
        self.Infernal_job_limit         = self.value_assignment(config, "Settings", "Infernal_job_limit", Infernal_job_limit_default)
        self.Barrnap_job_limit          = self.value_assignment(config, "Settings", "Barrnap_job_limit", Barrnap_job_limit_default)
        self.BWA_pp_job_limit           = self.value_assignment(config, "Settings", "BWA_pp_job_limit", BWA_pp_job_limit_default)
        self.BLAT_pp_job_limit          = self.value_assignment(config, "Settings", "BLAT_pp_job_limit", BLAT_pp_job_limit_default)
        self.DIAMOND_pp_job_limit       = self.value_assignment(config, "Settings", "DIAMOND_pp_job_limit", DIAMOND_pp_job_limit_default)
        self.GA_final_merge_job_limit   = self.value_assignment(config, "Settings", "GA_final_merge_job_limit", GA_final_merge_job_limit_default)
        self.TA_job_limit               = self.value_assignment(config, "Settings", "TA_job_limit", TA_job_limit_default)
        self.repop_job_limit            = self.value_assignment(config, "Settings", "repop_job_limit", repop_job_limit_default)
        self.EC_job_limit               = self.value_assignment(config, "Settings", "EC_job_limit", EC_job_limit_default)
        
        #------------------------------------------------------------------------
        
        self.Infernal_job_delay         = self.value_assignment(config, "Settings", "Infernal_job_delay", Infernal_job_delay_default)
        self.Barrnap_job_delay          = self.value_assignment(config, "Settings", "Barrnap_job_delay", Barrnap_job_delay_default)
        self.BWA_job_delay              = self.value_assignment(config, "Settings", "BWA_job_delay", BWA_job_delay_default)
        self.BLAT_job_delay             = self.value_assignment(config, "Settings", "BLAT_job_delay", BLAT_job_delay_default)
        self.DIAMOND_job_delay          = self.value_assignment(config, "Settings", "DIAMOND_job_delay", DIAMOND_job_delay_default)
        self.DETECT_job_delay           = self.value_assignment(config, "Settings", "DETECT_job_delay", DETECT_job_delay_default)
        self.BWA_pp_job_delay           = self.value_assignment(config, "Settings", "BWA_pp_job_delay", BWA_pp_job_delay_default)
        self.BLAT_pp_job_delay          = self.value_assignment(config, "Settings", "BLAT_pp_job_delay", BLAT_pp_job_delay_default)
        self.DIAMOND_pp_job_delay       = self.value_assignment(config, "Settings", "DIAMOND_pp_job_delay", DIAMOND_pp_job_delay_default)
        self.GA_final_merge_job_delay   = self.value_assignment(config, "Settings", "GA_final_merge_job_delay", GA_final_merge_job_delay_default)
        self.TA_job_delay               = self.value_assignment(config, "Settings", "TA_job_delay", TA_job_delay_default)
        self.repop_job_delay            = self.value_assignment(config, "Settings", "repop_job_delay", repop_job_delay_default)
        self.EC_job_delay               = self.value_assignment(config, "Settings", "EC_job_delay", EC_job_delay_default)

        #------------------------------------------------------------------------------------------------
        self.keep_all                   = self.value_assignment(config, "Settings", "keep_all", keep_all_default)
        self.keep_quality               = self.value_assignment(config, "Settings", "keep_quality", keep_quality_default)
        self.keep_host                  = self.value_assignment(config, "Settings", "keep_host", keep_host_default)
        self.keep_vector                = self.value_assignment(config, "Settings", "keep_vector", keep_vector_default)
        self.keep_rRNA                  = self.value_assignment(config, "Settings", "keep_rRNA", keep_rRNA_default)
        self.keep_repop                 = self.value_assignment(config, "Settings", "keep_repop", keep_repop_default)
        self.keep_assemble_contigs      = self.value_assignment(config, "Settings", "keep_assemble_contigs", keep_assemble_contigs_default)
        self.keep_GA_BWA                = self.value_assignment(config, "Settings", "keep_GA_BWA", keep_GA_BWA_default)
        self.keep_GA_BLAT               = self.value_assignment(config, "Settings", "keep_GA_BLAT", keep_GA_BLAT_default)
        self.keep_GA_DIAMOND            = self.value_assignment(config, "Settings", "keep_GA_DIAMOND", keep_GA_DIAMOND_default)
        self.keep_GA_final              = self.value_assignment(config, "Settings", "keep_GA_final", keep_GA_final_default)
        self.keep_TA                    = self.value_assignment(config, "Settings", "keep_TA", keep_TA_default)
        self.keep_EC                    = self.value_assignment(config, "Settings", "keep_EC", keep_EC_default)
        self.keep_outputs               = self.value_assignment(config, "Settings", "keep_outputs", keep_outputs_default)
        

        
        
        #--------------------------------------------------------------------------------------
        
        self.filter_stringency          = self.value_assignment(config, "Settings", "filter_stringency", filter_stringency_default)
        self.GA_chunksize               = self.value_assignment(config, "Settings", "GA_chunk_size", GA_chunksize_default)
        self.EC_chunksize               = self.value_assignment(config, "Settings", "EC_chunk_size", EC_chunksize_default)
        self.rRNA_chunksize             = self.value_assignment(config, "Settings", "rRNA_chunk_size", rRNA_chunksize_default)
        
     
        
    
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
        

        self.quality_filter_label                   = self.value_assignment(config, "Labels", "quality_filter",                     quality_filter_label_default)
        self.host_filter_label                      = self.value_assignment(config, "Labels", "host_filter",                        host_filter_label_default)
        self.vector_filter_label                    = self.value_assignment(config, "Labels", "vector_filter",                      vector_filter_label_default)
        self.rRNA_filter_label                      = self.value_assignment(config, "Labels", "rRNA_filter",                        rRNA_filter_label_default)
        self.rRNA_filter_split_label                = self.value_assignment(config, "Labels", "rRNA_filter_split",                  rRNA_filter_split_label_default)   
        self.rRNA_filter_convert_label              = self.value_assignment(config, "Labels", "rRNA_filter_convert",                rRNA_filter_convert_label_default)
        self.rRNA_filter_barrnap_label              = self.value_assignment(config, "Labels", "rRNA_filter_barrnap",                rRNA_filter_barrnap_label_default)
        self.rRNA_filter_barrnap_merge_label        = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_merge",          rRNA_filter_barrnap_merge_label_default)
        self.rRNA_filter_barrnap_pp_label           = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_pp",             rRNA_filter_barrnap_pp_label_default)
        self.rRNA_filter_infernal_label             = self.value_assignment(config, "Labels", "rRNA_filter_infernal",               rRNA_filter_infernal_label_default)
        self.rRNA_filter_infernal_prep_label        = self.value_assignment(config, "Labels", "rRNA_filter_infernal_prep",          rRNA_filter_infernal_prep_label_default)
        self.rRNA_filter_splitter_label             = self.value_assignment(config, "Labels", "rRNA_filter_splitter",               rRNA_filter_splitter_label_default)
        self.rRNA_filter_post_label                 = self.value_assignment(config, "Labels", "rRNA_filter_post",                   rRNA_filter_post_label_default)
        self.repop_label                            = self.value_assignment(config, "Labels", "repop",                              repop_label_default)
        self.assemble_contigs_label                 = self.value_assignment(config, "Labels", "assemble_contigs",                   assemble_contigs_label_default)
        self.destroy_contigs_label                  = self.value_assignment(config, "Labels", "destroy_contigs",                    destroy_contigs_label_default)
        self.GA_pre_scan_label                      = self.value_assignment(config, "Labels", "GA_pre_scan",                        GA_pre_scan_label_default)
        self.GA_split_label                         = self.value_assignment(config, "Labels", "GA_split",                           GA_split_label_default)
        self.GA_BWA_label                           = self.value_assignment(config, "Labels", "GA_BWA",                             GA_BWA_label_default)
        self.GA_BWA_pp_label                        = self.value_assignment(config, "Labels", "GA_BWA_pp",                          GA_BWA_pp_label_default)
        self.GA_BWA_merge_label                     = self.value_assignment(config, "Labels", "GA_BWA_merge",                       GA_BWA_merge_label_default)
        self.GA_BLAT_label                          = self.value_assignment(config, "Labels", "GA_BLAT",                            GA_BLAT_label_default)
        self.GA_BLAT_cleanup_label                  = self.value_assignment(config, "Labels", "GA_BLAT_cleanup",                    GA_BLAT_cleanup_label_default)
        self.GA_BLAT_cat_label                      = self.value_assignment(config, "Labels", "GA_BLAT_cat",                        GA_BLAT_cat_label_default)
        self.GA_BLAT_pp_label                       = self.value_assignment(config, "Labels", "GA_BLAT_pp",                         GA_BLAT_pp_label_default)
        self.GA_BLAT_merge_label                    = self.value_assignment(config, "Labels", "GA_BLAT_merge",                      GA_BLAT_merge_label_default)
        self.GA_DIAMOND_label                       = self.value_assignment(config, "Labels", "GA_DIAMOND",                         GA_DIAMOND_label_default)
        self.GA_DIAMOND_pp_label                    = self.value_assignment(config, "Labels", "GA_DIAMOND_pp",                      GA_DIAMOND_pp_label_default)
        self.GA_final_merge_label                   = self.value_assignment(config, "Labels", "GA_final_merge",                     GA_final_merge_label_default)
        self.ta_label                               = self.value_assignment(config, "Labels", "ta",                                 taxon_annotation_label_default)
        self.ec_label                               = self.value_assignment(config, "Labels", "ec",                                 ec_annotation_label_default)
        self.ec_detect_label                        = self.value_assignment(config, "Labels", "ec_detect",                          ec_annotation_detect_label_default)
        self.ec_priam_label                         = self.value_assignment(config, "Labels", "ec_priam",                           ec_annotation_priam_label_default)
        self.ec_priam_split_label                   = self.value_assignment(config, "Labels", "ec_priam_split",                     ec_annotation_priam_split_label_default)
        self.ec_priam_cat_label                     = self.value_assignment(config, "Labels", "ec_priam_cat",                       ec_annotation_priam_cat_label_default)
        self.ec_DIAMOND_label                       = self.value_assignment(config, "Labels", "ec_DIAMOND",                         ec_annotation_DIAMOND_label_default)
        self.ec_pp_label                            = self.value_assignment(config, "Labels", "ec_pp",                              ec_annotation_pp_label_default)
        self.output_label                           = self.value_assignment(config, "Labels", "outputs",                            output_label_default)
        self.output_copy_gene_map_label             = self.value_assignment(config, "Labels", "output_copy_gene_map",               output_copy_gene_map_label_default)
        self.output_clean_ec_label                  = self.value_assignment(config, "Labels", "output_clean_ec",                    output_clean_EC_label_default)
        self.output_copy_taxa_label                 = self.value_assignment(config, "Labels", "output_copy_taxa",                   output_copy_taxa_label_default)
        self.output_network_generation_label        = self.value_assignment(config, "Labels", "output_network_generation",          output_network_gen_label_default)
        self.output_unique_hosts_singletons_label   = self.value_assignment(config, "Labels", "output_unique_hosts_singletons",     output_unique_hosts_singletons_label_default)
        self.output_unique_hosts_pair_1_label       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_1",         output_unique_hosts_pair_1_label_default)
        self.output_unique_hosts_pair_2_label       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_2",         output_unique_hosts_pair_2_label_default)
        self.output_unique_vectors_singletons_label = self.value_assignment(config, "Labels", "output_unique_vectors_singletons",   output_unique_vectors_singletons_label_default)
        self.output_unique_vectors_pair_1_label     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_1",       output_unique_vectors_pair_1_label_default)
        self.output_unique_vectors_pair_2_label     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_2",       output_unique_vectors_pair_2_label_default)
        self.output_combine_hosts_label             = self.value_assignment(config, "Labels", "output_combine_hosts",               output_combine_hosts_label_default)
        self.output_per_read_scores_label           = self.value_assignment(config, "Labels", "output_per_read_scores",             output_per_read_scores_label_default)
        self.output_contig_stats_label              = self.value_assignment(config, "Labels", "output_contig_stats",                output_contig_stats_label_default)
        self.output_ec_heatmap_label                = self.value_assignment(config, "Labels", "output_ec_heatmap",                  output_ec_heatmap_label_default)
        self.output_taxa_groupby_label              = self.value_assignment(config, "Labels", "output_taxa_groupby",                output_taxa_groupby_label_default)
        self.output_read_count_label                = self.value_assignment(config, "Labels", "output_read_count",                  output_read_count_label_default)
        
        
        
        #----------------------------------------------------------
        # Reference Databases
        # Note: default host is Mouse CDS
        
        #if config:
        self.UniVec_Core        = self.value_assignment(config, "Databases", "UniVec_Core", os.path.join(database_path, "univec_core/UniVec_Core.fasta")) 
        self.Adapter            = self.value_assignment(config, "Databases", "Adapter", os.path.join(database_path, "Trimmomatic_adapters/TruSeq3-PE-2.fa"))
        self.Host               = self.value_assignment(config, "Databases", "Host",  os.path.join(database_path, "Mouse_cds/Mouse_cds.fasta"))
        self.Rfam               = self.value_assignment(config, "Databases", "Rfam", os.path.join(database_path, "Rfam/Rfam.cm"))
        self.DNA_DB             = self.value_assignment(config, "Databases", "DNA_DB", os.path.join(database_path, "ChocoPhlAn/ChocoPhlAn.fasta"))
        self.source_taxa_DB     = self.value_assignment(config, "Databases", "source_taxa_db", os.path.join(database_path, "family_llbs"))
        self.Prot_DB            = self.value_assignment(config, "Databases", "Prot_DB", os.path.join(database_path, "nr/nr"))
        self.Prot_DB_reads      = self.value_assignment(config, "Databases", "Prot_DB_reads", os.path.join(database_path, "nr/nr"))
        self.accession2taxid    = self.value_assignment(config, "Databases", "accession2taxid", os.path.join(database_path, "accession2taxid/accession2taxid"))
        self.nodes              = self.value_assignment(config, "Databases", "nodes", os.path.join(database_path, "WEVOTE_db", "nodes.dmp"))
        self.names              = self.value_assignment(config, "Databases", "names", os.path.join(database_path, "WEVOTE_db", "names.dmp"))
        self.Kaiju_db           = self.value_assignment(config, "Databases", "Kaiju_db", os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi"))
        self.Centrifuge_db      = self.value_assignment(config, "Databases", "Centrifuge_db", os.path.join(database_path, "centrifuge_db/nt"))
        self.SWISS_PROT         = self.value_assignment(config, "Databases", "SWISS_PROT", os.path.join(database_path, "swiss_prot_db/swiss_prot_db"))
        self.SWISS_PROT_map     = self.value_assignment(config, "Databases", "SWISS_PROT_map", os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv"))
        self.PriamDB            = self.value_assignment(config, "Databases", "PriamDB", os.path.join(database_path, "PRIAM_db/"))
        self.DetectDB           = self.value_assignment(config, "Databases", "DetectDB", os.path.join(database_path, "DETECTv2"))
        self.WEVOTEDB           = self.value_assignment(config, "Databases", "WEVOTEDB", os.path.join(database_path, "WEVOTE_db/"))
        self.EC_pathway         = self.value_assignment(config, "Databases", "EC_pathway", os.path.join(database_path, "EC_pathway.txt"))
        self.path_to_superpath  = self.value_assignment(config, "Databases", "path_to_superpath", os.path.join(custom_database_path, "pathway_to_superpathway.csv"))
        self.mgm_model          = self.value_assignment(config, "Databases", "MetaGeneMark_model", os.path.join(tool_path, "mgm/MetaGeneMark_v1.mod"))
        self.enzyme_db          = self.value_assignment(config, "Databases", "enzyme_db", os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt"))
        self.taxid_tree         = self.value_assignment(config, "Databases", "taxid_tree", os.path.join(custom_database_path, "taxid_trees", "family_tree.tsv"))
        self.kraken2_db         = self.value_assignment(config, "Databases", "kraken2_db", os.path.join(custom_database_path, "kraken2_db"))

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
        self.Python         = self.value_assignment(config, "Tools", "Python", "python3")
        self.Java           = self.value_assignment(config, "Tools", "Java", "java -jar")
        self.cdhit_dup      = self.value_assignment(config, "Tools", "cdhit_dup",  os.path.join(tool_path, "cdhit_dup/cd-hit-dup"))
        self.AdapterRemoval = self.value_assignment(config, "Tools", "AdapterRemoval", os.path.join(tool_path, "adapterremoval/AdapterRemoval"))
        self.vsearch        = self.value_assignment(config, "Tools", "vsearch", os.path.join(tool_path, "vsearch/vsearch"))
        self.BWA            = self.value_assignment(config, "Tools", "BWA", os.path.join(tool_path, "BWA/bwa"))
        self.SAMTOOLS       = self.value_assignment(config, "Tools", "SAMTOOLS", os.path.join(tool_path, "samtools/samtools"))
        self.BLAT           = self.value_assignment(config, "Tools", "BLAT", os.path.join(tool_path, "PBLAT/pblat"))
        self.DIAMOND        = self.value_assignment(config, "Tools", "DIAMOND", os.path.join(tool_path, "DIAMOND/diamond"))
        self.Blastp         = self.value_assignment(config, "Tools", "Blastp", os.path.join(tool_path, "BLAST_p/blastp"))
        self.Needle         = self.value_assignment(config, "Tools", "Needle", os.path.join(tool_path, "EMBOSS-6.6.0/emboss/stretcher"))
        self.Makeblastdb    = self.value_assignment(config, "Tools", "Makeblastdb", os.path.join(tool_path, "BLAST_p/makeblastdb"))
        self.Barrnap        = self.value_assignment(config, "Tools", "Barrnap", os.path.join(tool_path, "Barrnap/bin/barrnap"))
        self.Infernal       = self.value_assignment(config, "Tools", "Infernal", os.path.join(tool_path, "infernal/cmsearch"))
        self.Kaiju          = self.value_assignment(config, "Tools", "Kaiju", os.path.join(tool_path, "kaiju/kaiju"))
        self.Centrifuge     = self.value_assignment(config, "Tools", "Centrifuge", os.path.join(tool_path, "centrifuge/centrifuge"))
        self.Priam          = self.value_assignment(config, "Tools", "Priam", os.path.join(tool_path, "PRIAM_search/PRIAM_search.jar"))
        self.Detect         = self.value_assignment(config, "Tools", "Detect", os.path.join(script_path, "Detect_2.2.10.py"))
        self.BLAST_dir      = self.value_assignment(config, "Tools", "BLAST_dir", os.path.join(tool_path, "BLAST_p"))
        self.WEVOTE         = self.value_assignment(config, "Tools", "WEVOTE", os.path.join(tool_path, "WEVOTE/WEVOTE"))
        self.Spades         = self.value_assignment(config, "Tools", "Spades", os.path.join(tool_path, "SPAdes/bin/spades.py"))
        self.MetaGeneMark   = self.value_assignment(config, "Tools", "MetaGeneMark", os.path.join(tool_path, "mgm/gmhmmp"))
        self.kraken2        = self.value_assignment(config, "Tools", "kraken2", os.path.join(tool_path, "kraken2/kraken2"))
            
        #--------------------------------------------
        # Python scripts
        #if config:
        self.sam_trimmer                = self.value_assignment(config, "code", "sam_trimmer", os.path.join(script_path, "read_sam.py"))
        self.sort_reads                 = self.value_assignment(config, "code", "sort_reads", os.path.join(script_path, "read_sort.py"))
        self.duplicate_repopulate       = self.value_assignment(config, "code", "duplicate_repopulation", os.path.join(script_path, "read_repopulation.py"))
        self.orphaned_read_filter       = self.value_assignment(config, "code", "orphaned_read_filter", os.path.join(script_path, "read_orphan.py"))
        self.remove_tag                 = self.value_assignment(config, "code", "remove_tag", os.path.join(script_path, "read_remove_tag.py"))
        self.BLAT_Contaminant_Filter    = self.value_assignment(config, "code", "blat_contaminant_filter", os.path.join(script_path, "read_BLAT_filter_v3.py"))
        self.File_splitter              = self.value_assignment(config, "code", "file_splitter", os.path.join(script_path, "read_split.py"))
        self.barrnap_post               = self.value_assignment(config, "code", "barrnap_post", os.path.join(script_path, "read_rRNA_barrnap.py"))
        self.rRNA_filter                = self.value_assignment(config, "code", "rRNA_filter", os.path.join(script_path, "read_rRNA_infernal.py"))
        self.Map_contig                 = self.value_assignment(config, "code", "map_contig", os.path.join(script_path, "assembly_make_contig_map.py"))
        self.flush_bad_contigs          = self.value_assignment(config, "code", "flush_bad_contigs", os.path.join(script_path, "assembly_flush_bad_contigs.py"))
        self.contig_duplicate_remover   = self.value_assignment(config, "code", "contig_duplicate_remover", os.path.join(script_path, "assembly_deduplicate.py"))
        self.Map_reads_gene_BWA         = self.value_assignment(config, "code", "ga_bwa_pp", os.path.join(script_path, "ga_BWA_generic_v2.py"))
        self.Map_reads_gene_BLAT        = self.value_assignment(config, "code", "ga_blat_pp", os.path.join(script_path, "ga_BLAT_generic_v3.py"))
        self.Map_reads_prot_DMND        = self.value_assignment(config, "code", "ga_dmd_pp", os.path.join(script_path, "ga_Diamond_generic_v2.py"))
        self.GA_final_merge             = self.value_assignment(config, "code", "ga_final_merge", os.path.join(script_path, "ga_Final_merge_v4.py"))
        self.GA_merge_fasta             = self.value_assignment(config, "code", "ga_merge_fasta", os.path.join(script_path, "ga_merge_fasta.py"))
        self.GA_final_merge_fasta       = self.value_assignment(config, "code", "ga_final_merge_fasta", os.path.join(script_path, "ga_final_merge_fastq.py"))
        self.GA_final_merge_proteins    = self.value_assignment(config, "code", "ga_final_merge_proteins", os.path.join(script_path, "ga_final_merge_proteins.py"))
        self.GA_final_merge_maps        = self.value_assignment(config, "code", "ga_final_merge_maps", os.path.join(script_path, "ga_final_merge_map.py"))
        self.EC_Annotation_Post         = self.value_assignment(config, "code", "ec_combine", os.path.join(script_path, "ea_combine_v5.py"))
        self.Annotated_taxid            = self.value_assignment(config, "code", "ta_taxid", os.path.join(script_path, "ta_taxid_v3.py"))
        self.Constrain_classification   = self.value_assignment(config, "code", "ta_constrain", os.path.join(script_path, "ta_constrain_taxonomy_v2.py"))
        self.Classification_combine     = self.value_assignment(config, "code", "ta_combine", os.path.join(script_path, "ta_combine_v3.py"))
        self.Wevote_parser              = self.value_assignment(config, "code", "ta_wevote", os.path.join(script_path, "ta_wevote_parser.py"))
        self.taxa_table                 = self.value_assignment(config, "code", "output_taxa", os.path.join(script_path, "output_taxa_groupby.py"))
        self.RPKM                       = self.value_assignment(config, "code", "output_rpkm", os.path.join(script_path, "output_table_v3.py"))
        self.format_RPKM                = self.value_assignment(config, "code", "output_reformat", os.path.join(script_path, "output_reformat_rpkm_table.py"))
        self.read_count                 = self.value_assignment(config, "code", "output_read_count", os.path.join(script_path, "output_read_counts_v2.py"))
        self.read_quality_metrics       = self.value_assignment(config, "code", "output_qual", os.path.join(script_path, "output_read_quality_metrics.py"))
        self.contig_stats               = self.value_assignment(config, "code", "output_contig_stats", os.path.join(script_path, "output_contig_stats.py"))
        self.ec_heatmap                 = self.value_assignment(config, "code", "output_heatmap", os.path.join(script_path, "output_EC_metrics.py"))
        self.data_change_metrics        = self.value_assignment(config, "code", "output_data_change", os.path.join(script_path, "output_data_change_metrics.py"))
        self.get_unique_host_reads      = self.value_assignment(config, "code", "output_unique_hosts", os.path.join(script_path, "output_get_host_reads.py"))
        self.remove_gaps_in_fasta       = self.value_assignment(config, "code", "remove_fasta_gaps", os.path.join(script_path, "remove_gaps_in_fasta.py"))
        self.parse_sam                  = self.value_assignment(config, "code", "output_sam", os.path.join(script_path, "output_parse_sam.py"))
        self.are_you_in_a_contig        = self.value_assignment(config, "code", "output_contig_check", os.path.join(script_path, "output_are_you_in_a_contig.py"))
        self.convert_contig_segments    = self.value_assignment(config, "code", "output_convert_gene_map", os.path.join(script_path, "output_convert_gene_map_contig_segments.py"))
        self.output_filter_taxa         = self.value_assignment(config, "code", "output_filter_taxa", os.path.join(script_path, "output_filter_taxa.py"))
        self.output_filter_ECs          = self.value_assignment(config, "code", "output_filter_ec", os.path.join(script_path, "output_filter_ECs.py"))
        self.bwa_read_sorter            = self.value_assignment(config, "code", "bwa_read_sorter", os.path.join(script_path, "bwa_read_sorter.py"))
        self.ta_contig_name_convert     = self.value_assignment(config, "code", "ta_name_convert", os.path.join(script_path, "ta_contig_name_convert.py"))
        self.GA_pre_scan_get_lib        = self.value_assignment(config, "code", "ga_pre_scan_get_lib", os.path.join(script_path, "ga_pre_scan_get_libs.py"))
        self.GA_pre_scan_assemble_lib   = self.value_assignment(config, "code", "ga_pre_scan_assemble_lib", os.path.join(script_path, "ga_pre_scan_assemble_libs.py"))
        