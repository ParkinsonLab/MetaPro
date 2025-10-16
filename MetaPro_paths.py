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
import psutil as psu


        

class mpro_config:    
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
            value = value.strip("\"")
            list_str = value.split(",")
            value_list = list()
            for item in list_str:
                print("host id entry type:", type(item))
                value_list.append(item)
                print("added to list:", value_list)
                return value_list
        
        elif(filetype == "int"):
            value = int(value)

        elif((filetype == "float") or (filetype == "flt")):
            value = float(value)

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
        dir_name = os.path.dirname(self.config_dict["Prot_DB"])
        basename = os.path.basename(self.config_dict["Prot_DB"])
        print("dir name:", dir_name)
        file_name = basename.split(".")[0]
        print("file name:", file_name)#, "extension:", extension)
        
        
        print(self.config_dict["Prot_DB"])
        if not(os.path.exists(self.config_dict["Prot_DB"])):
            sys.exit("file does not exists")
        else:
            print(dt.today(), self.config_dict["Prot_DB"], "exists")
        dmd_index_path = os.path.join(dir_name, file_name + ".dmnd")
        db_size = os.path.getsize(self.config_dict["Prot_DB"])
        
        if(os.path.exists(dmd_index_path)):
            if(os.path.getsize(dmd_index_path) >= db_size * 0.9):
                print(dt.today(), "DMD index is ok")

        
    def check_BT2_valid(self, bt2_id_path, index_name):
        #just check if the bt2 indices exist. 
        bt2_loc = os.path.abspath(bt2_id_path)
        bt2_id_name = index_name

        file_list = os.listdir(bt2_loc)
        for item in file_list:
            full_path = os.path.join(bt2_loc, item)
            if item.endswith(".bt2"):
                split_name = item.split(".")
                first_name = split_name[0]
                number_part = str(split_name[1])
                if(number_part != "1"):
                    continue
                if(first_name == bt2_id_name):
                    if(os.path.getsize(full_path) > 0):
                        return True

        return False
    
    def get_config_dict(self):
        return self.config_dict
    
    

    def __init__ (self, config_path, output_path="None"):
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

        print(dt.today(), "checking onboard resources")
        self.config_dict = dict()
        
        self.config_dict["max_cpu"] = int(psu.cpu_count())
        self.config_dict["max_mem"] = int((psu.virtual_memory().total) * 0.8)

        

        script_path             = "/pipeline/Scripts"
        tool_path               = "/pipeline_tools/"
        database_path           = self.value_assignment("path", config, "Databases", "database_path", "None")
        
        custom_database_path    = "/pipeline/custom_databases/"
        if(output_path == "None"):
            output_folder_default = "metapro_" + dt.today().strftime("%m%d%Y_%H%M%S")
            output_path_default = os.path.join(os.getcwd(), output_folder_default)
            self.config_dict["out_dir"] = self.value_assignment("path", config, "Input", "out_dir",output_path_default)
            self.out_dir = self.config_dict["out_dir"]
        else:
            self.config_dict["out_dir"] = os.path.abspath(output_path)
            self.out_dir = os.path.abspath(output_path)

        if not(os.path.isabs(self.out_dir)):
            self.out_dir = os.path.abspath(self.out_dir)
            self.config_dict["out_dir"] = self.out_dir
            print(dt.today(), "output destination:", self.out_dir)
        
        #--------------------------------------------------
        # miscellaneous values
        BT2_mem_default = 50
        BLAT_mem_default = 10 #100MB
        DIAMOND_mem_default = 50 #60GB
        DETECT_mem_default = 50
        Infernal_mem_default = 50
        Barrnap_mem_default = 50
        repop_mem_default = 50
        TA_mem_threshold_default = 50
        BT2_pp_mem_default = 50
        BLAT_pp_mem_default = 50
        DIAMOND_pp_mem_default = 50
        GA_final_merge_mem_default = 5 
        EC_mem_threshold_default = 5
        
        cpu_default = int(math.ceil(os.cpu_count() * 0.75))
        rRNA_chunksize_default = 500000
        EC_chunksize_default = 50000
        GA_chunksize_default = 25000
        
        BT2_mem_footprint_default = 5
        BLAT_mem_footprint_default = 5
        DMD_mem_footprint_default = 10
        
        BT2_job_limit_default               = cpu_default
        BLAT_job_limit_default              = cpu_default
        DIAMOND_job_limit_default           = cpu_default
        DETECT_job_limit_default            = cpu_default
        Infernal_job_limit_default          = cpu_default
        Barrnap_job_limit_default           = cpu_default
        BWA_pp_job_limit_default            = cpu_default
        BLAT_pp_job_limit_default           = cpu_default
        DIAMOND_pp_job_limit_default        = cpu_default
        GA_final_merge_job_limit_default    = cpu_default
        repop_job_limit_default             = 1
        TA_job_limit_default                = cpu_default
        EC_job_limit_default                = cpu_default

        Barrnap_job_delay_default           = 5
        Infernal_job_delay_default          = 5
        BT2_job_delay_default               = 5
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
        dmd_speed_default                   = "fast"
        
        dmd_hit_count_default               = 1
        dmd_id_score_default                = 90
        dmd_query_cover_default             = 70
        dmd_min_score_default               = 80

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

        
        BT2_cigar_default = 90
        BLAT_identity_default = 85
        BLAT_length_default = 0.65
        BLAT_score_default = 60
        DIAMOND_identity_default = 85
        DIAMOND_length_default = 0.65
        DIAMOND_score_default = 60
        
        #identity_cutoff= 85
        #length_cutoff= 0.65
        #score_cutoff= 60
        
        self.config_dict["single"] = self.value_assignment("path", config, "Input", "singleton","None")
        self.config_dict["pair_1"] = self.value_assignment("path", config, "Input", "pair_1", "None")
        self.config_dict["pair_2"] = self.value_assignment("path", config, "Input", "pair_2", "None")
        self.config_dict["contig"] = self.value_assignment("path", config, "Input", "contig", "None")
        
        self.config_dict["read_mode"] = "p"
        if(self.config_dict["pair_1"] != "None"):
            self.config_dict["pair_1"] = os.path.abspath(self.config_dict["pair_1"])
            if(not os.path.exists(self.config_dict["pair_1"])):
                print(dt.today(), "pair 1 not a valid file. exiting")
                print(self.config_dict["pair_1"])
                sys.exit()
            

        if(self.config_dict["pair_2"] != "None"):
            self.config_dict["pair_2"] = os.path.abspath(self.config_dict["pair_2"])
            if(not os.path.exists(self.config_dict["pair_2"])):
                print(dt.today(), "pair 2 not a valid file. exiting")
                print(self.config_dict["pair_2"])
                sys.exit()
            else:
                print(dt.today(), "MetaPro operating in paired-mode")

        if(self.config_dict["single"] != "None"):
            self.config_dict["single"] = os.path.abspath(self.config_dict["single"])
            if(not os.path.exists(self.config_dict["single"])):
                print(dt.today(), "single not a valid file. exiting")
                print(self.config_dict["single"])
                sys.exit()
            else:
                self.config_dict["read_mode"] = "s"
                print(dt.today(), "MetaPro operating in SINGLE-end mode")
                time.sleep(10)
        if(self.config_dict["contig"] != "None"):
            self.config_dict["contig"] = os.path.abspath(self.config_dict["contig"])
            if(not os.path.exists(self.config_dict["contig"])):
                print(dt.today(), "contig not a valid file. exiting")
                print(self.config_dict["contig"])
                sys.exit()
            

        
         

        self.config_dict["bypass_log"] = os.path.abspath(os.path.join(self.config_dict["out_dir"], self.value_assignment("path", config, "Settings", "bypass_log", "bypass_log.txt")))
        self.config_dict["tutorial_keyword"] = self.value_assignment("str", config, "Settings", "tutorial_keyword", "None")
        self.config_dict["tutorial_mode"] = self.value_assignment("str", config, "Settings", "tutorial_mode", "None")
        self.config_dict["target_rank"]                 = self.value_assignment("str", config, "Settings", "target_rank", "genus")
        self.config_dict["adapterremoval_minlength"]    = self.value_assignment("str", config, "Settings", "AdapterRemoval_minlength", 30)
        self.config_dict["show_unclassified"]           = self.value_assignment("str", config, "Settings", "Show_unclassified", "No")
        #self.config_dict["bypass_log_name"]             = self.value_assignment("str", config, "Settings", "bypass_log_name", "bypass_log.txt")
        self.config_dict["debug_stop_flag"]             = self.value_assignment("str", config, "Settings", "debug_stop_flag", "none")
        self.config_dict["num_threads"]                 = self.value_assignment("int", config, "Settings", "num_threads", os.cpu_count())
        if(self.config_dict["num_threads"] == 0):
            self.config_dict["num_threads"] = 1
        self.config_dict["taxa_exist_cutoff"]           = self.value_assignment("float", config, "Settings", "taxa_existence_cutoff", 0.1)
        self.config_dict["DNA_DB_mode"]                 = self.value_assignment("str", config, "Settings", "DNA_DB_mode", "chocophlan") #used to indicate custom DB, or our grouped library
        self.config_dict["q_enc"]                       = self.value_assignment("str", config, "Settings", "fastq_encoding", "None")
        #other setting is "custom"
        
        
        
        self.config_dict["RPKM_cutoff"]                = self.value_assignment("float", config, "Settings", "RPKM_cutoff", 0.01)
        self.config_dict["BT2_cigar_cutoff"]           = self.value_assignment("int", config, "Settings", "BT2_cigar_cutoff", BT2_cigar_default)            
        self.config_dict["DMD_identity_cutoff"]    = self.value_assignment("int", config, "Settings", "DMD_identity_cutoff", DIAMOND_identity_default)
        self.config_dict["DMD_length_cutoff"]      = self.value_assignment("float", config, "Settings", "DMD_length_cutoff", DIAMOND_length_default)
        self.config_dict["DMD_score_cutoff"]       = self.value_assignment("int", config, "Settings", "DMD_score_cutoff", DIAMOND_score_default)
        #-----------------------------------------------------------------------------------------------   

        
        self.config_dict["BT2_mem_footprint"]      = self.value_assignment("int", config, "Settings", "BT2_mem_footprint", BT2_mem_footprint_default)
        self.config_dict["BLAT_mem_footprint"]     = self.value_assignment("int", config, "Settings", "BLAT_mem_footprint", BLAT_mem_footprint_default)
        self.config_dict["DMD_mem_footprint"]      = self.value_assignment("int", config, "Settings", "DMD_mem_footprint", DMD_mem_footprint_default)

        #----------------------------------------------------------------------------------------
        self.config_dict["DMD_speed"] = self.value_assignment("str", config, "Settings", "DMD_speed", dmd_speed_default)    
        self.config_dict["DMD_hit_count"] = self.value_assignment("int", config, "Settings", "DMD_hit_count", dmd_hit_count_default)
        self.config_dict["DMD_id_score"] = self.value_assignment("int", config, "Settings", "DMD_id_score", dmd_id_score_default)
        self.config_dict["DMD_query_cover"] = self.value_assignment("int", config, "Settings", "DMD_query_cover", dmd_query_cover_default)
        self.config_dict["DMD_min_score"] = self.value_assignment("int", config, "Settings", "DMD_min_score", dmd_min_score_default)

        #-------------------------------------------------------------------------------------------------
        self.config_dict["BT2_mem_threshold"]              = self.value_assignment("int", config, "Settings", "BT2_mem_threshold", BT2_mem_default)
        self.config_dict["BLAT_mem_threshold"]             = self.value_assignment("int", config, "Settings", "BLAT_mem_threshold", BLAT_mem_default)
        self.config_dict["DMD_mem_threshold"]          = self.value_assignment("int", config, "Settings", "DIAMOND_mem_threshold", DIAMOND_mem_default)
        self.config_dict["DETECT_mem_threshold"]           = self.value_assignment("int", config, "Settings", "DETECT_mem_threshold", DETECT_mem_default)
        self.config_dict["Infernal_mem_threshold"]         = self.value_assignment("int", config, "Settings", "Infernal_mem_threshold", Infernal_mem_default)
        self.config_dict["Barrnap_mem_threshold"]          = self.value_assignment("int", config, "Settings", "Barrnap_mem_threshold", Barrnap_mem_default)
        self.config_dict["BT2_pp_mem_threshold"]           = self.value_assignment("int", config, "Settings", "BT2_pp_mem_threshold", BT2_pp_mem_default)
        self.config_dict["BLAT_pp_mem_threshold"]          = self.value_assignment("int", config, "Settings", "BLAT_pp_mem_threshold", BLAT_pp_mem_default)
        self.config_dict["DIAMOND_pp_mem_threshold"]       = self.value_assignment("int", config, "Settings", "DIAMOND_pp_mem_threshold", DIAMOND_pp_mem_default)
        self.config_dict["GA_final_merge_mem_threshold"]   = self.value_assignment("int", config, "Settings", "GA_final_merge_mem_threshold", GA_final_merge_mem_default)
        self.config_dict["TA_mem_threshold"]               = self.value_assignment("int", config, "Settings", "TA_mem_threshold", TA_mem_threshold_default)
        self.config_dict["repop_mem_threshold"]            = self.value_assignment("int", config, "Settings", "repop_mem_threshold", repop_mem_default)
        self.config_dict["EC_mem_threshold"]               = self.value_assignment("int", config, "Settings", "EC_mem_threshold", EC_mem_threshold_default)
            
        #-----------------------------------------------------------------------------------------------    
        
        self.config_dict["BT2_job_limit"]              = self.value_assignment("int", config, "Settings", "BT2_job_limit", BT2_job_limit_default)
        self.config_dict["BLAT_job_limit"]             = self.value_assignment("int", config, "Settings", "BLAT_job_limit", BLAT_job_limit_default)
        self.config_dict["DMD_job_limit"]          = self.value_assignment("int", config, "Settings", "DIAMOND_job_limit", DIAMOND_job_limit_default)
        self.config_dict["DETECT_job_limit"]           = self.value_assignment("int", config, "Settings", "DETECT_job_limit", DETECT_job_limit_default)
        self.config_dict["Infernal_job_limit"]         = self.value_assignment("int", config, "Settings", "Infernal_job_limit", Infernal_job_limit_default)
        self.config_dict["Barrnap_job_limit"]          = self.value_assignment("int", config, "Settings", "Barrnap_job_limit", Barrnap_job_limit_default)
        self.config_dict["BWA_pp_job_limit"]           = self.value_assignment("int", config, "Settings", "BWA_pp_job_limit", BWA_pp_job_limit_default)
        self.config_dict["BLAT_pp_job_limit"]          = self.value_assignment("int", config, "Settings", "BLAT_pp_job_limit", BLAT_pp_job_limit_default)
        self.config_dict["DMD_pp_job_limit"]       = self.value_assignment("int", config, "Settings", "DIAMOND_pp_job_limit", DIAMOND_pp_job_limit_default)
        self.config_dict["GA_final_merge_job_limit"]   = self.value_assignment("int", config, "Settings", "GA_final_merge_job_limit", GA_final_merge_job_limit_default)
        self.config_dict["TA_job_limit"]               = self.value_assignment("int", config, "Settings", "TA_job_limit", TA_job_limit_default)
        self.config_dict["repop_job_limit"]            = self.value_assignment("int", config, "Settings", "repop_job_limit", repop_job_limit_default)
        self.config_dict["EC_job_limit"]               = self.value_assignment("int", config, "Settings", "EC_job_limit", EC_job_limit_default)
        
        #------------------------------------------------------------------------
        
        self.config_dict["Infernal_job_delay"]         = self.value_assignment("float", config, "Settings", "Infernal_job_delay", Infernal_job_delay_default)
        self.config_dict["Barrnap_job_delay"]          = self.value_assignment("float", config, "Settings", "Barrnap_job_delay", Barrnap_job_delay_default)
        self.config_dict["BT2_job_delay"]              = self.value_assignment("float", config, "Settings", "BT2_job_delay", BT2_job_delay_default)
        self.config_dict["BLAT_job_delay"]             = self.value_assignment("float", config, "Settings", "BLAT_job_delay", BLAT_job_delay_default)
        self.config_dict["DMD_job_delay"]          = self.value_assignment("float", config, "Settings", "DIAMOND_job_delay", DIAMOND_job_delay_default)
        self.config_dict["DETECT_job_delay"]           = self.value_assignment("float", config, "Settings", "DETECT_job_delay", DETECT_job_delay_default)
        self.config_dict["BWA_pp_job_delay"]           = self.value_assignment("float", config, "Settings", "BWA_pp_job_delay", BWA_pp_job_delay_default)
        self.config_dict["BLAT_pp_job_delay"]          = self.value_assignment("float", config, "Settings", "BLAT_pp_job_delay", BLAT_pp_job_delay_default)
        self.config_dict["DIAMOND_pp_job_delay"]       = self.value_assignment("float", config, "Settings", "DIAMOND_pp_job_delay", DIAMOND_pp_job_delay_default)
        self.config_dict["GA_final_merge_job_delay"]   = self.value_assignment("float", config, "Settings", "GA_final_merge_job_delay", GA_final_merge_job_delay_default)
        self.config_dict["TA_job_delay"]               = self.value_assignment("float", config, "Settings", "TA_job_delay", TA_job_delay_default)
        self.config_dict["repop_job_delay"]            = self.value_assignment("float", config, "Settings", "repop_job_delay", repop_job_delay_default)
        self.config_dict["EC_job_delay"]               = self.value_assignment("float", config, "Settings", "EC_job_delay", EC_job_delay_default)

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

        #--------------------------------------------------------------------------------------
        self.config_dict["deepec_cpu_count"]            = self.value_assignment("int", config, "Settings", "deepec_cpu_count", int(os.cpu_count()))
        self.config_dict["deepec_batch_size"]           = self.value_assignment("int", config, "Settings", "deepec_batch_size", 512)
        
     
        
    
        
        
        #----------------------------------------------------------
        # Reference Databases
        # Note: default host is Mouse CDS
        
        #if config:
        self.config_dict["vector_db"]             = self.value_assignment("path", config, "Databases", "vector_db", os.path.join(database_path, "univec"))
        self.config_dict["vector_ID"]          = self.value_assignment("str", config, "Databases", "vector_IDs","univec") 
        self.config_dict["Adapter"]             = self.value_assignment("path", config, "Databases", "Adapter", os.path.join(database_path, "Trimmomatic_adapters/TruSeq3-PE-2.fa"))
        self.config_dict["Host_db"]             = self.value_assignment("path", config, "Databases", "Host_db",  os.path.join(database_path, "Mouse_cds"))
        self.config_dict["Host_IDs"]            = self.value_assignment("list",  config, "Databases", "Host_IDs", "none")
        self.config_dict["Rfam"]                = self.value_assignment("path", config, "Databases", "Rfam", os.path.join(database_path, "Rfam/Rfam.cm"))
        self.config_dict["DNA_DB"]              = self.value_assignment("path", config, "Databases", "DNA_DB", os.path.join(database_path, "ChocoPhlAn/ChocoPhlAn.fasta"))
        self.config_dict["source_taxa_DB"]      = self.value_assignment("path", config, "Databases", "source_taxa_db", os.path.join(database_path, "family_llbs"))
        self.config_dict["Prot_DB"]             = self.value_assignment("path", config, "Databases", "Prot_DB", os.path.join(database_path, "nr/nr.fasta"))
        self.config_dict["Prot_DB_reads"]       = self.value_assignment("path", config, "Databases", "Prot_DB_reads", os.path.join(database_path, "nr/nr.fasta"))
        self.config_dict["accession2taxid"]     = self.value_assignment("path", config, "Databases", "accession2taxid", os.path.join(database_path, "accession2taxid/accession2taxid"))
        self.config_dict["nodes"]               = self.value_assignment("path", config, "Databases", "nodes", os.path.join(database_path, "WEVOTE_db", "nodes.dmp"))
        self.config_dict["names"]               = self.value_assignment("path", config, "Databases", "names", os.path.join(database_path, "WEVOTE_db", "names.dmp"))
        self.config_dict["Kaiju_db"]            = self.value_assignment("path", config, "Databases", "Kaiju_db", os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi"))
        self.config_dict["Centrifuge_db"]       = self.value_assignment("path", config, "Databases", "Centrifuge_db", os.path.join(database_path, "centrifuge_db/nt"))
        self.config_dict["SWISS_PROT"]          = self.value_assignment("path", config, "Databases", "SWISS_PROT", os.path.join(database_path, "swiss_prot_db/swiss_prot_db"))
        self.config_dict["SWISS_PROT_map"]      = self.value_assignment("path", config, "Databases", "SWISS_PROT_map", os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv"))
        self.config_dict["PriamDB"]             = self.value_assignment("path", config, "Databases", "PriamDB", os.path.join(database_path, "PRIAM_db/"))
        self.config_dict["DetectDB"]            = self.value_assignment("path", config, "Databases", "DetectDB", os.path.join(database_path, "DETECTv2"))
        self.config_dict["WEVOTEDB"]            = self.value_assignment("path", config, "Databases", "WEVOTEDB", os.path.join(database_path, "WEVOTE_db/"))
        self.config_dict["EC_pathway"]          = self.value_assignment("path", config, "Databases", "EC_pathway", os.path.join(database_path, "EC_pathway.txt"))
        self.config_dict["path_to_superpath"]   = self.value_assignment("path", config, "Databases", "path_to_superpath", os.path.join(custom_database_path, "pathway_to_superpathway.csv"))
        #self.config_dict["mgm_model"]           = self.value_assignment("path", config, "Databases", "MetaGeneMark_model", os.path.join(tool_path, "mgm/MetaGeneMark_v1.mod"))
        self.config_dict["mgm_model"]           = self.value_assignment("path", config, "Databases", "MetaGeneMark_model", os.path.join(tool_path, "mgm2/src/mgm2_11.mod"))
        self.config_dict["enzyme_db"]           = self.value_assignment("path", config, "Databases", "enzyme_db", os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt"))
        self.config_dict["taxid_tree"]          = self.value_assignment("path", config, "Databases", "taxid_tree", os.path.join(custom_database_path, "taxid_trees", "family_tree.tsv"))
        self.config_dict["kraken2_db"]          = self.value_assignment("path", config, "Databases", "kraken2_db", os.path.join(custom_database_path, "kraken2_db"))

        #-------------------------------------------------------
        # test DBs
        self.config_dict["custom_ga_lib_list"] = self.value_assignment("path", "config", "Databases", "custom_ga_lib_list", "none")
        self.check_dmd_valid()
        
        
        

        #----------------------------------------------------------
        # external tools
        
        #if config:
        self.config_dict["Python"]         = self.value_assignment("path", config, "Tools", "Python", "python3")
        self.config_dict["Java"]           = self.value_assignment("path", config, "Tools", "Java", "java -jar")
        self.config_dict["cdhit_dup"]      = self.value_assignment("path", config, "Tools", "cdhit_dup",  "cd-hit-dup")
        self.config_dict["AdapterRemoval"] = self.value_assignment("path", config, "Tools", "AdapterRemoval", os.path.join(tool_path, "adapterremoval/AdapterRemoval"))
        self.config_dict["vsearch"]        = self.value_assignment("path", config, "Tools", "vsearch", "vsearch")
        self.config_dict["BT2"]             = self.value_assignment("path", config, "Tools", "BT2", "bowtie2")
        self.config_dict["BT2_index"]       = self.value_assignment("path", config, "Tools", "BT2_index", "bowtie2-build")
        self.config_dict["samtools"]       = self.value_assignment("path", config, "Tools", "SAMTOOLS", "samtools")
        self.config_dict["DMD"]        = self.value_assignment("path", config, "Tools", "DIAMOND", "diamond")
        self.config_dict["Blastp"]         = self.value_assignment("path", config, "Tools", "Blastp", "blastp")
        self.config_dict["Makeblastdb"]    = self.value_assignment("path", config, "Tools", "Makeblastdb", "makeblastdb")
        self.config_dict["Barrnap"]        = self.value_assignment("path", config, "Tools", "Barrnap", "barrnap")
        self.config_dict["Infernal"]       = self.value_assignment("path", config, "Tools", "Infernal", "cmsearch")
        self.config_dict["BLAST_dir"]      = self.value_assignment("path", config, "Tools", "BLAST_dir", os.path.join(tool_path, "BLAST_p"))
        self.config_dict["Spades"]         = self.value_assignment("path", config, "Tools", "Spades", os.path.join(tool_path, "SPAdes-4.2.0-Linux/bin/spades.py"))
        self.config_dict["mgm2"]   = self.value_assignment("path", config, "Tools", "MetaGeneMark2", "gmhmmp2")
        self.config_dict["mgm1"]   = self.value_assignment("path", config, "Tools", "MetaGeneMark1", "gmhmmp")
        self.config_dict["kraken2"]        = self.value_assignment("path", config, "Tools", "kraken2", "kraken2")
        self.config_dict["DeepEC"]      = self.value_assignment("path", config, "Tools", "DeepEC", "run_deepectransformer.py")

        #--------------------------------------------
        # Python scripts
        #if config:
        self.config_dict["sam_trimmer"]                = self.value_assignment("path", config, "code", "sam_trimmer", os.path.join(script_path, "read_sam.py"))
        self.config_dict["sort_reads"]                 = self.value_assignment("path", config, "code", "sort_reads", os.path.join(script_path, "read_sort.py"))
        self.config_dict["duplicate_repopulate"]       = self.value_assignment("path", config, "code", "repop", os.path.join(script_path, "read_repopulation.py"))
        self.config_dict["orphaned_read_filter"]       = self.value_assignment("path", config, "code", "orphaned_read_filter", os.path.join(script_path, "read_orphan.py"))
        self.config_dict["remove_tag"]                 = self.value_assignment("path", config, "code", "remove_tag", os.path.join(script_path, "read_remove_tag.py"))
        self.config_dict["BLAT_Contaminant_Filter"]    = self.value_assignment("path", config, "code", "blat_contaminant_filter", os.path.join(script_path, "read_BLAT_filter_v3.py"))
        self.config_dict["File_splitter"]              = self.value_assignment("path", config, "code", "file_splitter", os.path.join(script_path, "read_split.py"))
        self.config_dict["barrnap_post"]               = self.value_assignment("path", config, "code", "bnap_post", os.path.join(script_path, "read_rRNA_barrnap.py"))
        self.config_dict["rRNA_inf_pp"]                = self.value_assignment("path", config, "code", "inf_post", os.path.join(script_path, "read_rRNA_infernal.py"))
        self.config_dict["Map_contig"]                 = self.value_assignment("path", config, "code", "map_contig", os.path.join(script_path, "assembly_make_contig_map.py"))
        self.config_dict["flush_bad_contigs"]          = self.value_assignment("path", config, "code", "flush_bad_contigs", os.path.join(script_path, "assembly_flush_bad_contigs.py"))
        self.config_dict["contig_duplicate_remover"]   = self.value_assignment("path", config, "code", "contig_duplicate_remover", os.path.join(script_path, "assembly_deduplicate.py"))
        self.config_dict["GA_BT2_pp"]                   = self.value_assignment("path", config, "code", "ga_bt2_pp", os.path.join(script_path, "GA_samfile.py"))
        self.config_dict["GA_dmd_pp"]        = self.value_assignment("path", config, "code", "ga_dmd_pp", os.path.join(script_path, "ga_dmd_pp.py"))
        self.config_dict["GA_final_merge"]             = self.value_assignment("path", config, "code", "ga_final_merge", os.path.join(script_path, "ga_Final_merge_v4.py"))
        self.config_dict["GA_merge_fasta"]             = self.value_assignment("path", config, "code", "ga_merge_fasta", os.path.join(script_path, "ga_merge_fasta.py"))
        self.config_dict["GA_final_merge_fasta"]       = self.value_assignment("path", config, "code", "ga_final_merge_fasta", os.path.join(script_path, "ga_final_merge_fasta.py"))
        self.config_dict["GA_final_merge_proteins"]    = self.value_assignment("path", config, "code", "ga_final_merge_proteins", os.path.join(script_path, "ga_final_merge_proteins.py"))
        self.config_dict["GA_final_merge_maps"]        = self.value_assignment("path", config, "code", "ga_final_merge_maps", os.path.join(script_path, "ga_final_merge_map.py"))
        self.config_dict["EC_Annotation_Post"]         = self.value_assignment("path", config, "code", "ec_combine", os.path.join(script_path, "ea_combine_v5.py"))
        self.config_dict["TA_apply_names"]   = self.value_assignment("path", config, "code", "ta_apply_names", os.path.join(script_path, "ta_apply_names.py"))
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
        self.config_dict["bt2_read_sorter"]            = self.value_assignment("path", config, "code", "bt2_read_sorter", os.path.join(script_path, "bwa_read_sorter.py"))
        self.config_dict["ta_contig_name_convert"]     = self.value_assignment("path", config, "code", "ta_name_convert", os.path.join(script_path, "ta_contig_name_convert.py"))
        self.config_dict["GA_pre_scan_get_lib"]        = self.value_assignment("path", config, "code", "ga_pre_scan_get_lib", os.path.join(script_path, "ga_pre_scan_get_libs.py"))
        


        
