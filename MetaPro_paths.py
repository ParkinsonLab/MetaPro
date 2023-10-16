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
from random import expovariate
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
        
def isint(num):
    try:
        int(num)
        return True
    except ValueError:
        return False
        

class tool_path_obj:
    
    def value_assignment(self, config, config_section, var_name, default, expected_type = "float"):
        value = ""
        #print("CONFIG:", config)
        if config:
            if(config_section in config):
                if(var_name in config[config_section]):
                    value = config[config_section][var_name]
                    if('"' in value):
                        value = value.strip('"')

                    if(expected_type == "int"):
                        if not isint(value):
                            exit_msg = var_name + ", wrong data, expecting int"
                            sys.exit(exit_msg)
                        else:
                            value = int(value)
                    elif(expected_type == "float"):
                        if not isfloat(value):
                            exit_msg = var_name + ", wrong data. expecting float"
                            sys.exit(exit_msg)
                        
                    

                    elif(expected_type == "string"):
                        wrong_value_flag = False
                        if isfloat(value):
                            wrong_value_flag = True
                        if isint(value):
                            wrong_value_flag = True
                        if(wrong_value_flag):
                            exit_msg = var_name + ", wrong data. expected a name"
                            sys.exit(exit_msg)

                    elif(expected_type == "path"):
                        if(os.path.exists(value)):
                            return value
                        else:
                            #Special cases. these vars don't need filenames.
                            if(var_name == "Centrifuge_db"):
                                temp_name = value + ".1.cf"
                                if(os.path.exists(temp_name)):
                                    return value
                            elif(var_name == "SWISS_PROT"):
                                temp_name = value + ".dmnd"
                                if(os.path.exists(temp_name)):
                                    return value
                            print(var_name, "does not exist:", value)
                            sys.exit("path not valid")
                    elif(expected_type == "dir"):
                        if(os.path.isdir(value)):
                            return value
                        else:
                            exit_msg = var_name + ", wrong data. expected a file directory"
                            sys.exit(exit_msg)


                    
                        #print("quotes cleaned")
                        #print(var_name, "found! using:", value)
                        #time.sleep(1)
                    print(var_name, "found! using:", value)
                    if(value == "0"):
                        print(var_name, "zero-setting detected: using default")
                else:
                    print(var_name, "no inner section found. using default:", default)
                    value = default
            else:
                print(config_section, "no config section found, using default:", default)
                value = default
        else:
            print("no config, using default:", default)
            value = default
            
        #print(var_name, type(value))
        return value
        
        
    def check_if_indexed(self, db_file):
        db_path = os.path.dirname(db_file)
        if(os.path.exists(os.path.join(db_path, self.index_complete_marker))):
            return True
        else:
            print("not indexed")
            return False
            
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

    def __init__ (self, config_path, output_path):
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
        self.output_path        = output_path
        print("output path:", self.output_path)
        
        

        #----------------------------------------------------------
        # Reference Databases
        # Note: default host is Mouse CDS
        
        #if config:
        self.Vector_DB          = self.value_assignment(config, "Databases", "Vector", os.path.join(database_path, "univec_core/UniVec_Core.fasta"), "path") 
        self.Adapter_DB         = self.value_assignment(config, "Databases", "Adapter", os.path.join(database_path, "Trimmomatic_adapters/TruSeq3-PE-2.fa"), "path")
        self.Host_DB            = self.value_assignment(config, "Databases", "Host",  os.path.join(database_path, "Mouse_cds/Mouse_cds.fasta"), "path")
        
        

        self.Rfam               = self.value_assignment(config, "Databases", "Rfam", os.path.join(database_path, "Rfam/Rfam.cm"), "path")
        self.DNA_DB             = self.value_assignment(config, "Databases", "DNA_DB", "None", "string")
        self.source_taxa_DB     = self.value_assignment(config, "Databases", "source_taxa_db", os.path.join(database_path, "family_llbs"), "dir")
        self.Prot_DB            = self.value_assignment(config, "Databases", "Prot_DB", os.path.join(database_path, "nr/nr"), "path")
        self.Prot_DB_reads      = self.value_assignment(config, "Databases", "Prot_DB_reads", os.path.join(database_path, "nr/nr"), "path")
        self.accession2taxid    = self.value_assignment(config, "Databases", "accession2taxid", os.path.join(database_path, "accession2taxid/accession2taxid"), "path")
        self.nodes              = self.value_assignment(config, "Databases", "nodes", os.path.join(database_path, "WEVOTE_db", "nodes.dmp"), "path")
        self.names              = self.value_assignment(config, "Databases", "names", os.path.join(database_path, "WEVOTE_db", "names.dmp"), "path")
        self.Kaiju_db           = self.value_assignment(config, "Databases", "Kaiju_db", os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi"), "path")
        self.Centrifuge_db      = self.value_assignment(config, "Databases", "Centrifuge_db", os.path.join(database_path, "centrifuge_db/nt"), "path")
        self.SWISS_PROT         = self.value_assignment(config, "Databases", "SWISS_PROT", os.path.join(database_path, "swiss_prot_db/swiss_prot_db"), "path")
        self.SWISS_PROT_map     = self.value_assignment(config, "Databases", "SWISS_PROT_map", os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv"), "path")
        self.PriamDB            = self.value_assignment(config, "Databases", "PriamDB", os.path.join(database_path, "PRIAM_db/"), "path")
        self.DetectDB           = self.value_assignment(config, "Databases", "DetectDB", os.path.join(database_path, "DETECTv2"), "path")
        self.WEVOTEDB           = self.value_assignment(config, "Databases", "WEVOTEDB", os.path.join(database_path, "WEVOTE_db/"), "path")
        self.EC_pathway         = self.value_assignment(config, "Databases", "EC_pathway", os.path.join(database_path, "EC_pathway.txt"), "path")
        self.path_to_superpath  = self.value_assignment(config, "Databases", "path_to_superpath", os.path.join(custom_database_path, "pathway_to_superpathway.csv"), "path")
        self.enzyme_db          = self.value_assignment(config, "Databases", "enzyme_db", os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt"), "path")
        self.taxid_tree         = self.value_assignment(config, "Databases", "taxid_tree", os.path.join(custom_database_path, "taxid_trees", "family_tree.tsv"), "path")
        self.kraken2_db         = self.value_assignment(config, "Databases", "kraken2_db", os.path.join(custom_database_path, "kraken2_db"), "path")
        self.taxa_lib_list      = self.value_assignment(config, "Databases", "taxa_lib_list", os.path.join(database_path, "taxa_lib_list.txt"), "path")

        
        #-------------------------------------------------------
        # test DBs
        self.GA_DB_mode = "multi" #by default for the new DB changes.
        self.check_dmd_valid()
        if(self.DNA_DB_mode == "custom"):
            self.check_bwa_valid(self.DNA_DB)
            self.check_blat_valid(self.DNA_DB)
        

        
        #--------------------------------------------------
        # Memory percentage limits for applications. in INT
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
        GA_merge_mem_default = 5 
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
        GA_merge_job_limit_default          = cpu_default
        repop_job_limit_default             = 1
        TA_job_limit_default                = cpu_default
        EC_job_limit_default                = cpu_default
        Centrifuge_job_limit_default        = 1

        Barrnap_job_delay_default           = 5
        Infernal_job_delay_default          = 5
        BWA_job_delay_default               = 5
        BLAT_job_delay_default              = 5
        DIAMOND_job_delay_default           = 5
        DETECT_job_delay_default            = 5
        BWA_pp_job_delay_default            = 5
        BLAT_pp_job_delay_default           = 5
        DIAMOND_pp_job_delay_default        = 5
        GA_merge_job_delay_default    = 5
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
        keep_GA_DMD_default = "no"
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

        
#---------------------------------------------------------
        #top-level marker names
        self.top_qc_marker              = "qc_marker"
        self.index_complete_marker      = "index_complete"
        self.top_host_rem_marker        = "host_rem"
        self.top_vec_marker             = "vec_rem"
        self.top_rRNA_marker            = "rRNA_rem"
        self.top_repop_marker           = "repop"
        self.contig_marker              = "assemble_contigs"

#-----------------------------------------------------------------
        
        self.target_rank                = self.value_assignment(config, "Settings", "target_rank", "genus", "string")
        self.adapterremoval_minlength   = self.value_assignment(config, "Settings", "AdapterRemoval_minlength", 30)
        self.show_unclassified          = self.value_assignment(config, "Settings", "Show_unclassified", "No", "string")
        self.bypass_log_name            = self.value_assignment(config, "Settings", "bypass_log_name", "bypass_log.txt", "string")
        self.debug_stop_flag            = self.value_assignment(config, "Settings", "debug_stop_flag", "none", "string")
        self.num_threads                = self.value_assignment(config, "Settings", "num_threads", os.cpu_count())
        if(self.num_threads == 0):
            self.num_threads = 1
        self.taxa_exist_cutoff          = self.value_assignment(config, "Settings", "taxa_existence_cutoff", 0.1, "float")
        self.DNA_DB_mode                = self.value_assignment(config, "Settings", "DNA_DB_mode", "chocophlan", "string") #used to indicate custom DB, or our grouped library
        #other setting is "custom"
        
        
        
        self.RPKM_cutoff                = self.value_assignment(config, "Settings", "RPKM_cutoff", 0.01, "float")
        self.BWA_cigar_cutoff           = self.value_assignment(config, "Settings", "BWA_cigar_cutoff", BWA_cigar_default)
        self.BLAT_identity_cutoff       = self.value_assignment(config, "Settings", "BLAT_identity_cutoff", BLAT_identity_default)
        self.BLAT_length_cutoff         = self.value_assignment(config, "Settings", "BLAT_length_cutoff", BLAT_length_default, "float")
        self.BLAT_score_cutoff          = self.value_assignment(config, "Settings", "BLAT_score_cutoff", BLAT_score_default)            
        self.DIAMOND_identity_cutoff    = self.value_assignment(config, "Settings", "DIAMOND_identity_cutoff", DIAMOND_identity_default)
        self.DIAMOND_length_cutoff      = self.value_assignment(config, "Settings", "DIAMOND_length_cutoff", DIAMOND_length_default, "float")
        self.DIAMOND_score_cutoff       = self.value_assignment(config, "Settings", "DIAMOND_score_cutoff", DIAMOND_score_default)
        #-----------------------------------------------------------------------------------------------   

        
        self.BWA_mem_footprint      = self.value_assignment(config, "Settings", "BWA_mem_footprint", BWA_mem_footprint_default)
        self.BLAT_mem_footprint     = self.value_assignment(config, "Settings", "BLAT_mem_footprint", BLAT_mem_footprint_default)
        self.DMD_mem_footprint      = self.value_assignment(config, "Settings", "DMD_mem_footprint", DMD_mem_footprint_default)
        

        #-------------------------------------------------------------------------------------------------
        self.BWA_mem_threshold              = self.value_assignment(config, "Settings", "BWA_mem_threshold", BWA_mem_default)
        self.BLAT_mem_threshold             = self.value_assignment(config, "Settings", "BLAT_mem_threshold", BLAT_mem_default)
        self.DIAMOND_mem_threshold          = self.value_assignment(config, "Settings", "DIAMOND_mem_threshold", DIAMOND_mem_default)
        self.DETECT_mem_threshold           = self.value_assignment(config, "Settings", "DETECT_mem_threshold", DETECT_mem_default)
        self.Infernal_mem_threshold         = self.value_assignment(config, "Settings", "Infernal_mem_threshold", Infernal_mem_default)
        self.Barrnap_mem_threshold          = self.value_assignment(config, "Settings", "Barrnap_mem_threshold", Barrnap_mem_default)
        self.BWA_pp_mem_threshold           = self.value_assignment(config, "Settings", "BWA_pp_mem_threshold", BWA_pp_mem_default)
        self.BLAT_pp_mem_threshold          = self.value_assignment(config, "Settings", "BLAT_pp_mem_threshold", BLAT_pp_mem_default)
        self.DIAMOND_pp_mem_threshold       = self.value_assignment(config, "Settings", "DIAMOND_pp_mem_threshold", DIAMOND_pp_mem_default)
        self.GA_merge_mem_threshold   = self.value_assignment(config, "Settings", "GA_merge_mem_threshold", GA_merge_mem_default)
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
        self.GA_merge_job_limit   = self.value_assignment(config, "Settings", "GA_merge_job_limit", GA_merge_job_limit_default)
        self.TA_job_limit               = self.value_assignment(config, "Settings", "TA_job_limit", TA_job_limit_default)
        self.repop_job_limit            = self.value_assignment(config, "Settings", "repop_job_limit", repop_job_limit_default)
        self.EC_job_limit               = self.value_assignment(config, "Settings", "EC_job_limit", EC_job_limit_default)
        self.Centrifuge_job_limit       = self.value_assignment(config, "Settings", "Centrifuge_job_limit", Centrifuge_job_limit_default)
        
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
        self.GA_merge_job_delay   = self.value_assignment(config, "Settings", "GA_merge_job_delay", GA_merge_job_delay_default)
        self.TA_job_delay               = self.value_assignment(config, "Settings", "TA_job_delay", TA_job_delay_default)
        self.repop_job_delay            = self.value_assignment(config, "Settings", "repop_job_delay", repop_job_delay_default)
        self.EC_job_delay               = self.value_assignment(config, "Settings", "EC_job_delay", EC_job_delay_default)

        #------------------------------------------------------------------------------------------------
        self.keep_all                   = self.value_assignment(config, "Settings", "keep_all", keep_all_default, "string")
        self.keep_quality               = self.value_assignment(config, "Settings", "keep_quality", keep_quality_default, "string")
        self.keep_host                  = self.value_assignment(config, "Settings", "keep_host", keep_host_default, "string")
        self.keep_vector                = self.value_assignment(config, "Settings", "keep_vector", keep_vector_default, "string")
        self.keep_rRNA                  = self.value_assignment(config, "Settings", "keep_rRNA", keep_rRNA_default, "string")
        self.keep_repop                 = self.value_assignment(config, "Settings", "keep_repop", keep_repop_default, "string")
        self.keep_assemble_contigs      = self.value_assignment(config, "Settings", "keep_assemble_contigs", keep_assemble_contigs_default, "string")
        self.keep_GA_BWA                = self.value_assignment(config, "Settings", "keep_GA_BWA", keep_GA_BWA_default, "string")
        self.keep_GA_BLAT               = self.value_assignment(config, "Settings", "keep_GA_BLAT", keep_GA_BLAT_default, "string")
        self.keep_GA_DMD            = self.value_assignment(config, "Settings", "keep_GA_DMD", keep_GA_DMD_default, "string")
        self.keep_GA_final              = self.value_assignment(config, "Settings", "keep_GA_final", keep_GA_final_default, "string")
        self.keep_TA                    = self.value_assignment(config, "Settings", "keep_TA", keep_TA_default, "string")
        self.keep_EC                    = self.value_assignment(config, "Settings", "keep_EC", keep_EC_default, "string")
        self.keep_outputs               = self.value_assignment(config, "Settings", "keep_outputs", keep_outputs_default, "string")
        

        
        
        #--------------------------------------------------------------------------------------
        
        self.filter_stringency          = self.value_assignment(config, "Settings", "filter_stringency", filter_stringency_default, "string")
        self.GA_chunksize               = self.value_assignment(config, "Settings", "GA_chunk_size", GA_chunksize_default)
        self.EC_chunksize               = self.value_assignment(config, "Settings", "EC_chunk_size", EC_chunksize_default)
        self.rRNA_chunksize             = self.value_assignment(config, "Settings", "rRNA_chunk_size", rRNA_chunksize_default)
        
     
        
    
        #--------------------------------------------------------------------------------------------
        # Labels.  
        # why? to change them during integration + new feature testing

        qc_label_default                                = "quality_filter"
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
        GA_DMD_label_default                            = "GA_DMD"
        GA_DMD_pp_label_default                         = "GA_DMD_pp"
        GA_merge_label_default                    = "GA_merge"
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
        

        self.qc_label                               = self.value_assignment(config, "Labels", "quality_filter",                     qc_label_default, "string")
        
        self.host_filter_label                      = self.value_assignment(config, "Labels", "host_filter",                        host_filter_label_default, "string")
        self.vector_filter_label                    = self.value_assignment(config, "Labels", "vector_filter",                      vector_filter_label_default, "string")
        self.rRNA_filter_label                      = self.value_assignment(config, "Labels", "rRNA_filter",                        rRNA_filter_label_default, "string")
        self.rRNA_filter_split_label                = self.value_assignment(config, "Labels", "rRNA_filter_split",                  rRNA_filter_split_label_default, "string")
        self.rRNA_filter_convert_label              = self.value_assignment(config, "Labels", "rRNA_filter_convert",                rRNA_filter_convert_label_default, "string")
        self.rRNA_filter_barrnap_label              = self.value_assignment(config, "Labels", "rRNA_filter_barrnap",                rRNA_filter_barrnap_label_default, "string")
        self.rRNA_filter_barrnap_merge_label        = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_merge",          rRNA_filter_barrnap_merge_label_default, "string")
        self.rRNA_filter_barrnap_pp_label           = self.value_assignment(config, "Labels", "rRNA_filter_barrnap_pp",             rRNA_filter_barrnap_pp_label_default, "string")
        self.rRNA_filter_infernal_label             = self.value_assignment(config, "Labels", "rRNA_filter_infernal",               rRNA_filter_infernal_label_default, "string")
        self.rRNA_filter_infernal_prep_label        = self.value_assignment(config, "Labels", "rRNA_filter_infernal_prep",          rRNA_filter_infernal_prep_label_default, "string")
        self.rRNA_filter_splitter_label             = self.value_assignment(config, "Labels", "rRNA_filter_splitter",               rRNA_filter_splitter_label_default, "string")
        self.rRNA_filter_post_label                 = self.value_assignment(config, "Labels", "rRNA_filter_post",                   rRNA_filter_post_label_default, "string")
        self.repop_label                            = self.value_assignment(config, "Labels", "repop",                              repop_label_default, "string")
        self.assemble_contigs_label                 = self.value_assignment(config, "Labels", "assemble_contigs",                   assemble_contigs_label_default, "string")
        self.destroy_contigs_label                  = self.value_assignment(config, "Labels", "destroy_contigs",                    destroy_contigs_label_default, "string")
        self.GA_pre_scan_label                      = self.value_assignment(config, "Labels", "GA_pre_scan",                        GA_pre_scan_label_default, "string")
        self.GA_split_label                         = self.value_assignment(config, "Labels", "GA_split",                           GA_split_label_default, "string")
        self.GA_BWA_label                           = self.value_assignment(config, "Labels", "GA_BWA",                             GA_BWA_label_default, "string")
        self.GA_BWA_pp_label                        = self.value_assignment(config, "Labels", "GA_BWA_pp",                          GA_BWA_pp_label_default, "string")
        self.GA_BWA_merge_label                     = self.value_assignment(config, "Labels", "GA_BWA_merge",                       GA_BWA_merge_label_default, "string")
        self.GA_BLAT_label                          = self.value_assignment(config, "Labels", "GA_BLAT",                            GA_BLAT_label_default, "string")
        self.GA_BLAT_cleanup_label                  = self.value_assignment(config, "Labels", "GA_BLAT_cleanup",                    GA_BLAT_cleanup_label_default, "string")
        self.GA_BLAT_cat_label                      = self.value_assignment(config, "Labels", "GA_BLAT_cat",                        GA_BLAT_cat_label_default, "string")
        self.GA_BLAT_pp_label                       = self.value_assignment(config, "Labels", "GA_BLAT_pp",                         GA_BLAT_pp_label_default, "string")
        self.GA_BLAT_merge_label                    = self.value_assignment(config, "Labels", "GA_BLAT_merge",                      GA_BLAT_merge_label_default, "string")
        self.GA_DMD_label                           = self.value_assignment(config, "Labels", "GA_DMD",                             GA_DMD_label_default, "string")
        self.GA_DMD_pp_label                        = self.value_assignment(config, "Labels", "GA_DMD_pp",                          GA_DMD_pp_label_default, "string")
        self.GA_merge_label                   = self.value_assignment(config, "Labels", "GA_merge",                     GA_merge_label_default, "string")
        self.ta_label                               = self.value_assignment(config, "Labels", "ta",                                 taxon_annotation_label_default, "string")
        self.EC_label                               = self.value_assignment(config, "Labels", "ec",                                 ec_annotation_label_default, "string")
        self.ec_detect_label                        = self.value_assignment(config, "Labels", "ec_detect",                          ec_annotation_detect_label_default, "string")
        self.ec_priam_label                         = self.value_assignment(config, "Labels", "ec_priam",                           ec_annotation_priam_label_default, "string")
        self.ec_priam_split_label                   = self.value_assignment(config, "Labels", "ec_priam_split",                     ec_annotation_priam_split_label_default, "string")
        self.ec_priam_cat_label                     = self.value_assignment(config, "Labels", "ec_priam_cat",                       ec_annotation_priam_cat_label_default, "string")
        self.ec_DIAMOND_label                       = self.value_assignment(config, "Labels", "ec_DIAMOND",                         ec_annotation_DIAMOND_label_default, "string")
        self.ec_pp_label                            = self.value_assignment(config, "Labels", "ec_pp",                              ec_annotation_pp_label_default, "string")
        self.output_label                           = self.value_assignment(config, "Labels", "outputs",                            output_label_default, "string")
        self.output_copy_gene_map_label             = self.value_assignment(config, "Labels", "output_copy_gene_map",               output_copy_gene_map_label_default, "string")
        self.output_clean_EC_label                  = self.value_assignment(config, "Labels", "output_clean_ec",                    output_clean_EC_label_default, "string")
        self.output_copy_taxa_label                 = self.value_assignment(config, "Labels", "output_copy_taxa",                   output_copy_taxa_label_default, "string")
        self.output_network_generation_label        = self.value_assignment(config, "Labels", "output_network_generation",          output_network_gen_label_default, "string")
        self.output_unique_hosts_singletons_label   = self.value_assignment(config, "Labels", "output_unique_hosts_singletons",     output_unique_hosts_singletons_label_default, "string")
        self.output_unique_hosts_pair_1_label       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_1",         output_unique_hosts_pair_1_label_default, "string")
        self.output_unique_hosts_pair_2_label       = self.value_assignment(config, "Labels", "output_unique_hosts_pair_2",         output_unique_hosts_pair_2_label_default, "string")
        self.output_unique_vectors_singletons_label = self.value_assignment(config, "Labels", "output_unique_vectors_singletons",   output_unique_vectors_singletons_label_default, "string")
        self.output_unique_vectors_pair_1_label     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_1",       output_unique_vectors_pair_1_label_default, "string")
        self.output_unique_vectors_pair_2_label     = self.value_assignment(config, "Labels", "output_unique_vectors_pair_2",       output_unique_vectors_pair_2_label_default, "string")
        self.output_combine_hosts_label             = self.value_assignment(config, "Labels", "output_combine_hosts",               output_combine_hosts_label_default, "string")
        self.output_per_read_scores_label           = self.value_assignment(config, "Labels", "output_per_read_scores",             output_per_read_scores_label_default, "string")
        self.output_contig_stats_label              = self.value_assignment(config, "Labels", "output_contig_stats",                output_contig_stats_label_default, "string")
        self.output_ec_heatmap_label                = self.value_assignment(config, "Labels", "output_ec_heatmap",                  output_ec_heatmap_label_default, "string")
        self.output_taxa_groupby_label              = self.value_assignment(config, "Labels", "output_taxa_groupby",                output_taxa_groupby_label_default, "string")
        self.output_read_count_label                = self.value_assignment(config, "Labels", "output_read_count",                  output_read_count_label_default, "string")
        
        
        
        
        #if config:
        self.Vector_DB_file     = self.value_assignment(config, "Databases", "Vector_DB", os.path.join(database_path, "univec_core/UniVec_Core.fasta"), "path") 
        if(self.check_if_indexed(os.path.dirname(self.Vector_DB_file))):
            print("Vector DB indexed!")
        else:
            print("Vector DB not indexed")

        sys.exit("paused")

        self.Host_DB            = self.value_assignment(config, "Databases", "Host_DB",  os.path.join(database_path, "Mouse_cds/Mouse_cds.fasta"), "path")
        self.Rfam               = self.value_assignment(config, "Databases", "Rfam", os.path.join(database_path, "Rfam/Rfam.cm"), "path")
        self.DNA_DB             = self.value_assignment(config, "Databases", "DNA_DB", "None", "string")
        self.source_taxa_DB     = self.value_assignment(config, "Databases", "source_taxa_db", os.path.join(database_path, "family_llbs"), "dir")
        self.Prot_DB            = self.value_assignment(config, "Databases", "Prot_DB", os.path.join(database_path, "nr/nr"), "path")
        self.Prot_DB_reads      = self.value_assignment(config, "Databases", "Prot_DB_reads", os.path.join(database_path, "nr/nr"), "path")
        self.accession2taxid    = self.value_assignment(config, "Databases", "accession2taxid", os.path.join(database_path, "accession2taxid/accession2taxid"), "path")
        self.nodes              = self.value_assignment(config, "Databases", "nodes", os.path.join(database_path, "WEVOTE_db", "nodes.dmp"), "path")
        self.names              = self.value_assignment(config, "Databases", "names", os.path.join(database_path, "WEVOTE_db", "names.dmp"), "path")
        self.Kaiju_db           = self.value_assignment(config, "Databases", "Kaiju_db", os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi"), "path")
        self.Centrifuge_db      = self.value_assignment(config, "Databases", "Centrifuge_db", os.path.join(database_path, "centrifuge_db/nt"), "path")
        self.SWISS_PROT         = self.value_assignment(config, "Databases", "SWISS_PROT", os.path.join(database_path, "swiss_prot_db/swiss_prot_db"), "path")
        self.SWISS_PROT_map     = self.value_assignment(config, "Databases", "SWISS_PROT_map", os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv"), "path")
        self.PriamDB            = self.value_assignment(config, "Databases", "PriamDB", os.path.join(database_path, "PRIAM_db/"), "path")
        self.DetectDB           = self.value_assignment(config, "Databases", "DetectDB", os.path.join(database_path, "DETECTv2"), "path")
        self.WEVOTEDB           = self.value_assignment(config, "Databases", "WEVOTEDB", os.path.join(database_path, "WEVOTE_db/"), "path")
        self.EC_pathway         = self.value_assignment(config, "Databases", "EC_pathway", os.path.join(database_path, "EC_pathway.txt"), "path")
        self.path_to_superpath  = self.value_assignment(config, "Databases", "path_to_superpath", os.path.join(custom_database_path, "pathway_to_superpathway.csv"), "path")
        self.enzyme_db          = self.value_assignment(config, "Databases", "enzyme_db", os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt"), "path")
        self.taxid_tree         = self.value_assignment(config, "Databases", "taxid_tree", os.path.join(custom_database_path, "taxid_trees", "family_tree.tsv"), "path")
        self.kraken2_db         = self.value_assignment(config, "Databases", "kraken2_db", os.path.join(custom_database_path, "kraken2_db"), "path")
        self.taxa_lib_list      = self.value_assignment(config, "Databases", "taxa_lib_list", os.path.join(database_path, "taxa_lib_list.txt"), "path")

        
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
        self.Python         = self.value_assignment(config, "Tools", "Python", "python3", "string")
        self.Java           = self.value_assignment(config, "Tools", "Java", "java -jar", "string")
        self.cdhit_dup      = self.value_assignment(config, "Tools", "cdhit_dup",  os.path.join(tool_path, "cdhit_dup/cd-hit-dup"), "path")
        self.AdapterRemoval = self.value_assignment(config, "Tools", "AdapterRemoval", os.path.join(tool_path, "adapterremoval/AdapterRemoval"), "path")
        self.vsearch        = self.value_assignment(config, "Tools", "vsearch", os.path.join(tool_path, "vsearch/vsearch"), "path")
        self.BWA            = self.value_assignment(config, "Tools", "BWA", os.path.join(tool_path, "BWA/bwa"), "path")
        self.SAMTOOLS       = self.value_assignment(config, "Tools", "SAMTOOLS", os.path.join(tool_path, "samtools/samtools"), "path")
        self.BLAT           = self.value_assignment(config, "Tools", "BLAT", os.path.join(tool_path, "PBLAT/pblat"), "path")
        self.DIAMOND        = self.value_assignment(config, "Tools", "DIAMOND", os.path.join(tool_path, "DIAMOND/diamond"), "path")
        self.Blastp         = self.value_assignment(config, "Tools", "Blastp", os.path.join(tool_path, "BLAST_p/blastp"), "path")
        self.Needle         = self.value_assignment(config, "Tools", "Needle", os.path.join(tool_path, "EMBOSS-6.6.0/emboss/stretcher"), "path")
        self.Makeblastdb    = self.value_assignment(config, "Tools", "Makeblastdb", os.path.join(tool_path, "BLAST_p/makeblastdb"), "path")
        self.Barrnap        = self.value_assignment(config, "Tools", "Barrnap", os.path.join(tool_path, "Barrnap/bin/barrnap"), "path")
        self.Infernal       = self.value_assignment(config, "Tools", "Infernal", os.path.join(tool_path, "infernal/cmsearch"), "path")
        self.Kaiju          = self.value_assignment(config, "Tools", "Kaiju", os.path.join(tool_path, "kaiju/kaiju"), "path")
        self.Centrifuge     = self.value_assignment(config, "Tools", "Centrifuge", os.path.join(tool_path, "centrifuge/centrifuge"), "path")
        self.Priam          = self.value_assignment(config, "Tools", "Priam", os.path.join(tool_path, "PRIAM_search/PRIAM_search.jar"), "path")
        self.Detect         = self.value_assignment(config, "Tools", "Detect", os.path.join(script_path, "Detect_2.2.10.py"), "path")
        self.BLAST_dir      = self.value_assignment(config, "Tools", "BLAST_dir", os.path.join(tool_path, "BLAST_p"), "dir")
        self.WEVOTE         = self.value_assignment(config, "Tools", "WEVOTE", os.path.join(tool_path, "WEVOTE/WEVOTE"), "path")
        self.Spades         = self.value_assignment(config, "Tools", "Spades", os.path.join(tool_path, "SPAdes/bin/spades.py"), "path")
        self.MetaGeneMark   = self.value_assignment(config, "Tools", "MetaGeneMark", os.path.join(tool_path, "mgm/gmhmmp"), "path")
        self.kraken2        = self.value_assignment(config, "Tools", "kraken2", os.path.join(tool_path, "kraken2/kraken2"), "path")
            
        #--------------------------------------------
        # Python scripts
        #if config:
        self.sam_trimmer                = self.value_assignment(config, "code", "sam_trimmer", os.path.join(script_path, "read_sam.py"), "path")
        self.sort_reads                 = self.value_assignment(config, "code", "sort_reads", os.path.join(script_path, "read_sort.py"), "path")
        self.duplicate_repopulate       = self.value_assignment(config, "code", "duplicate_repopulation", os.path.join(script_path, "read_repopulation.py"), "path")
        self.orphaned_read_filter       = self.value_assignment(config, "code", "orphaned_read_filter", os.path.join(script_path, "read_orphan.py"), "path")
        self.remove_tag                 = self.value_assignment(config, "code", "remove_tag", os.path.join(script_path, "read_remove_tag.py"), "path")
        self.BLAT_Contaminant_Filter    = self.value_assignment(config, "code", "blat_contaminant_filter", os.path.join(script_path, "read_BLAT_filter_v3.py"), "path")
        self.File_splitter              = self.value_assignment(config, "code", "file_splitter", os.path.join(script_path, "read_split.py"), "path")
        self.rRNA_barrnap_pp            = self.value_assignment(config, "code", "rRNA_barrnap_pp", os.path.join(script_path, "read_rRNA_barrnap_v2.py"), "path")
        self.rRNA_infernal_pp           = self.value_assignment(config, "code", "rRNA_infernal_pp", os.path.join(script_path, "read_rRNA_infernal_v2.py"), "path")
        self.read_split_convert         = self.value_assignment(config, "code", "read_split_convert", os.path.join(script_path, "read_split_and_convert.py"), "path")
        self.Map_contig                 = self.value_assignment(config, "code", "map_contig", os.path.join(script_path, "assembly_make_contig_map.py"), "path")
        self.flush_bad_contigs          = self.value_assignment(config, "code", "flush_bad_contigs", os.path.join(script_path, "assembly_flush_bad_contigs.py"), "path")
        self.contig_duplicate_remover   = self.value_assignment(config, "code", "contig_duplicate_remover", os.path.join(script_path, "assembly_deduplicate.py"), "path")
        self.Map_reads_gene_BWA         = self.value_assignment(config, "code", "ga_bwa_pp", os.path.join(script_path, "ga_BWA_generic_v2.py"), "path")
        self.Map_reads_gene_BLAT        = self.value_assignment(config, "code", "ga_blat_pp", os.path.join(script_path, "ga_BLAT_generic_v3.py"), "path")
        self.Map_reads_prot_DMND        = self.value_assignment(config, "code", "ga_dmd_pp", os.path.join(script_path, "GA_DMD_generic_v2.py"), "path")
        self.GA_merge                   = self.value_assignment(config, "code", "GA_merge", os.path.join(script_path, "GA_merge_v4.py"), "path")
        self.GA_merge_fasta             = self.value_assignment(config, "code", "ga_merge_fasta", os.path.join(script_path, "ga_merge_fasta.py"), "path")
        self.GA_merge_fasta             = self.value_assignment(config, "code", "GA_merge_fasta", os.path.join(script_path, "GA_merge_fastq.py"), "path")
        self.GA_merge_proteins          = self.value_assignment(config, "code", "GA_merge_proteins", os.path.join(script_path, "GA_merge_proteins.py"), "path")
        self.GA_merge_maps              = self.value_assignment(config, "code", "GA_merge_maps", os.path.join(script_path, "GA_merge_map.py"), "path")
        self.EC_Annotation_Post         = self.value_assignment(config, "code", "ec_combine", os.path.join(script_path, "ea_combine_v5.py"), "path")
        self.Annotated_taxid            = self.value_assignment(config, "code", "ta_taxid", os.path.join(script_path, "ta_taxid_v3.py"), "path")
        self.Constrain_classification   = self.value_assignment(config, "code", "ta_constrain", os.path.join(script_path, "ta_constrain_taxonomy_v2.py"), "path")
        self.Classification_combine     = self.value_assignment(config, "code", "ta_combine", os.path.join(script_path, "ta_combine_v3.py"), "path")
        self.Wevote_parser              = self.value_assignment(config, "code", "ta_wevote", os.path.join(script_path, "ta_wevote_parser.py"), "path")
        self.taxa_table                 = self.value_assignment(config, "code", "output_taxa", os.path.join(script_path, "output_taxa_groupby.py"), "path")
        self.RPKM                       = self.value_assignment(config, "code", "output_rpkm", os.path.join(script_path, "output_table_v3.py"), "path")
        self.format_RPKM                = self.value_assignment(config, "code", "output_reformat", os.path.join(script_path, "output_reformat_rpkm_table.py"), "path")
        self.read_count                 = self.value_assignment(config, "code", "output_read_count", os.path.join(script_path, "output_read_counts_v2.py"), "path")
        self.read_quality_metrics       = self.value_assignment(config, "code", "output_qual", os.path.join(script_path, "output_read_quality_metrics.py"), "path")
        self.contig_stats               = self.value_assignment(config, "code", "output_contig_stats", os.path.join(script_path, "output_contig_stats.py"), "path")
        self.ec_heatmap                 = self.value_assignment(config, "code", "output_heatmap", os.path.join(script_path, "output_EC_metrics.py"), "path")
        self.data_change_metrics        = self.value_assignment(config, "code", "output_data_change", os.path.join(script_path, "output_data_change_metrics.py"), "path")
        self.get_unique_host_reads      = self.value_assignment(config, "code", "output_unique_hosts", os.path.join(script_path, "output_get_host_reads.py"), "path")
        self.remove_gaps_in_fasta       = self.value_assignment(config, "code", "remove_fasta_gaps", os.path.join(script_path, "remove_gaps_in_fasta.py"), "path")
        self.parse_sam                  = self.value_assignment(config, "code", "output_sam", os.path.join(script_path, "output_parse_sam.py"), "path")
        self.are_you_in_a_contig        = self.value_assignment(config, "code", "output_contig_check", os.path.join(script_path, "output_are_you_in_a_contig.py"), "path")
        self.convert_contig_segments    = self.value_assignment(config, "code", "output_convert_gene_map", os.path.join(script_path, "output_convert_gene_map_contig_segments.py"), "path")
        self.output_filter_taxa         = self.value_assignment(config, "code", "output_filter_taxa", os.path.join(script_path, "output_filter_taxa.py"), "path")
        self.output_filter_ECs          = self.value_assignment(config, "code", "output_filter_ec", os.path.join(script_path, "output_filter_ECs.py"), "path")
        self.bwa_read_sorter            = self.value_assignment(config, "code", "bwa_read_sorter", os.path.join(script_path, "bwa_read_sorter.py"), "path")
        self.ta_contig_name_convert     = self.value_assignment(config, "code", "ta_name_convert", os.path.join(script_path, "ta_contig_name_convert.py"), "path")
        self.GA_pre_scan_get_lib        = self.value_assignment(config, "code", "ga_pre_scan_get_lib", os.path.join(script_path, "ga_pre_scan_get_libs.py"), "path")
        self.GA_pre_scan_assemble_lib   = self.value_assignment(config, "code", "ga_pre_scan_assemble_lib", os.path.join(script_path, "ga_pre_scan_assemble_libs.py"), "path")
        
        
         


        #---------------------------------------------------------------------------------------
        # Folder names + paths

        self.qc_top_path                = os.path.join(self.output_path, self.qc_label)
        self.qc_data_path               = os.path.join(self.qc_top_path, "data")
        self.qc_jobs_path               = os.path.join(self.qc_top_path, "jobs")
        self.qc_sort_path               = os.path.join(self.qc_data_path, "0_sorted_raw_input")
        self.qc_adapter_path            = os.path.join(self.qc_data_path, "1_adapter_removal")
        self.qc_merge_path              = os.path.join(self.qc_data_path, "2_vsearch_pair_merge")
        self.qc_filter_path             = os.path.join(self.qc_data_path, "3_quality_filter")
        self.qc_orphan_path             = os.path.join(self.qc_data_path, "4_orphan_read_filter")
        self.qc_cdhit_path              = os.path.join(self.qc_data_path, "5_remove_duplicates")
        self.qc_final_path              = os.path.join(self.qc_top_path, "final_results")

        self.host_top_path          = os.path.join(self.output_path, self.host_filter_label)
        self.host_data_path         = os.path.join(self.host_top_path, "data")
        self.host_jobs_path         = os.path.join(self.host_top_path, "jobs")
        self.host_bwa_path          = os.path.join(self.host_data_path, "0_remove_host")
        self.host_blat_path         = os.path.join(self.host_data_path, "1_blat_host")
        self.host_final_path        = os.path.join(self.host_top_path, "final_results")

        self.vector_top_path        = os.path.join(self.output_path, self.vector_filter_label)
        self.vector_data_path       = os.path.join(self.vector_top_path, "data")
        self.vector_jobs_path       = os.path.join(self.vector_top_path, "jobs")
        self.vector_bwa_path        = os.path.join(self.vector_data_path, "0_vector_removal")
        self.vector_blat_path       = os.path.join(self.vector_data_path, "1_blat_containment_vr")
        self.vector_final_path      = os.path.join(self.vector_top_path, "final_results")

        self.rRNA_top_path          = os.path.join(self.output_path, self.rRNA_filter_label)
        self.rRNA_data_path         = os.path.join(self.rRNA_top_path, "data")
        self.rRNA_jobs_path         = os.path.join(self.rRNA_top_path, "jobs")
        self.rRNA_exe_path          = os.path.join(self.rRNA_top_path, "exe")
        
        self.rRNA_p1_fa_path        = os.path.join(self.rRNA_data_path, "pair_1_fasta")
        self.rRNA_p1_bar_path       = os.path.join(self.rRNA_data_path, "pair_1_barrnap")
        self.rRNA_p1_inf_path       = os.path.join(self.rRNA_data_path, "pair_1_infernal")
        self.rRNA_p1_inf_mRNA_path  = os.path.join(self.rRNA_data_path, "pair_1_infernal_mRNA")
        self.rRNA_p1_inf_tRNA_path  = os.path.join(self.rRNA_data_path, "pair_1_infernal_other") 

        self.rRNA_p2_fa_path        = os.path.join(self.rRNA_data_path, "pair_2_fasta")
        self.rRNA_p2_bar_path       = os.path.join(self.rRNA_data_path, "pair_2_barrnap")
        self.rRNA_p2_inf_path       = os.path.join(self.rRNA_data_path, "pair_2_infernal")
        self.rRNA_p2_inf_mRNA_path  = os.path.join(self.rRNA_data_path, "pair_2_infernal_mRNA")
        self.rRNA_p2_inf_tRNA_path  = os.path.join(self.rRNA_data_path, "pair_2_infernal_other")

        self.rRNA_s_fa_path         = os.path.join(self.rRNA_data_path, "singletons_fasta")
        self.rRNA_s_bar_path        = os.path.join(self.rRNA_data_path, "singletons_barrnap")
        self.rRNA_s_inf_path        = os.path.join(self.rRNA_data_path, "singletons_infernal")
        self.rRNA_s_inf_mRNA_path   = os.path.join(self.rRNA_data_path, "singletons_infernal_mRNA")
        self.rRNA_s_inf_tRNA_path   = os.path.join(self.rRNA_data_path, "singletons_infernal_other")        
        
        self.rRNA_final_path        = os.path.join(self.rRNA_top_path, "final_results")
        self.rRNA_final_mRNA_path   = os.path.join(self.rRNA_final_path, "mRNA")
        self.rRNA_final_tRNA_path   = os.path.join(self.rRNA_final_path, "other")
        
        self.repop_top_path         = os.path.join(self.output_path, self.repop_label)
        self.repop_data_path        = os.path.join(self.repop_top_path, "data")
        self.repop_work_path        = os.path.join(self.repop_data_path, "0_repop")
        self.repop_final_path       = os.path.join(self.repop_top_path, "final_results")


        self.contigs_top_path           = os.path.join(self.output_path, self.assemble_contigs_label)
        self.contigs_data_path          = os.path.join(self.contigs_top_path, "data")
        self.contigs_spades_path        = os.path.join(self.contigs_data_path, "0_spades")
        self.contigs_mgm_path           = os.path.join(self.contigs_data_path, "1_mgm")
        self.contigs_bwa_path           = os.path.join(self.contigs_data_path, "2_bwa_align")
        self.contigs_map_path           = os.path.join(self.contigs_data_path, "3_mapped_reads")
        self.contigs_final_path         = os.path.join(self.contigs_top_path, "final_results")


        self.GA_pre_scan_top_path       = os.path.join(self.output_path, self.GA_pre_scan_label)
        self.GA_pre_scan_data_path      = os.path.join(self.GA_pre_scan_top_path, "data")
        self.GA_pre_scan_jobs_path      = os.path.join(self.GA_pre_scan_top_path, "jobs")
        self.GA_pre_scan_kraken_path    = os.path.join(self.GA_pre_scan_data_path, "1_kraken2")
        self.GA_pre_scan_centr_path     = os.path.join(self.GA_pre_scan_data_path, "2_centrifuge")
        self.GA_pre_scan_wevote_path    = os.path.join(self.GA_pre_scan_data_path, "3_wevote")
        self.GA_pre_scan_libs_path      = os.path.join(self.GA_pre_scan_data_path, "4_libs")
        self.GA_pre_scan_final_path     = os.path.join(self.GA_pre_scan_top_path, "final_results")

        self.GA_split_top_path      = os.path.join(self.output_path, self.GA_split_label)
        self.GA_split_data_path     = os.path.join(self.GA_split_top_path, "data")
        self.GA_split_jobs_path     = os.path.join(self.GA_split_top_path, "jobs")
        self.GA_split_p1_path       = os.path.join(self.GA_split_data_path, "pair_1")
        self.GA_split_p2_path       = os.path.join(self.GA_split_data_path, "pair_2")
        self.GA_split_c_path        = os.path.join(self.GA_split_data_path, "contigs")
        self.GA_split_s_path        = os.path.join(self.GA_split_data_path, "contigs")
        self.GA_split_final_path    = os.path.join(self.GA_split_top_path, "final_result")

        self.GA_BWA_top_path        = os.path.join(self.output_path, self.GA_BWA_label)
        self.GA_BWA_data_path       = os.path.join(self.GA_BWA_top_path, "data")
        self.GA_BWA_jobs_path       = os.path.join(self.GA_BWA_top_path, "jobs")
        self.GA_BWA_split_path      = os.path.join(self.GA_BWA_data_path, "0_read_split")
        self.GA_BWA_run_path        = os.path.join(self.GA_BWA_data_path, "1_bwa")
        self.GA_BWA_pp_path         = os.path.join(self.GA_BWA_data_path, "2_bwa_pp")
        self.GA_BWA_final_path      = os.path.join(self.GA_BWA_top_path, "final_results")

        self.GA_BLAT_top_path           = os.path.join(self.output_path, self.GA_BLAT_label)
        self.GA_BLAT_data_path      = os.path.join(self.GA_BLAT_top_path, "data")
        self.GA_BLAT_jobs_path      = os.path.join(self.GA_BLAT_top_path, "jobs")
        self.GA_BLAT_run_path       = os.path.join(self.GA_BLAT_data_path, "0_blat")
        self.GA_BLAT_pp_path        = os.path.join(self.GA_BLAT_data_path, "1_pp")
        self.GA_BLAT_final_path     = os.path.join(self.GA_BLAT_top_path, "final_results")

        self.GA_DMD_top_path        = os.path.join(self.output_path, self.GA_DMD_label)
        self.GA_DMD_data_path       = os.path.join(self.GA_DMD_top_path, "data")
        self.GA_DMD_jobs_path       = os.path.join(self.GA_DMD_top_path, "jobs")
        self.GA_DMD_tool_path       = os.path.join(self.GA_DMD_data_path, "0_dmd")
        self.GA_DMD_temp_path       = os.path.join(self.GA_DMD_tool_path, "temp")
        self.GA_DMD_final_path      = os.path.join(self.GA_DMD_top_path, "final_results")

        self.GA_merge_top_path      = os.path.join(self.output_path, self.GA_merge_label)
        self.GA_merge_data_path     = os.path.join(self.GA_merge_top_path, "data")
        self.GA_merge_jobs_path     = os.path.join(self.GA_merge_top_path, "jobs")
        self.GA_merge_final_path    = os.path.join(self.GA_merge_top_path, "final_results")

        self.TA_top_path            = os.path.join(self.output_path, self.ta_label)
        self.TA_data_path           = os.path.join(self.TA_top_path, "data")
        self.TA_jobs_path           = os.path.join(self.TA_top_path, "jobs")
        self.TA_ga_path             = os.path.join(self.TA_data_path, "0_gene_taxa")
        self.TA_kraken_path         = os.path.join(self.TA_data_path, "1_kraken2")
        self.TA_centr_path          = os.path.join(self.TA_data_path, "2_centrifuge")
        self.TA_wevote_path         = os.path.join(self.TA_data_path, "3_wevote")
        self.TA_final_path          = os.path.join(self.TA_top_path, "final_results")

        self.EC_top_path            = os.path.join(self.output_path, self.EC_label)
        self.EC_data_path           = os.path.join(self.EC_top_path, "data")
        self.EC_jobs_path           = os.path.join(self.EC_top_path, "jobs")
        self.EC_detect_path         = os.path.join(self.EC_data_path, "0_detect")
        self.EC_priam_path          = os.path.join(self.EC_data_path, "1_priam")
        self.EC_dmd_path            = os.path.join(self.EC_data_path, "2_diamond")
        self.EC_final_path          = os.path.join(self.EC_top_path, "final_results")

        self.reports_top_path       = os.path.join(self.output_path, self.output_label)
        self.reports_data_path      = os.path.join(self.reports_top_path, "data")
        self.reports_jobs_path      = os.path.join(self.reports_top_path, "jobs")
        self.reports_uhosts_path    = os.path.join(self.reports_data_path, "1_unique_hosts")
        self.reports_fhosts_path    = os.path.join(self.reports_data_path, "2_full_hosts")
        self.reports_uvecs_path     = os.path.join(self.reports_data_path, "3_unique_vectors")
        self.reports_fvecs_path     = os.path.join(self.reports_data_path, "4_full_vectors")
        self.reports_final_path     = os.path.join(self.reports_top_path, "final_results")




#---------------------------------------------------------
        #internal markers

        self.rRNA_barrnap_marker    = "barrnap_"
        self.rRNA_inf_marker        = "_inf"
        self.rRNA_split_marker      = "split_fasta"


#----------------------------------------------------------
        #special contig-bypasser logic vars
        self.spades_done_file = os.path.join(self.contigs_spades_path, "pipeline_state", "stage_7_terminate")
        self.spades_transcripts_file = os.path.join(self.contigs_spades_path, "transcripts.fasta")      