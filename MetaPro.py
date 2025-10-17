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



#!/usr/bin/env python
from curses import meta
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_paths as mpp
import MetaPro_stages as mps
import MetaPro_dir as mpd
import MetaPro_time as mpt
import MetaPro_files as mpf
import MetaPro_marker as mpm
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q
import time
def debug_stop_check(self, stop_flag, signal):
    if(stop_flag == signal):
        sys.exit("paused at:", stop_flag)
        


def main(config_dict, dir_obj, time_obj, file_obj):
    print("in main:", config_dict["pair_1"])
    print("in main:", config_dict["single"]) 
     
    metapro_stage_obj = mps.mp_stage(config_dict, dir_obj, time_obj, file_obj) 
    
    #obj, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, args_pack, tutorial_mode)

    # This is the format we use to launch each stage of the pipeline.
    # We start a multiprocess that starts a subprocess.
    # The subprocess is created from the commands object

    # The quality filter stage
    metapro_stage_obj.mp_quality_filter()
    
    

    # The host read filter stage
    metapro_stage_obj.mp_host_filter()
        
    # The vector contaminant filter stage
    metapro_stage_obj.mp_vector_filter()

    # rRNA removal stage
    metapro_stage_obj.mp_rRNA_filter()

    # Duplicate repopulation
    metapro_stage_obj.mp_repop()

    # Assemble contigs
    metapro_stage_obj.mp_assemble()  

    #metapro_stage_obj.mp_TA()
    if(config_dict["DNA_DB_mode"] == "chocophlan"):
        metapro_stage_obj.mp_GA_pre_scan()
        
    #sys.exit("paused")
    
    # GA split
    #metapro_stage_obj.mp_GA_split()

    # GA lib check
    #metapro_stage_obj.mp_GA_lib_check()
    
    # BWA gene annotation
    metapro_stage_obj.mp_GA_BT2()
    
    
    
    metapro_stage_obj.mp_GA_BT2_pp()
   
    
    #DIAMOND gene annotation
    metapro_stage_obj.mp_GA_dmd()
    metapro_stage_obj.mp_GA_dmd_pp()
    
    # final GA merge()
    metapro_stage_obj.mp_GA_final_merge()

    # Taxonomic annotation
    metapro_stage_obj.mp_TA()

    # Detect EC annotation
    metapro_stage_obj.mp_EC()
    
    # RPKM Table and Cytoscape Network
    metapro_stage_obj.mp_output()


def tutorial_main(config_dict, dir_dict, label_dict, time_obj):
    metapro_stage_obj = mps.mp_stage(config_dict, dir_dict, label_dict, time_obj)
    if(config_dict["tutorial_mode"] == "quality"):
         # The quality filter stage
        metapro_stage_obj.mp_quality_filter()

    elif(config_dict["tutorial_mode"] == "host"):
        # The host read filter stage
        metapro_stage_obj.mp_host_filter()
            
    elif(config_dict["tutorial_mode"] == "vector"):
        # The vector contaminant filter stage
        metapro_stage_obj.mp_vector_filter()

    elif(config_dict["tutorial_mode"] == "rRNA"):
        # rRNA removal stage
        metapro_stage_obj.mp_rRNA_filter()

    elif(config_dict["tutorial_mode"] == "repop"):
        # Duplicate repopulation
        metapro_stage_obj.mp_repop()

    elif(config_dict["tutorial_mode"] == "contigs"):
        # Assemble contigs
        metapro_stage_obj.mp_assemble()    

    elif(config_dict["tutorial_mode"] == "GA"):
        #check the contig state
        metapro_stage_obj.mp_contig_statecheck()
        # GA split
        metapro_stage_obj.mp_GA_split()

        # BWA gene annotation
        
        metapro_stage_obj.mp_GA_BWA()
        metapro_stage_obj.mp_GA_BWA_pp()
        metapro_stage_obj.mp_GA_BWA_merge()
        
        # BLAT gene annotation
        metapro_stage_obj.mp_GA_BLAT()
        metapro_stage_obj.mp_GA_BLAT_pp()
        metapro_stage_obj.mp_GA_BLAT_merge()
        
        #DIAMOND gene annotation
        metapro_stage_obj.mp_GA_dmd()
        metapro_stage_obj.mp_GA_dmd_pp()
        
        # final GA merge()
        metapro_stage_obj.mp_GA_final_merge()

    elif(config_dict["tutorial_mode"] == "TA"):
        # Taxonomic annotation
        metapro_stage_obj.mp_TA()

    elif(config_dict["tutorial_mode"] == "EC"):
        # Detect EC annotation
        metapro_stage_obj.mp_EC()
    elif(config_dict["tutorial_mode"] == "output"):
        #check the contig state
        metapro_stage_obj.mp_contig_statecheck()
        # RPKM Table and Cytoscape Network
        metapro_stage_obj.mp_output()

    
    
if __name__ == "__main__":
    output_folder_default = "metapro_" + dt.today().strftime("%m%d%Y_%H%M%S")
    output_path_default = os.path.join(os.getcwd(), output_folder_default)

    print("METAPRO metatranscriptomic analysis pipeline")
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline.  Version 3.0.2 © 2025")

    parser.add_argument("-c", "--config",   type=str,   help="Path to the configureation file")
    parser.add_argument("-1", "--pair1",    type=str,   help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2",    type=str,   help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single",   type=str,   help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-con", "--contig",   type=str,   help="Tutorial use only: Path to the file containing the contig reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, help="Path of the folder for the output of the pipeline")
    parser.add_argument("--nhost", "--no-host", action='store_true', help="Skip the host read removal step of the pipeline")
    parser.add_argument("--verbose_mode", type=str, help = "Decide how to handle the interim files, Compress them, or leave them alone.  Values are: keep, compress, quiet")
    parser.add_argument("--tutorial", type = str, help = "tutorial operating mode for MetaPro")
    
    args = parser.parse_args()
    
    config_file     = args.config if args.config else "None"
    contig          = args.contig if args.contig else "None"
    pair_1          = args.pair1 if args.pair1 else "None"
    pair_2          = args.pair2 if args.pair2 else "None"
    single          = args.single if args.single else "None"
    output_folder   = args.output_folder if args.output_folder else "None"
    no_host         = args.nhost if args.nhost else False
    verbose_mode    = args.verbose_mode if args.verbose_mode else "quiet"
    tutorial_mode   = args.tutorial if args.tutorial else "None"


    if(config_file == "None"):
        print(dt.today(), "METAPRO needs a config.  exiting")
        sys.exit()


    if(not os.path.isabs(config_file)):
        config_file = os.path.abspath(config_file)
        print("full path:", config_file)
    output_folder = os.path.abspath(output_folder)
    print("Outputing to:", output_folder)
    config_obj = mpp.mpro_config(config_file, output_folder)
    config_dict = config_obj.get_config_dict()
    
    
    time_obj = mpt.mpro_timing()


    
    

    config_dict["no_host"] = no_host
    config_dict["verbose_mode"] = verbose_mode

    if(pair_1 != "None"):
        print(dt.today(), "input Pair_1 overrides config")
        pair_1 = os.path.abspath(pair_1)
        config_dict["pair_1"] = pair_1
        
    if(pair_2 != "None"):
        print(dt.today(), "input Pair_2 overrides config")
        pair_2 = os.path.abspath(pair_2)
        config_dict["pair_2"] = pair_2
    if(single != "None"):
        print(dt.today(), "input single overrides config")
        single = os.path.abspath(single)
        config_dict["single"] = single
    
    
    


    if(config_dict["tutorial_mode"] == "None"):
        if(config_dict["pair_1"] == "None"):
            if(config_dict["single"] == "None"):
                print(dt.today(), "MetaPro needs input data.  Pair 1 and singletons are blank")
                sys.exit()
            else:
                config_dict["read_mode"] = "single"
        
        else:
            if(config_dict["single"] != "None"):
                print(dt.today(), "MetaPro needs Paired or single-ended data. Not both")
                sys.exit()
            else:
                config_dict["read_mode"] = "paired"

    

    if not (os.path.exists(output_folder)):
        print("output folder does not exist.  Now building directory.")
        os.makedirs(output_folder)
    
    if not(os.path.isabs(output_folder)):
        output_folder = os.path.abspath(output_folder)
        print(dt.today(), "output destination:", output_folder)

    #Check DB integrity
    if(config_dict["no_host"] is True):
        print(dt.today(), "pre-flight check: host-DB indexing")
        host_list = config_dict["host_list"]
        for item in host_list:
            pass_flag = config_obj.check_BT2_valid(config_dict["Host_db"], item)
            if(not pass_flag):
                print(dt.today(), "MetaPro shutting down.  missing host index:", item)
                print(dt.today(), "can't find Bowtie2 index:", item, "in:", config_dict["Host_db"])
                sys.exit()
    #Check vector lib integrity
    config_obj.check_BT2_valid(config_dict["vector_db"], "vectors")
    dir_obj = mpd.mpro_dir(config_file, config_dict)
    dir_dict = dir_obj.get_dir_dict()
    label_dict = dir_obj.get_label_dict()
    file_obj = mpf.mpro_file_handler(config_dict, dir_dict)
    file_dict = file_obj.get_file_dict()

    #if (tutorial_mode != "none"):
    #    print("working in tutorial mode:", tutorial_mode)
    #    tutorial_main(config_dict, dir_dict, label_dict, time_obj)
    
    #else:
    main(config_dict, dir_obj, time_obj, file_obj)
