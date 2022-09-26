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
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import MetaPro_utilities as mpu
import MetaPro_stages as mps
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q


def main(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode):

    metapro_stage_obj = mps.mp_stage(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode)

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

    # Taxonomic annotation
    metapro_stage_obj.mp_TA()

    # Detect EC annotation
    metapro_stage_obj.mp_EC()
    
    # RPKM Table and Cytoscape Network
    metapro_stage_obj.mp_output()


def tutorial_main(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode):
    metapro_stage_obj = mps.mp_stage(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode)
    if(tutorial_mode == "quality"):
         # The quality filter stage
        metapro_stage_obj.mp_quality_filter()

    elif(tutorial_mode == "host"):
        # The host read filter stage
        metapro_stage_obj.mp_host_filter()
            
    elif(tutorial_mode == "vector"):
        # The vector contaminant filter stage
        metapro_stage_obj.mp_vector_filter()

    elif(tutorial_mode == "rRNA"):
        # rRNA removal stage
        metapro_stage_obj.mp_rRNA_filter()

    elif(tutorial_mode == "repop"):
        # Duplicate repopulation
        metapro_stage_obj.mp_repop()

    elif(tutorial_mode == "contigs"):
        # Assemble contigs
        metapro_stage_obj.mp_assemble()    

    elif(tutorial_mode == "GA"):
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

    elif(tutorial_mode == "TA"):
        # Taxonomic annotation
        metapro_stage_obj.mp_TA()

    elif(tutorial_mode == "EC"):
        # Detect EC annotation
        metapro_stage_obj.mp_EC()
    elif(tutorial_mode == "output"):
        #check the contig state
        metapro_stage_obj.mp_contig_statecheck()
        # RPKM Table and Cytoscape Network
        metapro_stage_obj.mp_output()

    
    
if __name__ == "__main__":
    print("METAPRO metatranscriptomic analysis pipeline")
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline"
                                        "Version 2.0.2 © 2022")

    parser.add_argument("-c", "--config",   type=str,   help="Path to the configureation file")
    parser.add_argument("-1", "--pair1",    type=str,   help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2",    type=str,   help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single",   type=str,   help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-con", "--contig",   type=str,   help="Tutorial use only: Path to the file containing the contig reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, required=True, help="Path of the folder for the output of the pipeline")
    parser.add_argument("-t", "--num_threads", type=int, help="Maximum number of threads used by the pipeline")
    parser.add_argument("--nhost", "--no-host", action='store_true', help="Skip the host read removal step of the pipeline")
    parser.add_argument("--verbose_mode", type=str, help = "Decide how to handle the interim files, Compress them, or leave them alone.  Values are: keep, compress, quiet")
    parser.add_argument("--tutorial", type = str, help = "tutorial operating mode for MetaPro")
    
    args = parser.parse_args()
    
    config_file     = args.config if args.config else ""
    contig          = args.contig if args.contig else "None"
    pair_1          = args.pair1 if args.pair1 else ""
    pair_2          = args.pair2 if args.pair2 else ""
    single          = args.single if args.single else ""
    output_folder   = args.output_folder
    num_threads     = args.num_threads if args.num_threads else 0
    no_host         = args.nhost if args.nhost else False
    verbose_mode    = args.verbose_mode if args.verbose_mode else "quiet"
    tutorial_mode   = args.tutorial if args.tutorial else "none"

    if(tutorial_mode == "none"):
        if (args.pair1 and not args.pair2) or (args.pair2 and not args.pair1):
            print("You must specify both forward and reverse reads for a paired-end run")
            sys.exit()
        elif args.single and (args.pair1 or args.pair2):
            print("You cannot specify both paired-end and single-end reads in a single run.")
            sys.exit()

    if not (os.path.exists(output_folder)):
        print("output folder does not exist.  Now building directory.")
        os.makedirs(output_folder)
    os.chdir(output_folder)

    config = ConfigParser(interpolation = ExtendedInterpolation())
    if args.config:
        config.read(config_file)
        if not args.pair1 and not args.pair2 and not args.single:
            pair_1 = config["Sequences"]["pair1"] if config["Sequences"]["pair1"] else ""
            pair_2 = config["Sequences"]["pair2"] if config["Sequences"]["pair2"] else ""
            single = config["Sequences"]["single"] if config["Sequences"]["single"] else ""

    if pair_1 == "" and pair_2 == "" and single == "":
        print("You must specify paired-end or single-end reads as input for the pipeline.")
        sys.exit()

    args_pack = dict()
    args_pack["no_host"] = no_host
    args_pack["verbose_mode"] = verbose_mode
    
    print("=====================================")
    print("no-host:", no_host)
    print("verbose_mode:", verbose_mode)

    if (tutorial_mode != "none"):
        print("working in tutorial mode:", tutorial_mode)
        tutorial_main(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode)
    
    else:
        main(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode)
