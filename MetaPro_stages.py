#!/usr/bin/env python
from re import split
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import MetaPro_utilities as mpu
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q

#stores code for stage-launch.
#makes for a neat package/capsule

class mp_stage:
    def __init__ (self, config_path, pair_1_path="None", pair_2_path="None", single_path="None", contig_path="None", output_folder_path="None", args_pack="None", tutorial_mode_string = None, test_mode = "no"):
        #make our util obj
        #refresher: self -> instance var.  not self: class var (shared among class obj instances)
        
        #---------------------------------------------------------
        #Operational flags and state-recorders
        
        
        self.test_mode = test_mode
        self.tutorial_string = tutorial_mode_string
        
        self.output_folder_path = output_folder_path
        
        self.paths = mpp.tool_path_obj(config_path, self.output_folder_path)
        self.mp_util = mpu.mp_util(self.output_folder_path, self.paths)
        self.GA_DB_mode = self.paths.GA_DB_mode
        self.segmented_chocophlan_flag = True
        if(self.paths.DNA_DB.endswith(".fasta")):
            self.segmented_chocophlan_flag = False
        self.no_host = args_pack["no_host"]
        self.verbose_mode = args_pack["verbose_mode"]

        self.config_path = config_path
        
        self.pair_1_path = None if (pair_1_path == "None") else pair_1_path
        self.pair_2_path = None if (pair_2_path == "None") else pair_2_path
        self.single_path = None if (single_path == "None") else single_path
        self.contig_path = None if (contig_path == "None") else contig_path #tutorial/single-shot use
        self.quality_encoding = ""
        self.read_mode = "none"
        if not self.single_path == None:
            self.read_mode = "single"
            self.quality_encoding = self.mp_util.determine_encoding(single_path)
            print("ENCODING USED:", self.quality_encoding)
            print("OPERATING IN SINGLE-ENDED MODE")
        else:
            self.read_mode = "paired"
            self.quality_encoding = self.mp_util.determine_encoding(pair_1_path)
            print("ENCODING USED:", self.quality_encoding)
            print("OPERATING IN PAIRED-MODE")
        
        self.debug_stop_flag = self.paths.debug_stop_flag
        



        #timing vars
        self.start_time                     = time.time()
        self.end_time                       = 0
        self.quality_start                  = 0
        self.quality_end                    = 0
        self.cleanup_quality_start          = 0
        self.cleanup_quality_end            = 0
        
        self.host_start                     = 0
        self.host_end                       = 0
        self.cleanup_host_start             = 0
        self.cleanup_host_end               = 0
        
        self.vector_start                   = 0
        self.vector_end                     = 0
        self.cleanup_vector_start           = 0
        self.cleanup_vector_end             = 0
        
        self.rRNA_filter_start              = 0  
        self.rRNA_filter_end                = 0
        self.cleanup_rRNA_filter_start      = 0
        self.cleanup_rRNA_filter_end        = 0   
        
        self.repop_start                    = 0
        self.repop_end                      = 0
        self.cleanup_repop_start            = 0
        self.cleanup_repop_end              = 0
        
        self.assemble_contigs_start         = 0
        self.assemble_contigs_end           = 0
        self.cleanup_assemble_contigs_start = 0
        self.cleanup_assemble_contigs_end   = 0
        
        self.destroy_contigs_start          = 0
        self.destroy_contigs_end            = 0
        self.cleanup_destroy_contigs_start  = 0 
        self.cleanup_destroy_contigs_end    = 0
        
        self.GA_BWA_start                   = 0
        self.GA_BWA_end                     = 0
        self.cleanup_GA_BWA_start           = 0
        self.cleanup_GA_BWA_end             = 0
        
        self.GA_BLAT_start                  = 0
        self.GA_BLAT_end                    = 0
        self.cleanup_GA_BLAT_start          = 0
        self.cleanup_GA_BLAT_end            = 0
        
        self.GA_DIAMOND_start               = 0
        self.GA_DIAMOND_end                 = 0
        self.cleanup_GA_DIAMOND_start       = 0
        self.cleanup_GA_DIAMOND_end         = 0
        
        self.TA_start                       = 0
        self.TA_end                         = 0
        self.cleanup_TA_start               = 0
        self.cleanup_TA_end                 = 0
        
        self.EC_start                       = 0
        self.EC_end                         = 0

        self.EC_DETECT_start                = 0  
        self.EC_DETECT_end                  = 0
        
        self.EC_PRIAM_start                 = 0
        self.EC_PRIAM_end                   = 0
        
        self.EC_DIAMOND_start               = 0
        self.EC_DIAMOND_end                 = 0
        
        self.cleanup_EC_start               = 0
        self.cleanup_EC_end                 = 0
        
        self.Cytoscape_start                = 0
        self.Cytoscape_end                  = 0
        self.cleanup_cytoscape_start        = 0
        self.cleanup_cytoscape_end          = 0
        

        
            
             
                
        mp_store = []  # stores the multiprocessing processes

        # Creates our command object, for creating shellscripts.

        #if self.read_mode == "single":
        #    self.commands = mpcom.mt_pipe_commands(self.no_host, Config_path=config_path, Quality_score=self.quality_encoding, tutorial_keyword=self.tutorial_string, sequence_path_1=self.pair_1_path, sequence_path_2=self.pair_2_path, sequence_single=self.single_path, sequence_contigs = self.contig_path)
        #elif self.read_mode == "paired":
        self.commands = mpcom.mt_pipe_commands(self.no_host, self.paths, Quality_score=self.quality_encoding, tutorial_keyword=self.tutorial_string, sequence_path_1=self.pair_1_path, sequence_path_2=self.pair_2_path, sequence_single=self.single_path, sequence_contigs = self.contig_path)
    
        #special var to track contigs.
        self.contigs_present = True
        
        
    
    #--------------------------------------------------------------------------------
    # helpers
    def mp_contig_statecheck(self):
        
        if(os.path.exists(self.paths.spades_done_file)):
            if(os.path.exists(self.paths.spades_transcripts_file)):
                self.contigs_present = True
                print(dt.today(), "contigs present")
            else:
                self.contigs_present = False
                print(dt.today(), "contigs not present. but that's ok")
        else:
            sys.exit(dt.today(), "contig assembly incomplete. run the pipeline in --tutorial contigs to continue")
            
    def import_lib_names(self, names_file):
    #import the list of libs from GA_pre_scan
        files_list = []
        with open(names_file, "r") as names_in:
            for line in names_in:
                cleaned_line = line.strip("\n")
                if("can't" in line):
                    continue
                #print("line:", cleaned_line)
                src_path = cleaned_line.split("|")[3]
                if(os.path.exists(src_path)):
                
                #not complete. the exist path needs to point to the DNA_db folder
                    files_list.append(src_path)
                    #print("src:", src_path)
                    #time.sleep(1)
                    
        
        return files_list
    
    def debug_stop_check(self, stop_signal):
        if(self.debug_stop_flag == stop_signal):
            exit_string = "stoppped after: " + stop_signal
            sys.exit(exit_string)
        else:
            print(dt.today(), "continuing from:", stop_signal)
            

    #--------------------------------------------------------------------------------------------------------------
    # main calls
    def mp_quality_filter(self):
        self.quality_start = time.time()
        command_list = self.commands.create_quality_control_command()
        self.cleanup_quality_start, self.cleanup_quality_end = self.mp_util.launch_stage_simple(self.paths.qc_label, self.paths.qc_top_path, self.commands, command_list, self.paths.keep_all, self.paths.keep_quality)
        self.quality_end = time.time()
        print("quality filter:", '%1.1f' % (self.quality_end - self.quality_start - (self.cleanup_quality_end - self.cleanup_quality_start)), "s")
        print("quality filter cleanup:", '%1.1f' %(self.cleanup_quality_end - self.cleanup_quality_start), "s")
        self.debug_stop_check(self.paths.qc_label)
        

        
    def mp_host_filter(self):
        if not self.no_host:
            self.host_start = time.time()
            
            if not (self.paths.check_if_indexed(self.paths.Host_DB)):
                command_list = self.commands.create_host_index_command(self.paths.Host_DB_index_marker)
                

            else:
                print("skipping Host index")
            command_list = self.commands.create_host_filter_command(self.paths.host_final_marker)
            if(self.test_mode == "yes"):
                for item in command_list:
                    print(item)
            else:        
                self.cleanup_host_start, self.cleanup_host_end = self.mp_util.launch_stage_with_cleanup(self.commands, command_list, self.paths.host_final_marker, self.paths.host_data_path, self.paths.host_filter_label, self.paths.keep_all, self.paths.keep_host)
                self.host_end = time.time()
                print("host filter:", '%1.1f' % (self.host_end - self.host_start - (self.cleanup_host_end - self.cleanup_host_start)), "s")
                print("host filter cleanup:", '%1.1f' %(self.cleanup_host_end - self.cleanup_host_start),"s")
                self.debug_stop_check(self.paths.host_filter_label)
                
    def mp_vector_filter(self):
        self.vector_start = time.time()
        
        if not (self.paths.check_if_indexed(self.paths.Vector_DB)):
            command_list = self.commands.create_vector_filter_index_command(self.paths.Vector_DB_index_marker)
            self.mp_util.launch_only_with_mp_store(self.commands, command_list)
            
            #get dep args from quality filter
            #if not check_where_resume(vector_path, None, self.quality_path):
        command_list = self.commands.create_vector_filter_command(self.paths.vector_final_marker, self.no_host)
        self.cleanup_vector_start, self.cleanup_vector_end = self.mp_util.launch_stage_with_cleanup(self.commands, command_list, self.paths.vector_final_marker, self.paths.vector_data_path, self.paths.vector_filter_label, self.paths.keep_all, self.paths.keep_vector)

        
        self.vector_end = time.time()
        print("vector filter:", '%1.1f' % (self.vector_end - self.vector_start - (self.cleanup_vector_end - self.cleanup_vector_start)), "s")
        print("vector filter cleanup:", '%1.1f' % (self.cleanup_vector_end - self.cleanup_vector_start), "s")
        self.debug_stop_check(self.paths.vector_filter_label)

    def mp_rRNA_filter(self):
        self.rRNA_filter_start = time.time()
        self.rRNA_filter_end = time.time()
        self.cleanup_rRNA_filter_start = time.time()
        self.cleanup_rRNA_filter_end = time.time()

        #if not check_where_resume(self.rRNA_filter_path, None, self.vector_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.paths.rRNA_filter_label): 
            
            #split and convert fastq -> fasta
            marker_path = os.path.join(self.paths.rRNA_jobs_path, self.paths.rRNA_split_marker)
            command_list = self.commands.create_rRNA_filter_split_command(marker_path)
            if(not os.path.exists(marker_path)):
                self.mp_util.launch_only_with_mp_store(self.commands, command_list)
            else:
                print("skipping rRNA initial split")
            self.mp_util.wait_for_mp_store()
            

            #-------------------------------------------------------------------------------------------------
            # BARRNAP

            for split_fasta in os.listdir(self.paths.rRNA_p1_fa_path):
                
                file_name = split_fasta.split(".")[0]
                fasta_path = os.path.join(self.paths.rRNA_p1_fa_path, split_fasta)
                barrnap_out = os.path.join(self.paths.rRNA_p1_bar_path, file_name + ".barrnap_out")
                mRNA_out = os.path.join(self.paths.rRNA_p1_bar_path, file_name + ".fasta")
                marker_path = os.path.join(self.paths.rRNA_jobs_path, self.paths.rRNA_barrnap_marker + file_name)
                command_list = self.commands.create_rRNA_filter_barrnap_command(fasta_path, barrnap_out, mRNA_out, marker_path)
                #print(command_list)
                
                if not os.path.exists(marker_path):
                    print("running barrnap:", file_name)
                    self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                    #job_path = os.path.join(self.paths.rRNA_exe_path, file_name + "_barrnap.sh")
                    #self.mp_util.launch_and_create_v2(job_path, self.commands, command_list)
                else:
                    print("skipping barrnap:", file_name)
                    
            for split_fasta in os.listdir(self.paths.rRNA_p2_fa_path):
                file_name = split_fasta.split(".")[0]
                fasta_segment = os.path.join(self.paths.rRNA_p2_fa_path, split_fasta)
                barrnap_out = os.path.join(self.paths.rRNA_p2_bar_path, file_name + ".barrnap_out")
                mRNA_out = os.path.join(self.paths.rRNA_p2_bar_path, file_name + ".fasta")
                marker_path = os.path.join(self.paths.rRNA_jobs_path, self.paths.rRNA_barrnap_marker + file_name)
                command_list = self.commands.create_rRNA_filter_barrnap_command(fasta_segment, barrnap_out, mRNA_out, marker_path)
                if not os.path.exists(marker_path):
                    print("running barrnap:", file_name)
                    self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                    #job_path = os.path.join(self.paths.rRNA_exe_path, file_name + "_barrnap.sh")
                    #self.mp_util.launch_and_create_v2(job_path, self.commands, command_list)
                else:
                    print("skipping barrnap:", file_name)

            for split_fasta in os.listdir(self.paths.rRNA_s_fa_path):
                file_name = split_fasta.split(".")[0]
                barrnap_out = os.path.join(self.paths.rRNA_s_bar_path, file_name + ".barrnap_out")
                fasta_segment = os.path.join(self.paths.rRNA_s_fa_path, split_fasta)
                mRNA_out = os.path.join(self.paths.rRNA_s_bar_path, file_name + ".fasta")
                marker_path = os.path.join(self.paths.rRNA_jobs_path, self.paths.rRNA_barrnap_marker + file_name)
                command_list = self.commands.create_rRNA_filter_barrnap_command(fasta_segment, barrnap_out, mRNA_out, marker_path)
                if not os.path.exists(marker_path):
                    print("running barrnap:", file_name)
                    self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                else:
                    print("skipping barrnap:", file_name)
                
            self.mp_util.wait_for_mp_store()
            
                
            #----------------------------------------------------------------------------
            # INFERNAL + final splitting

            for split_fasta in os.listdir(self.paths.rRNA_p1_bar_path):
                if(split_fasta.endswith(".fasta")):
                    file_name = split_fasta.split(".")[0]
                    fasta_segment = os.path.join(self.paths.rRNA_p1_bar_path, split_fasta)
                    infernal_out_file = os.path.join(self.paths.rRNA_p1_inf_path, file_name + ".inf_out")
                    marker_path = os.path.join(self.paths.rRNA_jobs_path, file_name + self.paths.rRNA_inf_marker)
                    command_list = self.commands.create_rRNA_filter_infernal_command(fasta_segment, infernal_out_file)
                    if not os.path.exists(marker_path):
                        print("running inf:", file_name)
                        print("marker path:", marker_path)
                        time.sleep(1)
                        self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                    else:
                        print("skipping inf:", file_name) 

            for split_fasta in os.listdir(self.paths.rRNA_p2_bar_path):
                if(split_fasta.endswith(".fasta")):
                    file_name = split_fasta.split(".")[0]
                    fasta_segment = os.path.join(self.paths.rRNA_p2_bar_path, split_fasta)
                    infernal_out_file = os.path.join(self.paths.rRNA_p2_inf_path, file_name + ".inf_out")
                    marker_path = os.path.join(self.paths.rRNA_jobs_path, file_name + self.paths.rRNA_inf_marker)
                    command_list = self.commands.create_rRNA_filter_infernal_command(fasta_segment, infernal_out_file)
                    if not os.path.exists(marker_path):
                        print("running inf:", file_name)
                        print("marker path:", marker_path)
                        time.sleep(1)
                        self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                    else:
                        print("skipping inf:", file_name)
                    
            for split_fasta in os.listdir(self.paths.rRNA_s_bar_path):   
                if(split_fasta.endswith(".fasta")):
                    file_name = split_fasta.split(".")[0]
                    fasta_segment = os.path.join(self.paths.rRNA_s_bar_path, split_fasta)
                    infernal_out_file = os.path.join(self.paths.rRNA_s_inf_path, file_name + ".inf_out")
                    marker_path = os.path.join(self.paths.rRNA_jobs_path, file_name + self.paths.rRNA_inf_marker)
                    
                    command_list = self.commands.create_rRNA_filter_infernal_command(fasta_segment, infernal_out_file)
                    if not os.path.exists(marker_path):
                        print("running inf:", file_name)
                        print("marker path:", marker_path)
                        time.sleep(1)
                        self.mp_util.launch_only_with_mp_store(self.commands, command_list)
                    else:
                        print("skipping inf:", file_name)
            self.mp_util.wait_for_mp_store()

            #--------------------------------------------------------------------
            # rRNA cleanup
            command_list = self.commands.create_rRNA_cleanup_command("paired", self.paths.rRNA_p_final_marker)
            print("command list:", command_list)
            if not os.path.exists(self.paths.rRNA_p_final_marker):
                print("running rRNA-final: paired")
                self.mp_util.launch_only_with_mp_store(self.commands, command_list)
            else:
                print("skipping rRNA-final: paired" )

        
            
            command_list = self.commands.create_rRNA_cleanup_command("single", self.paths.rRNA_s_final_marker)
            print("command list:", command_list)
            if not os.path.exists(self.paths.rRNA_s_final_marker):
                print("running rRNA-final: single")
                self.mp_util.launch_only_with_mp_store(self.commands, command_list)
            else:
                print("skipping rRNA-final: single")

            self.mp_util.wait_for_mp_store()

            rRNA_complete_flag = False
            if(len(os.listdir(self.paths.rRNA_p1_inf_mRNA_path)) > 0) and (len(os.listdir(self.paths.rRNA_s_inf_mRNA_path))> 0):
                if(os.path.exists(self.paths.rRNA_p_final_marker) and (os.path.exists(self.paths.rRNA_s_final_marker))):
                    rRNA_complete_flag = True
            else:
                if(os.path.exists(self.paths.rRNA_s_final_marker)):
                    rRNA_complete_flag = True    

            if(rRNA_complete_flag):            
                self.debug_stop_check(self.paths.rRNA_filter_label)
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.rRNA_filter_label)
                self.cleanup_rRNA_filter_start = time.time()
                self.mp_util.clean_or_compress(self.paths.rRNA_data_path, self.paths.keep_all, self.paths.keep_rRNA)
                self.cleanup_rRNA_filter_end = time.time()


        print("rRNA filter:", '%1.1f' % (self.rRNA_filter_end - self.rRNA_filter_start - (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start)), "s")
        print("rRNA filter cleanup:", '%1.1f' % (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start), "s")
            

    def mp_repop(self):
        
        self.repop_start = time.time()
        #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.paths.repop_label):
            job_name = self.paths.repop_label
            command_list = self.commands.create_repop_command_v2_step_1(self.paths.repop_label, self.paths.qc_label, self.paths.rRNA_filter_label)
            self.mp_util.subdivide_and_launch(self.paths.repop_job_delay, self.paths.repop_mem_threshold, self.paths.repop_job_limit, self.paths.repop_label, job_name, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            
            if(self.read_mode == "paired"):
                job_name = self.paths.repop_label
                command_list = self.commands.create_repop_command_v2_step_2(self.paths.repop_label, self.paths.qc_label, self.paths.rRNA_filter_label)
                self.mp_util.subdivide_and_launch(self.paths.repop_job_delay, self.paths.repop_mem_threshold, self.paths.repop_job_limit, self.paths.repop_label, job_name, self.commands, command_list)
                self.mp_util.wait_for_mp_store()
            
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.repop_label)
            
            self.cleanup_repop_start = time.time()
            self.mp_util.clean_or_compress(self.paths.repop_data_path, self.paths.keep_all, self.paths.keep_repop)
            self.cleanup_repop_end = time.time()
            
        self.repop_end = time.time()
        print("repop:", '%1.1f' % (self.repop_end - self.repop_start - (self.cleanup_repop_end - self.cleanup_repop_start)), "s")
        print("repop cleanup:", '%1.1f' % (self.cleanup_repop_end - self.cleanup_repop_start), "s")
        self.debug_stop_check(self.paths.repop_label)

    def mp_assemble(self):
        self.assemble_contigs_start = time.time()
        
        #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
        mgm_gene_report = os.path.join(self.paths.contigs_data_path, "1_mgm", "gene_report.txt")
        
        spades_done_file = os.path.join(self.paths.contigs_data_path, "0_spades", "pipeline_state", "stage_7_terminate")
        spades_transcript_file = os.path.join(self.paths.contigs_data_path, "0_spades", "transcripts.fasta")

        mgm_fail_flag = True
        spades_fail_flag = True

        if self.mp_util.check_bypass_log(self.output_folder_path, self.paths.assemble_contigs_label):
            job_name = self.paths.assemble_contigs_label
            command_list = self.commands.create_assemble_contigs_command(self.paths.assemble_contigs_label, self.paths.repop_label)
            self.mp_util.launch_and_create_simple(self.paths.assemble_contigs_label, job_name, self.commands, command_list)
            
            if(os.path.exists(spades_done_file)):
                if(os.path.exists(spades_transcript_file)):
                    spades_fail_flag = False
                    print(dt.today(), "SPADes OK")
                else:
                    spades_fail_flag = True
                    print(dt.today(), "SPADes ran, but did not create contigs")
            else:
                err_msg = str(dt.today()) +  " SPADes did not run. this is not normal. Contact admin immediately"
                sys.exit(err_msg)

            if(os.path.exists(mgm_gene_report)):
                if(os.path.getsize(mgm_gene_report) > 0):
                    mgm_fail_flag = False
                    print(dt.today(), "MGM OK")
                else:
                    mgm_fail_flag = True
                    print(dt.today(), "MGM produced an empty gene report. this is not normal")
                    sys.exit("MGM_empty_report")
            else:
                mgm_fail_flag = True
                print(dt.today(), "MGM did not produce a report. likely it didn't run")

            if(spades_fail_flag and mgm_fail_flag):        
                print(dt.today(), "moving contig files to compensate")
                bypass_contig_map_path = os.path.join(self.paths.contigs_final_path, "contig_map.tsv")
                bypass_contig_path = os.path.join(self.paths.contigs_final_path, "contigs.fasta")
                s_src_path = os.path.join(self.paths.rRNA_final_mRNA_path, "singletons.fastq")
                p1_src_path = os.path.join(self.paths.rRNA_final_mRNA_path, "pair_1.fastq")
                p2_src_path = os.path.join(self.paths.rRNA_final_mRNA_path, "pair_2.fastq")
                s_dest_path = os.path.join(self.paths.contigs_final_path, "singletons.fastq")
                p1_dest_path = os.path.join(self.paths.contigs_final_path, "pair_1.fastq")
                p2_dest_path = os.path.join(self.paths.contigs_final_path, "pair_2.fastq")
                make_map = open(bypass_contig_map_path, "w")
                make_contig = open(bypass_contig_path, "w")
                shutil.copyfile(s_src_path, s_dest_path)
                shutil.copyfile(p1_src_path, p1_dest_path)
                shutil.copyfile(p2_src_path, p2_dest_path)
                self.contigs_present = False
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.assemble_contigs_label)

            elif(spades_fail_flag and (not mgm_fail_flag)):
                print(dt.today(), "SPADes fails but MGM runs anyway: this shouldn't happen. contact admin")
                sys.exit("SPADES_fail_MGM_ok")
            elif((not spades_fail_flag) and mgm_fail_flag):
                print(dt.today(), "SPADes ran fine, but MGM failed. Check your MetaGeneMark license")
                sys.exit("SPADES_ok_MGM_fail")
            else:
                self.contigs_present = True
                print(dt.today(), "Assemble-contigs pass")
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.assemble_contigs_label)
            
            self.cleanup_assemble_contigs_start = time.time()
            self.mp_util.clean_or_compress(self.paths.contigs_data_path, self.paths.keep_all, self.paths.keep_assemble_contigs)
            self.cleanup_assemble_contigs_end = time.time()
        
        else:
            if(os.path.exists(mgm_gene_report)):
                if(os.path.getsize(mgm_gene_report) == 0):
                    sys.exit("MGM did not run. gene report is 0")
                else:
                    print(dt.today(), "MGM OK. contigs present")   
                    self.contigs_present = True
            else:
                if(os.path.exists(spades_done_file)):
                    print(dt.today(), "No contigs were assembled.")
                    self.contigs_present = False
                
        self.assemble_contigs_end = time.time()
        print("assemble contigs:", '%1.1f' % (self.assemble_contigs_end - self.assemble_contigs_start - (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start)), "s")    
        print("assemble contigs cleanup:", '%1.1f' % (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start), "s")
        
        self.debug_stop_check(self.paths.assemble_contigs_label)
    
    def mp_GA_pre_scan(self):
        #scans the mRNA with a TA scanner to pick out a taxa trend.
        if self.mp_util.check_bypass_log(self.output_folder_path, self.paths.GA_pre_scan_label):
            marker_path_list = []
            #----------------------------------------------------------------------
            #kaiju on reads
            sections = ["singletons"]
            if self.read_mode == "paired":
                sections.extend(["paired"])
            if(self.contigs_present):
                sections.extend(["contigs"])    
            
            for section in sections:
                marker_file = "mp_ta_kraken2_" + section
                marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_GA_pre_scan_kraken2_command(section, marker_path)
                    self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)        
            
            

            #---------------------------------------------------------------------
            #centrifuge on reads
            
            sections_list = ["reads"]
            if(self.contigs_present):
                sections_list.extend(["contigs"])
                
            for read_cat in sections_list:
                marker_file = "mp_ta_centrifuge_" + read_cat
                marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
                print(marker_path)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_GA_pre_scan_centrifuge_command(read_cat, marker_path)
                    self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.Centrifuge_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
            
            self.mp_util.wait_for_mp_store()
            
            #------------------------------------------
            #merge kraken + centrifuge into single file
            
            marker_file = "TA_kraken2_pp"
            marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_kraken2_pp_command(marker_path)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            
            
            marker_path_list = []
            marker_file = "TA_centrifuge_pp"
            marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_centrifuge_pp_command(marker_path)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()

            
            
            #----------------------------------------------------
            #
            
            marker_file = "TA_wevote_combine"
            marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_wevote_combine_command(marker_file)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
                print(dt.today(), "running:", marker_file)
            self.mp_util.wait_for_mp_store()
            
            
            #-------------------------------------------------------
            #collect the wevote results, get the unique taxa, and match them to a class-level taxa
            
            marker_file = "ga_collect_db"
            marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_command(marker_path)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
                print(dt.today(), "running:", marker_file)
            self.mp_util.wait_for_mp_store()
            
            
            marker_file = "ga_assemble_db"
            marker_path = os.path.join(self.paths.GA_pre_scan_jobs_path, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_assemble_lib_command(marker_path)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.paths.GA_pre_scan_label, marker_file, self.commands, command_list)
                print(dt.today(), "running:", marker_file)
            self.mp_util.wait_for_mp_store()
            
            
            #---------------------------------------------------------
            
            
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.GA_pre_scan_label)
            
        self.debug_stop_check(self.paths.GA_pre_scan_label)
    
    def mp_GA_split(self):
        #separating GA split-data from GA_BWA for a few reasons:
        #1) so the pipe has the option to not split all the time
        #2) modularity
        if self.mp_util.check_bypass_log(self.output_folder_path, self.paths.GA_split_label):
            marker_path_list = []
            if(self.contigs_present):
                print(dt.today(), "splitting contigs")
                marker_file = "GA_split_fasta_contigs"
                marker_path = os.path.join(self.paths.GA_split_jobs_path, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping", marker_file)
                else:
                    job_name = "GA_prep_split_contigs"
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_split_ga_fasta_data_command("contigs", marker_path)
                    self.mp_util.launch_and_create_with_mp_store(self.paths.GA_split_label, job_name, self.commands, command_list)
            else:
                print(dt.today(), "no contigs present. skipping split")
            
            sections = ["singletons"]
            if(self.read_mode == "paired"):
                sections.extend(["pair_1", "pair_2"])
            for section in sections: 
                marker_file = "GA_split_fastq_" + section
                marker_path = os.path.join(self.paths.GA_split_jobs_path, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    job_name = "GA_prep_split_" + section
                    command_list = self.commands.create_split_ga_fastq_data_command(section, marker_path)
                    self.mp_util.launch_and_create_with_mp_store(self.paths.GA_split_label, job_name, self.commands, command_list)
            self.mp_util.wait_for_mp_store()

            
            #self.mp_util.write_to_bypass_log(self.output_folder_path, self.paths.GA_split_label)
            
        self.debug_stop_check(self.paths.GA_split_label)
        
    def mp_GA_lib_check(self):
        print(dt.today(), "Running GA lib check")
        
        if(self.paths.DNA_DB_mode == "chocophlan"):
            self.paths.DNA_DB = os.path.join(self.GA_pre_scan_path, "final_results")
        self.paths.check_bwa_valid(self.paths.DNA_DB)
        self.paths.check_blat_valid(self.paths.DNA_DB)
        
    def mp_GA_BWA(self):
        self.GA_BWA_start = time.time()
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BWA_label):
            marker_path_list = []
            
            
            if not self.mp_util.check_where_resume(self.GA_BWA_path, None, self.GA_split_path):
            
                #-------------------------------------------------------------------------
                sections = ["singletons"]
                if self.read_mode == "paired":
                    sections.extend(["pair_1", "pair_2"])
                if(self.contigs_present):
                    
                    sections.extend(["contigs"])
                
                for section in sections:
                    for split_sample in os.listdir(os.path.join(self.GA_split_path, "final_results", section)):
                        full_sample_path = os.path.join(os.path.join(self.GA_split_path, "final_results",section, split_sample))
                        print("split sample:", full_sample_path)
                        file_tag = os.path.basename(split_sample)
                        file_tag = os.path.splitext(file_tag)[0]
                        ref_path = self.GA_pre_scan_final_path #self.paths.DNA_DB
                            
                        command_list = ""
                        if (ref_path.endswith(".fasta")):
                            ref_tag = os.path.basename(ref_path)
                            ref_tag = ref_tag.strip(".fasta")
                        
                            file_tag = file_tag + "_" + ref_tag
                            job_name = "BWA" + "_" + file_tag
                            marker_file = file_tag + "_bwa"
                            marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
                        #this checker assumes that BWA only exports a file when it's finished running
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                            
                                #aug 10, 2021: new bigger chocophlan (from humann3) is in segments because we can't index it as a whole.  
                                #if the DB is still an old version, the tag should just say "chocophlan".  otherwise, it will say the chocophlan chunk name
                                
                                command_list = self.commands.create_BWA_annotate_command_v2(self.GA_BWA_label, ref_path, ref_tag, full_sample_path, marker_file)
                                #self.mp_util.launch_and_create_with_hold(self.BWA_mem_threshold, self.BWA_job_limit, self.BWA_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)
                                self.mp_util.launch_and_create_with_mem_footprint(self.BWA_mem_footprint, self.BWA_job_limit, self.GA_BWA_label, job_name, self.commands, command_list)
                                
                        else:
                            split_db = os.listdir(ref_path)
                            for db_segments in split_db:
                                if(db_segments.endswith(".fasta")):
                                    segment_ref_path = os.path.join(ref_path, db_segments)
                                    ref_tag = db_segments.strip(".fasta")
                                    segment_file_tag = file_tag + "_" + ref_tag
                                    job_name = "BWA" + "_" + segment_file_tag
                                    marker_file = segment_file_tag + "_bwa"
                                    marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
                                    
                                    if(os.path.exists(marker_path)):
                                        print(dt.today(), "skipping:", marker_file)
                                        continue
                                    else:
                                        marker_path_list.append(marker_path)
                                        command_list = self.commands.create_BWA_annotate_command_v2(self.GA_BWA_label, segment_ref_path, ref_tag, full_sample_path, marker_file)
                                        #footprint doesn't apply to a single-file BWA DB
                                        #self.mp_util.launch_and_create_with_hold(self.BWA_mem_threshold, self.BWA_job_limit, self.BWA_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)
                                        self.mp_util.launch_and_create_with_mem_footprint(self.BWA_mem_footprint, self.BWA_job_limit, self.GA_BWA_label, job_name, self.commands, command_list)
                                        

                print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
                self.mp_util.wait_for_mp_store()
                final_checklist = os.path.join(self.GA_BWA_path, "GA_BWA.txt")
                self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BWA_label)
    
        self.debug_stop_check(self.GA_BWA_label)
        
    def mp_GA_BWA_pp(self):                
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BWA_pp_label):
            marker_path_list = []
            sections = ["singletons"]
            if self.read_mode == "paired":
                sections.extend(["pair_1", "pair_2"])
                
            if(self.contigs_present):
                sections.extend(["contigs"])
                
                
            for section in sections:
                for split_sample in os.listdir(os.path.join(self.GA_split_path, "final_results", section)):
                    full_sample_path = os.path.join(os.path.join(self.GA_split_path, "final_results",section, split_sample))
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    
                    
                    ref_path = self.paths.DNA_DB
                    #chocophlan in many mutiple segments
                    if (ref_path.endswith(".fasta")):
                        ref_tag = os.path.basename(ref_path)
                        ref_tag = ref_tag.strip(".fasta")
                
                    
                        job_name = "BWA_pp" + "_" + file_tag + "_" + ref_tag
                        marker_file = file_tag + "_" + ref_tag +  "_bwa_pp"
                        marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
                        
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "skipping:", marker_file)
                            continue
                        else:
                            marker_path_list.append(marker_path)
                            command_list = self.commands.create_BWA_pp_command_v2(self.GA_BWA_label, self.paths.assemble_contigs_label, ref_tag, ref_path, full_sample_path, marker_file)
                            self.mp_util.launch_and_create_with_hold(self.BWA_pp_mem_threshold, self.BWA_pp_job_limit, self.BWA_pp_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)
                            
                    else:
                        #chocophlan in chunks
                        split_db = os.listdir(ref_path)
                        for db_segments in split_db:
                            if(db_segments.endswith(".fasta")):
                                segment_ref_path = os.path.join(ref_path, db_segments)
                                ref_tag = db_segments.strip(".fasta")
                                job_name = "BWA_pp" + "_" + file_tag + "_" + ref_tag
                                marker_file = file_tag + "_" + ref_tag + "_bwa_pp"
                                marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
                                
                                if(os.path.exists(marker_path)):
                                    print(dt.today(), "skipping:", marker_file)
                                    continue
                                else:
                                    marker_path_list.append(marker_path)
                                    command_list = self.commands.create_BWA_pp_command_v2(self.GA_BWA_label, self.paths.assemble_contigs_label, ref_tag, segment_ref_path, full_sample_path, marker_file)
                                    #print(dt.today(), "segmented BWA:", command_list)
                                    #time.sleep(2)
                                    self.mp_util.launch_and_create_with_hold(self.BWA_pp_mem_threshold, self.BWA_pp_job_limit, self.BWA_pp_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)

                            
            print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
            self.mp_util.wait_for_mp_store()
            marker_file = "BWA_copy_contig_map"
            marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:   
                marker_path_list.append(marker_path)
                command_list = self.commands.create_BWA_copy_contig_map_command(self.GA_BWA_label, self.paths.assemble_contigs_label, marker_file)
                self.mp_util.launch_and_create_simple(self.GA_BWA_label, self.GA_BWA_label + "_copy_contig_map", self.commands, command_list)

            
            final_checklist = os.path.join(self.GA_BWA_path, "GA_BWA_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BWA_pp_label)
        
        self.debug_stop_check("GA_BWA_pp")
        
    def mp_GA_BWA_merge(self):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BWA_merge_label):
            #merge 
            marker_path_list = []
            sections = ["singletons"]
            if self.read_mode == "paired":
                sections.extend(["pair_1", "pair_2"])
            if(self.contigs_present):
                sections.extend(["contigs"])
            
            for section in sections:
                for split_sample in os.listdir(os.path.join(self.GA_split_path, "final_results", section)):
                    full_sample_path = os.path.join(os.path.join(self.GA_split_path, "final_results",section, split_sample))
                    print("split sample:", full_sample_path)
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    ref_path = self.paths.DNA_DB

                    marker_file = file_tag + "_merge_fasta"
                    marker_path = os.path.join(self.GA_BWA_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        job_name = "BWA_fasta_merge_" + file_tag
                        command_list = self.commands.create_merge_BWA_fasta_command(self.GA_BWA_label, full_sample_path, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.BWA_pp_mem_threshold, self.BWA_pp_job_limit, self.BWA_pp_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)

            print(dt.today(), "All BWA merge jobs have launched. waiting for sync")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.GA_BWA_path, "GA_BWA_merge.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BWA_merge_label)
            
     
        self.cleanup_GA_BWA_start = time.time()
        self.mp_util.delete_folder_simple(self.GA_BWA_jobs_folder)
        self.mp_util.clean_or_compress(self.GA_BWA_path, self.keep_all, self.keep_GA_BWA)
        
        self.cleanup_GA_BWA_end = time.time()
        self.GA_BWA_end = time.time()
        print("GA BWA:", '%1.1f' % (self.GA_BWA_end - self.GA_BWA_start - (self.cleanup_GA_BWA_end - self.cleanup_GA_BWA_start)), "s")
        print("GA BWA cleanup:", '%1.1f' % (self.cleanup_GA_BWA_end - self.cleanup_GA_BWA_start), "s")
        self.debug_stop_check("GA_BWA_merge")

    def mp_GA_BLAT(self):    
        # ------------------------------------------------
        # BLAT gene annotation
        GA_BLAT_start = time.time()
        print("new DNA DB path:", self.paths.DNA_DB)
        
        self.mp_util.make_folder(self.GA_BLAT_path)
        self.mp_util.make_folder(self.GA_BLAT_data_folder)
        self.mp_util.make_folder(self.GA_BLAT_jobs_folder)
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BLAT_label):
            marker_path_list = []
            
            for split_sample in os.listdir(os.path.join(self.GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(self.GA_BWA_path, "final_results", split_sample))

                    delay_count = 0
                    for fasta_db in os.listdir(self.paths.DNA_DB):
                        if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                            job_name = "BLAT_" + file_tag + "_" + fasta_db
                            marker_file = file_tag + "_blat_" + fasta_db
                            marker_path = os.path.join(self.GA_BLAT_jobs_folder, marker_file)
                            blatout_path = os.path.join(self.GA_BLAT_path, "data", "0_blat", file_tag + "_"+fasta_db + ".blatout")
                            blat_queue_package = blatout_path+"|" + marker_file
                            #ref_db = os.path.join(self.paths.DNA_DB, fasta_db)
                            
                            #This checker assume BLAT only exports a file when it's finished running
                            if(os.path.exists(marker_path)):
                                if(os.path.exists(blatout_path)):
                                    #recover from a restart.  there will be files that have been missed.  thread would have deleted the file
                                    print(dt.today(), "file still exists. adding to merge thread:", blatout_path)
                                    #blat_file_queue.put(blat_queue_package)
                                    
                                else:
                                    print(dt.today(), "file doesn't exist anymore already merged", blatout_path)
                                    
                                print(dt.today(), "BLAT job ran already, skipping:", marker_file)
                                #time.sleep(1)
                                continue
                                
                            else:
                                print(dt.today(), "RUNNING:", marker_file)
                                
                                marker_path_list.append(marker_path)
                                command_list = self.commands.create_BLAT_annotate_command_v2(self.GA_BLAT_label, full_sample_path, self.paths.DNA_DB, fasta_db, marker_file)
                                #self.mp_util.launch_only_with_hold(self.BLAT_mem_threshold, self.BLAT_job_limit, self.BLAT_job_delay, job_name, self.commands, command_list)
                                self.mp_util.launch_and_create_with_mem_footprint(self.BLAT_mem_footprint, self.BLAT_job_limit, self.GA_BLAT_label, job_name, self.commands, command_list)
            #---------------------------------------------------------------------------

            
            print(dt.today(), "final BLAT job removal. now waiting for mp-store flush")
            #note: this wait is disabled because we now have a separate thread.  it will hang if we enable it.
            print(dt.today(), "flushing mp_store")
            #self.mp_util.mp_store[:] = []        
            self.mp_util.wait_for_mp_store()
            print(dt.today(), "moving onto BLAT PP")
            final_checklist = os.path.join(self.GA_BLAT_path, "GA_BLAT.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BLAT_label)
        
        self.debug_stop_check(self.GA_BLAT_label)
        
    def mp_GA_BLAT_pp(self):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BLAT_pp_label):
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(self.GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    
                    ref_path = self.paths.DNA_DB#_Split  #the chocophlan chunks
                    if (ref_path.endswith(".fasta")):
                        #single chocophlan mode
                        job_name = "BLAT_" + file_tag + "_pp"
                        full_sample_path = os.path.join(os.path.join(self.GA_BWA_path, "final_results", split_sample))
                        marker_file = file_tag + "_blat_pp"
                        marker_path = os.path.join(self.GA_BLAT_jobs_folder, marker_file)
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "skipping:", marker_file)
                            continue
                        else:
                            marker_path_list.append(marker_path)
                            command_list = self.commands.create_BLAT_pp_command_v2(self.GA_BLAT_label, full_sample_path, self.GA_BWA_label, ref_path, marker_file)
                            #self.mp_util.launch_and_create_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, self.GA_BLAT_label, job_name, self.commands, command_list)
                            self.mp_util.launch_only_with_hold(self.BLAT_pp_mem_threshold, self.BLAT_pp_job_limit, self.BLAT_pp_job_delay, job_name, self.commands, command_list)
                            
                    else:
                        for fasta_db in os.listdir(ref_path):
                            if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                                #split chocophlan mode
                                #decode the chocophlan chunk, and supply the appropriate one.
                                #print("file tag:", file_tag.split("chocophlan"))
                                choco_chunk = fasta_db.split(".fasta")[0]
                                ref_file = os.path.join(ref_path, fasta_db)
                                #print("BLAT file tag:", file_tag, "|chunk: ", ref_file)
                                
                                job_name = "BLAT_" + file_tag + "_" + choco_chunk + "_pp"
                                full_sample_path = os.path.join(os.path.join(self.GA_BWA_path, "final_results", split_sample))
                                #print("query file:", full_sample_path)
                                marker_file = file_tag + "_" + choco_chunk + "_blat_pp"
                                print("MARKER FILE:", marker_file)
                                marker_path = os.path.join(self.GA_BLAT_jobs_folder, marker_file)
                                
                                if(os.path.exists(marker_path)):
                                    print(dt.today(), "skipping:", marker_file)
                                    continue
                                else:
                                    marker_path_list.append(marker_path)
                                    command_list = self.commands.create_BLAT_pp_command_v3(self.GA_BLAT_label, full_sample_path, self.GA_BWA_label, ref_file, marker_file)
                                    #print("Command list:", command_list)
                                    #time.sleep(10)
                                    #self.mp_util.launch_and_create_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, self.GA_BLAT_label, job_name, self.commands, command_list)
                                    self.mp_util.launch_only_with_hold(self.BLAT_pp_mem_threshold, self.BLAT_pp_job_limit, self.BLAT_pp_job_delay, job_name, self.commands, command_list)
                                #time.sleep(10)

                    
            print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            
            job_name = "GA_BLAT_copy_contigs"
            marker_file = "blat_copy_contig_map"
            marker_path = os.path.join(self.GA_BLAT_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_BLAT_copy_contig_map_command(self.GA_BLAT_label, self.GA_BWA_label, marker_file)
                self.mp_util.launch_and_create_simple(self.GA_BLAT_label, job_name, self.commands, command_list)
            final_checklist = os.path.join(self.GA_BLAT_path, "GA_BLAT_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BLAT_pp_label)
        
        self.debug_stop_check(self.GA_BLAT_pp_label)


    def mp_GA_BLAT_merge(self):
        # GA BLAT merge    
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BLAT_merge_label):
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(self.GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]

                    marker_file = "BLAT_merge_" + file_tag
                    marker_path = os.path.join(self.GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping: ", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_BLAT_merge_fasta_command(self.GA_BLAT_label, file_tag, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.BLAT_pp_mem_threshold, self.BLAT_pp_job_limit, self.BLAT_pp_job_delay, self.GA_BLAT_label, marker_file, self.commands, command_list)

            print(dt.today(), "submitted all BLAT merge jobs. waiting for sync")
            self.mp_util.wait_for_mp_store()

            final_checklist = os.path.join(self.GA_BLAT_path, "GA_BLAT_merge.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BLAT_merge_label)

        #print(dt.today(), "stopping for a sanity check: BLAT merge")
        #sys.exit()

        self.cleanup_GA_BLAT_start = time.time()
        self.mp_util.delete_folder_simple(self.GA_BLAT_jobs_folder)
        self.mp_util.clean_or_compress(self.GA_BLAT_path, self.keep_all, self.keep_GA_BLAT)

        self.cleanup_GA_BLAT_end = time.time()
        GA_BLAT_end = time.time()
        print("GA BLAT:", '%1.1f' % (self.GA_BLAT_end - self.GA_BLAT_start - (self.cleanup_GA_BLAT_end - self.cleanup_GA_BLAT_start)), "s")
        print("GA BLAT cleanup:", '%1.1f' % (self.cleanup_GA_BLAT_end - self.cleanup_GA_BLAT_start), "s")
        
        self.debug_stop_check(self.GA_BLAT_merge_label)
    
    
    def mp_GA_dmd(self):
        # ------------------------------------------------------
        # Diamond gene annotation
        self.GA_DIAMOND_start = time.time()
        #GA_DIAMOND_tool_output_path = os.path.join(self.GA_DIAMOND_path, "data", "0_diamond")
        #if not check_where_resume(None, self.GA_DIAMOND_tool_output_path, self.GA_BLAT_path, file_check_bypass = True):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_DIAMOND_label):
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(self.GA_BLAT_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    job_name = "DIAMOND_" + file_tag
                    full_sample_path = os.path.join(os.path.join(self.GA_BLAT_path, "final_results", split_sample))
                    marker_file = file_tag + "_diamond"
                    marker_path = os.path.join(self.GA_DIAMOND_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_path)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_DIAMOND_annotate_command_v2(self.GA_DIAMOND_label, full_sample_path, marker_file)
                        #self.mp_util.launch_and_create_with_hold(self.DIAMOND_mem_threshold, self.DIAMOND_job_limit, self.DIAMOND_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                        self.mp_util.launch_and_create_with_mem_footprint(self.DMD_mem_footprint, self.DIAMOND_job_limit, self.GA_DIAMOND_label, job_name, self.commands, command_list)

            print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.GA_DIAMOND_path, "GA_DIAMOND.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_DIAMOND_label)
        
        self.debug_stop_check(self.GA_DIAMOND_label)
        
    def mp_GA_dmd_pp(self):        
        #if not check_where_resume(GA_DIAMOND_path, None, self.GA_DIAMOND_tool_output_path, file_check_bypass = True):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_DIAMOND_pp_label):
            #print(dt.today(), "DIAMOND PP threads used:", self.paths.num_threads/2)
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(self.GA_BLAT_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    job_name = "DIAMOND_pp_" + file_tag
                    full_sample_path = os.path.join(os.path.join(self.GA_BLAT_path, "final_results", split_sample))
                    marker_file = file_tag + "_diamond_pp"
                    marker_path = os.path.join(self.GA_DIAMOND_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_DIAMOND_pp_command_v2(self.GA_DIAMOND_label, self.GA_BLAT_label, full_sample_path, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.DIAMOND_pp_mem_threshold, self.DIAMOND_pp_job_limit, self.DIAMOND_pp_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                                        
            print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.GA_DIAMOND_path, "GA_DIAMOND_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_DIAMOND_pp_label)
        
            
        
            self.cleanup_GA_DIAMOND_start = time.time()
            self.mp_util.delete_folder_simple(self.GA_DIAMOND_jobs_folder)
            self.mp_util.clean_or_compress(self.GA_DIAMOND_path, self.keep_all, self.keep_GA_DIAMOND)
            self.cleanup_GA_DIAMOND_end = time.time()
        self.GA_DIAMOND_end = time.time()
        print("GA DIAMOND:", '%1.1f' % (self.GA_DIAMOND_end - self.GA_DIAMOND_start - (self.cleanup_GA_DIAMOND_end - self.cleanup_GA_DIAMOND_start)), "s")
        print("GA DIAMOND cleanup:", '%1.1f' % (self.cleanup_GA_DIAMOND_end - self.cleanup_GA_DIAMOND_start), "s")
        
        self.debug_stop_check(self.GA_DIAMOND_pp_label)
        
    def mp_GA_final_merge(self):
        self.GA_final_merge_start = time.time()
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_final_merge_label):
            marker_file = "GA_final_merge"
            marker_path_p = os.path.join(self.ga_final_merge_path, "data", "jobs", "GA_final_merge_proteins")
            marker_path_m = os.path.join(self.ga_final_merge_path, "data", "jobs", "GA_final_merge_maps")
            marker_path_f = os.path.join(self.ga_final_merge_path, "data", "jobs", "GA_final_merge_fastq")
            if(os.path.exists(marker_path_p) and os.path.exists(marker_path_m) and os.path.exists(marker_path_f)):
                print(dt.today(), "skipping: GA final merge")
            else:
                command_list = self.commands.create_GA_final_merge_command(self.GA_final_merge_label, self.paths.assemble_contigs_label, self.GA_BWA_label, self.GA_BLAT_label, self.GA_DIAMOND_label,  marker_file)
                job_name = "GA_final_merge"
                self.mp_util.subdivide_and_launch(self.GA_final_merge_job_delay, self.GA_final_merge_mem_threshold, self.GA_final_merge_job_limit, self.GA_final_merge_label, job_name, self.commands, command_list)
            
            #check if all_proteins.faa was generated
            all_proteins_path = os.path.join(self.output_folder_path, self.GA_final_merge_label, "final_results", "all_proteins.faa")
            if(os.path.exists(marker_path_p)):
                if(os.path.getsize(all_proteins_path) > 0):
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_final_merge_label)
                    print(dt.today(), "All_proteins.faa is OK.  Continuing")
                else:
                    sys.exit("GA final merge failed.  proteins weren't translated")
                
        self.GA_final_merge_end = time.time()
        print("GA final merge:", '%1.1f' % (self.GA_final_merge_end - self.GA_final_merge_start), "s")
        self.mp_util.clean_or_compress(self.ga_final_merge_path, self.keep_all, self.keep_GA_final)
        
        self.debug_stop_check(self.GA_final_merge_label)

    def mp_TA(self):
        self.TA_start = time.time()
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ta_label):
            #-----------------------------------------
            # stage 1
            marker_path_list = []
            #----------------------------------------------
            #centrifuge is too much of a RAM hog.  can't run more than 1 at a time
            sections = ["reads"]
            for section in sections:
                marker_file = "TA_centrifuge_" + section
                marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_centrifuge_command(self.ta_label, self.paths.rRNA_filter_label, self.paths.assemble_contigs_label, section, marker_file)
                    self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
                    
            sections = ["singletons"]
            if self.read_mode == "paired":
                sections.extend(["paired"])
            if(self.contigs_present):
                sections.extend(["contigs"])    
            
            for section in sections:
                marker_file = "TA_kraken2_" + section
                marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_kraken2_command(self.ta_label, self.paths.assemble_contigs_label, section, marker_file)
                    self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)        
            marker_file = "TA_taxon_pull"
            marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_taxon_pull_command(self.ta_label, self.GA_final_merge_label, marker_file)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            print(dt.today(), "waiting for TA stage 1")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.TA_path, "TA_stage_1.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
            #--------------------------------------------------
            # stage 2
            marker_path_list = []
            sections = ["contigs"]
            for section in sections:
                marker_file = "TA_centrifuge_" + section
                marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_centrifuge_command(self.ta_label, self.paths.rRNA_filter_label, self.paths.assemble_contigs_label, section, marker_file)
                    self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.Centrifuge_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            
            marker_file = "TA_kraken2_pp"
            marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_kraken2_pp_command(self.ta_label, marker_file)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.TA_path, "TA_stage_2.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            #------------------------------------------------------------------

            #-----------------------------------------------------------------
            # stage 3
            marker_path_list = []
            marker_file = "TA_centrifuge_pp"
            marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_centrifuge_pp_command(self.ta_label, marker_file)
                self.mp_util.launch_and_create_with_hold(self.paths.TA_mem_threshold, self.paths.TA_job_limit, self.paths.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.TA_path, "TA_stage_3.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            #-----------------------------------------------
            # stage 4
            marker_path_list = []
            
            marker_file = "TA_final"
            marker_path = os.path.join(self.paths.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_final_command(self.ta_label, self.paths.assemble_contigs_label, marker_file)
                self.mp_util.launch_and_create_simple(self.ta_label, marker_file, self.commands, command_list)
            final_checklist = os.path.join(self.TA_path, "TA_final.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
            if(os.path.exists(marker_path)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.ta_label)
                
        self.cleanup_TA_start = time.time()
        self.mp_util.clean_or_compress(self.TA_path, self.keep_all, self.keep_TA)
        self.cleanup_TA_end = time.time()
        self.TA_end = time.time()
        print("TA:", '%1.1f' % (self.TA_end - self.TA_start - (self.cleanup_TA_end - self.cleanup_TA_start)), "s")
        print("TA cleanup:", '%1.1f' % (self.cleanup_TA_end - self.cleanup_TA_start), "s")
        
        self.debug_stop_check(self.ta_label)

    def mp_EC(self):
        
        self.EC_start = time.time()
        #There's a 2-step check.  We don't want it ti re-run either DETECT, or PRIAM+DIAMOND because they're too slow
        #if not check_where_resume(ec_path, None, self.GA_DIAMOND_path):
        #if check_bypass_log(self.output_folder_path, self.ec_label):
        
        
        # --------------------------------------------------------------
        # Priam EC annotation.  Why isn't it parallel? computing restraints.  Not enough mem
        self.EC_PRIAM_start = time.time()

        #split the proteins
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_priam_split_label):
            marker_file = "ec_split"
            marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                command_list = self.commands.create_EC_PRIAM_split_command(self.ec_label, self.GA_final_merge_label, self.ec_split_path, marker_file)
                print(command_list)
                self.mp_util.launch_only_simple(self.commands, command_list)
                
                
        #if not check_where_resume(job_label = None, full_path = ec_priam_path, dep_job_path = GA_DIAMOND_path):
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_priam_label):
            split_count = 0
            
            print("ec split dir:", os.listdir(self.ec_split_path))
            
            for ec_split_file in os.listdir(self.ec_split_path):
                marker_file = "ec_priam_" + str(split_count)
                marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
                print(marker_path)
                #time.sleep(10)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    priam_file_path = os.path.join(self.ec_priam_path, "split_" + str(split_count))
                    split_path = os.path.join(self.ec_split_path, "protein_split_" + str(split_count) + ".fasta")
                    self.mp_util.make_folder(priam_file_path)    
                    command_list = self.commands.create_EC_PRIAM_command_v2(self.ec_label, self.GA_final_merge_label, priam_file_path, split_path, split_count, marker_file)
                    print(command_list)
                    if(os.path.exists(priam_file_path)):
                        print(dt.today(), "attempting PRIAM auto-resume")
                        
                    self.mp_util.launch_only_with_hold(self.EC_mem_threshold, self.EC_job_limit, self.EC_job_delay, marker_file, self.commands, command_list)
                split_count += 1
        
        
        
            #process.join()
        print(dt.today(), "Waiting for PRIAM jobs")
        self.mp_util.wait_for_mp_store()
        
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_priam_cat_label):
            marker_file = "ec_priam_cat"
            marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_path)
            else:
                command_list = self.commands.create_EC_PRIAM_cat_command(self.ec_label, marker_file)
                print("CAT:", command_list)
                self.mp_util.launch_only_simple(self.commands, command_list)
            
            
        self.EC_PRIAM_end = time.time()
        print("EC PRIAM:", '%1.1f' % (self.EC_PRIAM_end - self.EC_PRIAM_start), "s")
        
        #-----------------------------------------------------------------------
        # Detect
        self.EC_DETECT_start = time.time()
        #if not check_where_resume(job_label = None, full_path = ec_detect_path, dep_job_path = GA_DIAMOND_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_detect_label):
            marker_file = "ec_detect"
            marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                command_list = self.commands.create_EC_DETECT_command(self.ec_label, self.GA_final_merge_label, marker_file)
                self.mp_util.launch_and_create_with_mp_store(self.ec_label, marker_file, self.commands, command_list)
            
            
        self.EC_DETECT_end = time.time()
        print("EC DETECT:", '%1.1f' % (self.EC_DETECT_end - self.EC_DETECT_start), "s")
        
        
        
        # --------------------------------------------------------------
        # DIAMOND EC annotation 
        self.EC_DIAMOND_start = time.time()
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_DIAMOND_label):
            marker_file = "ec_diamond"
            marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                job_name = "ec_diamond"
                command_list = self.commands.create_EC_DIAMOND_command(self.ec_label, self.GA_final_merge_label, marker_file)
                self.mp_util.launch_and_create_with_mp_store(self.ec_label, job_name, self.commands, command_list)
            
        self.EC_DIAMOND_end = time.time()
        self.mp_util.wait_for_mp_store()
        
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_detect_label):
            if(os.path.exists(self.ec_detect_out)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.ec_detect_label)
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_priam_label):
            if(os.path.exists(self.ec_priam_out)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.ec_priam_label)
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_DIAMOND_label):
            if(os.path.exists(self.ec_diamond_out)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.ec_DIAMOND_label)
        
        #----------------------------------------------------------------------
        # EC post process
        self.EC_post_start = time.time()
        #if not (check_where_resume(ec_path, None, self.GA_DIAMOND_path)):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.ec_pp_label):
            
            marker_file = "ec_post"
            marker_path = os.path.join(self.ec_path, "data", "jobs", marker_file)
            
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                command_list = self.commands.create_EC_postprocess_command(self.ec_label, self.GA_final_merge_label, marker_file)
                self.mp_util.launch_and_create_simple(self.ec_label, marker_file, self.commands, command_list)
            
            if(os.path.exists(marker_path)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.ec_pp_label)
        
        self.cleanup_EC_start = time.time()
        self.mp_util.clean_or_compress(self.ec_path, self.keep_all, self.keep_EC)
        self.cleanup_EC_end = time.time()
        self.EC_post_end = time.time()
            
    
        self.EC_end = time.time()
        print("EC run:", '%1.1f' % (self.EC_end - self.EC_start), "s")
        print("EC cleanup:", '%1.1f' % (self.cleanup_EC_end - self.cleanup_EC_start), "s")
        
        self.debug_stop_check(self.ec_label)

    def mp_output(self):
        self.Cytoscape_start = time.time()
        #if not check_where_resume(network_path, None, self.ec_path):
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.output_label):
            
            #phase 1
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_copy_gene_map_label):
                job_name = self.output_copy_gene_map_label
                command_list = self.commands.create_output_copy_gene_map_command(self.output_label, self.GA_final_merge_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_copy_taxa_label):
                job_name = self.output_copy_taxa_label
                command_list = self.commands.create_output_copy_taxa_command(self.output_label, self.ta_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
            if(self.contigs_present): 
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_contig_stats_label):
                    job_name = self.output_contig_stats_label
                    command_list = self.commands.create_output_contig_stats_command(self.output_label, self.paths.assemble_contigs_label)
                    self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            
                
            if not(self.no_host):
                print(dt.today(), "repopulating hosts for output")
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_singletons_label):
                    job_name = self.output_unique_hosts_singletons_label
                    command_list = self.commands.create_output_unique_hosts_singletons_command(self.output_label, self.paths.qc_label, self.host_filter_label)
                    self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
                if(self.read_mode == "paired"):
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_pair_1_label):
                        job_name = self.output_unique_hosts_pair_1_label
                        command_list = self.commands.create_output_unique_hosts_pair_1_command(self.output_label, self.paths.qc_label, self.host_filter_label)
                        self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                        
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_pair_2_label):
                        job_name = self.output_unique_hosts_pair_2_label
                        command_list = self.commands.create_output_unique_hosts_pair_2_command(self.output_label, self.paths.qc_label, self.host_filter_label)
                        self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                        
                        
            #repop vectors
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_singletons_label):
                job_name = self.output_unique_vectors_singletons_label
                command_list = self.commands.create_output_unique_vectors_singletons_command(self.output_label, self.paths.qc_label, self.host_filter_label, self.vector_filter_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
            
            if(self.read_mode == "paired"):
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_pair_1_label):
                    job_name = self.output_unique_vectors_pair_1_label
                    command_list = self.commands.create_output_unique_vectors_pair_1_command(self.output_label, self.paths.qc_label, self.host_filter_label, self.vector_filter_label)
                    self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                    
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_pair_2_label):
                    job_name = self.output_unique_vectors_pair_2_label
                    command_list = self.commands.create_output_unique_vectors_pair_2_command(self.output_label, self.paths.qc_label, self.host_filter_label, self.vector_filter_label)
                    self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                    
            print(dt.today(), "output report phase 1 launched.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            
            self.mp_util.conditional_write_to_bypass_log(self.output_per_read_scores_label, "outputs/final_results", "input_per_seq_quality_report.csv")
            self.mp_util.conditional_write_to_bypass_log(self.output_copy_gene_map_label, "outputs/final_results", "final_gene_map.tsv")
            self.mp_util.conditional_write_to_bypass_log(self.output_copy_taxa_label, "outputs/final_results", "taxa_classifications.tsv")
            self.mp_util.conditional_write_to_bypass_log(self.output_contig_stats_label, "outputs/final_results", "contig_stats.txt")
            self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_singletons_label, "outputs/data/4_full_vectors", "singletons_full_vectors.fastq")
            if(self.read_mode == "paired"):
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_pair_1_label, "outputs/data/4_full_vectors", "pair_1_full_vectors.fastq")
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_pair_2_label, "outputs/data/4_full_vectors", "pair_2_full_vectors.fastq")
                
            if not (self.no_host):
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_singletons_label, "outputs/data/2_full_hosts", "singletons_full_hosts.fastq")
                if(self.read_mode == "paired"):
                    self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_pair_1_label, "outputs/data/2_full_hosts", "pair_1_full_hosts.fastq")
                    self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_pair_2_label, "outputs/data/2_full_hosts", "pair_2_full_hosts.fastq")
            #----------------------------------------------------------------------------
            #Phase 2
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_network_gen_label):
                command_list = self.commands.create_output_network_generation_command(self.output_label, self.GA_final_merge_label, self.ta_label, self.ec_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, self.output_network_gen_label, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_taxa_groupby_label):
                command_list = self.commands.create_output_taxa_groupby_command(self.output_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, self.output_taxa_groupby_label, self.commands, command_list)
        
            print(dt.today(), "output report phase 2 launched.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            self.mp_util.conditional_write_to_bypass_log(self.output_network_gen_label, "outputs/final_results", "RPKM_table.tsv")
            
            
            #-------------------------------------------------------------------
            #Phase 3
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_read_count_label):
                job_name = self.output_read_count_label
                command_list = self.commands.create_output_read_count_command(self.output_label, self.paths.qc_label, self.paths.repop_label, self.GA_final_merge_label, self.ec_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                                

            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_per_read_scores_label):
                job_name = self.output_per_read_scores_label
                command_list = self.commands.create_output_per_read_scores_command(self.output_label, self.paths.qc_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_ec_heatmap_label):
                job_name = self.output_ec_heatmap_label
                command_list = self.commands.create_output_EC_heatmap_command(self.output_label)
                self.mp_util.launch_and_create_with_mp_store(self.output_label, job_name, self.commands, command_list)    
            
            print(dt.today(), "output report phase 3 launched.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            self.mp_util.conditional_write_to_bypass_log(self.output_read_count_label, "outputs/final_results", "read_count.tsv")
            self.mp_util.conditional_write_to_bypass_log(self.output_ec_heatmap_label, "outputs/final_results", "EC_coverage.csv")
            self.mp_util.conditional_write_to_bypass_log(self.output_per_read_scores_label, "outputs/final_results", "quality_filter_hist.jpg")

            
        self.cleanup_cytoscape_start = time.time()
        self.mp_util.clean_or_compress(self.network_path, self.keep_all, self.keep_outputs)
        self.cleanup_cytoscape_end = time.time()
            
            
            
            
        self.Cytoscape_end = time.time()
        self.end_time = time.time()
        print("Outputs:", '%1.1f' % (self.Cytoscape_end - self.Cytoscape_start - (self.cleanup_cytoscape_end - self.cleanup_cytoscape_start)), "s")
        print("Outputs cleanup:", '%1.1f' % (self.cleanup_cytoscape_end - self.cleanup_cytoscape_start), "s")
        print("=============================================================================================")
        print("Final summary")
        print("--------------------------------------------------------")
        print("Total runtime:", '%1.1f' % (self.end_time - self.start_time), "s")
        print("quality filter:", '%1.1f' % (self.quality_end - self.quality_start - (self.cleanup_quality_end - self.cleanup_quality_start)), "s")
        print("quality filter cleanup:", '%1.1f' %(self.cleanup_quality_end - self.cleanup_quality_start), "s")
        if not self.no_host:
            print("host filter:", '%1.1f' % (self.host_end - self.host_start - (self.cleanup_host_end - self.cleanup_host_start)), "s")
            print("host filter cleanup:", '%1.1f' %(self.cleanup_host_end - self.cleanup_host_start),"s")
        print("vector filter:", '%1.1f' % (self.vector_end - self.vector_start - (self.cleanup_vector_end - self.cleanup_vector_start)), "s")
        print("vector filter cleanup:", '%1.1f' % (self.cleanup_vector_end - self.cleanup_vector_start), "s")
        print("rRNA filter:", '%1.1f' % (self.rRNA_filter_end - self.rRNA_filter_start - (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start)), "s")
        print("rRNA filter cleanup:", '%1.1f' % (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start), "s")
        print("repop:", '%1.1f' % (self.repop_end - self.repop_start - (self.cleanup_repop_end - self.cleanup_repop_start)), "s")
        print("repop cleanup:", '%1.1f' % (self.cleanup_repop_end - self.cleanup_repop_start), "s")
        print("assemble contigs:", '%1.1f' % (self.assemble_contigs_end - self.assemble_contigs_start - (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start)), "s")    
        print("assemble contigs cleanup:", '%1.1f' % (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start), "s")
        print("GA BWA:", '%1.1f' % (self.GA_BWA_end - self.GA_BWA_start - (self.cleanup_GA_BWA_end - self.cleanup_GA_BWA_start)), "s")
        print("GA BWA cleanup:", '%1.1f' % (self.cleanup_GA_BWA_end - self.cleanup_GA_BWA_start), "s")
        print("GA BLAT:", '%1.1f' % (self.GA_BLAT_end - self.GA_BLAT_start - (self.cleanup_GA_BLAT_end - self.cleanup_GA_BLAT_start)), "s")
        print("GA BLAT cleanup:", '%1.1f' % (self.cleanup_GA_BLAT_end - self.cleanup_GA_BLAT_start), "s")
        print("GA DIAMOND:", '%1.1f' % (self.GA_DIAMOND_end - self.GA_DIAMOND_start - (self.cleanup_GA_DIAMOND_end - self.cleanup_GA_DIAMOND_start)), "s")
        print("GA DIAMOND cleanup:", '%1.1f' % (self.cleanup_GA_DIAMOND_end - self.cleanup_GA_DIAMOND_start), "s")
        print("TA:", '%1.1f' % (self.TA_end - self.TA_start - (self.cleanup_TA_end - self.cleanup_TA_start)), "s")
        print("TA cleanup:", '%1.1f' % (self.cleanup_TA_end - self.cleanup_TA_start), "s")
        print("EC:", '%1.1f' % (self.EC_end - self.EC_start), "s")
        #print("---------------------------------------------")
        #print("Note: EC is in cloud-mode.  ignore individual timing")
        #print("EC DETECT:", '%1.1f' % (self.EC_DETECT_end - EC_DETECT_start), "s")
        #print("EC PRIAM:", '%1.1f' % (self.EC_PRIAM_end - EC_PRIAM_start), "s")
        #print("EC DIAMOND:", '%1.1f' % (self.EC_DIAMOND_end - EC_DIAMOND_start), "s")
        #print("EC cleanup:", '%1.1f' % (self.cleanup_EC_end - cleanup_EC_start), "s")
        #print("-------------------------------------------------")
        print("Outputs:", '%1.1f' % (self.Cytoscape_end - self.Cytoscape_start - (self.cleanup_cytoscape_end - self.cleanup_cytoscape_start)), "s")
        print("Outputs cleanup:", '%1.1f' % (self.cleanup_cytoscape_end - self.cleanup_cytoscape_start), "s")
        
#---------------------------------------------------------------------------------------------------------------------------------


    
    
    
