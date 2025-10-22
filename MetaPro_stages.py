#!/usr/bin/env python
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import MetaPro_utilities as mpu
import MetaPro_marker as mpm
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q
import shutil

#stores code for stage-launch.
#makes for a neat package/capsule

class mp_stage:
    def __init__ (self, config_dict, dir_control, time_control, file_control): #config_obj, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, args_pack, tutorial_mode_string = None):
        #make our util obj
        #refresher: self -> instance var.  not self: class var (shared among class obj instances)
        print("in stage:", config_dict["pair_1"])
        time.sleep(3)
        #---------------------------------------------------------
        #Operational flags and state-recorders
        self.dir_control = dir_control
        self.dir_dict = self.dir_control.get_dir_dict()
        self.label_dict = self.dir_control.get_label_dict()
        self.config_dict = config_dict
        print(dt.today(), "MP STAGE using:", self.dir_dict["main"])
        #time.sleep(5)
        #self.tutorial_string = tutorial_mode_string
        self.output_folder_path = self.dir_dict["main"] #output_folder_path
        self.mp_util = mpu.mp_util(self.config_dict, self.dir_dict)
        self.mp_seq_handler = mpu.mp_seq_handler()
        self.marker_control = mpm.mpro_marker(self.config_dict, self.dir_dict)
        self.marker_dict = self.marker_control.marker_dict
        self.m_name = self.marker_control.m_name
        #self.config_dict = config_dict
        
        self.process_list = self.mp_util.mp_store
        
        self.file_control = file_control
        self.file_dict = self.file_control.get_file_dict() #recall python passes by reference.
        self.time_control = time_control
        self.seq_handler = mpu.mp_seq_handler()
        

        #time.sleep(10)
        
        
        
        self.config_dict["no_host"] = config_dict["no_host"]
        self.verbose_mode = config_dict["verbose_mode"]
        self.rRNA_chunks = int(self.config_dict["rRNA_chunksize"])
        self.EC_chunksize = int(self.config_dict["EC_chunksize"])
        self.GA_chunksize = int(self.config_dict["GA_chunksize"])
        #self.config_path = config_path
        self.pair_1_path = self.config_dict["pair_1"]
        self.pair_2_path = self.config_dict["pair_2"]
        self.single_path = self.config_dict["single"]
        self.contig_path = self.config_dict["contig"] #_path  #tutorial/single-shot use
        
        self.read_mode = self.config_dict["read_mode"]
        
        
        
        self.debug_stop_flag = self.config_dict["debug_stop_flag"]
        
        mp_store = []  # stores the multiprocessing processes

        # Creates our command object, for creating shellscripts.

        self.commands = mpcom.mt_pipe_commands(self.config_dict, self.dir_dict, self.file_dict, self.marker_dict)
        
        #special contig-bypasser logic vars
        self.contigs_present = True  #for the contig/assembly bypasser
        self.spades_done_file = self.dir_dict["spades_done"]
        self.spades_transcripts_file = self.dir_dict["spades_transcripts"]
    
    #--------------------------------------------------------------------------------
    # helpers
    def mp_contig_statecheck(self):
        
        if(os.path.exists(self.spades_done_file)):
            if(os.path.exists(self.spades_transcripts_file)):
                self.contigs_present = True
                print(dt.today(), "contigs present")
            else:
                self.contigs_present = False
                print(dt.today(), "contigs not present. but that's ok")
        else:
            sys.exit(dt.today(), "contig assembly incomplete. run the pipeline in --tutorial contigs to continue")
            
    def import_lib_names(self, names_file):
    #import the list of libs from GA_pre_scan
    #lib list is also how custom DBs can enter.
    #also checks for validity.
        files_list = []
        with open(names_file, "r") as names_in:
            for line in names_in:
                cleaned_line = line.strip("\n")
                if("can't" in line):
                    continue
                #print("line:", cleaned_line)
                src_path = cleaned_line.split("|")[3]
                if(os.path.exists(src_path)):
                    index_sample = src_path + ".1.bt2"
                    if(os.path.exists(index_sample)):
                        files_list.append(src_path)
                    else:
                        print(dt.today(), "Ending pipe.  BT2 index is missing for:", src_path)
                        sys.exit()
    
                #not complete. the exist path needs to point to the DNA_db folder
                    #print("src:", src_path)
                    #time.sleep(1)
                    
        
        return files_list
    
    def import_ga_lib_list(self, lib_file):
        lib_list = list()
        with open(lib_file, "r") as libs_in:
            for line in libs_in:
                lib_line = line.strip("\n")
                lib_list.append(lib_line)
        return lib_list
    
    def debug_stop_check(self, stop_signal):
        if(self.debug_stop_flag == stop_signal):
            exit_string = "stopped after: " + stop_signal
            sys.exit(exit_string)
        else:
            print(dt.today(), "continuing from:", stop_signal)
            
    #--------------------------------------------------------------------------------------------------------------
    # main calls
    def mp_quality_filter(self):
        if(self.marker_control.check_marker("qf")):
            print(dt.today(), "running QF")
            if(self.config_dict["q_enc"] == "None"):
                if (self.single_path != "None"):
                    self.config_dict["q_enc"] = self.mp_util.determine_encoding(self.single_path)
                    print("ENCODING USED:", self.config_dict["q_enc"])
                else:
                    self.config_dict["q_enc"] = self.mp_util.determine_encoding(self.pair_1_path)
                    print("ENCODING USED:", self.config_dict["q_enc"])
            else:
                print(dt.today(), "skipping encode-finding. using config preset:", self.config_dict["q_enc"])

            self.time_control.measure_time("qf", "start")
            for item in self.dir_dict["qf_list"]:
                self.dir_control.make_dirs(self.dir_dict[item])

            command_list = self.commands.create_quality_control_command(self.marker_dict["qf"], self.config_dict["q_enc"])
            self.mp_util.launch_stage_simple(self.file_dict["qf_job"], command_list, self.config_dict["keep_all"], self.config_dict["keep_quality"])
            self.time_control.measure_time("qf", "end")
            self.marker_control.place_marker("qf")
            self.debug_stop_check(self.label_dict["qf"])
        else:
            print(dt.today(), "skipping: QF")
            

    def mp_host_filter(self):
        if not self.config_dict["no_host"]:
            print("host id list:", self.config_dict["Host_IDs"])
            if("none" in self.config_dict["Host_IDs"]):
                print(dt.today(), "no hosts specified. bypassing")
            
            else:
                print(dt.today(), "host elements detected")
                for host_id in self.config_dict["Host_IDs"]:
                    if((host_id != "none") or (host_id != "None")):
                        print(dt.today(), "using Host ID:", host_id)
                        print(dt.today(), "checking for marker:", host_id + "_host")
                        if(self.marker_control.check_marker(host_id + "_host")):
                            self.dir_control.make_dirs_from_list(host_id + "_dir_list")

                            print(dt.today(), "running:", host_id + "_host")
                            self.time_control.measure_time("host", "start")
                            host_count = 0
                            if(host_id == "none"):
                                print(dt.today(), "no host specified. skipping")
                                break
                            
                            host_mkr = self.marker_dict[host_id + "_host"]
                            host_job = os.path.join(self.dir_dict["host_jobs"], host_id + "_job.sh")
                            command_list = self.commands.create_host_filter_command(self.config_dict["Host_IDs"], host_count, host_mkr)
                            self.mp_util.launch_stage_simple(host_job, command_list, self.config_dict["keep_all"], self.config_dict["keep_host"])
                            self.time_control.measure_time("host", "end")
                            self.marker_control.place_marker(host_id + "_host")
                        else:
                            print(dt.today(), "skipping:", host_id + "_host")
                    else:
                        print(dt.today(), "no hosts to filter")
                        break
            
                #self.debug_stop_check(self.label_dict[host_id + "_host"])

                

    def mp_vector_filter(self):
        self.vector_start = time.time()
        if(self.marker_control.check_marker("vec")):
            #if self.config_dict["no_host"]:
                #get dep args from quality filter
                #if not check_where_resume(vector_path, None, self.quality_path):
            final_host = self.config_dict["Host_IDs"][-1]
            print("using final host:", final_host)

            for item in self.dir_dict["vec_list"]:
                print("making:", item, self.dir_dict[item])
                self.dir_control.make_dirs(self.dir_dict[item])
                

            command_list = self.commands.create_vector_filter_command(self.marker_dict["vec"], final_host)
            self.mp_util.launch_stage_simple(self.file_dict["vec_job"], command_list, self.config_dict["keep_all"], self.config_dict["keep_vector"])

        
            self.vector_end = time.time()
            self.debug_stop_check(self.label_dict["vec"])
            self.marker_control.place_marker("vec")

    def mp_rRNA_filter(self):
        #don't split bnap.
        #only split infernal.  but it can go over the 100k limit now.
        #bnap takes fasta, but it's not worth it to keep it in memory.


        self.rRNA_filter_start = time.time()
        
        #if not check_where_resume(self.rRNA_filter_path, None, self.vector_path):
        #if self.mp_util.check_bypass_log(self.output_folder_path, self.label_dict["rRNA"]): 
        if self.marker_control.check_marker("rRNA"):
        #run bnap

            for item in self.dir_dict["rRNA_list"]:
                print(dt.today(), "making:", item, self.dir_dict[item])
                self.dir_control.make_dirs(self.dir_dict[item])
            run_types = ["s"]
            #print("SELF.READ_MODE:", self.read_mode)
            if(self.read_mode == "paired"):
                run_types.extend(["p1", "p2"])

            print(dt.today(), "run types:", run_types)
            #time.sleep(10)
            for run_type in run_types:
                fa_in = self.file_dict["rRNA_"+str(run_type) + "_fa"]
                fq_in = self.file_dict["no_vec_" + str(run_type)]
                mRNA_out = self.file_dict["rRNA_bnap_mRNA_" + str(run_type)]
                rRNA_out = self.file_dict["rRNA_bnap_other_" + str(run_type)]
                bnap_out = self.file_dict["rRNA_bnap_all_" + str(run_type)]
                marker = self.marker_dict["rRNA_bnap_" + str(run_type)]
                self.mp_seq_handler.fastq_to_fasta(fq_in, fa_in)
                if not (os.path.exists(marker)):
                    print(dt.today(), "running:", marker)
                    command_list = self.commands.create_rRNA_filter_bnap_command(fa_in, fq_in, mRNA_out, rRNA_out, bnap_out, marker)
                    self.mp_util.launch_stage_simple(self.file_dict["rRNA_bnap_job_" + str(run_type)], command_list, self.config_dict["keep_all"], self.config_dict["keep_rRNA"])
                    #self.mp_util.run_subjob_with_mp_store()
                else:
                    print(dt.today(), "skipping:", marker)


            #wait for everything.
            self.mp_util.wait_for_mp_store()

                
            #----------------------------------------------------------------------------
            # INFERNAL
            #split the fasta mRNA here
            split_count_s = 0
            split_count_p1 = 0
            split_count_p2 = 0
            if not (os.path.exists(self.marker_dict["rRNA_split_s"])):
                split_count_s = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_s"], self.file_dict["rRNA_inf_in_s"], self.config_dict["rRNA_chunksize"], "fasta")
                self.marker_control.place_marker("rRNA_split_s")
            else:
                split_count_s = self.mp_seq_handler.get_split_count(self.dir_dict["rRNA_split"], "inf_s")

            if not (os.path.exists(self.marker_dict["rRNA_split_p1"])): 
                split_count_p1 = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_p1"], self.file_dict["rRNA_inf_in_p1"], self.config_dict["rRNA_chunksize"], "fasta")
                self.marker_control.place_marker("rRNA_split_p1")
            else:
                split_count_p1 = self.mp_seq_handler.get_split_count(self.dir_dict["rRNA_split"], "inf_p1")
            if not (os.path.exists(self.marker_dict["rRNA_split_p2"])):    
                split_count_p2 = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_p2"], self.file_dict["rRNA_inf_in_p2"], self.config_dict["rRNA_chunksize"], "fasta")
                self.marker_control.place_marker("rRNA_split_p2")
            else:
                split_count_p2 = self.mp_seq_handler.get_split_count(self.dir_dict["rRNA_split"], "inf_p2")
            
            for i in range(0, split_count_s):
                self.file_dict["rRNA_inf_in_s_" + str(i)] = self.file_dict["rRNA_inf_in_s"] + "_" + str(i) + ".fasta"
            for i in range(0, split_count_p1):
                self.file_dict["rRNA_inf_in_p1_" + str(i)] = self.file_dict["rRNA_inf_in_p1"] + "_" + str(i) + ".fasta"
            for i in range(0, split_count_p2):
                self.file_dict["rRNA_inf_in_p2_" + str(i)] = self.file_dict["rRNA_inf_in_p2"] + "_" + str(i) + ".fasta"

            #create inf out path
            self.file_control.create_split_filepaths("rRNA_inf_out_s", split_count_s, self.file_dict["rRNA_inf_out_s"], ".inf_out")
            self.file_control.create_split_filepaths("rRNA_inf_out_p1", split_count_p1, self.file_dict["rRNA_inf_out_p1"], ".inf_out")
            self.file_control.create_split_filepaths("rRNA_inf_out_p2", split_count_p2, self.file_dict["rRNA_inf_out_p2"], ".inf_out")
            self.file_control.create_split_filepaths("rRNA_inf_job_s", split_count_s, self.file_dict["rRNA_inf_job_s"], ".sh")
            self.file_control.create_split_filepaths("rRNA_inf_job_p1", split_count_p1, self.file_dict["rRNA_inf_job_p1"], ".sh")
            self.file_control.create_split_filepaths("rRNA_inf_job_p2", split_count_p2, self.file_dict["rRNA_inf_job_p2"], ".sh")
            print("split count [s, p1, p2]", split_count_s, split_count_p1, split_count_p2)
            self.marker_control.issue_split_markers(self.m_name["rRNA_inf_s"], split_count_s, "rRNA_mkrs")
            self.marker_control.issue_split_markers(self.m_name["rRNA_inf_p1"], split_count_p1, "rRNA_mkrs")
            self.marker_control.issue_split_markers(self.m_name["rRNA_inf_p2"], split_count_p2, "rRNA_mkrs")

            for i in range(0, split_count_s):
                s_seq = self.file_dict["rRNA_inf_in_s_" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_s_" + str(i)]
                self.file_dict["rRNA_inf_in_s_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_s_" + str(i)]
                job_path = self.file_dict["rRNA_inf_job_s_" + str(i)]
                if(self.marker_control.check_marker("rRNA_inf_s_" + str(i))):
                    command = self.commands.create_rRNA_filter_infernal_command(s_seq, inf_out, marker)
                    self.mp_util.run_subjob_no_final(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        job_path,
                        command)    
            self.mp_util.wait_for_mp_store()
            
            for i in range(0, split_count_p1):
                p1_seq = self.file_dict["rRNA_inf_in_p1_" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_p1_" + str(i)]
                self.file_dict["rRNA_inf_in_p1_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_p1_" + str(i)]
                job_path = self.file_dict["rRNA_inf_job_p1_" + str(i)]
                if(self.marker_control.check_marker("rRNA_inf_p1_" + str(i))):
                    command = self.commands.create_rRNA_filter_infernal_command(p1_seq, inf_out, marker)
                    self.mp_util.run_subjob_no_final(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        job_path,
                        command)    
            self.mp_util.wait_for_mp_store()

            for i in range(0, split_count_p2):
                p2_seq = self.file_dict["rRNA_inf_in_p2_" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_p2_" + str(i)]
                self.file_dict["rRNA_inf_in_p2_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_p2_" + str(i)]
                job_path = self.file_dict["rRNA_inf_job_p2_" + str(i)]
                if(self.marker_control.check_marker("rRNA_inf_p2_" + str(i))):
                    command = self.commands.create_rRNA_filter_infernal_command(p2_seq, inf_out, marker)
                    self.mp_util.run_subjob_no_final(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        job_path,
                        command)    

            #wait for everything.
            self.mp_util.wait_for_mp_store()
            #-----------------------------------------------
            #merge the inf and bnap files, then send it through inf pp.
            if(split_count_s > 0):
                with open(self.file_dict["rRNA_inf_all_s"], "wb") as out_file:
                    for i in range(0, split_count_s):
                        inf_out = self.file_dict["rRNA_inf_out_s_" + str(i)]
                        with open(inf_out, "rb") as in_file:
                            shutil.copyfileobj(in_file, out_file)
                        out_file.write(b"\n")
                

            if(split_count_p1 > 0):
                with open(self.file_dict["rRNA_inf_all_p1"], "wb") as out_file:
                    for i in range(0, split_count_p1):
                        inf_out = self.file_dict["rRNA_inf_out_p1_" + str(i)]
                        with open(inf_out, "rb") as in_file:
                            shutil.copyfileobj(in_file, out_file)
                        out_file.write(b"\n")
                

            if(split_count_p2 > 0):
                with open(self.file_dict["rRNA_inf_all_p2"], "wb") as out_file:
                    for i in range(0, split_count_p2):
                        inf_out = self.file_dict["rRNA_inf_out_p2_" + str(i)]
                        with open(inf_out, "rb") as in_file:
                            shutil.copyfileobj(in_file, out_file)
                        out_file.write(b"\n")
                

            


            #inf pp 3 sections
            if(split_count_p1 > 0):
                inf_p1 = self.file_dict["rRNA_inf_all_p1"]
                inf_p2 = self.file_dict["rRNA_inf_all_p2"]
                bnap_p1 = self.file_dict["rRNA_bnap_all_p1"]
                bnap_p2 = self.file_dict["rRNA_bnap_all_p2"]
                mRNA_p1 = self.file_dict["rRNA_mRNA_p1_fq"]
                mRNA_p2 = self.file_dict["rRNA_mRNA_p2_fq"]
                other_p1 = self.file_dict["rRNA_other_p1_fq"]
                other_p2 = self.file_dict["rRNA_other_p2_fq"]
                
                command = self.commands.create_rRNA_inf_pp_pair_command(inf_p1, inf_p2, 
                                                                        bnap_p1, bnap_p2,
                                                                        self.file_dict["no_vec_p1"],
                                                                        self.file_dict["no_vec_p2"],
                                                                        mRNA_p1, mRNA_p2,
                                                                        other_p1, other_p2,
                                                                        self.marker_dict["rRNA_inf_pp_paired"]
                                                                        )
                
                self.mp_util.run_subjob_no_final(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.file_dict["rRNA_pp_p_job"],
                        command) 


            command = self.commands.create_rRNA_inf_pp_s_command(self.file_dict["rRNA_inf_all_s"],
                                                                 self.file_dict["rRNA_bnap_all_s"],
                                                                 self.file_dict["no_vec_s"],
                                                                 self.file_dict["rRNA_mRNA_s_fq"],
                                                                 self.file_dict["rRNA_other_s_fq"],
                                                                 self.marker_dict["rRNA_inf_pp_s"]
                                                                 )
            self.mp_util.run_subjob_no_final(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.file_dict["rRNA_pp_s_job"],
                        command) 
            self.mp_util.wait_for_mp_store()
            self.marker_control.place_marker("rRNA")

#-----------------------------------------------------------------------------------------------------------------------        

    def mp_repop(self):
        
        self.repop_start = time.time()
        #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
        if self.marker_control.check_marker("repop"):
            for item in self.dir_dict["repop_list"]:
                print(dt.today(), "making:", item, self.dir_dict[item])
                self.dir_control.make_dirs(self.dir_dict[item])
            #job_name = self.repop_job_label
            
            command_list = self.commands.create_repop_command(self.marker_dict["repop"])
            self.mp_util.run_subjob_with_hold(  self.config_dict["repop_mem_threshold"], 
                                                self.config_dict["repop_job_limit"], 
                                                self.config_dict["repop_job_delay"], 
                                                self.file_dict["repop_job"], command_list)
            #self.mp_util.wait_for_mp_store()
            self.marker_control.place_marker("repop")
        self.repop_end = time.time()


    def mp_assemble(self):
        self.assemble_contigs_start = time.time()
        
        #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
        mgm_gene_report = self.file_dict["contigs_gene_report"]
        spades_done_file = self.file_dict["contigs_spades_done"]
        spades_transcript_file = self.file_dict["contigs_transcripts"]
        mgm_fail_flag = True
        spades_fail_flag = True

        if self.marker_control.check_marker("contigs"):
            self.dir_control.make_dirs_from_list("contigs_list")
            command_list = self.commands.create_assemble_contigs_command(self.marker_dict["contigs"])
            self.mp_util.run_subjob_with_hold(
                self.config_dict["repop_mem_threshold"], 
                self.config_dict["repop_job_limit"], 
                self.config_dict["repop_job_delay"],
                self.file_dict["contigs_job"], 
                command_list)
            if(os.path.exists(spades_done_file)):
                if(os.path.exists(spades_transcript_file)):
                    spades_fail_flag = False
                    print(dt.today(), "SPADes OK")
                else:
                    spades_fail_flag = True
                    print(dt.today(), "SPADes ran, but did not create contigs")
            else:
                exit_string = str(dt.today()) + " SPADes did not run. this is not normal. Contact admin immediately"
                sys.exit(exit_string)

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
                bypass_contig_map_path = self.file_dict["contigs_map"]                
                bypass_contig_path = self.file_dict["contigs_out_fa"]
                s_src_path = self.file_dict["repop_s"]
                p1_src_path = self.file_dict["repop_p1"]
                p2_src_path = self.file_dict["repop_p2"]
                s_dest_path = self.file_dict["contigs_s"]
                p1_dest_path = self.file_dict["contigs_p1"]
                p2_dest_path = self.file_dict["contigs_p2"]
                make_map = open(bypass_contig_map_path, "w")
                make_contig = open(bypass_contig_path, "w")
                shutil.copyfile(s_src_path, s_dest_path)
                shutil.copyfile(p1_src_path, p1_dest_path)
                shutil.copyfile(p2_src_path, p2_dest_path)
                self.contigs_present = False
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.assemble_contigs_label)

            elif(spades_fail_flag and (not mgm_fail_flag)):
                print(dt.today(), "SPADes fails but MGM runs anyway: this shouldn't happen. contact admin")
                sys.exit("SPADES_fail_MGM_ok")
            elif((not spades_fail_flag) and mgm_fail_flag):
                print(dt.today(), "SPADes ran fine, but MGM failed. Check your MetaGeneMark license")
                sys.exit("SPADES_ok_MGM_fail")
            else:
                self.contigs_present = True
                print(dt.today(), "Assemble-contigs pass")
                self.marker_control.place_marker("contigs")
                #self.mp_util.write_to_bypass_log(self.output_folder_path, self.assemble_contigs_label)
            
            
            #self.mp_util.clean_or_compress(self.assemble_contigs_path, self.config_dict["keep_all"], self.config_dict["keep_assemble_contigs"])
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

        #s
        #self.assemble_contigs_end = time.time()
        #print("assemble contigs:", '%1.1f' % (self.assemble_contigs_end - self.assemble_contigs_start - (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start)), "s")    
        #print("assemble contigs cleanup:", '%1.1f' % (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start), "s")
        
        #self.debug_stop_check(self.assemble_contigs_label)
        
    def mp_GA_pre_scan(self):
        #scans the mRNA with a TA scanner to pick out a taxa trend.
        if self.marker_control.check_marker("GA_ps"):
            self.dir_control.make_dirs_from_list("GA_ps_list")

            marker_path_list = []
            #----------------------------------------------------------------------
            #kraken2 on reads
            sections = ["s"]
            if self.read_mode == "paired":
                sections.append("p1")
            if(self.contigs_present):
                sections.append("c")    
            
            print("GA pre-scan:", sections)
            for section in sections:
                marker_tag = "ga_ps_" + section
                
                if(os.path.exists(self.marker_dict[marker_tag])):
                    print(dt.today(), "skipping:", marker_tag)
                else:
                    marker_path_list.append(self.marker_dict[marker_tag])
                    command_list = self.commands.create_GA_pre_scan_taxa_command(section, self.marker_dict[marker_tag])
                    self.mp_util.run_subjob_with_hold(
                        self.config_dict["TA_mem_threshold"], 
                        self.config_dict["TA_job_limit"], 
                        self.config_dict["TA_job_delay"], 
                        self.file_dict["ga_ps_k2_" + section + "_job"],
                        
                        command_list
                    )        
            
             
            
            #-------------------------------------------------------
            #use the kraken2 results to make the DB
            
            self.mp_util.wait_for_mp_store()
            #glue all k2 reports together
            k2_reports = ["ga_ps_k2_report_s", "ga_ps_k2_report_c", "ga_ps_k2_report_p"]
            with open(self.file_dict["ga_ps_k2_report_all"], "wb") as out_file:
                for item in k2_reports:
                    with open(self.file_dict[item], "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)


            if self.marker_control.check_marker("ga_ps_make"):
                command_list = self.commands.create_GA_pre_scan_command(self.marker_dict["ga_ps_make"])
                self.mp_util.run_subjob_with_hold(
                    self.config_dict["TA_mem_threshold"], 
                    self.config_dict["TA_job_limit"], 
                    self.config_dict["TA_job_delay"], 
                    self.file_dict["ga_ps_make_job"],
                    command_list)
                
                
            
            #---------------------------------------------------------
            
            self.marker_control.place_marker("GA_ps")  
            #self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_pre_scan_label)
        #sys.exit(dt.today(), "stop here")
        #self.debug_stop_check(self.GA_pre_scan_label)
    

        
    def mp_GA_BT2(self):
        #GA BT2 reads are not split
        self.GA_BT2_start = time.time()
        #if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BT2_label):
        if self.marker_control.check_marker("GA_BT2"):
            marker_path_list = []
            lib_list = list()
            if(os.path.exists(self.config_dict["custom_ga_lib_list"])):    
                lib_list = self.file_control.read_lib_list(self.config_dict["custom_ga_lib_list"])
            else:
                lib_list = self.file_control.read_lib_list(self.file_dict["ga_lib_list"])
                print(dt.today(), "custom ga lib list doesn't exist at:", self.config_dict["custom_ga_lib_list"])
                print(dt.today(), "DEFAULTING TO chocophlan subset with GA pre-scan")
                #time.sleep(10)

            for item in lib_list:
                print("lib list item:", item)
                
            #sys.exit("stop")    
            
            if (self.marker_control.check_marker("GA_BT2")):
                self.dir_control.make_dirs_from_list("GA_BT2_list")
                #-------------------------------------------------------------------------
                #no looping. just run s, p, and c manually
                #import the lib list
                #lib list should contain the full path
                

                p_mkr = ""
                p_job = ""
                for lib_entry in lib_list:
                    lib_basename = os.path.basename(lib_entry)
                    lib_tag = lib_basename
                    if(lib_basename.endswith(".fasta")):
                        lib_tag = lib_basename.strip(".fasta")
                    elif(lib_basename.endswith(".fa")):
                        lib_tag = lib_basename.strip(".fa")
                    
                    if((self.config_dict["read_mode"]=="p") or (self.config_dict["read_mode"] == "paired")):
                        p_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_p_" + lib_tag + "_job.sh")
                        p_mkr = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_p_" + lib_tag)
                        


                    if(os.path.exists(p_mkr)):
                        print(dt.today(), "skipping:", p_mkr)
                        
                    else:
                        marker_path_list.append(p_mkr)
                    
                        #aug 10, 2021: new bigger chocophlan (from humann3) is in segments because we can't index it as a whole.  
                        #if the DB is still an old version, the tag should just say "chocophlan".  otherwise, it will say the chocophlan chunk name
                        
                        command_list = self.commands.create_GA_BT2_command(
                            lib_entry, self.file_dict["contigs_p1"], self.file_dict["contigs_p2"], 
                            self.file_dict["ga_bt2_p_sam"], p_mkr, "p")
                        #self.mp_util.run_subjob_with_hold(self.BT2_mem_threshold, self.BT2_job_limit, self.BT2_job_delay, self.GA_BT2_label, job_name, self.commands, command_list)
                        self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"],  p_job, command_list)



                    s_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_s_" + lib_tag + "_job.sh")
                    s_mkr = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_s_" + lib_tag)
                    if(os.path.exists(s_mkr)):
                        print(dt.today(), "skipping:", s_mkr)
                    
                    else:
                        print(dt.today(), "running:", s_mkr)
                        marker_path_list.append(s_mkr)
                        command_list = self.commands.create_GA_BT2_command(lib_entry, self.file_dict["contigs_s"], "none", self.file_dict["ga_bt2_s_sam"], s_mkr, "s")
                        self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"],  s_job, command_list)


                    c_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_c_" + lib_tag + "_job.sh")
                    c_mkr = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_c_" + lib_tag)
                    if(os.path.exists(c_mkr)):
                        print(dt.today(), "skipping:", c_mkr)
                        
                    else:
                        marker_path_list.append(c_mkr)
                        print(dt.today(), "running:", c_mkr)
                        command_list = self.commands.create_GA_BT2_command(lib_entry, self.file_dict["contigs_out_fa"], "none", self.file_dict["ga_bt2_c_sam"], c_mkr, "s")
                        self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"], c_job, command_list)

                    print(dt.today(), "all BT2 jobs for", lib_entry, "have been launched. waiting")            
                    self.mp_util.wait_for_mp_store()
                
                #final_checklist = os.path.join(self.GA_BT2_path, "GA_BT2.txt")
                #self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                if(self.marker_control.check_marker_list(marker_path_list)):
                    self.marker_control.place_marker("GA_BT2")
    
        #self.debug_stop_check(self.GA_BT2_label)
        
    def mp_GA_BT2_pp(self):                
        #if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BT2_pp_label):
        if self.marker_control.check_marker("GA_BT2_PP"):
           
            marker_path_list = []
            sections = ["s"]
            if self.read_mode == "p":
                sections.extend(["p"])
                
            if(self.contigs_present):
                sections.extend(["c"])
                
            lib_list = list()
            if(os.path.exists(self.config_dict["custom_ga_lib_list"])):    
                lib_list = self.file_control.read_lib_list(self.config_dict["custom_ga_lib_list"])
            else:
                lib_list = self.file_control.read_lib_list(self.file_dict["ga_lib_list"])
                print(dt.today(), "custom ga lib list doesn't exist at:", self.config_dict["custom_ga_lib_list"])
                print(dt.today(), "DEFAULTING TO chocophlan subset with GA pre-scan")
                #time.sleep(10)

            for item in lib_list:
                print("lib list item:", item)

            
            for lib_entry in lib_list:
                lib_basename = os.path.basename(lib_entry)
                lib_tag = lib_basename
                if(lib_basename.endswith(".fasta")):
                    lib_tag = lib_basename.strip(".fasta")
                elif(lib_basename.endswith(".fa")):
                    lib_tag = lib_basename.strip(".fa")


                p_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_pp_p_" + lib_tag + "_job.sh")
                p_marker = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_pp_p_" + lib_tag)
                marker_path_list.append(p_marker)
                command_list = self.commands.create_GA_BT2_pp_command(
                    lib_entry, self.file_dict["gene_map_p"], self.file_dict["genes_p"], self.file_dict["bt2_prot_p"],
                    self.file_dict["contigs_p1"], self.file_dict["contigs_p2"], self.file_dict["ga_bt2_p_sam"],
                    self.file_dict["ga_rem_p1"], self.file_dict["ga_rem_p2"], "p", p_marker
                      )
                if(os.path.exists(p_marker)):
                    print(dt.today(), "skipping: ", p_marker)
                else:
                    self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"],  p_job, command_list)


                s_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_pp_s_" + lib_tag + "_job.sh")
                s_marker = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_pp_s_" + lib_tag)
                marker_path_list.append(s_marker)
                command_list = self.commands.create_GA_BT2_pp_command(
                    lib_entry, self.file_dict["gene_map_s"], self.file_dict["genes_s"], self.file_dict["bt2_prot_s"],
                    self.file_dict["contigs_s"], "None", self.file_dict["ga_bt2_s_sam"],
                    self.file_dict["ga_rem_s"], "None", "s", s_marker
                )
                if(os.path.exists(s_marker)):
                    print(dt.today(), "skipping: ", s_marker)
                else:
                    self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"],  s_job, command_list)
                
                c_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_pp_c_" + lib_tag + "_job.sh")
                c_marker = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_pp_c_" + lib_tag)
                marker_path_list.append(c_marker)
                command_list = self.commands.create_GA_BT2_pp_command(
                    lib_entry, self.file_dict["gene_map_c"], self.file_dict["genes_c"], self.file_dict["bt2_prot_c"],
                    self.file_dict["contigs_out_fa"], "None", self.file_dict["ga_bt2_c_sam"],
                    self.file_dict["ga_rem_c"], "None", "c", c_marker
                )
                if(os.path.exists(c_marker)):    
                    print(dt.today(), "skipping: ", c_marker)
                else:
                    self.mp_util.run_subjob_with_mem_footprint(self.config_dict["BT2_mem_footprint"], self.config_dict["BT2_job_limit"],  c_job, command_list)
                
                    print(dt.today(), "BT2 pp batch launched for: ", lib_entry)
                self.mp_util.wait_for_mp_store()

            print(dt.today(), "all BT2 PP jobs submitted.  waiting for sync")            
            self.mp_util.wait_for_mp_store()


            if(self.marker_control.check_marker_list(marker_path_list)):
                self.marker_control.place_marker("GA_BT2_PP")
            else:
                print(dt.today(), "not all markers finished")
                sys.exit()


            #self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            #self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BT2_pp_label)
        
        self.debug_stop_check("GA_BT2_pp")
        
        #there used to be a merge step, but it's no longer necessary with BT2 and the scaled-down parallelization
    
    def mp_GA_dmd(self):
        # ------------------------------------------------------
        # Diamond gene annotation
        self.GA_DIAMOND_start = time.time()
        #GA_DIAMOND_tool_output_path = os.path.join(self.GA_DIAMOND_path, "data", "0_diamond")
        #if not check_where_resume(None, self.GA_DIAMOND_tool_output_path, self.GA_BLAT_path, file_check_bypass = True):        
        #jul 21, 2025: it's been determined that sequential is the way to go, given that we can increase the block-size of the run
        #yes, it's true.  

        if(self.marker_control.check_marker("GA_DMD")):
            self.dir_control.make_dirs_from_list("GA_DMD_list")
            
            marker_path_list = []
            #remember: don't split the DMD jobs.  there's no consistency across mem use on input size.
            if(self.marker_control.check_marker("GA_DMD_p1")):
                block_size = self.mp_util.determine_dmd_mem_limit(self.file_dict["ga_rem_p1"])
                marker_path_list.append(self.marker_dict["GA_DMD_p1"])
                command_list = self.commands.create_GA_DMD_command(self.file_dict["ga_rem_p1"], self.file_dict["ga_dmd_p1_dmdout"], block_size, self.marker_dict["GA_DMD_p1"])
                self.mp_util.run_subjob_with_hold(self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], self.config_dict["DMD_job_delay"], self.file_dict["ga_dmd_p1_job"], command_list)
                

            if(self.marker_control.check_marker("GA_DMD_p2")):
                block_size = self.mp_util.determine_dmd_mem_limit(self.file_dict["ga_rem_p2"])
                marker_path_list.append(self.marker_dict["GA_DMD_p2"])
                command_list = self.commands.create_GA_DMD_command(self.file_dict["ga_rem_p2"], self.file_dict["ga_dmd_p2_dmdout"], block_size, self.marker_dict["GA_DMD_p2"])
                self.mp_util.run_subjob_with_hold(self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], self.config_dict["DMD_job_delay"], self.file_dict["ga_dmd_p2_job"], command_list)
                
            if(self.marker_control.check_marker("GA_DMD_c")):
                block_size = self.mp_util.determine_dmd_mem_limit(self.file_dict["ga_rem_c"])
                marker_path_list.append(self.marker_dict["GA_DMD_c"])
                command_list = self.commands.create_GA_DMD_command(self.file_dict["ga_rem_c"], self.file_dict["ga_dmd_c_dmdout"], block_size, self.marker_dict["GA_DMD_c"])
                self.mp_util.run_subjob_with_hold(self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], self.config_dict["DMD_job_delay"], self.file_dict["ga_dmd_c_job"], command_list)
                

            if(self.marker_control.check_marker("GA_DMD_s")):
                marker_path_list.append(self.marker_dict["GA_DMD_s"])
                block_size = self.mp_util.determine_dmd_mem_limit(self.file_dict["ga_rem_s"])
                command_list = self.commands.create_GA_DMD_command(self.file_dict["ga_rem_s"], self.file_dict["ga_dmd_s_dmdout"], block_size, self.marker_dict["GA_DMD_s"])
                self.mp_util.run_subjob_with_hold(self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], self.config_dict["DMD_job_delay"], self.file_dict["ga_dmd_s_job"], command_list)
                
            print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
            self.mp_util.wait_for_mp_store()
            if(self.marker_control.check_marker_list(marker_path_list)):
                self.marker_control.place_marker("GA_DMD")
            else:
                print(dt.today(), "DMD: not all markers finished")
                sys.exit()
        
        
    def mp_GA_dmd_pp(self):        
        #if not check_where_resume(GA_DIAMOND_path, None, self.GA_DIAMOND_tool_output_path, file_check_bypass = True):
        if(self.marker_control.check_marker("GA_DMD_pp")):
            #print(dt.today(), "DIAMOND PP threads used:", self.config_dict["num_threads/2)
            marker_path_list = []
            if(self.marker_control.check_marker("GA_DMD_pp_p1")):

                marker_path_list.append(self.marker_dict["GA_DMD_pp_p1"])
                command_list = self.commands.create_DIAMOND_pp_command_v2(
                    self.file_dict["ga_rem_p1"], self.file_dict["ga_dmd_p1_dmdout"], 
                    self.file_dict["ga_dmd_rem_p1"], self.file_dict["ga_dmd_prot_p1"],
                    self.file_dict["ga_dmd_prot_map_p1"], self.marker_dict["GA_DMD_pp_p1"]
                    )
                
                self.mp_util.run_subjob_with_mem_footprint(
                    self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], 
                    self.file_dict["ga_dmd_pp_p1_job"], command_list
                    )
                #self.mp_util.wait_for_mp_store()
                
            if(self.marker_control.check_marker("GA_DMD_pp_p2")):

                marker_path_list.append(self.marker_dict["GA_DMD_pp_p2"])
                command_list = self.commands.create_DIAMOND_pp_command_v2(
                    self.file_dict["ga_rem_p2"], self.file_dict["ga_dmd_p2_dmdout"], 
                    self.file_dict["ga_dmd_rem_p2"], self.file_dict["ga_dmd_prot_p2"],
                    self.file_dict["ga_dmd_prot_map_p2"], self.marker_dict["GA_DMD_pp_p2"]
                    )
                
                self.mp_util.run_subjob_with_mem_footprint(
                    self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], 
                    self.file_dict["ga_dmd_pp_p2_job"], command_list
                    )
                #self.mp_util.wait_for_mp_store()

            if(self.marker_control.check_marker("GA_DMD_pp_s")):

                marker_path_list.append(self.marker_dict["GA_DMD_pp_s"])
                command_list = self.commands.create_DIAMOND_pp_command_v2(
                    self.file_dict["ga_rem_s"], self.file_dict["ga_dmd_s_dmdout"], 
                    self.file_dict["ga_dmd_rem_s"], self.file_dict["ga_dmd_prot_s"],
                    self.file_dict["ga_dmd_prot_map_s"], self.marker_dict["GA_DMD_pp_s"]
                    )
                
                self.mp_util.run_subjob_with_mem_footprint(
                    self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], 
                    self.file_dict["ga_dmd_pp_s_job"], command_list
                    )    
                #self.mp_util.wait_for_mp_store()
                


            if(self.marker_control.check_marker("GA_DMD_pp_c")):

                marker_path_list.append(self.marker_dict["GA_DMD_pp_c"])
                command_list = self.commands.create_DIAMOND_pp_command_v2(
                    self.file_dict["ga_rem_c"], self.file_dict["ga_dmd_c_dmdout"], 
                    self.file_dict["ga_dmd_rem_c"], self.file_dict["ga_dmd_prot_c"], 
                    self.file_dict["ga_dmd_prot_map_c"], self.marker_dict["GA_DMD_pp_c"]
                    )
                
                self.mp_util.run_subjob_with_mem_footprint(
                    self.config_dict["DMD_mem_footprint"], self.config_dict["DMD_job_limit"], 
                    self.file_dict["ga_dmd_pp_c_job"], command_list
                    )   
                #self.mp_util.wait_for_mp_store()
                             
            print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
            self.mp_util.wait_for_mp_store()

            gene_map_list = [self.file_dict["ga_dmd_prot_map_s"], self.file_dict["ga_dmd_prot_map_c"], self.file_dict["ga_dmd_prot_map_p1"], self.file_dict["ga_dmd_prot_map_p2"]]
            #self.mp_util.concatenate_files_efficient(gene_map_list, self.file_dict["ga_dmd_prot_map"])

            prot_list = [self.file_dict["ga_dmd_prot_s"], self.file_dict["ga_dmd_prot_c"], self.file_dict["ga_dmd_prot_p1"], self.file_dict["ga_dmd_prot_p2"]]
            #self.mp_util.concatenate_files_efficient(prot_list, self.file_dict["ga_dmd_prot"])


            if(self.marker_control.check_marker_list(marker_path_list)):
                self.marker_control.place_marker("GA_DMD_pp")
            else:
                print(dt.today(), "Not all DMD pp markers finished")
                sys.exit()
            
            
        #self.debug_stop_check(self.GA_DIAMOND_pp_label)
        
    def mp_GA_final_merge(self):
        
        self.GA_final_merge_start = time.time()
        #if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_final_merge_label):
        
        if self.marker_control.check_marker("GA_FM"):
            self.dir_control.make_dirs_from_list("GA_FM_list")
            marker_path_list = [self.marker_dict["GA_FM_bt2"], self.marker_dict["GA_FM_dmd"],
                self.marker_dict["GA_FM_prot"], self.marker_dict["GA_FM_maps"]]
            command_list = self.commands.create_GA_final_merge_command(
                self.marker_dict["GA_FM_bt2"], self.marker_dict["GA_FM_dmd"],
                self.marker_dict["GA_FM_prot"], self.marker_dict["GA_FM_maps"]
                )
            job_name = "GA_final_merge"
            self.mp_util.subdivide_and_launch(
                self.config_dict["GA_final_merge_job_delay"], self.config_dict["GA_final_merge_mem_threshold"], 
                self.config_dict["GA_final_merge_job_limit"], self.file_dict["ga_fm_job"], command_list)
            if self.marker_control.check_marker_list(marker_path_list):
                self.marker_control.place_marker("GA_FM")

            self.mp_util.wait_for_mp_store()
     
        self.GA_final_merge_end = time.time()
        print("GA final merge:", '%1.1f' % (self.GA_final_merge_end - self.GA_final_merge_start), "s")
        #self.mp_util.clean_or_compress(self.ga_final_merge_path, self.config_dict["keep_all"], self.keep_GA_final)
        
        #self.debug_stop_check(self.GA_final_merge_label)

    def mp_TA(self):
        #this is a bit of a misnomer now.  It just repacks the kraken2 report.  We already did a TA scan to narrow down the library.
        self.TA_start = time.time()
        if self.marker_control.check_marker("TA"):
            self.dir_control.make_dirs_from_list("TA_list")
            command_list = self.commands.create_TA_repack_command(self.marker_dict["TA"])
            self.mp_util.run_subjob_with_mp_store(self.file_dict["ta_job"], command_list)
            self.mp_util.wait_for_mp_store()
     
        self.TA_end = time.time()
    

    def mp_EC(self):
        
        self.EC_start = time.time()
        if self.marker_control.check_marker("EC"):
            self.dir_control.make_dirs_from_list("EC_list")
            command_list = self.commands.create_deepec_command(self.marker_dict["EC"])
            self.mp_util.run_subjob_with_mp_store(self.file_dict["ec_job"], command_list)
            self.mp_util.wait_for_mp_store()
        self.cleanup_EC_start = time.time()
        #self.mp_util.clean_or_compress(self.ec_path, self.config_dict["keep_all"], self.keep_EC)
        self.cleanup_EC_end = time.time()
        self.EC_post_end = time.time()
            
    
        self.EC_end = time.time()
        print("EC run:", '%1.1f' % (self.EC_end - self.EC_start), "s")
        print("EC cleanup:", '%1.1f' % (self.cleanup_EC_end - self.cleanup_EC_start), "s")
        
        #self.debug_stop_check(self.ec_label)

    def mp_output(self):

        if self.marker_control.check_marker("Out"):
            self.dir_control.make_dirs_from_list("out_list")
            if self.marker_control.check_marker("Out_rpkm"):
                command_list = self.commands.create_rpkm_command(self.marker_dict["Out_rpkm"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_rpkm_job"], command_list)
            
            if self.marker_control.check_marker("Out_repop"):
                command_list = self.commands.create_output_unique_junk()
                self.mp_util.subdivide_and_launch(10, 50, int(os.cpu_count()), self.file_dict["out_repop_job"], command_list)
            
            if self.marker_control.check_marker("Out_per_read"):
                command_list = self.commands.create_output_per_read_scores_command(self.marker_dict["Out_per_read"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_per_read_job"], command_list)

            if self.marker_control.check_marker("Out_taxa_report"):
                command_list = self.commands.create_output_copy_taxa_command(self.marker_dict["Out_taxa_report"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_taxa_job"], command_list)
            
            if self.marker_control.check_marker("Out_contig_stats"):
                command_list = self.commands.create_output_contig_stats_command(self.marker_dict["Out_contig_stats"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_contig_stats_job"], command_list)

            
            
            
            if self.marker_control.check_marker("Out_taxa_groupby"):
                command_list = self.commands.create_output_taxa_groupby_command(self.marker_dict["Out_taxa_groupby"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_taxa_groupby_job"], command_list)

            self.mp_util.wait_for_mp_store()

            if self.marker_control.check_marker("Out_read_count"):
                command_list = self.commands.create_output_read_count_command(self.marker_dict["Out_read_count"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_read_count_job"], command_list)

            
            self.mp_util.wait_for_mp_store()

            if self.marker_control.check_marker("Out_ec_heatmap"):
                command_list = self.commands.create_output_EC_heatmap_command(self.marker_dict["Out_ec_heatmap"])
                self.mp_util.run_subjob_with_mp_store(self.file_dict["out_ec_heatmap_job"], command_list)

            self.mp_util.wait_for_mp_store()

            



            
    
