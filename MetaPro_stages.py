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
        self.seq_handler = mpu.mp_seq_handler(self.config_dict, self.dir_dict, self.file_dict)
        

        #time.sleep(10)
        
        
        
        self.config_dict["no_host"] = config_dict["no_host"]
        self.verbose_mode = config_dict["verbose_mode"]
        self.rRNA_chunks = int(self.config_dict["rRNA_chunksize"])
        self.EC_chunksize = int(self.config_dict["EC_chunksize"])
        self.GA_chunksize = int(self.config_dict["GA_chunksize"])
        #self.config_path = config_path
        self.pair_1_path = config_dict["pair_1"]
        self.pair_2_path = config_dict["pair_2"]
        self.single_path = config_dict["single"]
        self.contig_path = config_dict["contig"] #_path  #tutorial/single-shot use
        self.quality_encoding = ""
        self.read_mode = self.config_dict["read_mode"]
        if (self.single_path != "None"):
            self.config_dict["q_enc"] = self.mp_util.determine_encoding(self.single_path)
            print("ENCODING USED:", self.config_dict["q_enc"])
        else:
            self.config_dict["q_enc"] = self.mp_util.determine_encoding(self.pair_1_path)
            print("ENCODING USED:", self.config_dict["q_enc"])
        
        
        #self.debug_stop_flag = self.paths.debug_stop_flag
        
        mp_store = []  # stores the multiprocessing processes

        # Creates our command object, for creating shellscripts.

        self.commands = mpcom.mt_pipe_commands(self.config_dict, self.dir_dict, self.file_dict)
        
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
            exit_string = "stoppped after: " + stop_signal
            sys.exit(exit_string)
        else:
            print(dt.today(), "continuing from:", stop_signal)
            
    #--------------------------------------------------------------------------------------------------------------
    # main calls
    def mp_quality_filter(self):
        if(self.marker_control.check_marker(self.marker_dict["qf"])):
            self.time_control.measure_time("qc", "start")
            for item in self.dir_dict["qc_list"]:
                self.dir_control.make_dirs(self.dir_dict[item])

            command_list = self.commands.create_quality_control_command(self.marker_dict["qf"])
            self.mp_util.launch_stage_simple(self.label_dict["qc"], self.dir_dict["qc"], self.commands, command_list, self.keep_all, self.keep_quality)
            self.time_control.measure_time("qc", "end")
            
            self.debug_stop_check(self.quality_filter_label)
            

    def mp_host_filter(self):
        if not self.config_dict["no_host"]:
            if(self.marker_control.check_marker(self.marker_dict["host"])):
                self.time_control.measure_time("host", "start")
                for item in self.dir_dict["host_list"]:
                    self.dir_control.make_dirs(self.dir_dict[item])
                command_list = self.commands.create_host_filter_command(self.marker_dict["host"])
                self.mp_util.launch_stage_simple(self.host_filter_label, self.host_path, self.commands, command_list, self.keep_all, self.keep_host)
                self.time_control.measure_time("host", "end")

                self.debug_stop_check(self.host_filter_label)

    def mp_vector_filter(self):
        self.vector_start = time.time()
        if(self.marker_control.check_marker(self.marker_dict["vec"])):
            if self.config_dict["no_host"]:
                #get dep args from quality filter
                #if not check_where_resume(vector_path, None, self.quality_path):
                command_list = self.commands.create_vector_filter_command(self.vector_filter_label, self.quality_filter_label)
                self.cleanup_vector_start, self.cleanup_vector_end = self.mp_util.launch_stage_simple(self.vector_filter_label, self.vector_path, self.commands, command_list, self.keep_all, self.keep_vector)

            else:
                #get the dep args from host filter
                #if not check_where_resume(vector_path, None, self.host_path):
                command_list = self.commands.create_vector_filter_command(self.vector_filter_label, self.host_filter_label)
                self.cleanup_vector_start, self.cleanup_vector_end = self.mp_util.launch_stage_simple(self.vector_filter_label, self.vector_path, self.commands, command_list, self.keep_all, self.keep_vector)
                
            self.vector_end = time.time()
            print("vector filter:", '%1.1f' % (self.vector_end - self.vector_start - (self.cleanup_vector_end - self.cleanup_vector_start)), "s")
            print("vector filter cleanup:", '%1.1f' % (self.cleanup_vector_end - self.cleanup_vector_start), "s")
            self.debug_stop_check(self.vector_filter_label)

    def mp_rRNA_filter(self):
        #don't split bnap.
        #only split infernal.  but it can go over the 100k limit now.
        #bnap takes fasta, but it's not worth it to keep it in memory.


        self.rRNA_filter_start = time.time()
        
        #if not check_where_resume(self.rRNA_filter_path, None, self.vector_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.label_dict["rRNA"]): 
            #run bnap
            run_types = ["s"]
            if(self.read_mode == "p"):
                run_types.extend(["p1", "p2"])

            for run_type in run_types:
                fa_in = self.file_dict["rRNA_"+str(run_type) + "_fa"]
                fq_in = self.file_dict["no_vec_" + str(run_type)]
                mRNA_out = self.file_dict["rRNA_bnap_mRNA_" + str(run_type)]
                rRNA_out = self.file_dict["rRNA_bnap_rRNA_" + str(run_type)]
                bnap_out = self.file_dict["rRNA_bnap_all_" + str(run_type)]
                self.mp_seq_handler.fastq_to_fasta(fa_in, fq_in)
                command = self.commands.create_rRNA_filter_bnap_command(fa_in, fq_in, mRNA_out, rRNA_out, bnap_out)
                self.mp_util.run_subjob_with_mp_store()


            #wait for everything.
            self.mp_util.wait_for_mp_store()

                
            #----------------------------------------------------------------------------
            # INFERNAL
            #split the fasta mRNA here
            split_count_s = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_s"], self.file_dict["rRNA_inf_in_s"], self.config_dict["rRNA_chunksize"], "fasta")
            split_count_p1 = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_p1"], self.file_dict["rRNA_inf_in_p1"], self.config_dict["rRNA_chunksize"], "fasta")
            split_count_p2 = self.mp_seq_handler.split_fastq(self.file_dict["rRNA_bnap_mRNA_p2"], self.file_dict["rRNA_inf_in_p2"], self.config_dict["rRNA_chunksize"], "fasta")

            for i in range(0, split_count_s):
                self.file_dict["rRNA_inf_in_s_" + str(i)] = self.file_dict["rRNA_inf_in_s"] + "_" + str(i) + ".fasta"
            for i in range(0, split_count_p1):
                self.file_dict["rRNA_inf_in_p1_" + str(i)] = self.file_dict["rRNA_inf_in_p1"] + "_" + str(i) + ".fasta"
            for i in range(0, split_count_p2):
                self.file_dict["rRNA_inf_in_p2_" + str(i)] = self.file_dict["rRNA_inf_in_p2"] + "_" + str(i) + ".fasta"

            self.marker_control.issue_rRNA_markers(self.m_name["rRNA_inf_s"], split_count_s, "rRNA_mkrs")
            self.marker_control.issue_rRNA_markers(self.m_name["rRNA_inf_p1"], split_count_p1, "rRNA_mkrs")
            self.marker_control.issue_rRNA_markers(self.m_name["rRNA_inf_p2"], split_count_p2, "rRNA_mkrs")

            for i in range(0, split_count_s):
                s_seq = self.file_dict["rRNA_inf_in_s" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_s"] + "_" + str(i) + ".inf_out"
                self.file_dict["rRNA_inf_s_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_s_" + str(i)]
                job_name = "rRNA_inf_s" + "_" + str(i)
                if(self.marker_control.check_marker(marker)):
                    command = self.commands.create_rRNA_filter_infernal_command(s_seq, inf_out, marker)
                    self.mp_util.run_subjob_with_hold(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.dir_dict["rRNA_jobs"],
                        job_name,
                        self.commands,
                        command)    


            for i in range(0, split_count_p1):
                p1_seq = self.file_dict["rRNA_bnap_mRNA_p1_" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_p1"] + "_" + str(i) + ".inf_out"
                self.file_dict["rRNA_inf_p1_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_p1_" + str(i)]
                job_name = "rRNA_inf_p1" + "_" + str(i)
                if(self.marker_control.check_marker(marker)):
                    command = self.commands.create_rRNA_filter_infernal_command(p1_seq, inf_out, marker)
                    self.mp_util.run_subjob_with_hold(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.dir_dict["rRNA_jobs"],
                        job_name,
                        self.commands,
                        command)    

            for i in range(0, split_count_p2):
                p2_seq = self.file_dict["rRNA_bnap_mRNA_p2_" + str(i)]
                inf_out = self.file_dict["rRNA_inf_out_p2"] + "_" + str(i) + ".inf_out"
                self.file_dict["rRNA_inf_p2_" + str(i)] = inf_out
                marker = self.marker_dict["rRNA_inf_p2_" + str(i)]
                job_name = "rRNA_inf_p2" + "_" + str(i)
                if(self.marker_control.check_marker(marker)):
                    command = self.commands.create_rRNA_filter_infernal_command(p2_seq, inf_out, marker)
                    self.mp_util.run_subjob_with_hold(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.dir_dict["rRNA_jobs"],
                        job_name,
                        self.commands,
                        command)    

            #wait for everything.
            self.mp_util.wait_for_mp_store()
            #-----------------------------------------------
            #merge the inf and bnap files, then send it through inf pp.
            if(split_count_s > 0):
                with open(self.file_dict["rRNA_inf_all_s"], "wb") as out_file:
                    for i in range(0, split_count_s):
                        inf_out = self.file_dict["rRNA_inf_s_" + str(i)]
                        with open(inf_out, "rb") as in_file:
                            shutil.copyfileobj(in_file, out_file)
                        out_file.write(b"\n")
                

            if(split_count_p1 > 0):
                with open(self.file_dict["rRNA_inf_all_p1"], "wb") as out_file:
                    for i in range(0, split_count_p1):
                        inf_out = self.file_dict["rRNA_inf_p1_" + str(i)]
                        with open(inf_out, "rb") as in_file:
                            shutil.copyfileobj(in_file, out_file)
                        out_file.write(b"\n")
                

            if(split_count_p2 > 0):
                with open(self.file_dict["rRNA_inf_all_p2"], "wb") as out_file:
                    for i in range(0, split_count_p2):
                        inf_out = self.file_dict["rRNA_inf_p2_" + str(i)]
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
                
                self.mp_util.run_subjob_with_hold(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.dir_dict["rRNA_jobs"],
                        "rRNA_inf_pp_paired",
                        self.commands,
                        command) 


            command = self.commands.create_rRNA_inf_pp_s_command(self.file_dict["rRNA_inf_all_s"],
                                                                 self.file_dict["rRNA_bnap_all_s"],
                                                                 self.file_dict["no_vec_s"],
                                                                 self.file_dict["rRNA_mRNA_s_fq"],
                                                                 self.file_dict["rRNA_other_s_fq"],
                                                                 self.marker_dict["rRNA_inf_pp_s"]
                                                                 )
            self.mp_util.run_subjob_with_hold(
                        self.config_dict["Infernal_mem_threshold"],
                        self.config_dict["Infernal_job_limit"],
                        self.config_dict["Infernal_job_delay"],
                        self.dir_dict["rRNA_jobs"],
                        "rRNA_inf_pp_s",
                        self.commands,
                        command) 
            self.mp_util.wait_for_mp_store()

#-----------------------------------------------------------------------------------------------------------------------        

    def mp_repop(self):
        
        self.repop_start = time.time()
        #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.label_dict["repop"]):
            job_name = self.repop_job_label
            command_list = self.commands.create_repop_command(self.marker_dict["repop"])
            self.mp_util.subdivide_and_launch(self.repop_job_delay, self.repop_mem_threshold, self.repop_job_limit, self.repop_job_label, job_name, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            
        self.repop_end = time.time()


    def mp_assemble(self):
        self.assemble_contigs_start = time.time()
        
        #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
        mgm_gene_report = self.file_dict["contigs_gene_report"]
        spades_done_file = self.file_dict["contigs_spades_done"]
        spades_transcript_file = self.file_dict["contigs_transcripts"]
        mgm_fail_flag = True
        spades_fail_flag = True

        if self.mp_util.check_bypass_log(self.output_folder_path, self.label_dict["contigs"]):
            job_name = self.assemble_contigs_label
            command_list = self.commands.create_assemble_contigs_command(self.marker_dict["contigs"])
            self.mp_util.run_subjob_simple(self.assemble_contigs_label, job_name, self.commands, command_list)
            
            if(os.path.exists(spades_done_file)):
                if(os.path.exists(spades_transcript_file)):
                    spades_fail_flag = False
                    print(dt.today(), "SPADes OK")
                else:
                    spades_fail_flag = True
                    print(dt.today(), "SPADes ran, but did not create contigs")
            else:
                sys.exit(dt.today(), "SPADes did not run. this is not normal. Contact admin immediately")

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
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.assemble_contigs_label)
            
            self.cleanup_assemble_contigs_start = time.time()
            self.mp_util.clean_or_compress(self.assemble_contigs_path, self.keep_all, self.keep_assemble_contigs)
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
                
        #self.assemble_contigs_end = time.time()
        #print("assemble contigs:", '%1.1f' % (self.assemble_contigs_end - self.assemble_contigs_start - (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start)), "s")    
        #print("assemble contigs cleanup:", '%1.1f' % (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start), "s")
        
        #self.debug_stop_check(self.assemble_contigs_label)
    
    def mp_GA_pre_scan(self):
        #scans the mRNA with a TA scanner to pick out a taxa trend.
        if self.mp_util.check_bypass_log(self.output_folder_path, self.label_dict["GA_pre_scan"]):
            marker_path_list = []
            #----------------------------------------------------------------------
            #kraken2 on reads
            sections = ["s"]
            if self.read_mode == "p":
                sections.extend(["p"])
            if(self.contigs_present):
                sections.extend(["c"])    
            
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
                        self.dir_dict["GA_ps"], 
                        self.marker_dict[marker_tag], 
                        self.commands, 
                        command_list
                    )        
            
           
            
            #-------------------------------------------------------
            #use the kraken2 results to make the DB
            
            marker_file = "ga_collect_db"
            marker_path = os.path.join(self.GA_pre_scan_jobs_folder, marker_file)
            #glue all k2 reports together
            k2_reports = ["ga_ps_k2_report_s", "ga_ps_k2_report_c", "ga_ps_k2_report_p"]
            with open(self.file_dict["ga_ps_k2_report_all"], "wb") as out_file:
                for item in k2_reports:
                    with open(self.file_dict[item], "rb") as in_file:
                        shutil.copyfileobj(in_file, out_file)



            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_GA_pre_scan_command(marker_file)
                self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.GA_pre_scan_label, marker_file, self.commands, command_list)
                print(dt.today(), "running:", marker_file)
            self.mp_util.wait_for_mp_store()
            

            
            
            #---------------------------------------------------------
            
            
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_pre_scan_label)
            
        self.debug_stop_check(self.GA_pre_scan_label)
    

        
    def mp_GA_BT2(self):
        #GA BT2 reads are not split
        self.GA_BT2_start = time.time()
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BT2_label):
            marker_path_list = []
            
            
            if not self.mp_util.check_where_resume(self.GA_BT2_path, None, self.GA_split_path):
            
                #-------------------------------------------------------------------------
                #no looping. just run s, p, and c manually
                #import the lib list
                #lib list should contain the full path
                lib_list = list()
                if(self.config_dict["GA_DB_mode"] != "choco"):
                    if(os.path.exists(self.config_dict["custom_ga_lib_list"])):    
                        lib_list = self.import_ga_lib_list(self.config_dict["custom_ga_lib_list"])
                    else:
                        print(dt.today(), "custom ga lib list doesn't exist at:", self.config_dict["custom_ga_lib_list"])
                        print("must use absolute path")
                        sys.exit()
                else:
                    lib_list = self.import_ga_lib_list(self.file_dict["ga_lib_list"])

                
                for lib_entry in lib_list:
                    lib_basename = os.path.basename(lib_entry)
                    lib_tag = lib_basename
                    if(lib_basename.endswith(".fasta")):
                        lib_tag = lib_basename.strip(".fasta")
                    elif(lib_basename.endswith(".fa")):
                        lib_tag = lib_basename.strip(".fa")
                    
                    if(self.config_dict["op_mode"] == "paired")
                    p_job = os.path.join(self.dir_dict["GA_BT2_jobs"], "GA_BT2_p_" + lib_tag + "_job.sh")
                    p_marker = os.path.join(self.dir_dict["GA_BT2_mkrs"], "GA_BT2_p_" + lib_tag)
                    if(os.path.exists(p_marker)):
                        print(dt.today(), "skipping:", p_marker)
                        continue
                    else:
                        marker_path_list.append(p_marker)
                    
                        #aug 10, 2021: new bigger chocophlan (from humann3) is in segments because we can't index it as a whole.  
                        #if the DB is still an old version, the tag should just say "chocophlan".  otherwise, it will say the chocophlan chunk name
                        
                        command_list = self.commands.create_BT2_annotate_command_v2(lib_entry, self.file_dict["contigs_p1"], self.file_dict["contigs_p2"], self.file_dict["ga_bt2_p_sam"], p_marker, "p")
                        #self.mp_util.run_subjob_with_hold(self.BT2_mem_threshold, self.BT2_job_limit, self.BT2_job_delay, self.GA_BT2_label, job_name, self.commands, command_list)
                        self.mp_util.run_subjob_with_mem_footprint(self.BT2_mem_footprint, self.BT2_job_limit, self.GA_BT2_label, p_job, self.commands, command_list)
                        


                                    

                print(dt.today(), "all BT2 jobs have launched.  waiting for them to finish")            
                self.mp_util.wait_for_mp_store()
                final_checklist = os.path.join(self.GA_BT2_path, "GA_BT2.txt")
                self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BT2_label)
    
        self.debug_stop_check(self.GA_BT2_label)
        
    def mp_GA_BT2_pp(self):                
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BT2_pp_label):
            marker_path_list = []
            sections = ["s"]
            if self.read_mode == "p":
                sections.extend(["p"])
                
            if(self.contigs_present):
                sections.extend(["c"])
                
            lib_list = list()
            if(self.config_dict["GA_DB_mode"] != "choco"):
                if(os.path.exists(self.config_dict["custom_ga_lib_list"])):    
                    lib_list = self.import_ga_lib_list(self.config_dict["custom_ga_lib_list"])
                else:
                    print(dt.today(), "custom ga lib list doesn't exist at:", self.config_dict["custom_ga_lib_list"])
                    print("must use absolute path")
                    sys.exit()
            else:
                lib_list = self.import_ga_lib_list(self.file_dict["ga_lib_list"])

            
            for lib_entry in lib_list:
                lib_basename = os.path.basename(lib_entry)
                lib_tag = lib_basename
                if(lib_basename.endswith(".fasta")):
                    lib_tag = lib_basename.strip(".fasta")
                elif(lib_basename.endswith(".fa")):
                    lib_tag = lib_basename.strip(".fa")

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
                
                    
                        job_name = "BT2_pp" + "_" + file_tag + "_" + ref_tag
                        marker_file = file_tag + "_" + ref_tag +  "_BT2_pp"
                        marker_path = os.path.join(self.GA_BT2_jobs_folder, marker_file)
                        
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "skipping:", marker_file)
                            continue
                        else:
                            marker_path_list.append(marker_path)
                            command_list = self.commands.create_BT2_pp_command_v2(self.GA_BT2_label, self.assemble_contigs_label, ref_tag, ref_path, full_sample_path, marker_file)
                            self.mp_util.run_subjob_with_hold(self.BT2_pp_mem_threshold, self.BT2_pp_job_limit, self.BT2_pp_job_delay, self.GA_BT2_label, job_name, self.commands, command_list)
                            
                    else:
                        #chocophlan in chunks
                        split_db = os.listdir(ref_path)
                        for db_segments in split_db:
                            if(db_segments.endswith(".fasta")):
                                segment_ref_path = os.path.join(ref_path, db_segments)
                                ref_tag = db_segments.strip(".fasta")
                                job_name = "BT2_pp" + "_" + file_tag + "_" + ref_tag
                                marker_file = file_tag + "_" + ref_tag + "_BT2_pp"
                                marker_path = os.path.join(self.GA_BT2_jobs_folder, marker_file)
                                
                                if(os.path.exists(marker_path)):
                                    print(dt.today(), "skipping:", marker_file)
                                    continue
                                else:
                                    marker_path_list.append(marker_path)
                                    command_list = self.commands.create_BT2_pp_command_v2(self.GA_BT2_label, self.assemble_contigs_label, ref_tag, segment_ref_path, full_sample_path, marker_file)
                                    #print(dt.today(), "segmented BT2:", command_list)
                                    #time.sleep(2)
                                    self.mp_util.run_subjob_with_hold(self.BT2_pp_mem_threshold, self.BT2_pp_job_limit, self.BT2_pp_job_delay, self.GA_BT2_label, job_name, self.commands, command_list)

                            
            print(dt.today(), "all BT2 PP jobs submitted.  waiting for sync")            
            self.mp_util.wait_for_mp_store()
            marker_file = "BT2_copy_contig_map"
            marker_path = os.path.join(self.GA_BT2_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:   
                marker_path_list.append(marker_path)
                command_list = self.commands.create_BT2_copy_contig_map_command(self.GA_BT2_label, self.assemble_contigs_label, marker_file)
                self.mp_util.run_subjob_simple(self.GA_BT2_label, self.GA_BT2_label + "_copy_contig_map", self.commands, command_list)

            
            final_checklist = os.path.join(self.GA_BT2_path, "GA_BT2_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BT2_pp_label)
        
        self.debug_stop_check("GA_BT2_pp")
        
    def mp_GA_BT2_merge(self):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.GA_BT2_merge_label):
            #merge 
            marker_path_list = []
            sections = ["s"]
            if self.read_mode == "p":
                sections.extend(["p"])
            if(self.contigs_present):
                sections.extend(["c"])
            
            for section in sections:
                for split_sample in os.listdir(os.path.join(self.GA_split_path, "final_results", section)):
                    full_sample_path = os.path.join(os.path.join(self.GA_split_path, "final_results",section, split_sample))
                    print("split sample:", full_sample_path)
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    ref_path = self.paths.DNA_DB

                    marker_file = file_tag + "_merge_fasta"
                    marker_path = os.path.join(self.GA_BT2_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        job_name = "BT2_fasta_merge_" + file_tag
                        command_list = self.commands.create_merge_BT2_fasta_command(self.GA_BT2_label, full_sample_path, marker_file)
                        self.mp_util.run_subjob_with_hold(self.BT2_pp_mem_threshold, self.BT2_pp_job_limit, self.BT2_pp_job_delay, self.GA_BT2_label, job_name, self.commands, command_list)

            print(dt.today(), "All BT2 merge jobs have launched. waiting for sync")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.GA_BT2_path, "GA_BT2_merge.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.GA_BT2_merge_label)
            
     
        self.cleanup_GA_BT2_start = time.time()
        self.mp_util.delete_folder_simple(self.GA_BT2_jobs_folder)
        self.mp_util.clean_or_compress(self.GA_BT2_path, self.keep_all, self.keep_GA_BT2)
        
        self.cleanup_GA_BT2_end = time.time()
        self.GA_BT2_end = time.time()
        print("GA BT2:", '%1.1f' % (self.GA_BT2_end - self.GA_BT2_start - (self.cleanup_GA_BT2_end - self.cleanup_GA_BT2_start)), "s")
        print("GA BT2 cleanup:", '%1.1f' % (self.cleanup_GA_BT2_end - self.cleanup_GA_BT2_start), "s")
        self.debug_stop_check("GA_BT2_merge")

    
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
                        #self.mp_util.run_subjob_with_hold(self.DIAMOND_mem_threshold, self.DIAMOND_job_limit, self.DIAMOND_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                        self.mp_util.run_subjob_with_mem_footprint(self.DMD_mem_footprint, self.DIAMOND_job_limit, self.GA_DIAMOND_label, job_name, self.commands, command_list)

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
                        self.mp_util.run_subjob_with_hold(self.DIAMOND_pp_mem_threshold, self.DIAMOND_pp_job_limit, self.DIAMOND_pp_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                                        
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
                command_list = self.commands.create_GA_final_merge_command(self.GA_final_merge_label, self.assemble_contigs_label, self.GA_BT2_label, self.GA_BLAT_label, self.GA_DIAMOND_label,  marker_file)
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
                marker_path = os.path.join(self.TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_centrifuge_command(self.ta_label, self.rRNA_filter_label, self.assemble_contigs_label, section, marker_file)
                    self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
                    
            sections = ["singletons"]
            if self.read_mode == "p":
                sections.extend(["paired"])
            if(self.contigs_present):
                sections.extend(["contigs"])    
            
            for section in sections:
                marker_file = "TA_kraken2_" + section
                marker_path = os.path.join(self.TA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_kraken2_command(self.ta_label, self.assemble_contigs_label, section, marker_file)
                    self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)        
            marker_file = "TA_taxon_pull"
            marker_path = os.path.join(self.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_taxon_pull_command(self.ta_label, self.GA_final_merge_label, marker_file)
                self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
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
                marker_path = os.path.join(self.TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_TA_centrifuge_command(self.ta_label, self.rRNA_filter_label, self.assemble_contigs_label, section, marker_file)
                    self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            
            marker_file = "TA_kraken2_pp"
            marker_path = os.path.join(self.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_kraken2_pp_command(self.ta_label, marker_file)
                self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.TA_path, "TA_stage_2.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            #------------------------------------------------------------------

            #-----------------------------------------------------------------
            # stage 3
            marker_path_list = []
            marker_file = "TA_centrifuge_pp"
            marker_path = os.path.join(self.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_centrifuge_pp_command(self.ta_label, marker_file)
                self.mp_util.run_subjob_with_hold(self.TA_mem_threshold, self.TA_job_limit, self.TA_job_delay, self.ta_label, marker_file, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.TA_path, "TA_stage_3.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            #-----------------------------------------------
            # stage 4
            marker_path_list = []
            
            marker_file = "TA_final"
            marker_path = os.path.join(self.TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = self.commands.create_TA_final_command(self.ta_label, self.assemble_contigs_label, marker_file)
                self.mp_util.run_subjob_simple(self.ta_label, marker_file, self.commands, command_list)
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
                self.mp_util.run_subjob_with_mp_store(self.ec_label, marker_file, self.commands, command_list)
            
            
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
                self.mp_util.run_subjob_with_mp_store(self.ec_label, job_name, self.commands, command_list)
            
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
                self.mp_util.run_subjob_simple(self.ec_label, marker_file, self.commands, command_list)
            
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
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_copy_taxa_label):
                job_name = self.output_copy_taxa_label
                command_list = self.commands.create_output_copy_taxa_command(self.output_label, self.ta_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
            if(self.contigs_present): 
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_contig_stats_label):
                    job_name = self.output_contig_stats_label
                    command_list = self.commands.create_output_contig_stats_command(self.output_label, self.assemble_contigs_label)
                    self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            
                
            if not(self.config_dict["no_host"]):
                print(dt.today(), "repopulating hosts for output")
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_singletons_label):
                    job_name = self.output_unique_hosts_singletons_label
                    command_list = self.commands.create_output_unique_hosts_singletons_command(self.output_label, self.quality_filter_label, self.host_filter_label)
                    self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
                if(self.read_mode == "p"):
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_pair_1_label):
                        job_name = self.output_unique_hosts_pair_1_label
                        command_list = self.commands.create_output_unique_hosts_pair_1_command(self.output_label, self.quality_filter_label, self.host_filter_label)
                        self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                        
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_hosts_pair_2_label):
                        job_name = self.output_unique_hosts_pair_2_label
                        command_list = self.commands.create_output_unique_hosts_pair_2_command(self.output_label, self.quality_filter_label, self.host_filter_label)
                        self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                        
                        
            #repop vectors
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_singletons_label):
                job_name = self.output_unique_vectors_singletons_label
                command_list = self.commands.create_output_unique_vectors_singletons_command(self.output_label, self.quality_filter_label, self.host_filter_label, self.vector_filter_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
            
            if(self.read_mode == "p"):
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_pair_1_label):
                    job_name = self.output_unique_vectors_pair_1_label
                    command_list = self.commands.create_output_unique_vectors_pair_1_command(self.output_label, self.quality_filter_label, self.host_filter_label, self.vector_filter_label)
                    self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                    
                if self.mp_util.check_bypass_log(self.output_folder_path, self.output_unique_vectors_pair_2_label):
                    job_name = self.output_unique_vectors_pair_2_label
                    command_list = self.commands.create_output_unique_vectors_pair_2_command(self.output_label, self.quality_filter_label, self.host_filter_label, self.vector_filter_label)
                    self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                    
            print(dt.today(), "output report phase 1 launched.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            
            self.mp_util.conditional_write_to_bypass_log(self.output_per_read_scores_label, "outputs/final_results", "input_per_seq_quality_report.csv")
            self.mp_util.conditional_write_to_bypass_log(self.output_copy_gene_map_label, "outputs/final_results", "final_gene_map.tsv")
            self.mp_util.conditional_write_to_bypass_log(self.output_copy_taxa_label, "outputs/final_results", "taxa_classifications.tsv")
            self.mp_util.conditional_write_to_bypass_log(self.output_contig_stats_label, "outputs/final_results", "contig_stats.txt")
            self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_singletons_label, "outputs/data/4_full_vectors", "singletons_full_vectors.fastq")
            if(self.read_mode == "p"):
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_pair_1_label, "outputs/data/4_full_vectors", "pair_1_full_vectors.fastq")
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_vectors_pair_2_label, "outputs/data/4_full_vectors", "pair_2_full_vectors.fastq")
                
            if not (self.config_dict["no_host"]):
                self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_singletons_label, "outputs/data/2_full_hosts", "singletons_full_hosts.fastq")
                if(self.read_mode == "p"):
                    self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_pair_1_label, "outputs/data/2_full_hosts", "pair_1_full_hosts.fastq")
                    self.mp_util.conditional_write_to_bypass_log(self.output_unique_hosts_pair_2_label, "outputs/data/2_full_hosts", "pair_2_full_hosts.fastq")
            #----------------------------------------------------------------------------
            #Phase 2
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_network_gen_label):
                command_list = self.commands.create_output_network_generation_command(self.output_label, self.GA_final_merge_label, self.ta_label, self.ec_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, self.output_network_gen_label, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_taxa_groupby_label):
                command_list = self.commands.create_output_taxa_groupby_command(self.output_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, self.output_taxa_groupby_label, self.commands, command_list)
        
            print(dt.today(), "output report phase 2 launched.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            self.mp_util.conditional_write_to_bypass_log(self.output_network_gen_label, "outputs/final_results", "RPKM_table.tsv")
            
            
            #-------------------------------------------------------------------
            #Phase 3
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_read_count_label):
                job_name = self.output_read_count_label
                command_list = self.commands.create_output_read_count_command(self.output_label, self.quality_filter_label, self.repop_job_label, self.GA_final_merge_label, self.ec_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                                

            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_per_read_scores_label):
                job_name = self.output_per_read_scores_label
                command_list = self.commands.create_output_per_read_scores_command(self.output_label, self.quality_filter_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)
                
            if self.mp_util.check_bypass_log(self.output_folder_path, self.output_ec_heatmap_label):
                job_name = self.output_ec_heatmap_label
                command_list = self.commands.create_output_EC_heatmap_command(self.output_label)
                self.mp_util.run_subjob_with_mp_store(self.output_label, job_name, self.commands, command_list)    
            
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
        if not self.config_dict["no_host"]:
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
        print("GA BT2:", '%1.1f' % (self.GA_BT2_end - self.GA_BT2_start - (self.cleanup_GA_BT2_end - self.cleanup_GA_BT2_start)), "s")
        print("GA BT2 cleanup:", '%1.1f' % (self.cleanup_GA_BT2_end - self.cleanup_GA_BT2_start), "s")
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
        


    
    
    
