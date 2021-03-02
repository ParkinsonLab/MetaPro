#!/usr/bin/env python
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu

def mem_checker(threshold):
    #threshold is a percentage
    mem = psu.virtual_memory()
    available_mem = mem.available
    total_mem = mem.total
    
    available_pct = 100 * available_mem / total_mem
    
    if(float(available_pct) <= float(threshold)):
        return False
    else:
        return True
        
        
def make_folder(folder_path):
    if not (os.path.exists(folder_path)):
        os.makedirs(folder_path)
        
def delete_folder_simple(folder_path):
    if(os.path.exists(folder_path)):
        print(dt.today(), "Deleting:", folder_path)
        shutil.rmtree(folder_path)
        print(dt.today(), "finished deleting:", folder_path)

def delete_folder(folder_path):
    if (os.path.exists(os.path.join(folder_path, "data"))):
        print("deleting", os.path.join(folder_path, "data"))
        shutil.rmtree(os.path.join(folder_path, "data"))
    else:
        print("can't delete folder: doesn't exist:", folder_path)
        
def compress_folder(folder_path):
    zip_loc = os.path.join(folder_path, "data")
    z = zipfile.ZipFile(folder_path + "_data.zip", "a", zipfile.ZIP_DEFLATED)
    print("compressing interim files:", folder_path)
    for root, dirs, files in os.walk(zip_loc):
        #print("root:", root)
        #print("dirs:", dirs)
        #print("files:", files)
        #print("===============================")
        for file in files:
            z.write(os.path.join(root, file))
    z.close()
        
def write_to_bypass_log(folder_path, message):
    bypass_log_path = os.path.join(folder_path, "bypass_log.txt")
    with open(bypass_log_path, "a") as bypass_log:
        bypass_log.write("\n")
        new_message = message + "\n"
        bypass_log.write(new_message)
        


def check_bypass_log(folder_path, message):
    stop_message = "stop_" + message
    bypass_keys_list = list()
    bypass_log_path = os.path.join(folder_path, "bypass_log.txt")
    if(os.path.exists(bypass_log_path)):
        with open(bypass_log_path, "r") as bypass_log:
            for line in bypass_log:
                bypass_key = line.strip("\n")
                bypass_keys_list.append(bypass_key)
        
        if(stop_message in bypass_keys_list):
            print(dt.today(), "stopping at:", message)
            print("to continue, remove:", stop_message, "from the bypass_log")
            sys.exit("brakes engaged")
        
        elif(message in bypass_keys_list):
            print(dt.today(), "bypassing:", message)
            return False
        
        else:
            print(dt.today(), "running:", message) 
            return True
    else:
        open(bypass_log_path, "a").close()
        print(dt.today(), "no bypass log.  running:", message)
        return True
        
def conditional_write_to_bypass_log(label, stage_folder, file_name, output_folder_path): 
    #convenience for checking if a file exists, and writing to the bypass log
    if check_bypass_log (output_folder, label):
        file_path = os.path.join(output_folder, stage_folder, file_name)
        if(os.path.exists(file_path)):
            write_to_bypass_log(output_folder_path, label)
    

# Used to determine quality encoding of fastq sequences.
# Assumes Phred+64 unless there is a character within the first 10000 reads with encoding in the Phred+33 range.
def check_code(segment):
    encoding = 64
    for item in segment:
        if(ord(item) < 64):
            encoding = 33
            break
    return encoding

def determine_encoding(fastq):
    #import the first 10k lines, then check the quality scores.
    #if the quality score symbols are below 76, it's phred33.  
    fastq_df = pd.read_csv(fastq, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3, nrows=40000)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    quality_encoding = fastq_df["quality"].apply(lambda x: check_code(x)).mean() #condense into a single number.
    if(quality_encoding == 64): #all must be 64 or else it's 33
        quality_encoding = 64
    else:
        quality_encoding =  33
    return quality_encoding


# handles where to kill the pipeline, due to the prev step behaving badly
# logic is:  if the files inside the dep_path (or dep job label shortcut to the final_results)
#            are empty, then there's an error.  kill the pipeline 
def check_where_kill(dep_job_label=None, dep_path=None):
    if dep_job_label is None:
        if dep_path is None:
            return True
        else:
            dep_job_path = dep_path
    else:
        dep_job_path = os.path.join(dep_job_label, "final_results")

    file_list = os.listdir(dep_job_path)
    if len(file_list) > 0:
        for item in file_list:
            file_check_path = os.path.join(dep_job_path, item)
            if (os.path.getsize(file_check_path)) == 0:
                print("empty file detected: rerunning stage")
                sys.exit("bad dep")
        # run the job, silently
        return True
    else:
        print("stopping the pipeline.  dependencies don't exist")
        sys.exit("no dep")


# handles where to auto-resume the pipeline on a subsequent run
# label: used as a shorthand for paths we expect
# full path: a bypass for when we want to use it for detecting a location that doesn't fall into the normal format (final_results)
# dep: for checking if the job's dependencies are satisfied-> meant to point to the last stage's "final_results"
# logic is: if the full_path has no files (or the job label shortcut to final_results)
#           and the dependencies are ok, start the stage
#Aug 19, 2019: There's a tweak to this:  DIAMOND will generate zero-size files, due to no-matches
#it's allowable.
def check_where_resume(job_label=None, full_path=None, dep_job_path=None, file_check_bypass = False):
    if(not file_check_bypass):
        check_where_kill(dep_job_path)
    if job_label:
        job_path = os.path.join(job_label, "final_results")
    else:
        job_path = full_path

    print("looking at:", job_path)

    if os.path.exists(job_path):
        file_list = os.listdir(job_path)
        if(not file_check_bypass):
            if len(file_list) > 0:
                for item in file_list:
                    file_check_path = os.path.join(job_path, item)
                    if (os.path.getsize(file_check_path)) == 0:
                        print("empty file detected: rerunning stage")
                        return False
                print("bypassing!")
                return True
            else:
                print("no files detected: running")
                return False
        else:
            print("bypassing for special reasons")
            return True
    else:
        print("doesn't exist: running")
        return False

def launch_and_create_simple(job_location, job_label, command_obj, commands):
    #just launches a job.  no multi-process.
    process = mp.Process(
        target=command_obj.create_and_launch,
        args=(job_location, job_label, commands)
    )
    process.start()
    process.join()

def launch_and_create_with_mp_store(mp_store, job_location, job_label, command_obj, commands):
    #just launches a job.  no multi-process.
    process = mp.Process(
        target=command_obj.create_and_launch,
        args=(job_location, job_label, commands)
    )
    process.start()
    mp_store.append(process)

def launch_only_simple(command_obj, commands):
    process = mp.Process(
        target=command_obj.launch_only,
        args=(commands, len(commands))
    )
    process.start()
    process.join()
    
def subdivide_and_launch(mp_store, job_location, job_label, command_obj, commands):
    #just launches a job.  no multi-process.
    job_counter = 0
    for item in commands:
        new_job_label = job_label + "_" + str(job_counter)
        job_counter += 1
        process = mp.Process(
            target=command_obj.create_and_launch,
            args=(job_location, new_job_label, [item])
        )
        process.start()
        mp_store.append(process)
    
    
def launch_only_with_hold(mp_store, mem_threshold, job_limit, job_delay, job_name, command_obj, command):
    #launch a job in launch-only mode
    job_submitted = False
    while(not job_submitted):
        if(len(mp_store) < job_limit):
            if(mem_checker(mem_threshold)):
                process = mp.Process(
                    target = command_obj.launch_only,
                    args = (command, len(command))
                )
                process.start()
                mp_store.append(process)
                print(dt.today(), job_name, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                job_submitted = True
            else:
                print(dt.today(), job_name, "Pausing. mem limit reached:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                time.sleep(job_delay)
        else:
            print(dt.today(), "job limit reached.  waiting for queue to flush")
            for item in mp_store:
                item.join()
            mp_store[:] = []
    time.sleep(job_delay)
    

def launch_and_create_with_hold(mp_store, mem_threshold, job_limit, job_delay, job_location, job_name, command_obj, command):
    #launch a job in launch-with-create mode
    job_submitted = False
    while(not job_submitted):
        if(len(mp_store) < job_limit):
            if(mem_checker(mem_threshold)):
                process = mp.Process(
                    target = command_obj.create_and_launch,
                    args = (job_location, job_name, command)
                )
                process.start()
                mp_store.append(process)
                print(dt.today(), job_name, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                job_submitted = True
            else:
                print(dt.today(), job_name, "Pausing. mem limit reached:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                time.sleep(job_delay)
        else:
            print(dt.today(), "job limit reached.  waiting for queue to flush")
            for item in mp_store:
                item.join()
            mp_store[:] = []
    time.sleep(job_delay)

#check if all jobs ran
def check_all_job_markers(job_marker_list, final_folder_checklist):
    time.sleep(2)
    #if it's already been created, that means the job was killed.
    if(os.path.exists(final_folder_checklist)):
        print(dt.today(), final_folder_checklist, "exists: adding to it")
        #open it, import it.
        with open(final_folder_checklist, "r") as old_list:
            for line in old_list:
                cleaned_line = line.strip("\n")
                job_marker_list.append(cleaned_line)
        #then overwrite it
        with open(final_folder_checklist, "w") as checklist:
            for item in job_marker_list:
                checklist.write(item +"\n")
                
            for item in job_marker_list:
                if(not os.path.exists(item)):
                    print(dt.today(), item, "not found.  kill the pipe.  restart this stage")
                    sys.exit("not all jobs completed")
        
    else:
        
        with open(final_folder_checklist, "w") as checklist:
            for item in job_marker_list:
                checklist.write(item +"\n")
                
            for item in job_marker_list:
                if(not os.path.exists(item)):
                    print(dt.today(), item, "not found.  kill the pipe.  restart this stage")
                    sys.exit("not all jobs completed")


def main(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode):
    paths = mpp.tool_path_obj(config_path)
    no_host = args_pack["no_host"]
    verbose_mode = args_pack["verbose_mode"]
    rRNA_chunks = int(paths.rRNA_chunksize)
    GA_chunksize = int(paths.GA_chunksize)
    
    BWA_mem_threshold = int(paths.BWA_mem_threshold)
    BLAT_mem_threshold = int(paths.BLAT_mem_threshold)
    DIAMOND_mem_threshold = int(paths.DIAMOND_mem_threshold)
    BWA_pp_mem_threshold = int(paths.BWA_pp_mem_threshold)
    BLAT_pp_mem_threshold = int(paths.BLAT_pp_mem_threshold)
    DIAMOND_pp_mem_threshold = int(paths.DIAMOND_pp_mem_threshold)
    Infernal_mem_threshold = int(paths.Infernal_mem_threshold)
    Barrnap_mem_threshold = int(paths.Barrnap_mem_threshold)
    DETECT_mem_threshold = int(paths.DETECT_mem_threshold)
    TA_mem_threshold = int(paths.TA_mem_threshold)
    filter_stringency = paths.filter_stringency
    
    
    BWA_job_limit = int(paths.BWA_job_limit)
    BLAT_job_limit = int(paths.BLAT_job_limit)
    DIAMOND_job_limit = int(paths.DIAMOND_job_limit)
    BWA_pp_job_limit = int(paths.BWA_pp_job_limit)
    BLAT_pp_job_limit = int(paths.BLAT_pp_job_limit)
    DIAMOND_pp_job_limit = int(paths.DIAMOND_pp_job_limit)
    Infernal_job_limit = int(paths.Infernal_job_limit)
    Barrnap_job_limit = int(paths.Barrnap_job_limit)
    DETECT_job_limit = int(paths.DETECT_job_limit)
    TA_job_limit = int(paths.TA_job_limit)
    
    Infernal_job_delay      = float(paths.Infernal_job_delay)
    Barrnap_job_delay       = float(paths.Barrnap_job_delay)
    BWA_job_delay           = float(paths.BWA_job_delay)
    BLAT_job_delay          = float(paths.BLAT_job_delay)
    DIAMOND_job_delay       = float(paths.DIAMOND_job_delay)
    BWA_pp_job_delay        = float(paths.BWA_pp_job_delay)
    BLAT_pp_job_delay       = float(paths.BLAT_pp_job_delay)
    DIAMOND_pp_job_delay    = float(paths.DIAMOND_pp_job_delay)
    DETECT_job_delay        = float(paths.DETECT_job_delay)
    TA_job_delay            = float(paths.TA_job_delay)
    
    #-----------------------------------------------------
    keep_all                = paths.keep_all
    keep_quality            = paths.keep_quality
    keep_vector             = paths.keep_vector
    keep_host               = paths.keep_host
    keep_rRNA               = paths.keep_rRNA
    keep_repop              = paths.keep_repop
    keep_assemble_contigs   = paths.keep_assemble_contigs
    keep_GA_BWA             = paths.keep_GA_BWA
    keep_GA_BLAT            = paths.keep_GA_BLAT
    keep_GA_DIAMOND         = paths.keep_GA_DIAMOND
    keep_GA_final           = paths.keep_GA_final
    keep_TA                 = paths.keep_TA
    keep_EC                 = paths.keep_EC
    keep_outputs            = paths.keep_outputs
    
    #------------------------------------------------------------------------
    
    BWA_cigar_cutoff        = paths.BWA_cigar_cutoff
    BLAT_identity_cutoff    = paths.BLAT_identity_cutoff
    BLAT_length_cutoff      = paths.BLAT_length_cutoff
    BLAT_score_cutoff       = paths.BLAT_score_cutoff
    DIAMOND_identity_cutoff = paths.DIAMOND_identity_cutoff
    DIAMOND_length_cutoff   = paths.DIAMOND_length_cutoff
    DIAMOND_score_cutoff    = paths.DIAMOND_score_cutoff
    
    print("============================================================")
    print("data cleaner options:")
    print("keep all:", keep_all)
    print("keep quality:", keep_quality)
    print("keep vector:", keep_vector)
    print("keep host:", keep_host)
    print("keep rRNA:", keep_rRNA)
    print("keep repop:", keep_repop)
    print("keep assemble contigs:", keep_assemble_contigs)
    print("keep GA BWA:", keep_GA_BWA)
    print("keep GA BLAT:", keep_GA_BLAT)
    print("keep GA DIAMOND:", keep_GA_DIAMOND)
    print("keep GA final:", keep_GA_final)
    print("keep TA:", keep_TA)
    print("keep EC:", keep_EC)
    print("keep outputs:", keep_outputs)
    print("===============================================")
    print("Job delay options:")
    print("Infernal job delay:", Infernal_job_delay)
    print("Barrnap job delay:", Barrnap_job_delay)
    print("BWA job delay:", BWA_job_delay)
    print("BLAT job delay:", BLAT_job_delay)
    print("DIAMOND job delay:", DIAMOND_job_delay)
    print("BWA pp job delay:", BWA_pp_job_delay)
    print("BLAT pp job delay:", BLAT_pp_job_delay)
    print("DIAMOND pp job delay:", DIAMOND_pp_job_delay)
    print("DETECT job delay:", DETECT_job_delay)
    print("=================================================")
    print("memory thresholds")
    print("Barrnap mem threshold:", Barrnap_mem_threshold)
    print("Infernal mem threshold:", Infernal_mem_threshold)
    print("BWA mem threshold:", BWA_mem_threshold)
    print("BLAT mem threshold:", BLAT_mem_threshold)
    print("DIAMOND mem threshold:", DIAMOND_mem_threshold)
    print("BWA pp mem threshold:", BWA_pp_mem_threshold)
    print("BLAT pp mem threshold:", BLAT_pp_mem_threshold)
    print("DIAMOND pp mem threshold:", DIAMOND_pp_mem_threshold)
    print("DETECT mem threshold:", DETECT_mem_threshold)
    print("======================================================")
    print("Barrnap job limit:", Barrnap_job_limit)
    print("Infernal job limit:", Infernal_job_limit)
    print("BWA job limit:", BWA_job_limit)
    print("BLAT job limit:", BLAT_job_limit)
    print("DIAMOND job limit:", DIAMOND_job_limit)
    print("BWA pp job limit:", BWA_pp_job_limit)
    print("BLAT pp job limit:", BLAT_pp_job_limit)
    print("DIAMOND pp job limit:", DIAMOND_pp_job_limit)
    print("DETECT job limit:", DETECT_job_limit)
    print("===================================================")
    print("Filter stringency:", filter_stringency)
    print("rRNA filter Chunk size:", rRNA_chunks)
    print("GA chunk size:", GA_chunksize)
    print("===================================================")
    print("BWA cigar cutoff:", BWA_cigar_cutoff)
    print("BLAT identity cutoff:", BLAT_identity_cutoff)
    print("BLAT length cutoff:", BLAT_length_cutoff)
    print("BLAT score cutoff:", BLAT_score_cutoff)
    print("DIAMOND identity cutoff:", DIAMOND_identity_cutoff)
    print("DIAMOND length cutoff:", DIAMOND_length_cutoff)
    print("DIAMOND score cutoff:", DIAMOND_score_cutoff)
    print("---------------------------------")
    if not single_path == "":
        read_mode = "single"
        quality_encoding = determine_encoding(single_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN SINGLE-ENDED MODE")
    else:
        read_mode = "paired"
        quality_encoding = determine_encoding(pair_1_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN PAIRED-MODE")
        
    if threads == 0:
        real_thread_count = mp.cpu_count()
    else:
        real_thread_count = threads
       
    if(real_thread_count == 1):
        real_thread_count = 2
    print("number of threads used:", real_thread_count)         
            
    mp_store = []  # stores the multiprocessing processes

    # --------------------------------------------------
    # profiling vars
    # profiling vars are init here, in case a stage is skipped
    start_time = time.time()
    
    quality_start           = quality_end           = cleanup_quality_start             = cleanup_quality_end           = 0
    host_start              = host_end              = cleanup_host_start                = cleanup_host_end              = 0
    vector_start            = vector_end            = cleanup_vector_start              = cleanup_vector_end            = 0
    rRNA_filter_start       = rRNA_filter_end       = cleanup_rRNA_filter_start         = cleanup_rRNA_filter_end       = 0
    repop_start             = repop_end             = cleanup_repop_start               = cleanup_repop_end             = 0
    assemble_contigs_start  = assemble_contigs_end  = cleanup_assemble_contigs_start    = cleanup_assemble_contigs_end  = 0
    destroy_contigs_start   = destroy_contigs_end   = cleanup_destroy_contigs_start     = cleanup_destroy_contigs_end   = 0
    GA_BWA_start            = GA_BWA_end            = cleanup_GA_BWA_start              = cleanup_GA_BWA_end            = 0
    GA_BLAT_start           = GA_BLAT_end           = cleanup_GA_BLAT_start             = cleanup_GA_BLAT_end           = 0
    GA_DIAMOND_start        = GA_DIAMOND_end        = cleanup_GA_DIAMOND_start          = cleanup_GA_DIAMOND_end        = 0
    TA_start                = TA_end                = cleanup_TA_start                  = cleanup_TA_end                = 0
    EC_start                = EC_end                                                                                    = 0
    EC_DETECT_start         = EC_DETECT_end                                                                             = 0
    EC_PRIAM_start          = EC_PRIAM_end                                                                              = 0
    EC_DIAMOND_start        = EC_DIAMOND_end                                                                            = 0
    cleanup_EC_start        = cleanup_EC_end                                                                            = 0
    Cytoscape_start         = Cytoscape_end         = cleanup_cytoscape_start           = cleanup_cytoscape_end         = 0
    
    # the pipeline stages are all labelled.  This is for multiple reasons:  to keep the interim files organized properly
    # and to perform the auto-resume/kill features

    quality_filter_label                    = "quality_filter"
    host_filter_label                       = "host_read_filter"
    vector_filter_label                     = "vector_read_filter"
    rRNA_filter_label                       = "rRNA_filter"
    rRNA_filter_split_label                 = "rRNA_filter_split"
    rRNA_filter_convert_label               = "rRNA_filter_convert"
    rRNA_filter_barrnap_label               = "rRNA_filter_barrnap"
    rRNA_filter_barrnap_merge_label         = "rRNA_filter_barrnap_merge"
    rRNA_filter_barrnap_pp_label            = "rRNA_filter_barrnap_pp"
    rRNA_filter_infernal_label              = "rRNA_filter_infernal"
    rRNA_filter_infernal_prep_label         = "rRNA_filter_infernal_prep"
    rRNA_filter_splitter_label              = "rRNA_filter_splitter"
    rRNA_filter_post_label                  = "rRNA_filter_post"
    repop_job_label                         = "duplicate_repopulation"
    assemble_contigs_label                  = "assemble_contigs"
    destroy_contigs_label                   = "destroy_contigs"
    GA_BWA_label                            = "GA_BWA"
    GA_BWA_pp_label                         = "GA_BWA_pp"
    GA_BLAT_label                           = "GA_BLAT"
    GA_BLAT_cleanup_label                   = "GA_BLAT_cleanup"
    GA_BLAT_cat_label                       = "GA_BLAT_cat"
    GA_BLAT_pp_label                        = "GA_BLAT_pp"
    GA_DIAMOND_label                        = "GA_DIAMOND"
    GA_DIAMOND_pp_label                     = "GA_DIAMOND_pp"
    GA_final_merge_label                    = "GA_FINAL_MERGE"
    taxon_annotation_label                  = "taxonomic_annotation"
    ec_annotation_label                     = "enzyme_annotation"
    ec_annotation_detect_label              = "enzyme_annotation_detect"
    ec_annotation_priam_label               = "enzyme_annotation_priam"
    ec_annotation_DIAMOND_label             = "enzyme_annotation_DIAMOND"
    ec_annotation_pp_label                  = "enzyme_annotation_pp"
    output_label                            = "outputs"
    output_copy_gene_map_label              = "output_copy_gene_map"
    output_clean_EC_label                   = "output_clean_ec"
    output_copy_taxa_label                  = "output_copy_taxa"
    output_network_gen_label                = "output_network_generation"
    output_unique_hosts_singletons_label    = "output_unique_hosts_singletons"
    output_unique_hosts_pair_1_label        = "output_unique_hosts_pair_1"
    output_unique_hosts_pair_2_label        = "output_unique_hosts_pair_2"
    output_unique_vectors_singletons_label  = "output_unique_vectors_singletons"
    output_unique_vectors_pair_1_label      = "output_unique_vectors_pair_1"
    output_unique_vectors_pair_2_label      = "output_unique_vectors_pair_2"
    output_combine_hosts_label              = "output_combine_hosts"
    output_per_read_scores_label            = "output_per_read_scores"
    output_contig_stats_label               = "output_contig_stats"
    output_ec_heatmap_label                 = "output_ec_heatmap"
    output_taxa_groupby_label               = "output_taxa_groupby"
    output_read_count_label                 = "output_read_count"
    
    # Creates our command object, for creating shellscripts.
    if read_mode == "single":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = None, sequence_path_1=None, sequence_path_2=None, sequence_single=single_path)
    elif read_mode == "paired":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = None, sequence_path_1=pair_1_path, sequence_path_2=pair_2_path, sequence_single=None)
    

    # This is the format we use to launch each stage of the pipeline.
    # We start a multiprocess that starts a subprocess.
    # The subprocess is created from the commands object

    # The quality filter stage
    #------------------------------------------------------------------------------
    quality_start = time.time()
    quality_path = os.path.join(output_folder_path, quality_filter_label)
    #if not check_where_resume(quality_path):
    if check_bypass_log(output_folder, quality_filter_label):
        command_list = commands.create_quality_control_command(quality_filter_label)
        job_name = quality_filter_label
        launch_and_create_simple(quality_filter_label, job_name, commands, command_list) 
        write_to_bypass_log(output_folder_path, quality_filter_label)
        cleanup_quality_start = time.time()
        if(keep_all == "no" and keep_quality == "no"):
            delete_folder(quality_path)
        elif(keep_all == "compress" or keep_quality == "compress"):
            compress_folder(quality_path)
            delete_folder(quality_path)
        cleanup_quality_end = time.time()
    quality_end = time.time()
    print("quality filter:", '%1.1f' % (quality_end - quality_start - (cleanup_quality_end - cleanup_quality_start)), "s")
    print("quality filter cleanup:", '%1.1f' %(cleanup_quality_end - cleanup_quality_start), "s")

    
    # The host read filter stage
    #-------------------------------------------------------------------------
    if not no_host:
        host_start = time.time()
        host_path = os.path.join(output_folder_path, host_filter_label)
        #if not check_where_resume(host_path, None, quality_path):
        if check_bypass_log(output_folder_path, host_filter_label):
            job_name = host_filter_label
            command_list = commands.create_host_filter_command(host_filter_label, quality_filter_label)
            launch_and_create_simple(host_filter_label, job_name, commands, command_list)
            
            write_to_bypass_log(output_folder_path, host_filter_label)
            cleanup_host_start = time.time()
            if(keep_all == "no" and keep_host == "no"):
                delete_folder(host_path)
            elif(keep_all == "compress" or keep_host == "compress"):
                compress_folder(host_path)
                delete_folder(host_path)
            cleanup_host_end = time.time()
                
        host_end = time.time()
        print("host filter:", '%1.1f' % (host_end - host_start - (cleanup_host_end - cleanup_host_start)), "s")
        print("host filter cleanup:", '%1.1f' %(cleanup_host_end - cleanup_host_start),"s")
        
    #-----------------------------------------------------------------
    # The vector contaminant filter stage
    vector_start = time.time()
    vector_path = os.path.join(output_folder_path, vector_filter_label)
    if no_host:
        #get dep args from quality filter
        #if not check_where_resume(vector_path, None, quality_path):
        if check_bypass_log(output_folder_path, vector_filter_label):        
            job_name = vector_filter_label
            command_list = commands.create_vector_filter_command(vector_filter_label, quality_filter_label)
            launch_and_create_simple(vector_filter_label, job_name, commands, command_list)
            write_to_bypass_log(output_folder_path, vector_filter_label)
            cleanup_vector_start = time.time()
            if(keep_all == "no" and keep_vector == "no"):
                delete_folder(vector_path)
            elif(keep_all == "compress" or  keep_vector == "compress"):
                compress_folder(vector_path)
                delete_folder(vector_path)
            cleanup_vector_end = time.time()
    else:
        #get the dep args from host filter
        #if not check_where_resume(vector_path, None, host_path):
        if check_bypass_log(output_folder_path, vector_filter_label):
            job_name = vector_filter_label
            command_list = commands.create_vector_filter_command(vector_filter_label, host_filter_label)
            launch_and_create_simple(vector_filter_label, job_name, commands, command_list)
            
            write_to_bypass_log(output_folder_path, vector_filter_label)
            cleanup_vector_start = time.time()
            if(keep_all == "no" and keep_vector == "no"):
                delete_folder(vector_path)
            elif(keep_all == "compress" or keep_vector == "compress"):
                compress_folder(vector_path)
                delete_folder(vector_path)
            cleanup_vector_end = time.time()
    vector_end = time.time()
    print("vector filter:", '%1.1f' % (vector_end - vector_start - (cleanup_vector_end - cleanup_vector_start)), "s")
    print("vector filter cleanup:", '%1.1f' % (cleanup_vector_end - cleanup_vector_start), "s")
    

    # ----------------------------------------------
    # rRNA removal stage
    rRNA_filter_start = time.time()
    rRNA_filter_path = os.path.join(output_folder_path, rRNA_filter_label)
    rRNA_filter_jobs_folder = os.path.join(rRNA_filter_path, "data", "jobs")
    #if not check_where_resume(rRNA_filter_path, None, vector_path):
    if check_bypass_log(output_folder_path, rRNA_filter_label): 
        marker_path_list = []
        sections = ["singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in reversed(sections):  #we go backwards due to a request by Ana.  pairs first, if applicable, then singletons
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            split_path = os.path.join(rRNA_filter_path, "data", section + "_fastq")
            barrnap_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section, section + "_barrnap")
            infernal_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section, section + "_infernal") 
            marker_file = "rRNA_filter_prep_" + section
            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
            #if not check_where_resume(job_label = None, full_path = second_split_path, dep_job_path = vector_path):
            if check_bypass_log(output_folder_path, rRNA_filter_split_label + "_" + section):
                print(dt.today(), "splitting:", section, " for rRNA filtration")
                job_name = "rRNA_filter_prep_" + section
                marker_path_list.append(marker_path)
                command_list = commands.create_rRNA_filter_prep_command_v3(rRNA_filter_label, section, vector_filter_label, marker_file)
                launch_and_create_with_mp_store(mp_store, rRNA_filter_label, job_name, commands, command_list)
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_prep.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        for element in sections:
            if check_bypass_log(output_folder_path, rRNA_filter_split_label + "_" + element):
                write_to_bypass_log(output_folder_path, rRNA_filter_split_label + "_" + element)
            
        #-------------------------------------------------------------------------------------------------
        # Convert fastq segments to fasta
        
        for section in reversed(sections):
            split_path = os.path.join(rRNA_filter_path, "data", section + "_fastq")
            if check_bypass_log(output_folder_path, rRNA_filter_convert_label + "_" + section):
                marker_path_list = []
                for item in os.listdir(split_path):
                    root_name = item.split(".")[0]
                    fasta_path = os.path.join(rRNA_filter_path, "data", section + "_fasta")
                    fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                    marker_file = root_name + "_convert_fasta"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    
                    fasta_out_size = os.stat(fasta_file).st_size if (os.path.exists(fasta_file)) else 0
                    if(fasta_out_size > 0) or (os.path.exists(marker_path)):
                        print(dt.today(), item, "already converted to fasta.  skipping")
                        continue
                    else:
                        job_name = root_name + "_convert_to_fasta"
                        marker_path_list.append(marker_path)
                        command_list = commands.create_rRNA_filter_convert_fastq_command("rRNA_filter", section, root_name+".fastq", marker_file)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                        
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_convert_" + section + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                write_to_bypass_log(output_folder_path, rRNA_filter_convert_label + "_" + section)
            
                    
        #-------------------------------------------------------------------------------------------------
        # BARRNAP
        for section in reversed(sections):  
            #convert data to fasta, then run barrnap separately, then cat the barrnap, then run barrnap PP
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            split_path      = os.path.join(rRNA_filter_path, "data", section + "_fastq")
            fasta_path      = os.path.join(rRNA_filter_path, "data", section + "_fasta")
            barrnap_path    = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap")
            infernal_path   = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal") 
            
            mRNA_path       = os.path.join(rRNA_filter_path, "data", section + "_mRNA")
            
            #if not check_where_resume(job_label = None, full_path = barrnap_path, dep_job_path = vector_path):
            if check_bypass_log(output_folder_path, rRNA_filter_barrnap_label + "_" + section):
                concurrent_job_count = 0
                batch_count = 0
                marker_path_list = []
                
                for item in os.listdir(fasta_path):
                    root_name = item.split(".")[0]
                    final_marker_file = root_name + "_barrnap_concat"
                    final_marker_path = os.path.join(rRNA_filter_jobs_folder, final_marker_file)
                    barrnap_arc_out_file = os.path.join(barrnap_path, root_name + "_arc.barrnap_out")
                    barrnap_bac_out_file = os.path.join(barrnap_path, root_name + "_bac.barrnap_out")
                    barrnap_euk_out_file = os.path.join(barrnap_path, root_name + "_euk.barrnap_out")
                    barrnap_mit_out_file = os.path.join(barrnap_path, root_name + "_mit.barrnap_out")
                    final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                    fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                    fastq_file = os.path.join(split_path, root_name + ".fastq")
                    barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                    marker_file_arc = root_name + "_barrnap_arc"
                    marker_file_bac = root_name + "_barrnap_bac"
                    marker_file_euk = root_name + "_barrnap_euk"
                    marker_file_mit = root_name + "_barrnap_mit"
                    
                    marker_path_arc = os.path.join(rRNA_filter_jobs_folder, marker_file_arc)
                    marker_path_bac = os.path.join(rRNA_filter_jobs_folder, marker_file_bac)
                    marker_path_euk = os.path.join(rRNA_filter_jobs_folder, marker_file_euk)
                    marker_path_mit = os.path.join(rRNA_filter_jobs_folder, marker_file_mit)
                    
                    barrnap_arc_out_size    = os.stat(barrnap_arc_out_file).st_size if (os.path.exists(barrnap_arc_out_file)) else 0
                    barrnap_bac_out_size    = os.stat(barrnap_bac_out_file).st_size if (os.path.exists(barrnap_bac_out_file)) else 0
                    barrnap_euk_out_size    = os.stat(barrnap_euk_out_file).st_size if (os.path.exists(barrnap_euk_out_file)) else 0
                    barrnap_mit_out_size    = os.stat(barrnap_mit_out_file).st_size if (os.path.exists(barrnap_mit_out_file)) else 0
                    
                    
                    if(os.path.exists(final_marker_path)):
                        print(dt.today(), "skipping barrnap.  data already merged", final_marker_path)
                        continue
                    else:
                        if((barrnap_arc_out_size > 0) and (os.path.exists(marker_path_arc))):
                            print(dt.today(), "barrnap arc already run.  skipping:", item) 
                            continue
                        else:
                            job_name = root_name + "_barrnap_arc"
                            marker_path_list.append(marker_path_arc)
                            command_list = commands.create_rRNA_filter_barrnap_arc_command("rRNA_filter", section, root_name, marker_file_arc)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                            
                            
                        if((barrnap_bac_out_size > 0) and (os.path.exists(marker_path_bac))):
                            print(dt.today(), "barrnap bac already run.  skipping:", item) 
                            continue
                        else:
                            job_name = root_name + "_barrnap_bac"
                            marker_path_list.append(marker_path_bac)
                            command_list = commands.create_rRNA_filter_barrnap_bac_command("rRNA_filter", section, root_name, marker_file_bac)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                            
                        if((barrnap_euk_out_size > 0) and (os.path.join(marker_path_euk))):
                            print(dt.today(), "barrnap euk already run.  skipping:", item) 
                            continue
                        else:
                            job_name = root_name + "_barrnap_euk"
                            marker_path_list.append(marker_path_euk)
                            command_list = commands.create_rRNA_filter_barrnap_euk_command("rRNA_filter", section, root_name, marker_file_euk)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                            
                        if((barrnap_mit_out_size > 0) and (os.path.join(marker_path_mit))):
                            print(dt.today(), "barrnap mit already run.  skipping:", item) 
                            continue
                        else:
                            job_name = root_name + "_barrnap_mit"
                            marker_path_list.append(marker_path_mit)
                            command_list = commands.create_rRNA_filter_barrnap_mit_command("rRNA_filter", section, root_name, marker_file_mit)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                print(dt.today(), "waiting for Barrnap jobs to finish")
                for item in mp_store:
                    item.join()
                mp_store[:] = []
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_" + section + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                
                #------------------------------------------------------
                #merge the barrnap data
                if check_bypass_log(output_folder_path, rRNA_filter_barrnap_merge_label + "_" + section):
                    marker_path_list = []
                    for item in os.listdir(fasta_path):
                        root_name = item.split(".")[0]
                        marker_file = root_name + "_barrnap_cat"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                        final_barrnap_out_size  = os.stat(final_barrnap_out).st_size if (os.path.exists(final_barrnap_out)) else 0
                        
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "barrnap already merged. skipping:", item)
                            continue
                        else:
                            job_name = root_name + "_barrnap_cat"
                            marker_path_list.append(marker_path)
                            command_list = commands.create_rRNA_filter_barrnap_cat_command("rRNA_filter", section, root_name, marker_file)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                print(dt.today(), "waiting for Barrnap pp to finish")
                for item in mp_store:
                    item.join()
                mp_store[:] = []
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_cat_" + section  + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                write_to_bypass_log(output_folder_path, rRNA_filter_barrnap_merge_label + "_" + section)

                #-----------------------------------------------------
                #run the barrnap PP
                if check_bypass_log(output_folder_path, rRNA_filter_barrnap_pp_label + "_" + section):
                    marker_path_list = []
                    for item in os.listdir(fasta_path):
                        root_name = item.split(".")[0]
                        barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                        barrnap_mRNA_out_size   = os.stat(barrnap_mrna_file).st_size if (os.path.exists(barrnap_mrna_file)) else 0
                        marker_file = root_name + "_barrnap_pp"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "barrnap pp already run.  skipping:", item)
                            continue
                        else:
                            job_name = root_name + "_barrnap_pp"
                            marker_path_list.append(marker_path)
                            command_list = commands.create_rRNA_filter_barrnap_pp_command("rRNA_filter", section, root_name + ".fastq", marker_file)
                            launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                    
                    print(dt.today(), "waiting for Barrnap pp to finish")
                    for item in mp_store:
                        item.join()
                    mp_store[:] = []
                    final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_" + section +  ".txt")
                    check_all_job_markers(marker_path_list, final_checklist)
                    write_to_bypass_log(output_folder_path, rRNA_filter_barrnap_label + "_" + section)
            
        #----------------------------------------------------------------------------
        # INFERNAL
        for section in reversed(sections):  
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            barrnap_mRNA_fastq_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap_mRNA")
            infernal_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal") 
            barrnap_mRNA_fasta_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap_mRNA_fasta")
            splitter_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal_mRNA")
        
            if check_bypass_log(output_folder_path, rRNA_filter_infernal_prep_label + "_" + section):
                concurrent_job_count = 0
                batch_count = 0
                #these jobs now have to be launched in segments
                for item in os.listdir(barrnap_mRNA_fastq_path):
                   
                    if(item.endswith("_barrnap_mRNA.fastq")):
                        root_name = item.split(".")[0]
                        marker_file = root_name + "_infernal_prep"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        infernal_prep_out_file = os.path.join(barrnap_mRNA_fasta_path, root_name + ".fasta")
                        infernal_prep_file_size = os.stat(infernal_prep_out_file).st_size if (os.path.exists(infernal_prep_out_file)) else 0
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "Infernal prep already ran on this sample.  skipping", item)
                            continue
                        
                        else:
                            marker_path_list.append(marker_path)
                            job_name = "rRNA_filter_infernal_prep_" + root_name
                            command_list = commands.create_rRNA_filter_infernal_prep_command("rRNA_filter", section, item, root_name, marker_file)
                            launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
                            
                print(dt.today(), "final batch: infernal prep")
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []  # clear the list
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_prep_" + section + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                write_to_bypass_log(output_folder_path, rRNA_filter_infernal_prep_label + "_" + section)
            

            if check_bypass_log(output_folder_path, rRNA_filter_infernal_label + "_" + section):
                marker_path_list = []
                for item in os.listdir(barrnap_mRNA_fasta_path):
                    #using a job marker is ineffective.  The marker will still write 
                    root_name = item.split("_barrnap_mRNA")[0]
                    marker_file = root_name + "_infernal"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "infernal already run. skipping:", root_name + "_infernal")
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        inf_command = commands.create_rRNA_filter_infernal_command("rRNA_filter", section, root_name, marker_file)
                        job_name = "rRNA_filter_infernal_" + root_name
                        #launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, inf_command)
                        launch_and_create_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, rRNA_filter_label, job_name, commands, inf_command)
                        
                        
                print(dt.today(), "final batch: infernal")
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_" + section + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                write_to_bypass_log(output_folder_path, rRNA_filter_infernal_label + "_" + section)
            
            if (section != "pair_2"):
                if check_bypass_log(output_folder_path, rRNA_filter_splitter_label + "_" + section):
                    marker_path_list = []
                    for item in os.listdir(barrnap_mRNA_fasta_path):
                        root_name = item.split("_barrnap_mRNA")[0]
                        splitter_out_file = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal_mRNA", root_name + "_mRNA.fastq")
                        splitter_out_file_size = os.stat(splitter_out_file).st_size if os.path.exists(splitter_out_file) else 0
                        marker_file = root_name + "_infernal_pp"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "infernal mRNA splitter already run. skipping:", marker_file)
                            print("file size:", splitter_out_file_size, "file:", splitter_out_file)
                            continue
                        else:
                            job_name = "rRNA_filter_infernal_splitter_" + root_name
                            marker_path_list.append(marker_path)
                            command_list = commands.create_rRNA_filter_splitter_command("rRNA_filter", section, root_name, marker_file)
                            print(command_list)
                            launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
                            
                    print(dt.today(), "final batch: infernal splitter")
                    for p_item in mp_store:
                        p_item.join()
                    mp_store[:] = []
                    final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_splitter_" + section + ".txt")
                    check_all_job_markers(marker_path_list, final_checklist)
                    write_to_bypass_log(output_folder_path, rRNA_filter_splitter_label + "_" + section)
            else:
                print(dt.today(), "not calling Infernal rRNA splitter on pair 2.  data handled by pair 1 as a combination")
        
                    
                    
        marker_path_list = []
        for section in reversed(sections):
            if check_bypass_log(output_folder_path, rRNA_filter_post_label + "_" + section):
                print(dt.today(), "now running rRNA filter post:", section)
                marker_file = section + "_rRNA_packup"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                job_name = "rRNA_post_cat"
                marker_path_list.append(marker_path)
                command_list = commands.create_rRNA_filter_final_cat_command("rRNA_filter", section, marker_file)
                print("command list:", command_list)
                launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
                
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_final_cat.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        for section in reversed(sections):
            write_to_bypass_log(output_folder_path, rRNA_filter_splitter_label + "_" + section)
        
        write_to_bypass_log(output_folder_path, rRNA_filter_label)
        cleanup_rRNA_filter_start = time.time()
        delete_folder_simple(rRNA_filter_jobs_folder)
        if(keep_all == "no" and keep_rRNA == "no"):
            delete_folder(rRNA_filter_path)
        elif(keep_all == "compress" or keep_rRNA == "compress"):
            compress_folder(rRNA_filter_path)
            delete_folder(rRNA_filter_path)
        cleanup_rRNA_filter_end = time.time()
    rRNA_filter_end = time.time()
    
    print("rRNA filter:", '%1.1f' % (rRNA_filter_end - rRNA_filter_start - (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start)), "s")
    print("rRNA filter cleanup:", '%1.1f' % (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start), "s")
    

    #---------------------------------------------------------------------------------------------------------------------
    # Duplicate repopulation
    repop_start = time.time()
    repop_path = os.path.join(output_folder_path, repop_job_label)
    #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
    if check_bypass_log(output_folder_path, repop_job_label):
        job_name = repop_job_label
        command_list = commands.create_repop_command_v2_step_1(repop_job_label, quality_filter_label, rRNA_filter_label)
        subdivide_and_launch(mp_store, repop_job_label, job_name, commands, command_list)
    
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        
        job_name = repop_job_label
        command_list = commands.create_repop_command_v2_step_2(repop_job_label, quality_filter_label, rRNA_filter_label)
        subdivide_and_launch(mp_store, repop_job_label, job_name, commands, command_list)
    
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        
    
        write_to_bypass_log(output_folder_path, repop_job_label)
        
        cleanup_repop_start = time.time()
        if(keep_all == "no" and keep_repop == "no"):
            delete_folder(repop_path)
        elif(keep_all == "compress" or keep_repop == "compress"):
            compress_folder(repop_job_path)
            delete_folder(repop_path)
        cleanup_repop_end = time.time()
        
    repop_end = time.time()
    print("repop:", '%1.1f' % (repop_end - repop_start - (cleanup_repop_end - cleanup_repop_start)), "s")
    print("repop cleanup:", '%1.1f' % (cleanup_repop_end - cleanup_repop_start), "s")

    # -------------------------------------------------------------
    # Assemble contigs
    assemble_contigs_start = time.time()
    assemble_contigs_path = os.path.join(output_folder_path, assemble_contigs_label)
    
    
    #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
    
    if check_bypass_log(output_folder_path, assemble_contigs_label):
        job_name = assemble_contigs_label
        command_list = commands.create_assemble_contigs_command(assemble_contigs_label, repop_job_label)
        launch_and_create_simple(assemble_contigs_label, job_name, commands, command_list)
        mgm_file = os.path.join(assemble_contigs_path, "data", "1_mgm", "gene_report.txt")
        if(os.path.exists(mgm_file)):
            write_to_bypass_log(output_folder_path, assemble_contigs_label)
        else:
            sys.exit("mgm did not run.  look into it.  pipeline stopping here")
        
        cleanup_assemble_contigs_start = time.time()
        
        if(keep_all == "no" and keep_assemble_contigs == "no"):
            delete_folder(assemble_contigs_path)
        elif(keep_all == "compress" or keep_assemble_contigs == "compress"):
            compress_folder(assemble_contigs_path)
            delete_folder(assemble_contigs_path)
        cleanup_assemble_contigs_end = time.time()
        
        
    assemble_contigs_end = time.time()
    print("assemble contigs:", '%1.1f' % (assemble_contigs_end - assemble_contigs_start - (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start)), "s")    
    print("assemble contigs cleanup:", '%1.1f' % (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start), "s")


    # ----------------------------------------------
    # BWA gene annotation
    
    
    GA_BWA_start = time.time()
    GA_BWA_path = os.path.join(output_folder_path, GA_BWA_label)
    GA_BWA_jobs_folder = os.path.join(GA_BWA_path, "data", "jobs")
    #if not check_where_resume(GA_BWA_path, None, assemble_contigs_path):
    if check_bypass_log(output_folder_path, GA_BWA_label):
        marker_path_list = []
        marker_file = "GA_split_fasta_contigs"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping", marker_file)
        else:
            job_name = "GA_prep_split_contigs"
            marker_path_list.append(marker_path)
            command_list = commands.create_split_ga_fasta_data_command(GA_BWA_label, assemble_contigs_label, "contigs", marker_file)
            launch_and_create_with_mp_store(mp_store, GA_BWA_label, job_name, commands, command_list)
        
        
        sections = ["singletons"]
        if(read_mode == "paired"):
            sections.extend(["pair_1", "pair_2"])
        for section in sections: 
            marker_file = "GA_split_fastq_" + section
            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping", marker_file)
            else:
                marker_path_list.append(marker_path)
                job_name = "GA_prep_split_" + section
                command_list = commands.create_split_ga_fastq_data_command(GA_BWA_label, assemble_contigs_label, section, marker_file)
                launch_and_create_with_mp_store(mp_store, GA_BWA_label, job_name, commands, command_list)
        
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_prep.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        
        #-------------------------------------------------------------------------
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                job_submitted = False
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA" + "_" + file_tag
                marker_file = file_tag + "_bwa"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                #this checker assumes that BWA only exports a file when it's finished running
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BWA_annotate_command_v2(GA_BWA_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, BWA_mem_threshold, BWA_job_limit, BWA_job_delay, GA_BWA_label, job_name, commands, command_list)

        print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_BWA_label)
            
    if check_bypass_log(output_folder_path, GA_BWA_pp_label):
        marker_path_list = []
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA_pp" + "_" + file_tag
                marker_file = file_tag + "_bwa_pp"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BWA_pp_command_v2(GA_BWA_label, assemble_contigs_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, BWA_pp_mem_threshold, BWA_pp_job_limit, BWA_pp_job_delay, GA_BWA_label, job_name, commands, command_list)
                        
        print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
        for item in mp_store:
            item.join()
        mp_store[:] = []
        marker_file = "BWA_copy_contig_map"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:   
            marker_path_list.append(marker_path)
            command_list = commands.create_BWA_copy_contig_map_command(GA_BWA_label, assemble_contigs_label, marker_file)
            launch_and_create_simple(GA_BWA_label, GA_BWA_label + "_copy_contig_map", commands, command_list)
        
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_pp.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_BWA_pp_label)
        
        cleanup_GA_BWA_start = time.time()
        delete_folder_simple(GA_BWA_jobs_folder)
        if(keep_all == "no" and keep_GA_BWA == "no"):
            delete_folder(GA_BWA_path)
        elif(keep_all == "compress" or keep_GA_BWA == "compress"):
            compress_folder(GA_BWA_path)
            delete_folder(GA_BWA_path)
        cleanup_GA_BWA_end = time.time()
    GA_BWA_end = time.time()
    print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
    print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
    
    # ------------------------------------------------
    # BLAT gene annotation
    GA_BLAT_start = time.time()
    GA_BLAT_path = os.path.join(output_folder_path, GA_BLAT_label)
    GA_BLAT_jobs_folder = os.path.join(GA_BLAT_path, "data", "jobs")

    if check_bypass_log(output_folder_path, GA_BLAT_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                for fasta_db in os.listdir(paths.DNA_DB_Split):
                    if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                        job_name = "BLAT_" + file_tag + "_" + fasta_db
                        marker_file = file_tag + "_blat_" + fasta_db
                        marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                        #This checker assume BLAT only exports a file when it's finished running
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "BLAT job ran already, skipping:", marker_file)
                            continue
                        else:
                            marker_path_list.append(marker_path)
                            command_list = commands.create_BLAT_annotate_command_v2(GA_BLAT_label, full_sample_path, fasta_db, marker_file)
                            launch_only_with_hold(mp_store, BLAT_mem_threshold, BLAT_job_limit, BLAT_job_delay, job_name, commands, command_list)
                            
                                
        print(dt.today(), "final BLAT job removal")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_BLAT_label)
        

        
    if check_bypass_log(output_folder_path, GA_BLAT_cat_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                job_name = file_tag + "_cat"
                
                marker_file = file_tag + "_blat_cat"
                marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BLAT_cat_command_v2(GA_BLAT_label, full_sample_path, marker_file)
                    launch_only_with_hold(mp_store, BLAT_mem_threshold, BLAT_job_limit, BLAT_job_delay, job_name, commands, command_list)
                
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_cat.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_BLAT_cat_label)
        
    
    
    if check_bypass_log(output_folder_path, GA_BLAT_pp_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BLAT_" + file_tag + "_pp"
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                marker_file = file_tag + "_blat_pp"
                marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BLAT_pp_command_v2(GA_BLAT_label, full_sample_path, GA_BWA_label, marker_file)
                    launch_and_create_with_hold(mp_store, BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, GA_BLAT_label, job_name, commands, command_list)
                
        print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        
        job_name = "GA_BLAT_copy_contigs"
        marker_file = "blat_copy_contig_map"
        marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_BLAT_copy_contig_map_command(GA_BLAT_label, GA_BWA_label, marker_file)
            launch_and_create_simple(GA_BLAT_label, job_name, commands, command_list)
        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_pp.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_BLAT_pp_label)
        
        
    cleanup_GA_BLAT_start = time.time()
    delete_folder_simple(GA_BLAT_jobs_folder)
    if(keep_all == "no" and keep_GA_BLAT == "no"):
        delete_folder(GA_BLAT_path)
    elif(keep_all == "compress" or keep_GA_BLAT == "compress"):
        compress_folder(GA_BLAT_path)
        delete_folder(GA_BLAT_path)
    cleanup_GA_BLAT_end = time.time()
    GA_BLAT_end = time.time()
    print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
    print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
    
    # ------------------------------------------------------
    # Diamond gene annotation
    GA_DIAMOND_start = time.time()
    GA_DIAMOND_path = os.path.join(output_folder_path, GA_DIAMOND_label)
    GA_DIAMOND_tool_output_path = os.path.join(GA_DIAMOND_path, "data", "0_diamond")
    GA_DIAMOND_jobs_folder = os.path.join(GA_DIAMOND_path, "data", "jobs")
    #if not check_where_resume(None, GA_DIAMOND_tool_output_path, GA_BLAT_path, file_check_bypass = True):
    if check_bypass_log(output_folder_path, GA_DIAMOND_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_path)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_annotate_command_v2(GA_DIAMOND_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, DIAMOND_mem_threshold, DIAMOND_job_limit, DIAMOND_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                
        print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        write_to_bypass_log(output_folder_path, GA_DIAMOND_label)
        
        
    #if not check_where_resume(GA_DIAMOND_path, None, GA_DIAMOND_tool_output_path, file_check_bypass = True):
    if check_bypass_log(output_folder_path, GA_DIAMOND_pp_label):
        print(dt.today(), "DIAMOND PP threads used:", real_thread_count/2)
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_pp_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond_pp"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_pp_command_v2(GA_DIAMOND_label, GA_BLAT_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, DIAMOND_pp_mem_threshold, DIAMOND_pp_job_limit, DIAMOND_pp_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                                    
        print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND_pp.txt")
        check_all_job_markers(marker_path_list, final_checklist)
            
        write_to_bypass_log(output_folder_path, GA_DIAMOND_pp_label)
    
        
    
        cleanup_GA_DIAMOND_start = time.time()
        delete_folder_simple(GA_DIAMOND_jobs_folder)
        if(keep_all == "no" and keep_GA_DIAMOND == "no"):
            delete_folder(GA_DIAMOND_path)
        elif(keep_all == "compress" or keep_GA_DIAMOND == "compress"):
            compress_folder(GA_DIAMOND_path)
            delete_folder(GA_DIAMOND_path)
        cleanup_GA_DIAMOND_end = time.time()
    GA_DIAMOND_end = time.time()
    print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
    print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
    
    
    GA_final_merge_start = time.time()
    GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
    if check_bypass_log(output_folder_path, GA_final_merge_label):
        marker_file = "GA_final_merge"
        marker_path = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge")
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping: GA final merge")
        else:
            command_list = commands.create_GA_final_merge_command(GA_final_merge_label, GA_BWA_label, GA_BLAT_label, GA_DIAMOND_label, assemble_contigs_label, marker_file)
            job_name = "GA_final_merge"
            launch_and_create_simple(GA_final_merge_label, job_name, commands, command_list)
        
        #check if all_proteins.faa was generated
        all_proteins_path = os.path.join(output_folder_path, GA_final_merge_label, "final_results", "all_proteins.faa")
        if(os.path.getsize(all_proteins_path) > 0):
            write_to_bypass_log(output_folder_path, GA_final_merge_label)
            print(dt.today(), "All_proteins.faa is OK.  Continuing")
        else:
            sys.exit("GA final merge failed.  proteins weren't translated")
            
    GA_final_merge_end = time.time()
    print("GA final merge:", '%1.1f' % (GA_final_merge_end - GA_final_merge_start), "s")
    if(keep_all == "no" and keep_GA_final == "no"):
        delete_folder(GA_FINAL_MERGE_path)
    elif(keep_all == "compress" or keep_GA_final == "compress"):
        compress_folder(GA_FINAL_MERGE_path)
        delete_folder(GA_FINAL_MERGE_path)
    
    # ------------------------------------------------------
    # Taxonomic annotation
    TA_start = time.time()
    TA_path = os.path.join(output_folder_path, taxon_annotation_label)
    TA_jobs_folder = os.path.join(TA_path, "data", "jobs")
    if check_bypass_log(output_folder_path, taxon_annotation_label):
        #-----------------------------------------
        # stage 1
        marker_path_list = []
        #----------------------------------------------
        #centrifuge is too much of a RAM hog.  can't run more than 1 at a time
        sections = ["reads"]
        for section in sections:
            marker_file = "TA_centrifuge_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
                
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["paired"])
            
        for section in sections:
            marker_file = "TA_kaiju_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_kaiju_command(taxon_annotation_label, assemble_contigs_label, section, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)        
        marker_file = "TA_taxon_pull"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_taxon_pull_command(taxon_annotation_label, GA_final_merge_label, marker_file)
            launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        print(dt.today(), "waiting for TA stage 1")
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        final_checklist = os.path.join(TA_path, "TA_stage_1.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        
        #--------------------------------------------------
        # stage 2
        marker_path_list = []
        sections = ["contigs"]
        for section in sections:
            marker_file = "TA_centrifuge_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        
        marker_file = "TA_kaiju_pp"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_file)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_kaiju_pp_command(taxon_annotation_label, marker_file)
            launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        final_checklist = os.path.join(TA_path, "TA_stage_2.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        #------------------------------------------------------------------

        #-----------------------------------------------------------------
        # stage 3
        marker_path_list = []
        marker_file = "TA_centrifuge_pp"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_centrifuge_pp_command(taxon_annotation_label, marker_file)
            launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        final_checklist = os.path.join(TA_path, "TA_stage_3.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        #-----------------------------------------------
        # stage 4
        marker_path_list = []
        
        marker_file = "TA_final"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_final_command(taxon_annotation_label, assemble_contigs_label, marker_file)
            launch_and_create_simple(taxon_annotation_label, marker_file, commands, command_list)
        final_checklist = os.path.join(TA_path, "TA_final.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        
        if(os.path.exists(marker_path)):
            write_to_bypass_log(output_folder_path, taxon_annotation_label)
            
    cleanup_TA_start = time.time()
    if(keep_all == "no" and keep_TA == "no"):
        delete_folder(TA_path)
    elif(keep_all == "compress" or keep_TA == "compress"):
        compress_folder(TA_path)
        delete_folder(TA_path)
    cleanup_TA_end = time.time()
    TA_end = time.time()
    print("TA:", '%1.1f' % (TA_end - TA_start - (cleanup_TA_end - cleanup_TA_start)), "s")
    print("TA cleanup:", '%1.1f' % (cleanup_TA_end - cleanup_TA_start), "s")
    
    
    # ------------------------------------------------------
    # Detect EC annotation
    ec_annotation_path = os.path.join(output_folder_path, ec_annotation_label)
    EC_start = time.time()
    #There's a 2-step check.  We don't want it ti re-run either DETECT, or PRIAM+DIAMOND because they're too slow
    #if not check_where_resume(ec_annotation_path, None, GA_DIAMOND_path):
    #if check_bypass_log(output_folder_path, ec_annotation_label):
    EC_DETECT_start = time.time()
    ec_detect_path = os.path.join(ec_annotation_path, "data", "0_detect")
    #if not check_where_resume(job_label = None, full_path = ec_detect_path, dep_job_path = GA_DIAMOND_path):
    if check_bypass_log(output_folder_path, ec_annotation_detect_label):
        marker_file = "ec_detect"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            command_list = commands.create_EC_DETECT_command(ec_annotation_label, GA_final_merge_label, marker_file)
            launch_and_create_with_mp_store(mp_store, ec_annotation_label, marker_file, commands, command_list)
        
        
    EC_DETECT_end = time.time()
    print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    
    # --------------------------------------------------------------
    # Priam EC annotation.  Why isn't it parallel? computing restraints.  Not enough mem
    EC_PRIAM_start = time.time()
    
    ec_priam_path = os.path.join(ec_annotation_path, "data", "1_priam")
    #if not check_where_resume(job_label = None, full_path = ec_priam_path, dep_job_path = GA_DIAMOND_path):
    if check_bypass_log(output_folder_path, ec_annotation_priam_label):
        marker_file = "ec_priam"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            if(os.path.exists(ec_priam_path)):
                print(dt.today(), "starting with a fresh PRIAM run")
                shutil.rmtree(ec_priam_path)
                
            command_list = commands.create_EC_PRIAM_command(ec_annotation_label, GA_final_merge_label, marker_file)
            launch_and_create_with_mp_store(mp_store, ec_annotation_label, marker_file, commands, command_list)
        
      
        #process.join()
    EC_PRIAM_end = time.time()
    print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
    # --------------------------------------------------------------
    # DIAMOND EC annotation 
    EC_DIAMOND_start = time.time()
    ec_diamond_path = os.path.join(ec_annotation_path, "data", "2_diamond")
    if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        marker_file = "ec_diamond"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            job_name = "ec_diamond"
            command_list = commands.create_EC_DIAMOND_command(ec_annotation_label, GA_final_merge_label, marker_file)
            launch_and_create_with_mp_store(mp_store, ec_annotation_label, job_name, commands, command_list)
        
    EC_DIAMOND_end = time.time()
    
    for item in mp_store:
        item.join()
    mp_store[:] = []
    
    ec_detect_out   = os.path.join(ec_annotation_path, "data", "jobs", "ec_detect")
    ec_priam_out    = os.path.join(ec_annotation_path, "data", "jobs", "ec_priam")
    ec_diamond_out  = os.path.join(ec_annotation_path, "data", "jobs", "ec_diamond")
    if check_bypass_log(output_folder_path, ec_annotation_detect_label):
        if(os.path.exists(ec_detect_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_detect_label)
    if check_bypass_log(output_folder_path, ec_annotation_priam_label):
        if(os.path.exists(ec_priam_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_priam_label)
    if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        if(os.path.exists(ec_diamond_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_DIAMOND_label)
    
    #----------------------------------------------------------------------
    # EC post process
    EC_post_start = time.time()
    #if not (check_where_resume(ec_annotation_path, None, GA_DIAMOND_path)):
    if check_bypass_log(output_folder_path, ec_annotation_pp_label):
        
        marker_file = "ec_post"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            command_list = commands.create_EC_postprocess_command(ec_annotation_label, GA_final_merge_label, marker_file)
            launch_and_create_simple(ec_annotation_label, marker_file, commands, command_list)
        
        if(os.path.exists(marker_path)):
            write_to_bypass_log(output_folder_path, ec_annotation_pp_label)
    
    cleanup_EC_start = time.time()
    if(keep_all == "no" and keep_EC == "no"):
        delete_folder(ec_annotation_path)
    elif(keep_all == "compress" or keep_EC == "compress"):
        compress_folder(ec_annotation_path)
        delete_folder(ec_annotation_path)
    cleanup_EC_end = time.time()
    EC_post_end = time.time()
        
   
    EC_end = time.time()
    print("EC run:", '%1.1f' % (EC_end - EC_start), "s")
    print("EC cleanup:", '%1.1f' % (cleanup_EC_end - cleanup_EC_start), "s")
    
    # ------------------------------------------------------
    # RPKM Table and Cytoscape Network
    Cytoscape_start = time.time()
    network_path = os.path.join(output_folder_path, output_label)
    #if not check_where_resume(network_path, None, ec_annotation_path):
    
    if check_bypass_log(output_folder, output_label):
        
        #phase 1
        if check_bypass_log(output_folder, output_copy_gene_map_label):
            job_name = output_copy_gene_map_label
            command_list = commands.create_output_copy_gene_map_command(output_label, GA_final_merge_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
        if check_bypass_log(output_folder, output_copy_taxa_label):
            job_name = output_copy_taxa_label
            command_list = commands.create_output_copy_taxa_command(output_label, taxon_annotation_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
        
        if check_bypass_log(output_folder, output_contig_stats_label):
            job_name = output_contig_stats_label
            command_list = commands.create_output_contig_stats_command(output_label, assemble_contigs_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
        
            
        if not(no_host):
            if check_bypass_log(output_folder, output_unique_hosts_singletons_label):
                job_name = output_unique_hosts_singletons_label
                command_list = commands.create_output_unique_hosts_singletons_command(output_label, quality_filter_label, host_filter_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
            if(read_mode == "paired"):
                if check_bypass_log(output_folder, output_unique_hosts_pair_1_label):
                    job_name = output_unique_hosts_pair_1_label
                    command_list = commands.create_output_unique_hosts_pair_1_command(output_label, quality_filter_label, host_filter_label)
                    launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                    
                if check_bypass_log(output_folder, output_unique_hosts_pair_2_label):
                    job_name = output_unique_hosts_pair_2_label
                    command_list = commands.create_output_unique_hosts_pair_2_command(output_label, quality_filter_label, host_filter_label)
                    launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                    
                    
        #repop vectors
        if check_bypass_log(output_folder, output_unique_vectors_singletons_label):
            job_name = output_unique_vectors_singletons_label
            command_list = commands.create_output_unique_vectors_singletons_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
        
        if(read_mode == "paired"):
            if check_bypass_log(output_folder, output_unique_vectors_pair_1_label):
                job_name = output_unique_vectors_pair_1_label
                command_list = commands.create_output_unique_vectors_pair_1_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                
            if check_bypass_log(output_folder, output_unique_vectors_pair_2_label):
                job_name = output_unique_vectors_pair_2_label
                command_list = commands.create_output_unique_vectors_pair_2_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                
        print(dt.today(), "output report phase 1 launched.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        
        
        conditional_write_to_bypass_log(output_per_read_scores_label, "outputs/final_results", "input_per_seq_quality_report.csv", output_folder_path)
        conditional_write_to_bypass_log(output_copy_gene_map_label, "outputs/final_results", "final_gene_map.tsv", output_folder_path)
        conditional_write_to_bypass_log(output_copy_taxa_label, "outputs/final_results", "taxa_classifications.tsv", output_folder_path)
        conditional_write_to_bypass_log(output_contig_stats_label, "outputs/final_results", "contig_stats.txt", output_folder_path)
        conditional_write_to_bypass_log(output_unique_vectors_singletons_label, "outputs/data/4_full_vectors", "singletons_full_vectors.fastq", output_folder_path)
        if(read_mode == "paired"):
            conditional_write_to_bypass_log(output_unique_vectors_pair_1_label, "outputs/data/4_full_vectors", "pair_1_full_vectors.fastq", output_folder_path)
            conditional_write_to_bypass_log(output_unique_vectors_pair_2_label, "outputs/data/4_full_vectors", "pair_2_full_vectors.fastq", output_folder_path)
            
        if not (no_host):
            conditional_write_to_bypass_log(output_unique_hosts_singletons_label, "outputs/data/2_full_hosts", "singletons_full_hosts.fastq", output_folder_path)
            if(read_mode == "paired"):
                conditional_write_to_bypass_log(output_unique_hosts_pair_1_label, "outputs/data/2_full_hosts", "pair_1_full_hosts.fastq", output_folder_path)
                conditional_write_to_bypass_log(output_unique_hosts_pair_2_label, "outputs/data/2_full_hosts", "pair_2_full_hosts.fastq", output_folder_path)
        #----------------------------------------------------------------------------
        #Phase 2
        if check_bypass_log(output_folder, output_network_gen_label):
            command_list = commands.create_output_network_generation_command(output_label, GA_final_merge_label, taxon_annotation_label, ec_annotation_label)
            launch_and_create_with_mp_store(mp_store, output_label, output_network_gen_label, commands, command_list)
            
        if check_bypass_log(output_folder, output_taxa_groupby_label):
            command_list = commands.create_output_taxa_groupby_command(output_label)
            launch_and_create_with_mp_store(mp_store, output_label, output_taxa_groupby_label, commands, command_list)
       
        print(dt.today(), "output report phase 2 launched.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        conditional_write_to_bypass_log(output_network_gen_label, "outputs/final_results", "RPKM_table.tsv", output_folder_path)
        
        
        #-------------------------------------------------------------------
        #Phase 3
        if check_bypass_log(output_folder, output_read_count_label):
            job_name = output_read_count_label
            command_list = commands.create_output_read_count_command(output_label, quality_filter_label, repop_job_label, GA_final_merge_label, ec_annotation_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                               

        if check_bypass_log(output_folder, output_per_read_scores_label):
            job_name = output_per_read_scores_label
            command_list = commands.create_output_per_read_scores_command(output_label, quality_filter_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
        if check_bypass_log(output_folder, output_ec_heatmap_label):
            job_name = output_ec_heatmap_label
            command_list = commands.create_output_EC_heatmap_command(output_label)
            launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)    
        
        print(dt.today(), "output report phase 3 launched.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        conditional_write_to_bypass_log(output_read_count_label, "outputs/final_results", "read_count.tsv", output_folder_path)
        conditional_write_to_bypass_log(output_ec_heatmap_label, "outputs/final_results", "EC_coverage.csv", output_folder_path)

        
    cleanup_cytoscape_start = time.time()
    if(keep_all == "no" and keep_outputs == "no"):
        delete_folder(network_path)
    elif(keep_all == "compress" or keep_outputs == "compress"):
        compress_folder(network_path)
        delete_folder(network_path)
    cleanup_cytoscape_end = time.time()
        
        
        
        
    Cytoscape_end = time.time()
    end_time = time.time()
    print("Outputs:", '%1.1f' % (Cytoscape_end - Cytoscape_start - (cleanup_cytoscape_end - cleanup_cytoscape_start)), "s")
    print("Outputs cleanup:", '%1.1f' % (cleanup_cytoscape_end - cleanup_cytoscape_start), "s")
    print("=============================================================================================")
    print("Final summary")
    print("--------------------------------------------------------")
    print("Total runtime:", '%1.1f' % (end_time - start_time), "s")
    print("quality filter:", '%1.1f' % (quality_end - quality_start - (cleanup_quality_end - cleanup_quality_start)), "s")
    print("quality filter cleanup:", '%1.1f' %(cleanup_quality_end - cleanup_quality_start), "s")
    if not no_host:
        print("host filter:", '%1.1f' % (host_end - host_start - (cleanup_host_end - cleanup_host_start)), "s")
        print("host filter cleanup:", '%1.1f' %(cleanup_host_end - cleanup_host_start),"s")
    print("vector filter:", '%1.1f' % (vector_end - vector_start - (cleanup_vector_end - cleanup_vector_start)), "s")
    print("vector filter cleanup:", '%1.1f' % (cleanup_vector_end - cleanup_vector_start), "s")
    print("rRNA filter:", '%1.1f' % (rRNA_filter_end - rRNA_filter_start - (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start)), "s")
    print("rRNA filter cleanup:", '%1.1f' % (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start), "s")
    print("repop:", '%1.1f' % (repop_end - repop_start - (cleanup_repop_end - cleanup_repop_start)), "s")
    print("repop cleanup:", '%1.1f' % (cleanup_repop_end - cleanup_repop_start), "s")
    print("assemble contigs:", '%1.1f' % (assemble_contigs_end - assemble_contigs_start - (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start)), "s")    
    print("assemble contigs cleanup:", '%1.1f' % (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start), "s")
    print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
    print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
    print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
    print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
    print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
    print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
    print("TA:", '%1.1f' % (TA_end - TA_start - (cleanup_TA_end - cleanup_TA_start)), "s")
    print("TA cleanup:", '%1.1f' % (cleanup_TA_end - cleanup_TA_start), "s")
    print("EC:", '%1.1f' % (EC_end - EC_start), "s")
    #print("---------------------------------------------")
    #print("Note: EC is in cloud-mode.  ignore individual timing")
    #print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    #print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
    #print("EC DIAMOND:", '%1.1f' % (EC_DIAMOND_end - EC_DIAMOND_start), "s")
    #print("EC cleanup:", '%1.1f' % (cleanup_EC_end - cleanup_EC_start), "s")
    #print("-------------------------------------------------")
    print("Outputs:", '%1.1f' % (Cytoscape_end - Cytoscape_start - (cleanup_cytoscape_end - cleanup_cytoscape_start)), "s")
    print("Outputs cleanup:", '%1.1f' % (cleanup_cytoscape_end - cleanup_cytoscape_start), "s")
    

def tutorial_main(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode_string): 
    paths = mpp.tool_path_obj(config_path)
    no_host = args_pack["no_host"]
    verbose_mode = args_pack["verbose_mode"]
    rRNA_chunks = int(paths.rRNA_chunksize)
    GA_chunksize = int(paths.GA_chunksize)
    
    BWA_mem_threshold = int(paths.BWA_mem_threshold)
    BLAT_mem_threshold = int(paths.BLAT_mem_threshold)
    DIAMOND_mem_threshold = int(paths.DIAMOND_mem_threshold)
    BWA_pp_mem_threshold = int(paths.BWA_pp_mem_threshold)
    BLAT_pp_mem_threshold = int(paths.BLAT_pp_mem_threshold)
    DIAMOND_pp_mem_threshold = int(paths.DIAMOND_pp_mem_threshold)
    Infernal_mem_threshold = int(paths.Infernal_mem_threshold)
    Barrnap_mem_threshold = int(paths.Barrnap_mem_threshold)
    DETECT_mem_threshold = int(paths.DETECT_mem_threshold)
    TA_mem_threshold = int(paths.TA_mem_threshold)
    filter_stringency = paths.filter_stringency
    
    
    BWA_job_limit = int(paths.BWA_job_limit)
    BLAT_job_limit = int(paths.BLAT_job_limit)
    DIAMOND_job_limit = int(paths.DIAMOND_job_limit)
    BWA_pp_job_limit = int(paths.BWA_pp_job_limit)
    BLAT_pp_job_limit = int(paths.BLAT_pp_job_limit)
    DIAMOND_pp_job_limit = int(paths.DIAMOND_pp_job_limit)
    Infernal_job_limit = int(paths.Infernal_job_limit)
    Barrnap_job_limit = int(paths.Barrnap_job_limit)
    DETECT_job_limit = int(paths.DETECT_job_limit)
    TA_job_limit = int(paths.TA_job_limit)
    
    Infernal_job_delay      = float(paths.Infernal_job_delay)
    Barrnap_job_delay       = float(paths.Barrnap_job_delay)
    BWA_job_delay           = float(paths.BWA_job_delay)
    BLAT_job_delay          = float(paths.BLAT_job_delay)
    DIAMOND_job_delay       = float(paths.DIAMOND_job_delay)
    BWA_pp_job_delay        = float(paths.BWA_pp_job_delay)
    BLAT_pp_job_delay       = float(paths.BLAT_pp_job_delay)
    DIAMOND_pp_job_delay    = float(paths.DIAMOND_pp_job_delay)
    DETECT_job_delay        = float(paths.DETECT_job_delay)
    TA_job_delay            = float(paths.TA_job_delay)
    
    #-----------------------------------------------------
    keep_all                = paths.keep_all
    keep_quality            = paths.keep_quality
    keep_vector             = paths.keep_vector
    keep_host               = paths.keep_host
    keep_rRNA               = paths.keep_rRNA
    keep_repop              = paths.keep_repop
    keep_assemble_contigs   = paths.keep_assemble_contigs
    keep_GA_BWA             = paths.keep_GA_BWA
    keep_GA_BLAT            = paths.keep_GA_BLAT
    keep_GA_DIAMOND         = paths.keep_GA_DIAMOND
    keep_GA_final           = paths.keep_GA_final
    keep_TA                 = paths.keep_TA
    keep_EC                 = paths.keep_EC
    keep_outputs            = paths.keep_outputs
    
    #------------------------------------------------------------------------
    
    
    print("============================================================")
    print("data cleaner options:")
    print("keep all:", keep_all)
    print("keep quality:", keep_quality)
    print("keep vector:", keep_vector)
    print("keep host:", keep_host)
    print("keep rRNA:", keep_rRNA)
    print("keep repop:", keep_repop)
    print("keep assemble contigs:", keep_assemble_contigs)
    print("keep GA BWA:", keep_GA_BWA)
    print("keep GA BLAT:", keep_GA_BLAT)
    print("keep GA DIAMOND:", keep_GA_DIAMOND)
    print("keep GA final:", keep_GA_final)
    print("keep TA:", keep_TA)
    print("keep EC:", keep_EC)
    print("keep outputs:", keep_outputs)
    print("===============================================")
    print("Job delay options:")
    print("Infernal job delay:", Infernal_job_delay)
    print("Barrnap job delay:", Barrnap_job_delay)
    print("BWA job delay:", BWA_job_delay)
    print("BLAT job delay:", BLAT_job_delay)
    print("DIAMOND job delay:", DIAMOND_job_delay)
    print("BWA pp job delay:", BWA_pp_job_delay)
    print("BLAT pp job delay:", BLAT_pp_job_delay)
    print("DIAMOND pp job delay:", DIAMOND_pp_job_delay)
    print("DETECT job delay:", DETECT_job_delay)
    print("=================================================")
    print("memory thresholds")
    print("Barrnap mem threshold:", Barrnap_mem_threshold)
    print("Infernal mem threshold:", Infernal_mem_threshold)
    print("BWA mem threshold:", BWA_mem_threshold)
    print("BLAT mem threshold:", BLAT_mem_threshold)
    print("DIAMOND mem threshold:", DIAMOND_mem_threshold)
    print("BWA pp mem threshold:", BWA_pp_mem_threshold)
    print("BLAT pp mem threshold:", BLAT_pp_mem_threshold)
    print("DIAMOND pp mem threshold:", DIAMOND_pp_mem_threshold)
    print("DETECT mem threshold:", DETECT_mem_threshold)
    print("======================================================")
    print("Barrnap job limit:", Barrnap_job_limit)
    print("Infernal job limit:", Infernal_job_limit)
    print("BWA job limit:", BWA_job_limit)
    print("BLAT job limit:", BLAT_job_limit)
    print("DIAMOND job limit:", DIAMOND_job_limit)
    print("BWA pp job limit:", BWA_pp_job_limit)
    print("BLAT pp job limit:", BLAT_pp_job_limit)
    print("DIAMOND pp job limit:", DIAMOND_pp_job_limit)
    print("DETECT job limit:", DETECT_job_limit)
    print("===================================================")
    print("Filter stringency:", filter_stringency)
    print("rRNA filter Chunk size:", rRNA_chunks)
    print("GA chunk size:", GA_chunksize)
    print("---------------------------------")
    if not single_path == "":
        read_mode = "single"
        quality_encoding = determine_encoding(single_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN SINGLE-ENDED MODE")
    else:
        read_mode = "paired"
        quality_encoding = determine_encoding(pair_1_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN PAIRED-MODE")
        
    if threads == 0:
        real_thread_count = mp.cpu_count()
    else:
        real_thread_count = threads
       
    if(real_thread_count == 1):
        real_thread_count = 2
    print("number of threads used:", real_thread_count)         
            
    
    quality_filter_label                    = "quality_filter"
    host_filter_label                       = "host_read_filter"
    vector_filter_label                     = "vector_read_filter"
    rRNA_filter_label                       = "rRNA_filter"
    rRNA_filter_split_label                 = "rRNA_filter_split"
    rRNA_filter_convert_label               = "rRNA_filter_convert"
    rRNA_filter_barrnap_label               = "rRNA_filter_barrnap"
    rRNA_filter_barrnap_merge_label         = "rRNA_filter_barrnap_merge"
    rRNA_filter_barrnap_pp_label            = "rRNA_filter_barrnap_pp"
    rRNA_filter_infernal_label              = "rRNA_filter_infernal"
    rRNA_filter_infernal_prep_label         = "rRNA_filter_infernal_prep"
    rRNA_filter_splitter_label              = "rRNA_filter_splitter"
    rRNA_filter_post_label                  = "rRNA_filter_post"
    repop_job_label                         = "duplicate_repopulation"
    assemble_contigs_label                  = "assemble_contigs"
    destroy_contigs_label                   = "destroy_contigs"
    GA_BWA_label                            = "GA_BWA"
    GA_BWA_pp_label                         = "GA_BWA_pp"
    GA_BLAT_label                           = "GA_BLAT"
    GA_BLAT_cleanup_label                   = "GA_BLAT_cleanup"
    GA_BLAT_cat_label                       = "GA_BLAT_cat"
    GA_BLAT_pp_label                        = "GA_BLAT_pp"
    GA_DIAMOND_label                        = "GA_DIAMOND"
    GA_DIAMOND_pp_label                     = "GA_DIAMOND_pp"
    GA_final_merge_label                    = "GA_FINAL_MERGE"
    taxon_annotation_label                  = "taxonomic_annotation"
    ec_annotation_label                     = "enzyme_annotation"
    ec_annotation_detect_label              = "enzyme_annotation_detect"
    ec_annotation_priam_label               = "enzyme_annotation_priam"
    ec_annotation_DIAMOND_label             = "enzyme_annotation_DIAMOND"
    ec_annotation_pp_label                  = "enzyme_annotation_pp"
    output_label                            = "outputs"
    output_copy_gene_map_label              = "output_copy_gene_map"
    output_clean_EC_label                   = "output_clean_ec"
    output_copy_taxa_label                  = "output_copy_taxa"
    output_network_gen_label                = "output_network_generation"
    output_unique_hosts_singletons_label    = "output_unique_hosts_singletons"
    output_unique_hosts_pair_1_label        = "output_unique_hosts_pair_1"
    output_unique_hosts_pair_2_label        = "output_unique_hosts_pair_2"
    output_unique_vectors_singletons_label  = "output_unique_vectors_singletons"
    output_unique_vectors_pair_1_label      = "output_unique_vectors_pair_1"
    output_unique_vectors_pair_2_label      = "output_unique_vectors_pair_2"
    output_combine_hosts_label              = "output_combine_hosts"
    output_per_read_scores_label            = "output_per_read_scores"
    output_contig_stats_label               = "output_contig_stats"
    output_ec_heatmap_label                 = "output_ec_heatmap"
    output_taxa_groupby_label               = "output_taxa_groupby"
    output_read_count_label                 = "output_read_count"
            
            
    mp_store = []  # stores the multiprocessing processes    
    
    # Creates our command object, for creating shellscripts.
    if read_mode == "single":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = tutorial_mode_string, sequence_path_1=None, sequence_path_2=None, sequence_single=single_path, sequence_contigs = contig_path)
    elif read_mode == "paired":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = tutorial_mode_string, sequence_path_1=pair_1_path, sequence_path_2=pair_2_path, sequence_single=single_path, sequence_contigs = contig_path)
        
    if(tutorial_mode_string == "quality"):
        print(dt.today(), "working on:", tutorial_mode_string)
        command_list = commands.create_quality_control_command(quality_filter_label)
        job_name = quality_filter_label
        launch_and_create_simple(quality_filter_label, job_name, commands, command_list)
        
    elif(tutorial_mode_string == "host"):
        print(dt.today(), "working on:", tutorial_mode_string)
        job_name = host_filter_label
        command_list = commands.create_host_filter_command(host_filter_label, quality_filter_label)
        launch_and_create_simple(host_filter_label, job_name, commands, command_list)
        
    elif(tutorial_mode_string == "vector"):
        print(dt.today(), "working on:", tutorial_mode_string)
        job_name = vector_filter_label
        command_list = commands.create_vector_filter_command(vector_filter_label, quality_filter_label)
        launch_and_create_simple(vector_filter_label, job_name, commands, command_list)
        
    elif(tutorial_mode_string == "rRNA"):
        print(dt.today(), "working on:", tutorial_mode_string)
        rRNA_filter_path = os.path.join(output_folder_path, rRNA_filter_label)
        rRNA_filter_jobs_folder = os.path.join(rRNA_filter_path, "data", "jobs")
        marker_path_list = []
        sections = ["singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        job_name = "rRNA_filter_prep_tutorial"
        command_list = commands.create_rRNA_filter_prep_command_v3(rRNA_filter_label, "tutorial", vector_filter_label, "tutorial")
        launch_and_create_simple(rRNA_filter_label, job_name, commands, command_list)
        
        
        for section in reversed(sections):
            split_path = os.path.join(rRNA_filter_path, "data", section + "_fastq")
            if check_bypass_log(output_folder_path, rRNA_filter_convert_label + "_" + section):
                marker_path_list = []
                for item in os.listdir(split_path):
                    root_name = item.split(".")[0]
                    fasta_path = os.path.join(rRNA_filter_path, "data", section + "_fasta")
                    fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                    marker_file = root_name + "_convert_fasta"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    
                    fasta_out_size = os.stat(fasta_file).st_size if (os.path.exists(fasta_file)) else 0
                    if(os.path.exists(marker_path)):
                        print(dt.today(), item, "already converted to fasta.  skipping")
                        continue
                    else:
                        job_name = root_name + "_convert_to_fasta"
                        marker_path_list.append(marker_path)
                        command_list = commands.create_rRNA_filter_convert_fastq_command("rRNA_filter", section, root_name+".fastq", marker_file)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                        
                        
        # BARRNAP
        for section in reversed(sections):  
            #convert data to fasta, then run barrnap separately, then cat the barrnap, then run barrnap PP
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            split_path      = os.path.join(rRNA_filter_path, "data", section + "_fastq")
            fasta_path      = os.path.join(rRNA_filter_path, "data", section + "_fasta")
            barrnap_path    = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap")
            infernal_path   = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal") 
            
            mRNA_path       = os.path.join(rRNA_filter_path, "data", section + "_mRNA")
            
            concurrent_job_count = 0
            batch_count = 0
            marker_path_list = []
            
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                final_marker_file = root_name + "_barrnap_cat"
                final_marker_path = os.path.join(rRNA_filter_jobs_folder, final_marker_file)
                barrnap_arc_out_file = os.path.join(barrnap_path, root_name + "_arc.barrnap_out")
                barrnap_bac_out_file = os.path.join(barrnap_path, root_name + "_bac.barrnap_out")
                barrnap_euk_out_file = os.path.join(barrnap_path, root_name + "_euk.barrnap_out")
                barrnap_mit_out_file = os.path.join(barrnap_path, root_name + "_mit.barrnap_out")
                final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                fastq_file = os.path.join(split_path, root_name + ".fastq")
                barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                marker_file_arc = root_name + "_barrnap_arc"
                marker_file_bac = root_name + "_barrnap_bac"
                marker_file_euk = root_name + "_barrnap_euk"
                marker_file_mit = root_name + "_barrnap_mit"
                
                marker_path_arc = os.path.join(rRNA_filter_jobs_folder, marker_file_arc)
                marker_path_bac = os.path.join(rRNA_filter_jobs_folder, marker_file_bac)
                marker_path_euk = os.path.join(rRNA_filter_jobs_folder, marker_file_euk)
                marker_path_mit = os.path.join(rRNA_filter_jobs_folder, marker_file_mit)
                
                barrnap_arc_out_size    = os.stat(barrnap_arc_out_file).st_size if (os.path.exists(barrnap_arc_out_file)) else 0
                barrnap_bac_out_size    = os.stat(barrnap_bac_out_file).st_size if (os.path.exists(barrnap_bac_out_file)) else 0
                barrnap_euk_out_size    = os.stat(barrnap_euk_out_file).st_size if (os.path.exists(barrnap_euk_out_file)) else 0
                barrnap_mit_out_size    = os.stat(barrnap_mit_out_file).st_size if (os.path.exists(barrnap_mit_out_file)) else 0
                
                if(os.path.exists(final_marker_path)):
                    print(dt.today(), "Job already run. skipping:", final_marker_file, final_marker_path)
                    time.sleep(5)
                    continue
                else:
                    if((barrnap_arc_out_size > 0) and (os.path.exists(marker_path_arc))):
                        print(dt.today(), "barrnap arc already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_arc"
                        marker_path_list.append(marker_path_arc)
                        command_list = commands.create_rRNA_filter_barrnap_arc_command("rRNA_filter", section, root_name, marker_file_arc)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                        
                        
                    if((barrnap_bac_out_size > 0) and (os.path.exists(marker_path_bac))):
                        print(dt.today(), "barrnap bac already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_bac"
                        marker_path_list.append(marker_path_bac)
                        command_list = commands.create_rRNA_filter_barrnap_bac_command("rRNA_filter", section, root_name, marker_file_bac)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                        
                    if((barrnap_euk_out_size > 0) and (os.path.join(marker_path_euk))):
                        print(dt.today(), "barrnap euk already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_euk"
                        marker_path_list.append(marker_path_euk)
                        command_list = commands.create_rRNA_filter_barrnap_euk_command("rRNA_filter", section, root_name, marker_file_euk)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
                        
                    if((barrnap_mit_out_size > 0) and (os.path.join(marker_path_mit))):
                        print(dt.today(), "barrnap mit already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_mit"
                        marker_path_list.append(marker_path_mit)
                        command_list = commands.create_rRNA_filter_barrnap_mit_command("rRNA_filter", section, root_name, marker_file_mit)
                        launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
            print(dt.today(), "waiting for Barrnap jobs to finish")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_" + section + ".txt")
            check_all_job_markers(marker_path_list, final_checklist)
        
            #------------------------------------------------------
            #merge the barrnap data
            marker_path_list = []
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                marker_file = root_name + "_barrnap_cat"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                final_barrnap_out_size  = os.stat(final_barrnap_out).st_size if (os.path.exists(final_barrnap_out)) else 0
                
                if((final_barrnap_out_size > 0) and (os.path.exists(marker_path))):
                    print(dt.today(), "barrnap already merged. skipping:", item)
                    continue
                else:
                    job_name = root_name + "_barrnap_cat"
                    marker_path_list.append(marker_path)
                    command_list = commands.create_rRNA_filter_barrnap_cat_command("rRNA_filter", section, root_name, marker_file)
                    launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
            print(dt.today(), "waiting for Barrnap pp to finish")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_cat_" + section  + ".txt")
            check_all_job_markers(marker_path_list, final_checklist)
            
            #-----------------------------------------------------
            #run the barrnap PP
            marker_path_list = []
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                barrnap_mRNA_out_size   = os.stat(barrnap_mrna_file).st_size if (os.path.exists(barrnap_mrna_file)) else 0
                marker_file = root_name + "_barrnap_pp"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                if(barrnap_mRNA_out_size > 0):
                    print(dt.today(), "barrnap pp already run.  skipping:", item)
                    continue
                else:
                    job_name = root_name + "_barrnap_pp"
                    marker_path_list.append(marker_path)
                    command_list = commands.create_rRNA_filter_barrnap_pp_command("rRNA_filter", section, root_name + ".fastq", marker_file)
                    launch_only_with_hold(mp_store, Barrnap_mem_threshold, Barrnap_job_limit, Barrnap_job_delay, job_name, commands, command_list)
            
            print(dt.today(), "waiting for Barrnap pp to finish")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_barrnap_" + section +  ".txt")
            check_all_job_markers(marker_path_list, final_checklist)
                
        #----------------------------------------------------------------------------
        # INFERNAL
        for section in reversed(sections):  
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            barrnap_mRNA_fastq_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap_mRNA")
            infernal_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal") 
            barrnap_mRNA_fasta_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_barrnap_mRNA_fasta")
            splitter_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal_mRNA")
        
            concurrent_job_count = 0
            batch_count = 0
            #these jobs now have to be launched in segments
            for item in os.listdir(barrnap_mRNA_fastq_path):
               
                if(item.endswith("_barrnap_mRNA.fastq")):
                    root_name = item.split(".")[0]
                    marker_file = root_name + "_infernal_prep"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    infernal_prep_out_file = os.path.join(barrnap_mRNA_fasta_path, root_name + ".fasta")
                    infernal_prep_file_size = os.stat(infernal_prep_out_file).st_size if (os.path.exists(infernal_prep_out_file)) else 0
                    if((infernal_prep_file_size > 0) and (os.path.exists(marker_path))):
                        print(dt.today(), "Infernal prep already ran on this sample.  skipping", item)
                        continue
                    
                    else:
                        marker_path_list.append(marker_path)
                        job_name = "rRNA_filter_infernal_prep_" + root_name
                        command_list = commands.create_rRNA_filter_infernal_prep_command("rRNA_filter", section, item, root_name, marker_file)
                        launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
                        
            print(dt.today(), "final batch: infernal prep")
            for p_item in mp_store:
                p_item.join()
            mp_store[:] = []  # clear the list
            final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_prep_" + section + ".txt")
            check_all_job_markers(marker_path_list, final_checklist)
            

            marker_path_list = []
            for item in os.listdir(barrnap_mRNA_fasta_path):
                #using a job marker is ineffective.  The marker will still write 
                root_name = item.split("_barrnap_mRNA")[0]
                marker_file = root_name + "_infernal"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "infernal already run. skipping:", root_name + "_infernal")
                    continue
                else:
                    marker_path_list.append(marker_path)
                    inf_command = commands.create_rRNA_filter_infernal_command("rRNA_filter", section, root_name, marker_file)
                    job_name = "rRNA_filter_infernal_" + root_name
                    #launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, inf_command)
                    launch_and_create_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, rRNA_filter_label, job_name, commands, inf_command)
                    
                    
            print(dt.today(), "final batch: infernal")
            for p_item in mp_store:
                p_item.join()
            mp_store[:] = []
            final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_" + section + ".txt")
            check_all_job_markers(marker_path_list, final_checklist)
            
            if (section != "pair_2"):
                marker_path_list = []
                for item in os.listdir(barrnap_mRNA_fasta_path):
                    root_name = item.split("_barrnap_mRNA")[0]
                    splitter_out_file = os.path.join(output_folder_path, rRNA_filter_label, "data", section + "_infernal_mRNA", root_name + "_mRNA.fastq")
                    splitter_out_file_size = os.stat(splitter_out_file).st_size if os.path.exists(splitter_out_file) else 0
                    marker_file = root_name + "_infernal_pp"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    if((splitter_out_file_size > 0) and (os.path.exists(marker_path))):
                        print(dt.today(), "infernal mRNA splitter already run. skipping:", marker_file)
                        print("file size:", splitter_out_file_size, "file:", splitter_out_file)
                        continue
                    else:
                        job_name = "rRNA_filter_infernal_splitter_" + root_name
                        marker_path_list.append(marker_path)
                        command_list = commands.create_rRNA_filter_splitter_command("rRNA_filter", section, root_name, marker_file)
                        print(command_list)
                        launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
                        
                print(dt.today(), "final batch: infernal splitter")
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []
                final_checklist = os.path.join(rRNA_filter_path, "rRNA_filter_infernal_splitter_" + section + ".txt")
                check_all_job_markers(marker_path_list, final_checklist)
                
            else:
                print(dt.today(), "not calling Infernal rRNA splitter on pair 2.  data handled by pair 1 as a combination")
        
                    
                    
        marker_path_list = []
        for section in reversed(sections):
            print(dt.today(), "now running rRNA filter post:", section)
            marker_file = section + "_rRNA_packup"
            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
            job_name = "rRNA_post_cat"
            marker_path_list.append(marker_path)
            command_list = commands.create_rRNA_filter_final_cat_command("rRNA_filter", section, marker_file)
            print("command list:", command_list)
            launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, commands, command_list)
            
        for p_item in mp_store:
            p_item.join()
        mp_store[:] = []
        
    elif(tutorial_mode_string == "repop"):
        print(dt.today(), "working on:", tutorial_mode_string)
        quality_path = os.path.join(output_folder_path, quality_filter_label)
        if(os.path.exists(quality_path)):
            job_name = repop_job_label
            command_list = commands.create_repop_command(repop_job_label, quality_filter_label, rRNA_filter_label)
            launch_and_create_simple(repop_job_label, job_name, commands, command_list)
        
            
        else:
            print(dt.today(), "no quality filter was run.  please run MetaPro with option  <--tutorial quality> to continue")
        
    elif(tutorial_mode_string == "contigs"):
        print(dt.today(), "working on:", tutorial_mode_string)
        assemble_contigs_path = os.path.join(output_folder_path, assemble_contigs_label)
        job_name = assemble_contigs_label
        command_list = commands.create_assemble_contigs_command(assemble_contigs_label, repop_job_label)
        launch_and_create_simple(assemble_contigs_label, job_name, commands, command_list)
        mgm_file = os.path.join(assemble_contigs_path, "data", "1_mgm", "gene_report.txt")
        if(os.path.exists(mgm_file)):
            print(dt.today(), "MetaGeneMark ran successfully.  Everything is great")
        else:
            sys.exit("mgm did not run.  look into it.  Proabably a license issue. pipeline stopping here")
        
    elif(tutorial_mode_string == "GA"):
    
        GA_BWA_start = time.time()
        GA_BWA_path = os.path.join(output_folder_path, GA_BWA_label)
        GA_BWA_jobs_folder = os.path.join(GA_BWA_path, "data", "jobs")
        #if not check_where_resume(GA_BWA_path, None, assemble_contigs_path):
        marker_path_list = []
        marker_file = "GA_split_fasta_contigs"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping", marker_file)
        else:
            job_name = "GA_prep_split_contigs"
            marker_path_list.append(marker_path)
            command_list = commands.create_split_ga_fasta_data_command(GA_BWA_label, assemble_contigs_label, "contigs", marker_file)
            launch_and_create_with_mp_store(mp_store, GA_BWA_label, job_name, commands, command_list)
        
        
        sections = ["singletons"]
        if(read_mode == "paired"):
            sections.extend(["pair_1", "pair_2"])
        for section in sections: 
            marker_file = "GA_split_fastq_" + section
            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping", marker_file)
            else:
                marker_path_list.append(marker_path)
                job_name = "GA_prep_split_" + section
                command_list = commands.create_split_ga_fastq_data_command(GA_BWA_label, assemble_contigs_label, section, marker_file)
                launch_and_create_with_mp_store(mp_store, GA_BWA_label, job_name, commands, command_list)
        
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_prep.txt")

            
        #-------------------------------------------------------------------------
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                job_submitted = False
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA" + "_" + file_tag
                marker_file = file_tag + "_bwa"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                #this checker assumes that BWA only exports a file when it's finished running
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BWA_annotate_command_v2(GA_BWA_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, BWA_mem_threshold, BWA_job_limit, BWA_job_delay, GA_BWA_label, job_name, commands, command_list)

        print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA.txt")
        check_all_job_markers(marker_path_list, final_checklist)
            
        marker_path_list = []
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA_pp" + "_" + file_tag
                marker_file = file_tag + "_bwa_pp"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BWA_pp_command_v2(GA_BWA_label, assemble_contigs_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, BWA_pp_mem_threshold, BWA_pp_job_limit, BWA_pp_job_delay, GA_BWA_label, job_name, commands, command_list)
                        
        print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
        for item in mp_store:
            item.join()
        mp_store[:] = []
        marker_file = "BWA_copy_contig_map"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:   
            marker_path_list.append(marker_path)
            command_list = commands.create_BWA_copy_contig_map_command(GA_BWA_label, assemble_contigs_label, marker_file)
            launch_and_create_simple(GA_BWA_label, GA_BWA_label + "_copy_contig_map", commands, command_list)
        
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_pp.txt")
        check_all_job_markers(marker_path_list, final_checklist)
        
        cleanup_GA_BWA_start = time.time()
        
        GA_BWA_data_path = os.path.join(GA_BWA_path, "data")
        if(keep_all == "no" and keep_GA_BWA == "no"):
            delete_folder(GA_BWA_data_path)
        elif(keep_all == "compress" or keep_GA_BWA == "compress"):
            compress_folder(GA_BWA_path)
            delete_folder(GA_BWA_data_path)
        cleanup_GA_BWA_end = time.time()
        GA_BWA_end = time.time()
        print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
        print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
        
        # ------------------------------------------------
        # BLAT gene annotation
        GA_BLAT_start = time.time()
        GA_BLAT_path = os.path.join(output_folder_path, GA_BLAT_label)
        GA_BLAT_jobs_folder = os.path.join(GA_BLAT_path, "data", "jobs")
        GA_BLAT_final_job_marker = os.path.join(GA_BLAT_path, "all_BLAT")
        if (os.path.exists(GA_BLAT_final_job_marker)):
            print(dt.today(), "BLAT was run, skipping")
        else:
            
        
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    for fasta_db in os.listdir(paths.DNA_DB_Split):
                        if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                            job_name = "BLAT_" + file_tag + "_" + fasta_db
                            marker_file = file_tag + "_blat_" + fasta_db
                            marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                            #This checker assume BLAT only exports a file when it's finished running
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "BLAT job ran already, skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                                command_list = commands.create_BLAT_annotate_command_v2(GA_BLAT_label, full_sample_path, fasta_db, marker_file)
                                launch_only_with_hold(mp_store, BLAT_mem_threshold, BLAT_job_limit, BLAT_job_delay, job_name, commands, command_list)
                                
                                    
            print(dt.today(), "final BLAT job removal")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            

                
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    job_name = file_tag + "_cat"
                    
                    marker_file = file_tag + "_blat_cat"
                    marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = commands.create_BLAT_cat_command_v2(GA_BLAT_label, full_sample_path, marker_file)
                        launch_only_with_hold(mp_store, BLAT_mem_threshold, BLAT_job_limit, BLAT_job_delay, job_name, commands, command_list)
                    
            for item in mp_store:
                item.join()
            mp_store[:] = []
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_cat.txt")
            check_all_job_markers(marker_path_list, final_checklist)
                
            
            
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    job_name = "BLAT_" + file_tag + "_pp"
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    marker_file = file_tag + "_blat_pp"
                    marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = commands.create_BLAT_pp_command_v2(GA_BLAT_label, full_sample_path, GA_BWA_label, marker_file)
                        launch_and_create_with_hold(mp_store, BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, GA_BLAT_label, job_name, commands, command_list)
                    
            print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            
            job_name = "GA_BLAT_copy_contigs"
            marker_file = "blat_copy_contig_map"
            marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_BLAT_copy_contig_map_command(GA_BLAT_label, GA_BWA_label, marker_file)
                launch_and_create_simple(GA_BLAT_label, job_name, commands, command_list)
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_pp.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            delete_folder_simple(GA_BLAT_jobs_folder)
            open(GA_BLAT_final_job_marker, "a").close()
            
        cleanup_GA_BLAT_start = time.time()
        if(keep_all == "no" and keep_GA_BLAT == "no"):
            delete_folder(GA_BLAT_path)
        elif(keep_all == "compress" or keep_GA_BLAT == "compress"):
            compress_folder(GA_BLAT_path)
            delete_folder(GA_BLAT_path)
        cleanup_GA_BLAT_end = time.time()
        GA_BLAT_end = time.time()
        print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
        print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
        
        # ------------------------------------------------------
        # Diamond gene annotation
        GA_DIAMOND_start = time.time()
        GA_DIAMOND_path = os.path.join(output_folder_path, GA_DIAMOND_label)
        GA_DIAMOND_tool_output_path = os.path.join(GA_DIAMOND_path, "data", "0_diamond")
        GA_DIAMOND_jobs_folder = os.path.join(GA_DIAMOND_path, "data", "jobs")
        #if not check_where_resume(None, GA_DIAMOND_tool_output_path, GA_BLAT_path, file_check_bypass = True):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_path)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_annotate_command_v2(GA_DIAMOND_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, DIAMOND_mem_threshold, DIAMOND_job_limit, DIAMOND_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                
        print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND.txt")
        check_all_job_markers(marker_path_list, final_checklist)
            
            
        #if not check_where_resume(GA_DIAMOND_path, None, GA_DIAMOND_tool_output_path, file_check_bypass = True):
        print(dt.today(), "DIAMOND PP threads used:", real_thread_count/2)
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_pp_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond_pp"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_pp_command_v2(GA_DIAMOND_label, GA_BLAT_label, full_sample_path, marker_file)
                    launch_and_create_with_hold(mp_store, DIAMOND_pp_mem_threshold, DIAMOND_pp_job_limit, DIAMOND_pp_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                                    
        print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND_pp.txt")
        check_all_job_markers(marker_path_list, final_checklist)
   
            
        
        cleanup_GA_DIAMOND_start = time.time()
        if(keep_all == "no" and keep_GA_DIAMOND == "no"):
            delete_folder(GA_DIAMOND_path)
        elif(keep_all == "compress" or keep_GA_DIAMOND == "compress"):
            compress_folder(GA_DIAMOND_path)
            delete_folder(GA_DIAMOND_path)
        cleanup_GA_DIAMOND_end = time.time()
        GA_DIAMOND_end = time.time()
        print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
        print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
        
        
        GA_final_merge_start = time.time()
        GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
        marker_file = "GA_final_merge"
        marker_path = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge")
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping: GA final merge")
        else:
            command_list = commands.create_GA_final_merge_command(GA_final_merge_label, GA_BWA_label, GA_BLAT_label, GA_DIAMOND_label, assemble_contigs_label, marker_file)
            job_name = "GA_final_merge"
            launch_and_create_simple(GA_final_merge_label, job_name, commands, command_list)
        
        #check if all_proteins.faa was generated
        all_proteins_path = os.path.join(output_folder_path, GA_final_merge_label, "final_results", "all_proteins.faa")
        if(os.path.getsize(all_proteins_path) > 0):
            print(dt.today(), "All_proteins.faa is OK.  Continuing")
        else:
            sys.exit("GA final merge failed.  proteins weren't translated")
            
        GA_final_merge_end = time.time()
        print("GA final merge:", '%1.1f' % (GA_final_merge_end - GA_final_merge_start), "s")
        if(keep_all == "no" and keep_GA_final == "no"):
            delete_folder(GA_FINAL_MERGE_path)
        elif(keep_all == "compress" or keep_GA_final == "compress"):
            compress_folder(GA_FINAL_MERGE_path)
            delete_folder(GA_FINAL_MERGE_path)
          
    elif(tutorial_mode_string == "TA"):
        print(dt.today(), "working on:", tutorial_mode_string)
        GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
        if not(os.path.exists(GA_FINAL_MERGE_path)):
            print(dt.today(), "MetaPro's Taxonomic annotation relies the gene map created by GA.  please run MetaPro again with --tutorial GA")
        else:
            # ------------------------------------------------------
            # Taxonomic annotation
            TA_start = time.time()
            TA_path = os.path.join(output_folder_path, taxon_annotation_label)
            TA_jobs_folder = os.path.join(TA_path, "data", "jobs")
            #-----------------------------------------
            # stage 1
            marker_path_list = []
            #----------------------------------------------
            #centrifuge is too much of a RAM hog.  can't run more than 1 at a time
            sections = ["reads"]
            for section in sections:
                marker_file = "TA_centrifuge_" + section
                marker_path = os.path.join(TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                    launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
                    
            sections = ["contigs", "singletons"]
            if read_mode == "paired":
                sections.extend(["paired"])
                
            for section in sections:
                marker_file = "TA_kaiju_" + section
                marker_path = os.path.join(TA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_TA_kaiju_command(taxon_annotation_label, assemble_contigs_label, section, marker_file)
                    launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)        
            marker_file = "TA_taxon_pull"
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_taxon_pull_command(taxon_annotation_label, GA_final_merge_label, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
            print(dt.today(), "waiting for TA stage 1")
            for p_item in mp_store:
                p_item.join()
            mp_store[:] = []
            final_checklist = os.path.join(TA_path, "TA_stage_1.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            
            #--------------------------------------------------
            # stage 2
            marker_path_list = []
            sections = ["contigs"]
            for section in sections:
                marker_file = "TA_centrifuge_" + section
                marker_path = os.path.join(TA_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                    launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
            
            marker_file = "TA_kaiju_pp"
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_file)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_kaiju_pp_command(taxon_annotation_label, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
            for p_item in mp_store:
                p_item.join()
            mp_store[:] = []
            final_checklist = os.path.join(TA_path, "TA_stage_2.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            #------------------------------------------------------------------

            #-----------------------------------------------------------------
            # stage 3
            marker_path_list = []
            marker_file = "TA_centrifuge_pp"
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_centrifuge_pp_command(taxon_annotation_label, marker_file)
                launch_and_create_with_hold(mp_store, TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
            for p_item in mp_store:
                p_item.join()
            mp_store[:] = []
            final_checklist = os.path.join(TA_path, "TA_stage_3.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            #-----------------------------------------------
            # stage 4
            marker_path_list = []
            
            marker_file = "TA_final"
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_final_command(taxon_annotation_label, assemble_contigs_label, marker_file)
                launch_and_create_simple(taxon_annotation_label, marker_file, commands, command_list)
            final_checklist = os.path.join(TA_path, "TA_final.txt")
            check_all_job_markers(marker_path_list, final_checklist)
            

                

    elif(tutorial_mode_string == "EC"):
        print(dt.today(), "working on:", tutorial_mode_string)
        # Detect EC annotation
        GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
        if not(os.path.exists(GA_FINAL_MERGE_path)):
            print(dt.today(), "MetaPro Enzyme Annotation depends on the Gene Annotation phase.  Please run MetaPro with the option --tutorial GA")
        else:
            ec_annotation_path = os.path.join(output_folder_path, ec_annotation_label)
            EC_start = time.time()
            #There's a 2-step check.  We don't want it ti re-run either DETECT, or PRIAM+DIAMOND because they're too slow
            #if not check_where_resume(ec_annotation_path, None, GA_DIAMOND_path):
            #if check_bypass_log(output_folder_path, ec_annotation_label):
            EC_DETECT_start = time.time()
            ec_detect_path = os.path.join(ec_annotation_path, "data", "0_detect")
            #if not check_where_resume(job_label = None, full_path = ec_detect_path, dep_job_path = GA_DIAMOND_path):
            if check_bypass_log(output_folder_path, ec_annotation_detect_label):
                marker_file = "ec_detect"
                marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    command_list = commands.create_EC_DETECT_command(ec_annotation_label, GA_final_merge_label, marker_file)
                    launch_and_create_with_mp_store(mp_store, ec_annotation_label, marker_file, commands, command_list)
                
                
            EC_DETECT_end = time.time()
            print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
            
            # --------------------------------------------------------------
            # Priam EC annotation.  Why isn't it parallel? computing restraints.  Not enough mem
            EC_PRIAM_start = time.time()
            
            ec_priam_path = os.path.join(ec_annotation_path, "data", "1_priam")
            #if not check_where_resume(job_label = None, full_path = ec_priam_path, dep_job_path = GA_DIAMOND_path):
            if check_bypass_log(output_folder_path, ec_annotation_priam_label):
                marker_file = "ec_priam"
                marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    if(os.path.exists(ec_priam_path)):
                        print(dt.today(), "starting with a fresh PRIAM run")
                        shutil.rmtree(ec_priam_path)
                        
                    command_list = commands.create_EC_PRIAM_command(ec_annotation_label, GA_final_merge_label, marker_file)
                    launch_and_create_with_mp_store(mp_store, ec_annotation_label, marker_file, commands, command_list)
                
              
                #process.join()
            EC_PRIAM_end = time.time()
            print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
            # --------------------------------------------------------------
            # DIAMOND EC annotation 
            EC_DIAMOND_start = time.time()
            ec_diamond_path = os.path.join(ec_annotation_path, "data", "2_diamond")
            if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
                marker_file = "ec_diamond"
                marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    job_name = "ec_diamond"
                    command_list = commands.create_EC_DIAMOND_command(ec_annotation_label, GA_final_merge_label, marker_file)
                    launch_and_create_with_mp_store(mp_store, ec_annotation_label, job_name, commands, command_list)
                
            EC_DIAMOND_end = time.time()
            
            for item in mp_store:
                item.join()
            mp_store[:] = []
            
            ec_detect_out   = os.path.join(ec_annotation_path, "data", "jobs", "ec_detect")
            ec_priam_out    = os.path.join(ec_annotation_path, "data", "jobs", "ec_priam")
            ec_diamond_out  = os.path.join(ec_annotation_path, "data", "jobs", "ec_diamond")
            if check_bypass_log(output_folder_path, ec_annotation_detect_label):
                if(os.path.exists(ec_detect_out)):
                    write_to_bypass_log(output_folder_path, ec_annotation_detect_label)
            if check_bypass_log(output_folder_path, ec_annotation_priam_label):
                if(os.path.exists(ec_priam_out)):
                    write_to_bypass_log(output_folder_path, ec_annotation_priam_label)
            if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
                if(os.path.exists(ec_diamond_out)):
                    write_to_bypass_log(output_folder_path, ec_annotation_DIAMOND_label)
            
            #----------------------------------------------------------------------
            # EC post process
            EC_post_start = time.time()
            #if not (check_where_resume(ec_annotation_path, None, GA_DIAMOND_path)):
            if check_bypass_log(output_folder_path, ec_annotation_pp_label):
                
                marker_file = "ec_post"
                marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    command_list = commands.create_EC_postprocess_command(ec_annotation_label, GA_final_merge_label, marker_file)
                    launch_and_create_simple(ec_annotation_label, marker_file, commands, command_list)
                
                if(os.path.exists(marker_path)):
                    write_to_bypass_log(output_folder_path, ec_annotation_pp_label)

    
    elif(tutorial_mode_string == "output"):
        network_path = os.path.join(output_folder_path, output_label)
        #if not check_where_resume(network_path, None, ec_annotation_path):
        GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
        TA_path = os.path.join(output_folder_path, taxon_annotation_label)
        EC_path = os.path.join(output_folder_path, ec_annotation_label)
        all_good_flag = False
        if not(os.path.exists(GA_FINAL_MERGE_path)):
            print(dt.today(), "MetaPro's final outputs rely on the gene annotation.  please run MetaPro with --tutorial GA")
            sys.exit("missing GA")
        if not(os.path.exists(TA_path)):
            print(dt.today(), "MetaPro's final outputs rely on the Taxonomic annotation.  please run MetaPro with --tutorial TA")
            sys.exit("missing TA")
        if not(os.path.exists(EC_path)):
            print(dt.today(), "MetaPRo's final outputs rely on the Enzyme Annotation. Please run MetaPro with --tutorial EC")
            sys.exit("missing EC")
            
        all_good_flag = True    
        if(all_good_flag):
            #phase 1
            if check_bypass_log(output_folder, output_copy_gene_map_label):
                job_name = output_copy_gene_map_label
                command_list = commands.create_output_copy_gene_map_command(output_label, GA_final_merge_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                
            if check_bypass_log(output_folder, output_copy_taxa_label):
                job_name = output_copy_taxa_label
                command_list = commands.create_output_copy_taxa_command(output_label, taxon_annotation_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
            if check_bypass_log(output_folder, output_contig_stats_label):
                job_name = output_contig_stats_label
                command_list = commands.create_output_contig_stats_command(output_label, assemble_contigs_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
            #repop hosts
            if not(no_host):
                if check_bypass_log(output_folder, output_unique_hosts_singletons_label):
                    job_name = output_unique_hosts_singletons_label
                    command_list = commands.create_output_unique_hosts_singletons_command(output_label, quality_filter_label, host_filter_label)
                    launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                
                if(read_mode == "paired"):
                    if check_bypass_log(output_folder, output_unique_hosts_pair_1_label):
                        job_name = output_unique_hosts_pair_1_label
                        command_list = commands.create_output_unique_hosts_pair_1_command(output_label, quality_filter_label, host_filter_label)
                        launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                        
                    if check_bypass_log(output_folder, output_unique_hosts_pair_2_label):
                        job_name = output_unique_hosts_pair_2_label
                        command_list = commands.create_output_unique_hosts_pair_2_command(output_label, quality_filter_label, host_filter_label)
                        launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)

            #repop vectors
            if check_bypass_log(output_folder, output_unique_vectors_singletons_label):
                job_name = output_unique_vectors_singletons_label
                command_list = commands.create_output_unique_vectors_singletons_command(output_label, quality_filter_label, vector_filter_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
            
            if(read_mode == "paired"):
                if check_bypass_log(output_folder, output_unique_vectors_pair_1_label):
                    job_name = output_unique_vectors_pair_1_label
                    command_list = commands.create_output_unique_vectors_pair_1_command(output_label, quality_filter_label, vector_filter_label)
                    launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                    
                if check_bypass_log(output_folder, output_unique_vectors_pair_2_label):
                    job_name = output_unique_vectors_pair_2_label
                    command_list = commands.create_output_unique_vectors_pair_2_command(output_label, quality_filter_label, vector_filter_label)
                    launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                        
                        
            print(dt.today(), "output report phase 1 launched.  waiting for sync")
            for item in mp_store:
                item.join()
            mp_store[:] = []
            
            #----------------------------------------------------------------------------
            #Phase 2
            if check_bypass_log(output_folder, output_network_gen_label):
                command_list = commands.create_output_network_generation_command(output_label, GA_final_merge_label, taxon_annotation_label, ec_annotation_label)
                launch_and_create_with_mp_store(mp_store, output_label, output_network_gen_label, commands, command_list)
                
            if check_bypass_log(output_folder, output_taxa_groupby_label):
                command_list = commands.create_output_taxa_groupby_command(output_label)
                launch_and_create_with_mp_store(mp_store, output_label, output_taxa_groupby_label, commands, command_list)
           
            print(dt.today(), "output report phase 2 launched.  waiting for sync")
            for item in mp_store:
                item.join()
            mp_store[:] = []
         
            
            #-------------------------------------------------------------------
            #Phase 3
            if check_bypass_log(output_folder, output_read_count_label):
                job_name = output_read_count_label
                command_list = commands.create_output_read_count_command(output_label, quality_filter_label, repop_job_label, GA_final_merge_label, ec_annotation_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                                   

            if check_bypass_log(output_folder, output_per_read_scores_label):
                job_name = output_per_read_scores_label
                command_list = commands.create_output_per_read_scores_command(output_label, quality_filter_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)
                
            if check_bypass_log(output_folder, output_ec_heatmap_label):
                job_name = output_ec_heatmap_label
                command_list = commands.create_output_EC_heatmap_command(output_label)
                launch_and_create_with_mp_store(mp_store, output_label, job_name, commands, command_list)    
            
            print(dt.today(), "output report phase 3 launched.  waiting for sync")
            for item in mp_store:
                item.join()
            mp_store[:] = []
     
    else:
        print(dt.today(), "tutorial mode not recognized:", tutorial_mode_string)
    
    
    
    
    
if __name__ == "__main__":
    print("METAPRO metatranscriptomic analysis pipeline")
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline"
                                        "Version 1.0 © 2018")

    parser.add_argument("-c", "--config",   type=str,   help="Path to the configureation file")
    parser.add_argument("-1", "--pair1",    type=str,   help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2",    type=str,   help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single",   type=str,   help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-con", "--contig",   type=str,   help="Tutorial use only: Path to the file containing the contig reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, required=True, help="Path of the folder for the output of the pipeline")
    parser.add_argument("-t", "--num_threads", type=int, help="Maximum number of threads used by the pipeline")
    parser.add_argument("--nhost", action='store_true', help="Skip the host read removal step of the pipeline")
    parser.add_argument("--verbose_mode", type=str, help = "Decide how to handle the interim files, Compress them, or leave them alone.  Values are: keep, compress, quiet")
    parser.add_argument("--tutorial", type = str, help = "tutorial operating mode for MetaPro")
    
    args = parser.parse_args()

    if (args.pair1 and not args.pair2) or (args.pair2 and not args.pair1):
        print("You must specify both forward and reverse reads for a paired-end run")
        sys.exit()
    elif args.single and (args.pair1 or args.pair2):
        print("You cannot specify both paired-end and single-end reads in a single run.")
        sys.exit()

    config_file     = args.config if args.config else ""
    contig          = args.contig if args.contig else ""
    pair_1          = args.pair1 if args.pair1 else ""
    pair_2          = args.pair2 if args.pair2 else ""
    single          = args.single if args.single else ""
    output_folder   = args.output_folder
    num_threads     = args.num_threads if args.num_threads else 0
    no_host         = args.nhost if args.nhost else False
    verbose_mode    = args.verbose_mode if args.verbose_mode else "quiet"
    tutorial_mode   = args.tutorial if args.tutorial else "none"

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
