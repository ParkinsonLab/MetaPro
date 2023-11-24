# The functions here generate the pipeline commands.
# The functions here generate the pipeline commands.
# Each command module is made up of sub stages that are used to get the final result.

import os
import os.path
from re import split

from numpy.core.arrayprint import _make_options_dict
import MetaPro_paths as mpp
import subprocess as sp
import time
import sys
from datetime import datetime as dt

class mt_pipe_commands:
    # --------------------------------------------------------------------
    # constructor:
    # there should only be one of these objects used for an entire pipeline.
    def __init__(self, no_host, path_obj, Quality_score=33, tutorial_keyword = None, sequence_path_1=None, sequence_path_2=None, sequence_single=None, sequence_contigs = None):

        self.no_host_flag = no_host
        self.path_obj = path_obj
        # path to the genome sequence file

        
        if(tutorial_keyword is None):
            print("MetaPro operating in auto-mode")
            self.tutorial_keyword = None
            
            if sequence_single is not None:
                self.sequence_single = sequence_single
                self.sequence_path_1 = ""
                self.sequence_path_2 = ""
                print("Reads:", self.sequence_single)
                self.read_mode = "single"
                self.sequence_contigs = ""
            else:
                self.sequence_single = ""
                self.sequence_path_1 = sequence_path_1
                self.sequence_path_2 = sequence_path_2
                print("Forward Reads:", self.sequence_path_1)
                print("Reverse Reads:", self.sequence_path_2)
                self.read_mode = "paired"
                self.sequence_contigs = ""
            
        else:
            print("MetaPro is in TUTORIAL MODE:", tutorial_keyword)
            self.tutorial_keyword = tutorial_keyword
            if sequence_path_1 is None:
                self.sequence_single = sequence_single
                self.sequence_path_1 = ""
                self.sequence_path_2 = ""
                self.sequence_contigs = ""
                if(sequence_contigs is not None):
                    self.sequence_contigs = sequence_contigs
                print("Reads:", self.sequence_single)
                print("potential contigs:", self.sequence_contigs)
                self.read_mode = "single"
            else:
                self.sequence_single = ""
                if(sequence_single is not None):
                    self.sequence_single = sequence_single
                self.sequence_contigs = ""
                if(sequence_contigs is not None):
                    self.sequence_contigs = sequence_contigs
                    
                self.sequence_path_1 = sequence_path_1
                self.sequence_path_2 = sequence_path_2
                print("Forward Reads:", self.sequence_path_1)
                print("Reverse Reads:", self.sequence_path_2)
                print("potential singletons:", self.sequence_single)
                print("potential contigs:", self.sequence_contigs)
                self.read_mode = "paired"
            
                
        self.Qual_str = str(Quality_score)
        self.Output_Path = os.getcwd()
        self.threads_str = str(self.path_obj.num_threads)
        self.thread_count = self.path_obj.num_threads
        self.DNA_DB = self.path_obj.DNA_DB
        
        
        print("Output filepath:", self.Output_Path)

    # -----------------------------------------------------------
    # support functions
    def make_folder(self, folder_path):
        if not (os.path.exists(folder_path)):
            os.makedirs(folder_path)
            
    def conditional_insert_job_into_queue(self, job_location, queue, command):
        if(not os.path.exists(job_location)):
            print(dt.today(), job_location, "doesn't exist. adding")
            if(not queue):
                queue = [command]
            else:
                queue.append(command)
                
            print(queue)
        return queue

    def create_and_launch(self, job_folder, inner_name, command_list):
        # create the pbs job, and launch items
        # job name: string tag for export file name
        # command list:  list of command statements for writing
        # mode: selection of which pbs template to use: default -> low memory
        # dependency_list: if not empty, will append wait args to sbatch subprocess call. it's polymorphic
        # returns back the job ID given from sbatch

        # docker mode: single cpu
        # no ID, no sbatch.  just run the command
        
        shell_script_full_path = os.path.join(self.Output_Path, job_folder, inner_name + ".sh")

        with open(shell_script_full_path, "w") as PBS_script_out:
            for item in command_list:
                PBS_script_out.write(item + "\n")
            PBS_script_out.close()
        #if not work_in_background:
        output = ""
        try:
            sp.check_output(["sh", shell_script_full_path])#, stderr = sp.STDOUT)
        except sp.CalledProcessError as e:
            return_code = e.returncode
            if return_code != 1:
                raise
                
    def create_and_launch_v2(self, job_path, command_list):
        # create the pbs job, and launch items
        # job name: string tag for export file name
        # command list:  list of command statements for writing
        # mode: selection of which pbs template to use: default -> low memory
        # dependency_list: if not empty, will append wait args to sbatch subprocess call. it's polymorphic
        # returns back the job ID given from sbatch

        # docker mode: single cpu
        # no ID, no sbatch.  just run the command
        
        #shell_script_full_path = os.path.join(self.Output_Path, job_folder, inner_name + ".sh")

        with open(job_path, "w") as PBS_script_out:
            for item in command_list:
                PBS_script_out.write(item + "\n")
            PBS_script_out.close()
        #if not work_in_background:
        output = ""
        try:
            sp.check_output(["sh", job_path])#, stderr = sp.STDOUT)
        except sp.CalledProcessError as e:
            return_code = e.returncode
            if return_code != 1:
                raise                
                
    def launch_only(self, command_list, command_list_length):
        #just launch the job.  Don't make a script file.
        #print(dt.today(), "inside launch_only:", len(command_list))
        
        if(command_list_length == 1):
            #print("0th item:", command_list[0])
            try:
                os.system(command_list[0])
            except sp.CalledProcessError as e:
                return_code = e.returncode
                if return_code != 1:
                    raise
            #else:
            #    sys.exit("something bad happened")
        else:
        
            for command_item in command_list:
                try:
                    os.system(command_item)
                except sp.CalledProcessError as e:
                    return_code = e.returncode
                    if return_code != 1:
                        raise
                
                
    def create_quality_control_command(self):
        
        self.make_folder(self.path_obj.qc_top_path)
        self.make_folder(self.path_obj.qc_data_path)
        self.make_folder(self.path_obj.qc_sort_path)
        self.make_folder(self.path_obj.qc_adapter_path)
        self.make_folder(self.path_obj.qc_merge_path)
        self.make_folder(self.path_obj.qc_filter_path)
        self.make_folder(self.path_obj.qc_orphan_path)
        self.make_folder(self.path_obj.qc_cdhit_path)
        self.make_folder(self.path_obj.qc_final_path)

        
        sort_pair_1 = ">&2 echo Sorting pair 1 | "
        sort_pair_1 += self.path_obj.Python + " "
        sort_pair_1 += self.path_obj.sort_reads + " "
        sort_pair_1 += self.sequence_path_1 + " "
        sort_pair_1 += os.path.join(self.path_obj.qc_sort_path, "pair_1_sorted.fastq") + " "
        sort_pair_1 += "forward"

        sort_pair_2 = ">&2 echo Sorting pair 2 | "
        sort_pair_2 += self.path_obj.Python + " "
        sort_pair_2 += self.path_obj.sort_reads + " "
        sort_pair_2 += self.sequence_path_2 + " "
        sort_pair_2 += os.path.join(self.path_obj.qc_sort_path, "pair_2_sorted.fastq") + " "
        sort_pair_2 += "reverse"

        adapter_removal_line = ">&2 echo Removing adapters | "
        adapter_removal_line += self.path_obj.AdapterRemoval
        if self.read_mode == "single":
            adapter_removal_line += " --file1 " + self.sequence_single
        elif self.read_mode == "paired":
            adapter_removal_line += " --file1 " + os.path.join(self.path_obj.qc_sort_path, "pair_1_sorted.fastq")
            adapter_removal_line += " --file2 " + os.path.join(self.path_obj.qc_sort_path, "pair_2_sorted.fastq")
        adapter_removal_line += " --qualitybase " + str(self.Qual_str)
        if(self.Qual_str == "33"):
            adapter_removal_line += " --qualitymax 75"
        adapter_removal_line += " --threads " + self.threads_str
        adapter_removal_line += " --minlength " + str(self.path_obj.adapterremoval_minlength)
        adapter_removal_line += " --basename " + self.path_obj.qc_adapter_path
        adapter_removal_line += "_AdapterRemoval"
        adapter_removal_line += " --trimqualities "
        if self.read_mode == "single":
            adapter_removal_line += " --output1 " + os.path.join(self.path_obj.qc_adapter_path, "singletons_adptr_rem.fastq")
        elif self.read_mode == "paired":
            adapter_removal_line += " --output1 " + os.path.join(self.path_obj.qc_adapter_path, "pair_1_adptr_rem.fastq")
            adapter_removal_line += " --output2 " + os.path.join(self.path_obj.qc_adapter_path, "pair_2_adptr_rem.fastq")
            adapter_removal_line += " --singleton " + os.path.join(self.path_obj.qc_adapter_path, "singletons_adptr_rem.fastq")

        #Sort-reads introduces tags at the read-level of the 
        #tag_remove_pair_1 = ">&2 echo Remove tags pair 1 | "
        #tag_remove_pair_1 += self.path_obj.Python + " "
        #tag_remove_pair_1 += self.path_obj.remove_tag + " "
        #tag_remove_pair_1 += os.path.join(self.path_obj.qc_adapter_path, "pair_1_adptr_rem.fastq") + " "
        #tag_remove_pair_1 += os.path.join(self.path_obj.qc_tag_path, "pair_1_no_tags.fastq")
        
        #tag_remove_pair_2 = ">&2 echo Remove tags pair 2 | "
        #tag_remove_pair_2 += self.path_obj.Python + " "
        #tag_remove_pair_2 += self.path_obj.remove_tag + " "
        #tag_remove_pair_2 += os.path.join(self.path_obj.qc_adapter_path, "pair_2_adptr_rem.fastq") + " "
        #tag_remove_pair_2 += os.path.join(self.path_obj.qc_tag_path, "pair_2_no_tags.fastq")

        #tag_remove_singletons =  ">&2 echo Remove tags singletons | " 
        #tag_remove_singletons += self.path_obj.Python + " "
        #tag_remove_singletons += self.path_obj.remove_tag + " "
        #tag_remove_singletons += os.path.join(self.path_obj.qc_adapter_path, "singletons_adptr_rem.fastq") + " "
        #tag_remove_singletons += os.path.join(self.path_obj.qc_tag_path, "singletons_no_tags.fastq")
        # tries to merge the cleaned pairs
        # rejects get sent out
        vsearch_merge = ">&2 echo " + "Vsearch Merge pairs | "
        vsearch_merge += self.path_obj.vsearch
        vsearch_merge += " --fastq_mergepairs " + os.path.join(self.path_obj.qc_adapter_path, "pair_1_adptr_rem.fastq")
        vsearch_merge += " --reverse " + os.path.join(self.path_obj.qc_adapter_path, "pair_2_adptr_rem.fastq")
        vsearch_merge += " --fastq_ascii " + str(self.Qual_str)
        vsearch_merge += " --fastqout " + os.path.join(self.path_obj.qc_merge_path, "merge_success.fastq")
        vsearch_merge += " --fastqout_notmerged_fwd " + os.path.join(self.path_obj.qc_merge_path, "pair_1_merge_reject.fastq")
        vsearch_merge += " --fastqout_notmerged_rev " + os.path.join(self.path_obj.qc_merge_path, "pair_2_merge_reject.fastq")

        # concatenate the merge overlaps with the singletons
        cat_glue = ">&2 echo concatenating singletons | "
        cat_glue += "cat "
        cat_glue += os.path.join(self.path_obj.qc_merge_path, "merge_success.fastq") + " "
        cat_glue += os.path.join(self.path_obj.qc_adapter_path, "singletons_adptr_rem.fastq")
        cat_glue += " > " + os.path.join(self.path_obj.qc_merge_path, "singletons.fastq")

        # Filter out low-quality reads
        # start with the singles / merged sections
        
        vsearch_filter_0 = ">&2 echo low-quality filter on singletons | "
        vsearch_filter_0 += self.path_obj.vsearch
        if self.read_mode == "single":
            vsearch_filter_0 += " --fastq_filter " + os.path.join(self.path_obj.qc_adapter_path, "singletons_adptr_rem.fastq")
        elif self.read_mode == "paired":
            vsearch_filter_0 += " --fastq_filter " + os.path.join(self.path_obj.qc_merge_path, "singletons.fastq")
        vsearch_filter_0 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_0 += " --fastq_maxee " + "2.0"
        vsearch_filter_0 += " --fastqout " + os.path.join(self.path_obj.qc_filter_path, "singletons_hq.fastq")

        # then move onto the standalones in pair 1
        vsearch_filter_1 = ">&2 echo low-quality filter on pair 1 | "
        vsearch_filter_1 += self.path_obj.vsearch
        vsearch_filter_1 += " --fastq_filter " + os.path.join(self.path_obj.qc_merge_path, "pair_1_merge_reject.fastq")
        vsearch_filter_1 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_1 += " --fastq_maxee " + "2.0"
        vsearch_filter_1 += " --fastqout " + os.path.join(self.path_obj.qc_filter_path, "pair_1_hq.fastq")

        vsearch_filter_2 = ">&2 echo low-quality filter on pair 2 | "
        vsearch_filter_2 += self.path_obj.vsearch
        vsearch_filter_2 += " --fastq_filter " + os.path.join(self.path_obj.qc_merge_path, "pair_2_merge_reject.fastq")
        vsearch_filter_2 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_2 += " --fastq_maxee " + "2.0"
        vsearch_filter_2 += " --fastqout " + os.path.join(self.path_obj.qc_filter_path, "pair_2_hq.fastq")

        # redistribute data into singletons, or paired-reads
        orphan_read_filter = ">&2 echo moving newly orphaned reads | "
        orphan_read_filter += self.path_obj.Python + " "
        orphan_read_filter += self.path_obj.orphaned_read_filter + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_filter_path, "pair_1_hq.fastq") + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_filter_path, "pair_2_hq.fastq") + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_filter_path, "singletons_hq.fastq") + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_orphan_path, "pair_1_match.fastq") + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_orphan_path, "pair_2_match.fastq") + " "
        orphan_read_filter += os.path.join(self.path_obj.qc_orphan_path, "singletons_with_duplicates.fastq")

        # remove duplicates (to shrink the data size)
        cdhit_singletons = ">&2 echo removing singleton duplicates | "
        cdhit_singletons += self.path_obj.cdhit_dup + " -i "
        if self.read_mode == "single":
            cdhit_singletons += os.path.join(self.path_obj.qc_filter_path, "singletons_hq.fastq")
        elif self.read_mode == "paired":
            cdhit_singletons += os.path.join(self.path_obj.qc_orphan_path, "singletons_with_duplicates.fastq")
        cdhit_singletons += " -o " + os.path.join(self.path_obj.qc_cdhit_path, "singletons_unique.fastq")

        # remove duplicates in the pairs
        cdhit_paired = ">&2 echo remove duplicates from paired | "
        cdhit_paired += self.path_obj.cdhit_dup + " "
        cdhit_paired += "-i"    + " " + os.path.join(self.path_obj.qc_orphan_path, "pair_1_match.fastq") + " "
        cdhit_paired += "-i2"   + " " + os.path.join(self.path_obj.qc_orphan_path, "pair_2_match.fastq") + " "
        cdhit_paired += "-o"    + " " + os.path.join(self.path_obj.qc_cdhit_path, "pair_1_unique.fastq") + " "
        cdhit_paired += "-o2"   + " " + os.path.join(self.path_obj.qc_cdhit_path, "pair_2_unique.fastq")

        #move data to appropriate places
        mv_singletons = "mv " + os.path.join(self.path_obj.qc_cdhit_path, "singletons_unique.fastq") + " "
        mv_singletons += os.path.join(self.path_obj.qc_final_path, "singletons.fastq")

        mv_pair_1 = "mv " + os.path.join(self.path_obj.qc_cdhit_path, "pair_1_unique.fastq") + " "
        mv_pair_1 += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq")

        mv_pair_2 = "mv " + os.path.join(self.path_obj.qc_cdhit_path, "pair_2_unique.fastq") + " "
        mv_pair_2 += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq")
        
        # move these particular files to final_folder because they'll be needed by another stage.
        mv_duplicate_singletons = "mv "
        if(self.read_mode == "single"):
            mv_duplicate_singletons += os.path.join(self.path_obj.qc_filter_path, "singletons_hq.fastq") + " "
            mv_duplicate_singletons += os.path.join(self.path_obj.qc_final_path, "singletons_hq.fastq")
        else:
            mv_duplicate_singletons += os.path.join(self.path_obj.qc_orphan_path, "singletons_with_duplicates.fastq") + " "
            mv_duplicate_singletons += os.path.join(self.path_obj.qc_final_path, "singletons_with_duplicates.fastq")

        mv_pair_1_match = "mv " + os.path.join(self.path_obj.qc_orphan_path, "pair_1_match.fastq") + " "
        mv_pair_1_match += os.path.join(self.path_obj.qc_final_path, "pair_1_match.fastq")

        mv_pair_2_match = "mv " + os.path.join(self.path_obj.qc_orphan_path, "pair_2_match.fastq") + " "
        mv_pair_2_match += os.path.join(self.path_obj.qc_final_path, "pair_2_match.fastq")

        mv_singletons_cluster = "mv " + os.path.join(self.path_obj.qc_cdhit_path, "singletons_unique.fastq.clstr") + " "
        mv_singletons_cluster += os.path.join(self.path_obj.qc_final_path, "singletons_unique.fastq.clstr")

        mv_paired_cluster = "mv " + os.path.join(self.path_obj.qc_cdhit_path, "pair_1_unique.fastq.clstr") + " "
        mv_paired_cluster += os.path.join(self.path_obj.qc_final_path, "pair_1_unique.fastq.clstr")

        if self.read_mode == "single":
            COMMANDS_qual = [
                adapter_removal_line,
                vsearch_filter_0,
                cdhit_singletons,
                mv_singletons,
                mv_duplicate_singletons,
                mv_singletons_cluster
            ]
        elif self.read_mode == "paired":
            COMMANDS_qual = [
                sort_pair_1,
                sort_pair_2,
                adapter_removal_line,
                vsearch_merge,
                cat_glue,
                vsearch_filter_0,
                vsearch_filter_1,
                vsearch_filter_2,
                orphan_read_filter,
                cdhit_singletons,
                cdhit_paired,
                mv_singletons,
                mv_pair_1,
                mv_pair_2,
                mv_duplicate_singletons,
                mv_singletons_cluster,
                mv_pair_1_match,
                mv_paired_cluster,
                mv_pair_2_match
            ]

        return COMMANDS_qual

    def create_host_index_command(self, marker_path):
        self.make_folder(self.path_obj.host_top_path)
        self.make_folder(self.path_obj.host_data_path)
        self.make_folder(self.path_obj.host_bwa_path)
        self.make_folder(self.path_obj.host_blat_path)
        self.make_folder(self.path_obj.host_final_path)
        #indexing takes long.  It's broken up so that we have control over it. 

        # craft a BWA index for the host sequences
        bwa_hr_prep = ">&2 echo make host contaminants index for BWA | "
        bwa_hr_prep += self.path_obj.BWA + " index -a bwtsw " + self.path_obj.Host_DB

        samtools_hr_prep = ">&2 echo SAMTOOLS host contaminant prep | "
        samtools_hr_prep += self.path_obj.SAMTOOLS + " faidx " + self.path_obj.Host_DB

        make_marker = "touch " + marker_path
        command = [
                    bwa_hr_prep,
                    samtools_hr_prep + " && " + make_marker
        ]

        return command


    def create_host_filter_command(self, marker_path):
        self.make_folder(self.path_obj.host_top_path)
        self.make_folder(self.path_obj.host_data_path)
        self.make_folder(self.path_obj.host_jobs_path)
        self.make_folder(self.path_obj.host_bwa_path)
        self.make_folder(self.path_obj.host_blat_path)
        self.make_folder(self.path_obj.host_final_path)

        #copy_host = ">&2 echo Copy the host file over | "
        #copy_host += "cp " + self.path_obj.Host + " " + self.path_obj.Host_DB

        

        # host removal on unique singletons
        bwa_hr_singletons = ">&2 echo BWA host remove on singletons | "
        bwa_hr_singletons += self.path_obj.BWA + " mem -t "
        bwa_hr_singletons += self.threads_str + " "
        bwa_hr_singletons += self.path_obj.Host_DB + " "
        bwa_hr_singletons += os.path.join(self.path_obj.qc_final_path, "singletons.fastq") + " " 
        bwa_hr_singletons += ">" + " "
        bwa_hr_singletons += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.sam")
        
        #Tutorial-use only.  
        bwa_hr_tut_singletons = ">&2 echo BWA host remove on singletons | "
        bwa_hr_tut_singletons += self.path_obj.BWA + " mem -t "
        bwa_hr_tut_singletons += self.threads_str + " "
        bwa_hr_tut_singletons += self.path_obj.Host_DB + " "
        bwa_hr_tut_singletons += self.sequence_single 
        bwa_hr_tut_singletons += " > " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.sam")

        # annoying type conversion pt 1
        samtools_hr_singletons_sam_to_bam = ">&2 echo convert singletons host reads | "
        samtools_hr_singletons_sam_to_bam += self.path_obj.SAMTOOLS
        samtools_hr_singletons_sam_to_bam += " view -bS " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.sam")
        samtools_hr_singletons_sam_to_bam += " > " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.bam")
        # annoying type conversion pt 2
        samtools_no_host_singletons_bam_to_fastq = self.path_obj.SAMTOOLS + " fastq -n -f 4" + " -0 "
        samtools_no_host_singletons_bam_to_fastq += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.fastq") + " "
        samtools_no_host_singletons_bam_to_fastq += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.bam")

        # apparently, we're to keep the host separation
        samtools_host_singletons_bam_to_fastq = self.path_obj.SAMTOOLS + " fastq -n -F 4" + " -0 "
        samtools_host_singletons_bam_to_fastq += os.path.join(self.path_obj.host_bwa_path, "singletons_host_only.fastq") + " "
        samtools_host_singletons_bam_to_fastq += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.bam")


        
        bwa_hr_paired = ">&2 echo bwa host-removal on paired | " 
        bwa_hr_paired += self.path_obj.BWA + " "
        bwa_hr_paired += "mem" + " "  + "-t" + " " + self.threads_str + " "
        bwa_hr_paired += self.path_obj.Host_DB + " "
        bwa_hr_paired += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
        bwa_hr_paired += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
        bwa_hr_paired += ">" + " "
        bwa_hr_paired += os.path.join(self.path_obj.host_bwa_path, "paired_on_host.sam")
        
        #Tutorial-use only
        bwa_hr_tut_paired = ">&2 echo bwa host-removal on paired | " 
        bwa_hr_tut_paired += self.path_obj.BWA + " "
        bwa_hr_tut_paired += "mem" + " "  + "-t" + " " + self.threads_str + " "
        bwa_hr_tut_paired += self.path_obj.Host_DB + " "
        bwa_hr_tut_paired += self.sequence_path_1 + " "
        bwa_hr_tut_paired += self.sequence_path_2 + " "
        bwa_hr_tut_paired += ">" + " "
        bwa_hr_tut_paired += os.path.join(self.path_obj.host_bwa_path, "paired_on_host.sam")
        
        
        bwa_hr_filter_paired = ">&2 echo BWA host-removal PP on paired | "
        bwa_hr_filter_paired += self.path_obj.Python + " "
        bwa_hr_filter_paired += self.path_obj.bwa_read_sorter + " "
        bwa_hr_filter_paired += "paired" + " "
        bwa_hr_filter_paired += self.path_obj.filter_stringency + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.host_bwa_path, "paired_on_host.sam") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fastq") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fastq") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_host_only.fastq") + " "
        bwa_hr_filter_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_host_only.fastq")
        
        bwa_hr_filter_tut_paired = ">&2 echo BWA host-removal PP on paired | "
        bwa_hr_filter_tut_paired += self.path_obj.Python + " "
        bwa_hr_filter_tut_paired += self.path_obj.bwa_read_sorter + " "
        bwa_hr_filter_tut_paired += "paired" + " "
        bwa_hr_filter_tut_paired += self.path_obj.filter_stringency + " "
        bwa_hr_filter_tut_paired += os.path.join(self.path_obj.host_bwa_path, "paired_on_host.sam") + " "
        bwa_hr_filter_tut_paired += self.sequence_path_1 + " "
        bwa_hr_filter_tut_paired += self.sequence_path_2 + " "
        bwa_hr_filter_tut_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fastq") + " "
        bwa_hr_filter_tut_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fastq") + " "
        bwa_hr_filter_tut_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_host_only.fastq") + " "
        bwa_hr_filter_tut_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_host_only.fastq")

        # blat prep
        make_blast_db_host = ">&2 echo Make BLAST db for host contaminants | "
        make_blast_db_host += self.path_obj.Makeblastdb + " -in " + self.path_obj.Host_DB + " -dbtype nucl"

        vsearch_filter_s = ">&2 echo Convert singletons for BLAT | "
        vsearch_filter_s += self.path_obj.vsearch
        vsearch_filter_s += " --fastq_filter " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.fastq")
        vsearch_filter_s += " --fastq_ascii " + self.Qual_str
        vsearch_filter_s += " --fastaout " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.fasta")

        vsearch_filter_p1 = ">&2 echo Convert pair 1 for BLAT | "
        vsearch_filter_p1 += self.path_obj.vsearch
        vsearch_filter_p1 += " --fastq_filter " + os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fastq")
        vsearch_filter_p1 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_p1 += " --fastaout " + os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fasta")

        vsearch_filter_p2 = ">&2 echo Convert pair 2 for BLAT | "
        vsearch_filter_p2 += self.path_obj.vsearch
        vsearch_filter_p2 += " --fastq_filter " + os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fastq")
        vsearch_filter_p2 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_p2 += " --fastaout " + os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fasta")

        blat_hr_singletons = ">&2 echo BLAT host singletons | "
        blat_hr_singletons += self.path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_hr_singletons += self.path_obj.Host_BLAT_DB + " "
        blat_hr_singletons += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.fasta")
        blat_hr_singletons += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str
        blat_hr_singletons += " " + os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.blatout")

        blat_hr_pair_1 = ">&2 echo BLAT host pair 1 | "
        blat_hr_pair_1 += self.path_obj.BLAT
        blat_hr_pair_1 += " -noHead -minIdentity=90 -minScore=65 " + self.path_obj.Host_BLAT_DB + " "
        blat_hr_pair_1 += os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fasta")
        blat_hr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str
        blat_hr_pair_1 += " " + os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.blatout")

        blat_hr_pair_2 = ">&2 echo BLAT host pair 2 | "
        blat_hr_pair_2 += self.path_obj.BLAT
        blat_hr_pair_2 += " -noHead -minIdentity=90 -minScore=65 " + self.path_obj.Host_BLAT_DB + " "
        blat_hr_pair_2 += os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fasta")
        blat_hr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str
        blat_hr_pair_2 += " " + os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.blatout")

        # HR BLAT
        hr_singletons = ">&2 echo BLAT contaminant singletons | "
        hr_singletons += self.path_obj.Python + " " + self.path_obj.BLAT_Contaminant_Filter + " "
        hr_singletons += "single" +  " "
        hr_singletons += self.path_obj.filter_stringency + " "
        hr_singletons += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.fastq") + " "  # in
        hr_singletons += os.path.join(self.path_obj.host_bwa_path, "singletons_no_host.blatout") + " "  # in
        hr_singletons += os.path.join(self.path_obj.host_blat_path, "singletons_no_host.fastq") + " "  # out
        hr_singletons += os.path.join(self.path_obj.host_blat_path, "singletons_host_only.fastq")  # out

        hr_paired = ">&2 echo BLAT contaminant paired | "
        hr_paired += self.path_obj.Python + " "
        hr_paired += self.path_obj.BLAT_Contaminant_Filter + " "
        hr_paired += "paired" + " "
        hr_paired += self.path_obj.filter_stringency + " "
        hr_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.fastq") + " "
        hr_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.fastq") + " "
        hr_paired += os.path.join(self.path_obj.host_bwa_path, "pair_1_no_host.blatout") + " "
        hr_paired += os.path.join(self.path_obj.host_bwa_path, "pair_2_no_host.blatout") + " "
        hr_paired += os.path.join(self.path_obj.host_blat_path, "pair_1_no_host.fastq") + " "
        hr_paired += os.path.join(self.path_obj.host_blat_path, "pair_2_no_host.fastq") + " "
        hr_paired += os.path.join(self.path_obj.host_blat_path, "pair_1_host_only.fastq") + " "
        hr_paired += os.path.join(self.path_obj.host_blat_path, "pair_2_host_only.fastq")
        
        

        
        #-----------------------------
        #orphan correction
        
        

        mv_singletons = "mv " + os.path.join(self.path_obj.host_blat_path, "singletons_no_host.fastq") + " "
        mv_singletons += os.path.join(self.path_obj.host_final_path, "singletons.fastq")

        mv_pair_1 = "mv " + os.path.join(self.path_obj.host_blat_path, "pair_1_no_host.fastq") + " "
        mv_pair_1 += os.path.join(self.path_obj.host_final_path, "pair_1.fastq")

        mv_pair_2 = "mv " + os.path.join(self.path_obj.host_blat_path, "pair_2_no_host.fastq") + " "
        mv_pair_2 += os.path.join(self.path_obj.host_final_path, "pair_2.fastq")
        
        #----------------------------------
        #final check against false-positive bypass-log writes
        make_marker = "touch" + " "
        make_marker += marker_path
        
        if(self.tutorial_keyword is None):
            if self.read_mode == "single":
                COMMANDS_host = [
                    bwa_hr_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    make_blast_db_host,
                    vsearch_filter_s,
                    blat_hr_singletons,
                    hr_singletons,
                    mv_singletons + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_host = [
                    bwa_hr_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    bwa_hr_paired,
                    bwa_hr_filter_paired,
                    make_blast_db_host,
                    vsearch_filter_s,
                    vsearch_filter_p1,
                    vsearch_filter_p2,
                    blat_hr_singletons,
                    blat_hr_pair_1,
                    blat_hr_pair_2,
                    hr_singletons,
                    hr_paired,
                    mv_singletons,
                    mv_pair_1,
                    mv_pair_2 + " && " + make_marker
                ]
        else:
            if self.read_mode == "single":
                print(dt.today(), "Host filter operating in tutorial-mode: SINGLE")
                COMMANDS_host = [
                    
                    bwa_hr_tut_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    make_blast_db_host,
                    vsearch_filter_s,
                    blat_hr_singletons,
                    hr_singletons,
                    mv_singletons + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                print(dt.today(), "Host filter operating in tutorial-mode: PAIRED")
                if(self.sequence_single == ""):
                    COMMANDS_host = [
                        
                        bwa_hr_tut_paired,
                        bwa_hr_filter_tut_paired,
                        make_blast_db_host,
                        vsearch_filter_p1,
                        vsearch_filter_p2,
                        blat_hr_pair_1,
                        blat_hr_pair_2,
                        hr_paired,
                        mv_pair_1,
                        mv_pair_2 + " && " + make_marker
                    ]

                else:
                    COMMANDS_host = [
                        bwa_hr_tut_singletons,
                        samtools_hr_singletons_sam_to_bam,
                        samtools_no_host_singletons_bam_to_fastq,
                        samtools_host_singletons_bam_to_fastq,
                        bwa_hr_tut_paired,
                        bwa_hr_filter_paired,
                        make_blast_db_host,
                        vsearch_filter_s,
                        vsearch_filter_p1,
                        vsearch_filter_p2,
                        blat_hr_singletons,
                        blat_hr_pair_1,
                        blat_hr_pair_2,
                        hr_singletons,
                        hr_paired,
                        mv_singletons,
                        mv_pair_1,
                        mv_pair_2 + " && " + make_marker
                    ]
                
            

                
        return COMMANDS_host

    def create_vector_filter_index_command(self, marker_path):
        
        #univec core is tiny. we can spend time indexing. 
        bwa_vr_prep = ">&2 echo BWA vector prep | "
        bwa_vr_prep += self.path_obj.BWA + " index -a bwtsw " + self.path_obj.Vector_DB
        samtools_vr_prep = ">&2 echo samtools vector prep | "
        samtools_vr_prep += self.path_obj.SAMTOOLS + " faidx " + self.path_obj.Vector_DB

        make_marker = "touch" + " "
        make_marker += marker_path
        command = [
                    bwa_vr_prep, 
                    samtools_vr_prep +  " && " + make_marker
        ]
        return command




    def create_vector_filter_command(self, marker_path, no_host_flag = False):

        self.make_folder(self.path_obj.vector_top_path)
        self.make_folder(self.path_obj.vector_jobs_path)
        self.make_folder(self.path_obj.vector_data_path)
        self.make_folder(self.path_obj.vector_bwa_path)
        self.make_folder(self.path_obj.vector_blat_path)
        self.make_folder(self.path_obj.vector_final_path)

        #Vector_Contaminants = os.path.join(self.path_obj.vector_bwa_path, "vector_contaminants_seq.fasta")
        Vector_Contaminants = self.path_obj.Vector_DB
        #copy_vector = ">&2 echo copy vector prep | "
        #copy_vector += "cp " + self.path_obj.UniVec_Core + " " + Vector_Contaminants

        

        

        bwa_vr_singletons = ">&2 echo BWA vector oprhans | "
        bwa_vr_singletons += self.path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_vr_singletons += Vector_Contaminants + " "
        bwa_vr_singletons += os.path.join(self.path_obj.qc_final_path, "singletons.fastq")
        bwa_vr_singletons += " > " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.sam")
        
        
        bwa_vr_tut_singletons = ">&2 echo BWA vector oprhans TUTORIAL MODE | "
        bwa_vr_tut_singletons += self.path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_vr_tut_singletons += Vector_Contaminants + " "
        bwa_vr_tut_singletons += self.sequence_single
        bwa_vr_tut_singletons += " > " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.sam")

        samtools_no_vector_singletons_convert = ">&2 echo samtools vector oprhans pt 1 | "
        samtools_no_vector_singletons_convert += self.path_obj.SAMTOOLS + " view -bS "
        samtools_no_vector_singletons_convert += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.sam")
        samtools_no_vector_singletons_convert += " > " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.bam")

        samtools_no_vector_singletons_export = ">&2 echo samtools vector singletons pt 2 | "
        samtools_no_vector_singletons_export += self.path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_vector_singletons_export += " -0 " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.fastq") + " "
        samtools_no_vector_singletons_export += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.bam")

        samtools_vector_singletons_export = ">&2 echo samtools vector singletons pt 3 | "
        samtools_vector_singletons_export += self.path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_vector_singletons_export += " -0 " + os.path.join(self.path_obj.vector_bwa_path, "singletons_vectors_only.fastq") + " "
        samtools_vector_singletons_export += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.bam")

        bwa_vr_paired = ">&2 echo bwa vector paired | "
        bwa_vr_paired += self.path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_vr_paired += Vector_Contaminants + " "
        if(no_host_flag):

            bwa_vr_paired += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
            bwa_vr_paired += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
        else:
            bwa_vr_paired += os.path.join(self.path_obj.host_final_path, "pair_1.fastq") + " "
            bwa_vr_paired += os.path.join(self.path_obj.host_final_path, "pair_2.fastq") + " "
            
        
        bwa_vr_paired += " > " + os.path.join(self.path_obj.vector_bwa_path, "paired_on_vectors.sam")

        bwa_vr_tut_paired = ">&2 echo bwa vector paired TUTORIAL MODE | "
        bwa_vr_tut_paired += self.path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_vr_tut_paired += Vector_Contaminants + " "
        bwa_vr_tut_paired += self.sequence_path_1 + " "
        bwa_vr_tut_paired += self.sequence_path_2 + " "
        bwa_vr_tut_paired += " > " + os.path.join(self.path_obj.vector_bwa_path, "paired_on_vectors.sam")
        
        bwa_vr_filter_paired = ">&2 echo BWA vector filter on paired | "
        bwa_vr_filter_paired += self.path_obj.Python + " "
        bwa_vr_filter_paired += self.path_obj.bwa_read_sorter + " "
        bwa_vr_filter_paired += "paired" + " "
        bwa_vr_filter_paired += self.path_obj.filter_stringency + " "
        bwa_vr_filter_paired += os.path.join(self.path_obj.vector_bwa_path, "paired_on_vectors.sam") + " "
        if(no_host_flag):
            bwa_vr_filter_paired += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
            bwa_vr_filter_paired += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
        else:
            bwa_vr_filter_paired += os.path.join(self.path_obj.host_final_path, "pair_1.fastq") + " "
            bwa_vr_filter_paired += os.path.join(self.path_obj.host_final_path, "pair_2.fastq") + " "

        bwa_vr_filter_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.fastq") + " "
        bwa_vr_filter_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.fastq") + " "
        bwa_vr_filter_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_1_vectors_only.fastq") + " "
        bwa_vr_filter_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_2_vectors_only.fastq")


        make_blast_db_vector = ">&2 echo BLAST make db vectors | "
        make_blast_db_vector += self.path_obj.Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl"

        vsearch_filter_6 = ">&2 echo convert vector singletons for BLAT | "
        vsearch_filter_6 += self.path_obj.vsearch
        vsearch_filter_6 += " --fastq_filter " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.fastq")
        vsearch_filter_6 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_6 += " --fastaout " + os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.fasta")

        vsearch_filter_7 = ">&2 echo convert vector pair 1 for BLAT | "
        vsearch_filter_7 += self.path_obj.vsearch
        vsearch_filter_7 += " --fastq_filter " + os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.fastq")
        vsearch_filter_7 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_7 += " --fastaout " + os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.fasta")

        vsearch_filter_8 = ">&2 echo convert vector pair 2 for BLAT | "
        vsearch_filter_8 += self.path_obj.vsearch
        vsearch_filter_8 += " --fastq_filter " + os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.fastq")
        vsearch_filter_8 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_8 += " --fastaout " + os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.fasta")

        blat_vr_singletons = ">&2 echo BLAT vector singletons | "
        blat_vr_singletons += self.path_obj.BLAT
        blat_vr_singletons += " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_singletons += Vector_Contaminants + " "
        blat_vr_singletons += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.fasta")
        blat_vr_singletons += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str + " "
        blat_vr_singletons += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.blatout")

        blat_vr_pair_1 = ">&2 echo BLAT vector pair 1 | "
        blat_vr_pair_1 += self.path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_1 += Vector_Contaminants + " "
        blat_vr_pair_1 += os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.fasta")
        blat_vr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str + " "
        blat_vr_pair_1 += os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.blatout")

        blat_vr_pair_2 = ">&2 echo BLAT vector pair 2 | "
        blat_vr_pair_2 += self.path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_2 += Vector_Contaminants + " "
        blat_vr_pair_2 += os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.fasta")
        blat_vr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.threads_str + " "
        blat_vr_pair_2 += os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.blatout")

        blat_filter_vector_singletons = ">&2 echo BLAT contaminant singletons | "
        blat_filter_vector_singletons += self.path_obj.Python + " " + self.path_obj.BLAT_Contaminant_Filter + " "
        blat_filter_vector_singletons += "single" + " "
        blat_filter_vector_singletons += self.path_obj.filter_stringency + " "
        blat_filter_vector_singletons += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.fastq") + " "  # in
        blat_filter_vector_singletons += os.path.join(self.path_obj.vector_bwa_path, "singletons_no_vectors.blatout") + " "  # in
        blat_filter_vector_singletons += os.path.join(self.path_obj.vector_blat_path, "singletons_no_vectors.fastq") + " "  # out
        blat_filter_vector_singletons += os.path.join(self.path_obj.vector_blat_path, "singletons_vectors_only.fastq")  # out

        blat_filter_vector_paired = ">&2 echo BLAT contaminant pair 1 | "
        blat_filter_vector_paired += self.path_obj.Python + " " + self.path_obj.BLAT_Contaminant_Filter + " "
        blat_filter_vector_paired += "paired" + " "
        blat_filter_vector_paired += self.path_obj.filter_stringency + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_1_no_vectors.blatout") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_bwa_path, "pair_2_no_vectors.blatout") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_blat_path, "pair_1_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_blat_path, "pair_2_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_blat_path, "pair_1_vectors_only.fastq") + " "
        blat_filter_vector_paired += os.path.join(self.path_obj.vector_blat_path, "pair_2_vectors_only.fastq")

        mv_singletons = "cp " + os.path.join(self.path_obj.vector_blat_path, "singletons_no_vectors.fastq") + " "
        mv_singletons += os.path.join(self.path_obj.vector_final_path, "singletons.fastq")

        mv_pair_1 = "cp " + os.path.join(self.path_obj.vector_blat_path, "pair_1_no_vectors.fastq") + " "
        mv_pair_1 += os.path.join(self.path_obj.vector_final_path, "pair_1.fastq")

        mv_pair_2 = "cp " + os.path.join(self.path_obj.vector_blat_path, "pair_2_no_vectors.fastq") + " "
        mv_pair_2 += os.path.join(self.path_obj.vector_final_path, "pair_2.fastq")

        make_marker = "touch" + " "
        make_marker += marker_path
        
        if(self.tutorial_keyword == "vectors" or self.tutorial_keyword == "vector"):
            if self.read_mode == "single":
                COMMANDS_vector = [
                    
                    #bwa_vr_prep,
                    #samtools_vr_prep,
                    bwa_vr_tut_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    make_blast_db_vector,
                    vsearch_filter_6,
                    blat_vr_singletons,
                    blat_filter_vector_singletons,
                    mv_singletons + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    #bwa_vr_prep,
                    #samtools_vr_prep,
                    bwa_vr_tut_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    bwa_vr_tut_paired,
                    bwa_vr_filter_paired, 
                    make_blast_db_vector,
                    vsearch_filter_6,
                    vsearch_filter_7,
                    vsearch_filter_8,
                    blat_vr_singletons,
                    blat_vr_pair_1,
                    blat_vr_pair_2,
                    blat_filter_vector_singletons,
                    blat_filter_vector_paired,
                    mv_singletons,
                    mv_pair_1,
                    mv_pair_2 + " && " + make_marker
                ]
        else:    
            if self.read_mode == "single":
                COMMANDS_vector = [
                    #bwa_vr_prep,
                    #samtools_vr_prep,
                    bwa_vr_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    make_blast_db_vector,
                    vsearch_filter_6,
                    blat_vr_singletons,
                    blat_filter_vector_singletons,
                    mv_singletons + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    #bwa_vr_prep,
                    #samtools_vr_prep,
                    bwa_vr_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    bwa_vr_paired,
                    bwa_vr_filter_paired, 
                    make_blast_db_vector,
                    vsearch_filter_6,
                    vsearch_filter_7,
                    vsearch_filter_8,
                    blat_vr_singletons,
                    blat_vr_pair_1,
                    blat_vr_pair_2,
                    blat_filter_vector_singletons,
                    blat_filter_vector_paired,
                    mv_singletons,
                    mv_pair_1,
                    mv_pair_2 + " && " + make_marker
                ]    

        return COMMANDS_vector
        
    def create_rRNA_filter_split_command(self, marker_path):
        #using a new splitter that doesn't need an external tool to convert fastq -> fasta
        #so we split and convert in 1 go.
        self.make_folder(self.path_obj.rRNA_top_path)
        self.make_folder(self.path_obj.rRNA_data_path)
        self.make_folder(self.path_obj.rRNA_jobs_path)
        self.make_folder(self.path_obj.rRNA_exe_path)
        self.make_folder(self.path_obj.rRNA_s_fa_path)
        self.make_folder(self.path_obj.rRNA_p1_fa_path)
        self.make_folder(self.path_obj.rRNA_p2_fa_path)

        split_s_fastq = self.path_obj.Python + " "
        split_s_fastq += self.path_obj.read_split_convert + " "
        split_s_fastq += os.path.join(self.path_obj.vector_final_path, "singletons.fastq") + " "
        split_s_fastq += self.path_obj.rRNA_s_fa_path + " "
        split_s_fastq += self.path_obj.rRNA_chunksize

        split_p1_fastq = self.path_obj.Python + " "
        split_p1_fastq += self.path_obj.read_split_convert + " "
        split_p1_fastq += os.path.join(self.path_obj.vector_final_path, "pair_1.fastq") + " "
        split_p1_fastq += self.path_obj.rRNA_p1_fa_path + " "  
        split_p1_fastq += self.path_obj.rRNA_chunksize 

        split_p2_fastq = self.path_obj.Python + " "
        split_p2_fastq += self.path_obj.read_split_convert  + " "
        split_p2_fastq += os.path.join(self.path_obj.vector_final_path, "pair_2.fastq") + " "
        split_p2_fastq += self.path_obj.rRNA_p2_fa_path + " "
        split_p2_fastq += self.path_obj.rRNA_chunksize 

        split_s_tut = self.path_obj.Python + " "
        split_s_tut += self.path_obj.read_split_convert + " " 
        split_s_tut += self.sequence_single + " "
        split_s_tut += self.path_obj.rRNA_s_fa_path + " "
        split_s_tut += self.path_obj.rRNA_chunksize

        split_p1_tut = self.path_obj.Python + " "
        split_p1_tut += self.path_obj.read_split_convert + " "
        split_p1_tut += self.sequence_path_1 + " "
        split_p1_tut += self.path_obj.rRNA_p1_fa_path + " "
        split_p1_tut += self.path_obj.rRNA_chunksize

        split_p2_tut = self.path_obj.Python + " "
        split_p2_tut += self.path_obj.read_split_convert + " "
        split_p2_tut += self.sequence_path_2 + " "
        split_p2_tut += self.path_obj.rRNA_p2_fa_path + " "
        split_p2_tut += self.path_obj.rRNA_chunksize
        
        make_marker = "touch" + " " + marker_path

        command = []
        if(self.tutorial_keyword == "rRNA"):
            if(self.read_mode == "single"):
                command = [split_s_tut]
            else:
                command = [split_p1_tut, split_p2_tut + " && " + make_marker]
        else:
            command = [split_s_fastq, split_p1_fastq, split_p2_fastq + " && " + make_marker]
        return command

    def create_rRNA_filter_barrnap_command(self, fasta_segment, barrnap_out_file, mRNA_file, marker_path):
        self.make_folder(self.path_obj.rRNA_s_bar_path)
        self.make_folder(self.path_obj.rRNA_p1_bar_path)
        self.make_folder(self.path_obj.rRNA_p2_bar_path)


        barrnap_a = self.path_obj.Barrnap + " " 
        barrnap_a += "--quiet --reject 0.01 --kingdom arc --threads " + self.threads_str + " "
        barrnap_a += fasta_segment + " "
        barrnap_a += ">> " + barrnap_out_file


        barrnap_b = self.path_obj.Barrnap + " "
        barrnap_b += "--quiet --reject 0.01 --kingdom bac --threads " + self.threads_str + " "
        barrnap_b += fasta_segment + " "
        barrnap_b += ">> " + barrnap_out_file

        barrnap_e = self.path_obj.Barrnap + " "
        barrnap_e += "--quiet --reject 0.01 --kingdom euk --threads " + self.threads_str + " "
        barrnap_e += fasta_segment + " "
        barrnap_e += ">> " + barrnap_out_file

        barrnap_m = self.path_obj.Barrnap + " "
        barrnap_m += "--quiet --reject 0.01 --kingdom mit --threads " + self.threads_str + " "
        barrnap_m += fasta_segment + " "
        barrnap_m += ">> " + barrnap_out_file


        barrnap_pp = self.path_obj.Python + " "
        barrnap_pp += self.path_obj.rRNA_barrnap_pp + " "
        barrnap_pp += barrnap_out_file + " "
        barrnap_pp += fasta_segment + " "
        barrnap_pp += mRNA_file

        make_marker = "touch" + " " + marker_path

        command = [barrnap_a, barrnap_b, barrnap_e, barrnap_m, barrnap_pp + " && " + make_marker]
        return command
              
    


    def create_rRNA_filter_infernal_command(self, fasta_segment, infernal_out_file):
        
        fasta_basename = os.path.basename(fasta_segment).split(".")[0]
        self.make_folder(self.path_obj.rRNA_s_inf_path)
        self.make_folder(self.path_obj.rRNA_p1_inf_path)
        self.make_folder(self.path_obj.rRNA_p2_inf_path)
        self.make_folder(self.path_obj.rRNA_final_path)
        self.make_folder(self.path_obj.rRNA_final_mRNA_path)
        self.make_folder(self.path_obj.rRNA_final_tRNA_path)
        
       
        infernal_command = self.path_obj.Infernal
        infernal_command += " -o /dev/null --tblout"        + " "
        infernal_command += infernal_out_file               + " "
        #infernal's parallelism is only good up to 4 CPUs
        if (int(self.threads_str) < 4):
            infernal_command += "--cpu 1"                   + " "
        else:
            infernal_command += "--cpu 4"                   + " "
        infernal_command += "--anytrunc --rfam -E 0.001"    + " "
        infernal_command += self.path_obj.Rfam              + " "
        infernal_command += fasta_segment

        marker_path = os.path.join(self.path_obj.rRNA_jobs_path, fasta_basename + "_inf")
        make_marker = "touch" + " " + marker_path
        return [infernal_command + " && " + make_marker]
          
    def create_rRNA_cleanup_command(self, data_style, marker_path):
        infernal_pp = self.path_obj.Python + " "
        infernal_pp += self.path_obj.rRNA_infernal_pp + " "
        infernal_pp += self.path_obj.filter_stringency + " "
        infernal_pp += data_style + " "
            
        if(data_style == "single"):
            infernal_pp += os.path.join(self.path_obj.rRNA_s_inf_path)    + " "
            infernal_pp += os.path.join(self.path_obj.vector_final_path,  "singletons.fastq")            + " " 
            infernal_pp += os.path.join(self.path_obj.rRNA_final_mRNA_path, "singletons_mRNA.fastq") + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_final_tRNA_path, "singletons_other.fastq")

        else:
            infernal_pp += os.path.join(self.path_obj.rRNA_p1_inf_path)       + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_p2_inf_path)       + " "
            infernal_pp += os.path.join(self.path_obj.vector_final_path, "pair_1.fastq")               + " "
            infernal_pp += os.path.join(self.path_obj.vector_final_path, "pair_2.fastq")               + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_final_mRNA_path, "pair_1_mRNA.fastq")    + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_final_mRNA_path, "pair_2_mRNA.fastq")    + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_final_tRNA_path, "pair_1_other.fastq")   + " "
            infernal_pp += os.path.join(self.path_obj.rRNA_final_tRNA_path, "pair_2_other.fastq") 
            
        make_marker = "touch" + " " + marker_path

        return [infernal_pp +  " && " +  make_marker]

    def create_repop_command(self, stage_name, preprocess_stage_name, dependency_stage_name):
        # This stage reintroduces the duplicate reads into the data.  We need it to count towards things.
        # Due to time, and hierarchical importance, we're leaving this stage alone.
        # Leaving it alone in a tangled state
        # But the issue is that by leaving it alone, we violate the design plan
        # The fix? We have to detect if preprocess has been run.  If so, pull the missing data there
        # if not,
        # What has to happen here:
        # -> detect if we've run the preprocess stage.
        # -> if it's run, grab data
        # -> if not, run our own custom preprocess up to what we need
        dep_loc                 = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        subfolder               = os.path.join(self.Output_Path, stage_name)
        data_folder             = os.path.join(subfolder, "data")
        repop_folder            = os.path.join(data_folder, "0_repop")
        final_folder            = os.path.join(subfolder, "final_results")
        preprocess_subfolder    = os.path.join(self.Output_Path, preprocess_stage_name)
        
        tut_keyword = "repop"

        # we ran a previous preprocess.  grab files
        # need 3, 5(clstr only), and mRNA from the 2nd stage.
        hq_path                 = os.path.join(preprocess_subfolder, "final_results")
        cluster_path            = os.path.join(preprocess_subfolder, "final_results")
        singleton_path          = os.path.join(preprocess_subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(repop_folder)
        self.make_folder(final_folder)

        repop_singletons = ">&2 echo " + str(dt.today()) + " Duplication repopulation singletons mRNA| "
        repop_singletons += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        #the reference data to be drawn from 
        if self.read_mode == "single":
            repop_singletons += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        print("TIP_TOP")
        repop_singletons += os.path.join(dep_loc, "mRNA", "singletons_mRNA.fastq") + " "  # in -> rRNA filtration output
        repop_singletons += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output

        
        if self.read_mode == "single":
            repop_singletons += os.path.join(final_folder, "singletons.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(repop_folder, "singletons.fastq")  # out
            
            

        repop_singletons_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulations singletons rRNA | "
        repop_singletons_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        if self.read_mode == "single":
            repop_singletons_rRNA += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons_rRNA += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        repop_singletons_rRNA += os.path.join(dep_loc, "other", "singletons_other.fastq") + " "  # in -> rRNA filtration output
        repop_singletons_rRNA += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output
        if self.read_mode == "single":
            repop_singletons_rRNA += os.path.join(final_folder, "singletons_rRNA.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons_rRNA += os.path.join(repop_folder, "singletons_rRNA.fastq")  # out

        repop_pair_1 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 mRNA | "
        repop_pair_1 += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_1 += os.path.join(hq_path, "pair_1_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_1 += self.sequence_path_1 + " "
        else:
            repop_pair_1 += os.path.join(dep_loc, "mRNA", "pair_1_mRNA.fastq") + " "
        repop_pair_1 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1 += os.path.join(repop_folder, "pair_1.fastq")

        repop_pair_1_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 rRNA | "
        repop_pair_1_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_1_rRNA += os.path.join(hq_path, "pair_1_match.fastq") + " "
        repop_pair_1_rRNA += os.path.join(dep_loc, "other", "pair_1_other.fastq") + " "
        repop_pair_1_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq")

        repop_pair_2 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2 += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_2 += os.path.join(hq_path, "pair_2_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_2 += self.sequence_path_2 + " "
        else:
            repop_pair_2 += os.path.join(dep_loc, "mRNA", "pair_2_mRNA.fastq") + " "
        repop_pair_2 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2 += os.path.join(repop_folder, "pair_2.fastq")

        repop_pair_2_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_2_rRNA += os.path.join(hq_path, "pair_2_match.fastq") + " "
        repop_pair_2_rRNA += os.path.join(dep_loc, "other", "pair_2_other.fastq") + " "
        repop_pair_2_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq")

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.path_obj.Python + " "
        singleton_repop_filter += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.path_obj.Python + " "
        singleton_repop_filter_rRNA += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "singletons_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "singletons_rRNA.fastq")
        
        if(self.tutorial_keyword == tut_keyword):
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_pair_1,
                    repop_pair_2,
                    singleton_repop_filter
                ]
        
        else:
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_singletons_rRNA
                ]
            elif self.read_mode == "paired":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_singletons_rRNA,
                    repop_pair_1,
                    repop_pair_1_rRNA,
                    repop_pair_2,
                    repop_pair_2_rRNA,
                    singleton_repop_filter,
                    singleton_repop_filter_rRNA
                ]

        return COMMANDS_Repopulate
        
    def create_repop_command_v2_step_1(self, stage_name, preprocess_stage_name, dependency_stage_name):
        # This stage reintroduces the duplicate reads into the data.  We need it to count towards things.
        # Due to time, and hierarchical importance, we're leaving this stage alone.
        # Leaving it alone in a tangled state
        # But the issue is that by leaving it alone, we violate the design plan
        # The fix? We have to detect if preprocess has been run.  If so, pull the missing data there
        # if not,
        # What has to happen here:
        # -> detect if we've run the preprocess stage.
        # -> if it's run, grab data
        # -> if not, run our own custom preprocess up to what we need
        
        dep_loc                 = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        subfolder               = os.path.join(self.Output_Path, stage_name)
        data_folder             = os.path.join(subfolder, "data")
        repop_folder            = os.path.join(data_folder, "0_repop")
        final_folder            = os.path.join(subfolder, "final_results")
        preprocess_subfolder    = os.path.join(self.Output_Path, preprocess_stage_name)
        
        tut_keyword = "repop"

        # we ran a previous preprocess.  grab files
        # need 3, 5(clstr only), and mRNA from the 2nd stage.
        hq_path                 = os.path.join(preprocess_subfolder, "final_results")
        cluster_path            = os.path.join(preprocess_subfolder, "final_results")
        singleton_path          = os.path.join(preprocess_subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(repop_folder)
        self.make_folder(final_folder)

        repop_singletons = ">&2 echo " + str(dt.today()) + " Duplication repopulation singletons mRNA| "
        repop_singletons += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        #the reference data to be drawn from 
        if self.read_mode == "single":
            repop_singletons += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        
        repop_singletons += os.path.join(dep_loc, "mRNA", "singletons_mRNA.fastq") + " "  # in -> rRNA filtration output
        repop_singletons += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output

        
        if self.read_mode == "single":
            repop_singletons += os.path.join(final_folder, "singletons.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(repop_folder, "singletons.fastq")  # out
            
            

        repop_singletons_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulations singletons rRNA | "
        repop_singletons_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        if self.read_mode == "single":
            repop_singletons_rRNA += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons_rRNA += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        repop_singletons_rRNA += os.path.join(dep_loc, "other", "singletons_other.fastq") + " "  # in -> rRNA filtration output
        repop_singletons_rRNA += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output
        if self.read_mode == "single":
            repop_singletons_rRNA += os.path.join(final_folder, "singletons_rRNA.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons_rRNA += os.path.join(repop_folder, "singletons_rRNA.fastq")  # out

        repop_pair_1 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 mRNA | "
        repop_pair_1 += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_1 += os.path.join(hq_path, "pair_1_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_1 += self.sequence_path_1 + " "
        else:
            repop_pair_1 += os.path.join(dep_loc, "mRNA", "pair_1_mRNA.fastq") + " "
        repop_pair_1 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1 += os.path.join(repop_folder, "pair_1.fastq")

        repop_pair_1_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 rRNA | "
        repop_pair_1_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_1_rRNA += os.path.join(hq_path, "pair_1_match.fastq") + " "
        repop_pair_1_rRNA += os.path.join(dep_loc, "other", "pair_1_other.fastq") + " "
        repop_pair_1_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq")

        repop_pair_2 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2 += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_2 += os.path.join(hq_path, "pair_2_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_2 += self.sequence_path_2 + " "
        else:
            repop_pair_2 += os.path.join(dep_loc, "mRNA", "pair_2_mRNA.fastq") + " "
        repop_pair_2 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2 += os.path.join(repop_folder, "pair_2.fastq")

        repop_pair_2_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2_rRNA += self.path_obj.Python + " " + self.path_obj.duplicate_repopulate + " "
        repop_pair_2_rRNA += os.path.join(hq_path, "pair_2_match.fastq") + " "
        repop_pair_2_rRNA += os.path.join(dep_loc, "other", "pair_2_other.fastq") + " "
        repop_pair_2_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq")

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.path_obj.Python + " "
        singleton_repop_filter += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.path_obj.Python + " "
        singleton_repop_filter_rRNA += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "singletons_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "singletons_rRNA.fastq")
        
        if(self.tutorial_keyword == tut_keyword):
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_pair_1,
                    repop_pair_2,
                    singleton_repop_filter
                ]
        
        else:
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_singletons_rRNA
                ]
            elif self.read_mode == "paired":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_singletons_rRNA,
                    repop_pair_1,
                    repop_pair_1_rRNA,
                    repop_pair_2,
                    repop_pair_2_rRNA#,
                    #singleton_repop_filter,
                    #singleton_repop_filter_rRNA
                ]

        return COMMANDS_Repopulate        
        
    def create_repop_command_v2_step_2(self, stage_name, preprocess_stage_name, dependency_stage_name):
        # This stage reintroduces the duplicate reads into the data.  We need it to count towards things.
        # Due to time, and hierarchical importance, we're leaving this stage alone.
        # Leaving it alone in a tangled state
        # But the issue is that by leaving it alone, we violate the design plan
        # The fix? We have to detect if preprocess has been run.  If so, pull the missing data there
        # if not,
        # What has to happen here:
        # -> detect if we've run the preprocess stage.
        # -> if it's run, grab data
        # -> if not, run our own custom preprocess up to what we need
        dep_loc                 = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        subfolder               = os.path.join(self.Output_Path, stage_name)
        data_folder             = os.path.join(subfolder, "data")
        repop_folder            = os.path.join(data_folder, "0_repop")
        final_folder            = os.path.join(subfolder, "final_results")
        preprocess_subfolder    = os.path.join(self.Output_Path, preprocess_stage_name)
        
        tut_keyword = "repop"

        # we ran a previous preprocess.  grab files
        # need 3, 5(clstr only), and mRNA from the 2nd stage.
        hq_path                 = os.path.join(preprocess_subfolder, "final_results")
        cluster_path            = os.path.join(preprocess_subfolder, "final_results")
        singleton_path          = os.path.join(preprocess_subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(repop_folder)
        self.make_folder(final_folder)

        

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.path_obj.Python + " "
        singleton_repop_filter += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.path_obj.Python + " "
        singleton_repop_filter_rRNA += self.path_obj.orphaned_read_filter + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(repop_folder, "singletons_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_1_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "pair_2_rRNA.fastq") + " "
        singleton_repop_filter_rRNA += os.path.join(final_folder, "singletons_rRNA.fastq")
        
        if(not self.tutorial_keyword == tut_keyword):
            if self.read_mode == "paired":
                COMMANDS_Repopulate = [
                    singleton_repop_filter,
                    singleton_repop_filter_rRNA
                ]

        return COMMANDS_Repopulate 

    def create_assemble_contigs_command(self, stage_name, dependency_stage_name):
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        dep_loc             = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        spades_folder       = os.path.join(data_folder, "0_spades")
        mgm_folder          = os.path.join(data_folder, "1_mgm")
        bwa_folder          = os.path.join(data_folder, "2_bwa_align")
        mapped_reads_folder = os.path.join(data_folder, "3_mapped_reads")
        final_folder        = os.path.join(subfolder, "final_results")
        

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(spades_folder)
        self.make_folder(bwa_folder)
        self.make_folder(mgm_folder)
        self.make_folder(final_folder)
        
        tut_keyword = "assembly"
        
        # this assembles contigs
        spades = ">&2 echo Spades Contig assembly | "
        spades += self.path_obj.Python + " "
        spades += self.path_obj.Spades + " --rna"
        if(self.tutorial_keyword == tut_keyword):
            if self.read_mode == "paired":
                spades += " -1 " + self.sequence_path_1  # in1 (pair 1)
                spades += " -2 " + self.sequence_path_2  # in2 (pair 2)
            spades += " -s " + self.sequence_single  # in_single (singletons)
        else:
            if self.read_mode == "paired":
                spades += " -1 " + os.path.join(dep_loc, "pair_1.fastq")  # in1 (pair 1)
                spades += " -2 " + os.path.join(dep_loc, "pair_2.fastq")  # in2 (pair 2)
            spades += " -s " + os.path.join(dep_loc, "singletons.fastq")  # in_single (singletons)
        spades += " -o " + spades_folder  # out

        #if there is no output, bypass contigs. -> But this is a v2 upgrade.  
        spades_rename = "cp " + os.path.join(spades_folder, "transcripts.fasta") + " " + os.path.join(spades_folder, "contigs.fasta")  # rename output
        
        original_contigs    = os.path.join(spades_folder, "contigs.fasta") 
        post_mgm_contig     = os.path.join(mgm_folder, "disassembled_contigs.fasta")
        mgm_report          = os.path.join(mgm_folder, "gene_report.txt")
        final_contigs       = os.path.join(mgm_folder, "contigs.fasta")
        contig_map          = os.path.join(final_folder, "contig_map.tsv")
        
        #-------------------------------------------------------
        #spades does too good of a job sometimes.  Disassemble it into genes.
        disassemble_contigs = ">&2 echo Disassembling contigs | "
        disassemble_contigs += self.path_obj.MetaGeneMark + " -o " + mgm_report + " "
        disassemble_contigs += "-D " + post_mgm_contig + " "
        disassemble_contigs += "-m " + self.path_obj.mgm_model + " "
        disassemble_contigs += os.path.join(spades_folder, "contigs.fasta")
        
        remove_whitespace = ">&2 echo Removing whitespace from fasta | " 
        remove_whitespace += self.path_obj.Python + " " + self.path_obj.remove_gaps_in_fasta + " "
        remove_whitespace += post_mgm_contig + " "
        remove_whitespace += final_contigs
        
        #BWA-ing against the final contigs gives us a proper contig-segment -> read map. 
        bwa_index = self.path_obj.BWA + " index -a bwtsw " + final_contigs
        
        
        # Build a report of what was consumed by contig transmutation (assemble/disassemble)
        bwa_paired_contigs = ">&2 echo BWA pair contigs | "
        bwa_paired_contigs += self.path_obj.BWA + " mem -t " + self.threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_paired_contigs += final_contigs + " "
        bwa_paired_contigs += os.path.join(dep_loc, "pair_1.fastq") + " "
        bwa_paired_contigs += os.path.join(dep_loc, "pair_2.fastq") + " "
        bwa_paired_contigs += ">" + " " 
        bwa_paired_contigs += os.path.join(bwa_folder, "paired_on_contigs.sam")

        bwa_singletons_contigs = ">&2 echo BWA singleton contigs | "
        bwa_singletons_contigs += self.path_obj.BWA + " mem -t " + self.threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_singletons_contigs += final_contigs + " "
        bwa_singletons_contigs += os.path.join(dep_loc, "singletons.fastq")
        bwa_singletons_contigs += " > " + os.path.join(bwa_folder, "singletons_on_contigs.sam")
        
        make_contig_map = ">&2 echo Making contig map | " 
        make_contig_map += self.path_obj.Python + " "
        make_contig_map += self.path_obj.Map_contig + " "
        make_contig_map += self.read_mode + " "
        make_contig_map += dep_loc + " "
        make_contig_map += final_folder + " "
        make_contig_map += os.path.join(bwa_folder, "singletons_on_contigs.sam") + " "
        if(self.read_mode == "paired"):
            make_contig_map += os.path.join(bwa_folder, "paired_on_contigs.sam")
            
            
        flush_bad_contigs = ">&2 echo flush bad contigs | " 
        flush_bad_contigs += self.path_obj.Python + " "
        flush_bad_contigs += self.path_obj.flush_bad_contigs + " "
        flush_bad_contigs += contig_map + " "
        flush_bad_contigs += final_contigs + " "
        flush_bad_contigs += os.path.join(final_folder, "contigs.fasta")
        
            
        move_gene_report = ">&2 echo moving gene report | "
        move_gene_report += "cp" + " "
        move_gene_report += os.path.join(mgm_folder, "gene_report.txt") + " "
        move_gene_report += os.path.join(final_folder, "gene_report.txt") 
            

        if self.read_mode == "single":
            COMMANDS_Assemble = [
                spades + " && " + 
                spades_rename + " && " +
                disassemble_contigs + " && " +
                remove_whitespace + " && " +
                bwa_index + " && " +
                bwa_singletons_contigs + " && " +
                make_contig_map + " && " +
                flush_bad_contigs + " && " +
                move_gene_report
            ]
        elif self.read_mode == "paired":
            COMMANDS_Assemble = [
                spades + " && " +
                spades_rename + " && " +
                disassemble_contigs + " && " +
                remove_whitespace + " && " +
                bwa_index + " && " +
                bwa_paired_contigs + " && " +
                bwa_singletons_contigs + " && " +
                make_contig_map + " && " + 
                flush_bad_contigs + " && " +
                move_gene_report
            ]

        return COMMANDS_Assemble

    def create_GA_pre_scan_command(self, marker_path):
        
        self.make_folder(self.path_obj.GA_pre_scan_libs_path)
        
        ga_get_lib = ">&2 echo GA pre-scan get libs | "
        ga_get_lib += self.path_obj.Python + " "
        ga_get_lib += self.path_obj.GA_pre_scan_get_lib + " "
        ga_get_lib += os.path.join(self.path_obj.GA_pre_scan_wevote_path, "taxonomic_classifications.tsv") + " "
        ga_get_lib += self.path_obj.taxid_tree + " "
        ga_get_lib += self.path_obj.nodes + " "
        ga_get_lib += os.path.join(self.path_obj.GA_pre_scan_libs_path, "lib_list.txt") + " "
        ga_get_lib += os.path.join(self.path_obj.GA_pre_scan_libs_path, "lib_reject.txt") + " "
        ga_get_lib += self.path_obj.source_taxa_DB + " "
        ga_get_lib += str(self.path_obj.taxa_exist_cutoff)
        
        make_marker = "touch" + " " + marker_path
        
        
        return [ga_get_lib + " && " + make_marker]
        
        
    def create_GA_pre_scan_assemble_lib_command(self, marker_path):
        
        self.make_folder(self.path_obj.GA_pre_scan_final_path)
    
        assemble_lib = ">&2 echo GA assemble libs | " 
        assemble_lib += self.path_obj.Python + " " 
        assemble_lib += self.path_obj.GA_pre_scan_assemble_lib + " "
        if(os.path.exists(self.path_obj.taxa_lib_list)):
            assemble_lib += self.path_obj.taxa_lib_list + " "
        else:
            assemble_lib += os.path.join(self.path_obj.GA_pre_scan_libs_path, "lib_list.txt") + " " 
            
        assemble_lib += self.path_obj.source_taxa_DB +  " "
        assemble_lib += self.path_obj.GA_pre_scan_final_path +  " " 
        assemble_lib += "all"
        
        #index_lib = "for i in $(ls " + final_folder + ");" + " "
        #index_lib += "do " + self.path_obj.BWA + " index" + " "
        #index_lib += final_folder + "/$i; done" 
        
        make_marker = "touch" + " " + marker_path
        
        self.path_obj.DNA_DB = self.path_obj.GA_pre_scan_final_path
        
        #return [assemble_lib + " && " + index_lib + " && " + make_marker]
        return [assemble_lib + " && " + make_marker]
        
 
    def create_split_ga_fastq_data_command(self, category, marker_path):
        
        self.make_folder(self.path_obj.GA_split_top_path)
        self.make_folder(self.path_obj.GA_split_data_path)
        self.make_folder(self.path_obj.GA_split_jobs_path)
        #self.make_folder(self.path_obj.GA_split_p1_path)
        #self.make_folder(self.path_obj.GA_split_p2_path)
        #self.make_folder(self.path_obj.GA_split_s_path)
        #self.make_folder(self.path_obj.GA_split_c_path)
        self.make_folder(self.path_obj.GA_split_final_path)

        
        
        
        if(self.tutorial_keyword == "GA"):
            if(category == "pair_1"):
            
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_path_1 + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(self.path_obj.GA_split_p1_path, category + "_")
                
                make_marker = "touch" + " " + marker_path
                
                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
                
            elif(category == "pair_2"):
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_path_2 + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(self.path_obj.GA_split_p2_path, category + "_")
                
                make_marker = "touch" + " " + marker_path
                
                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
            elif(category == "singletons"):
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_single + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(self.path_obj.GA_split_s_path, category + "_")
                
                make_marker = "touch" + " " + marker_path

                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
     
        else:
            split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
            split_fastq += "split -l " + str(int(self.path_obj.GA_chunksize) * 4) + " "        
            split_fastq += os.path.join(self.path_obj.contigs_final_path, category + ".fastq") + " "
            split_fastq += "--additional-suffix .fastq" + " "
            split_fastq += "-d" + " "
            split_fastq += os.path.join(self.path_obj.GA_split_final_path, category + "_")
            #if(category == "contigs"):
            #    split_fastq += os.path.join(self.path_obj.GA_split_c_path, category + "_")
            #elif(category == "pair_1"):
            #    split_fastq += os.path.join(self.path_obj.GA_split_p1_path, category + "_")
            #elif(category == "pair_2"):
            #    split_fastq += os.path.join(self.path_obj.GA_split_p2_path, category + "_")
            #elif(category == "singletons"):
            #    split_fastq += os.path.join(self.path_obj.GA_split_s_path, category + "_")

            make_marker = "touch" + " " + marker_path
            
            COMMANDS_GA_prep_fastq = [
                split_fastq + " && " + make_marker
            ]
            
        return COMMANDS_GA_prep_fastq

    def create_split_ga_fasta_data_command(self, category, marker_path):
        
        self.make_folder(self.path_obj.GA_split_top_path)
        self.make_folder(self.path_obj.GA_split_data_path)
        self.make_folder(self.path_obj.GA_split_final_path)
        #self.make_folder(self.path_obj.GA_split_p1_path)
        #self.make_folder(self.path_obj.GA_split_p2_path)
        #self.make_folder(self.path_obj.GA_split_c_path)
        #self.make_folder(self.path_obj.GA_split_s_path)
        self.make_folder(self.path_obj.GA_split_jobs_path)
        
        
        if(self.tutorial_keyword == "GA"):
            if(category == "singletons"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.path_obj.Python + " "    
                split_fasta += self.path_obj.File_splitter + " "
                split_fasta += self.sequence_single + " "
                split_fasta += os.path.join(self.path_obj.GA_split_s_path, category) + " "
                split_fasta += str(self.path_obj.GA_chunksize)
                
                make_marker = "touch" + " " + marker_path
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
                
            elif(category == "contigs"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.path_obj.Python + " "    
                split_fasta += self.path_obj.File_splitter + " "
                split_fasta += self.sequence_contigs + " "
                split_fasta += os.path.join(self.path_obj.GA_split_c_path, category) + " "
                print("split_fasta", split_fasta)

                time.sleep(10)
                split_fasta += str(self.path_obj.GA_chunksize)
                
                make_marker = "touch" + " " + marker_path
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
            elif(category == "pair_1"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.path_obj.Python + " "    
                split_fasta += self.path_obj.File_splitter + " "
                split_fasta += self.sequence_path_1 + " "
                split_fasta += os.path.join(self.path_obj.GA_split_p1_path, category) + " "
                split_fasta += str(self.path_obj.GA_chunksize)
                
                make_marker = "touch" + " " + marker_path   
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
            elif(category == "pair_2"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.path_obj.Python + " "    
                split_fasta += self.path_obj.File_splitter + " "
                split_fasta += self.sequence_path_2 + " "
                split_fasta += os.path.join(self.path_obj.GA_split_p2_path, category) + " "
                split_fasta += str(self.path_obj.GA_chunksize)
                
                make_marker = "touch" + " " + marker_path
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
                
        else:
            split_fasta = ">&2 echo splitting fasta for " + category + " | "
            split_fasta += self.path_obj.Python + " "    
            split_fasta += self.path_obj.File_splitter + " "
            split_fasta += os.path.join(self.path_obj.contigs_final_path, category +".fasta") + " "
            split_fasta += os.path.join(self.path_obj.GA_split_final_path, category) + " "
            #if(category == "contigs"):
            #    split_fasta += os.path.join(self.path_obj.GA_split_c_path, "contigs") + " "
            #elif(category == "pair_1"):
            #    split_fasta += os.path.join(self.path_obj.GA_split_p1_path, category) + " "
            #elif(category == "pair_2"):
            #    split_fasta += os.path.join(self.path_obj.GA_split_p2_path, category) + " "
            #elif(category == "singletons"):
            #    split_fasta += os.path.join(self.path_obj.GA_split_s_path, category) + " "
            split_fasta += str(self.path_obj.GA_chunksize)
            
            make_marker = "touch" + " " + marker_path
            
            COMMANDS_GA_prep_fasta = [
                split_fasta + " && " + make_marker
            ]
        
        return COMMANDS_GA_prep_fasta


    def create_BWA_annotate_command_v2(self, ref_path, query_path, marker_path):
        # meant to be called multiple times: query file is a split file
        # aug 10, 2021: changed ref path to accomodate new split-chocophlan

        self.make_folder(self.path_obj.GA_BWA_top_path)
        self.make_folder(self.path_obj.GA_BWA_data_path)
        self.make_folder(self.path_obj.GA_BWA_jobs_path)
        self.make_folder(self.path_obj.GA_BWA_run_path)

        file_tag = os.path.basename(query_path)
        file_tag = os.path.splitext(file_tag)[0]

        ref_tag = os.path.basename(ref_path)
        ref_tag = os.path.splitext(ref_tag)[0]
        
        bwa_job = ">&2 echo " + str(dt.today()) + " BWA on " + file_tag + " | "
        bwa_job += self.path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_job += ref_path + " "
        #bwa_job += os.path.join(dep_loc, section_file) + " | "
        bwa_job += query_path + " | "
        bwa_job += self.path_obj.SAMTOOLS + " view "
        bwa_job += "> " + os.path.join(self.path_obj.GA_BWA_run_path, file_tag +"_" + ref_tag + ".sam")
        
        #make_marker = ">&2 echo marking BWA job complete: " + file_tag + " | "
        make_marker = "touch" + " " + marker_path

        COMMANDS_BWA = [
            bwa_job + " && " + make_marker
        ]

        return COMMANDS_BWA
        
        
    def create_BWA_pp_command_v2(self, ref_path, query_file, marker_path):
        self.make_folder(self.path_obj.GA_BWA_final_path)
        self.make_folder(self.path_obj.GA_BWA_pp_path)
        self.make_folder(self.path_obj.GA_BWA_u_path)

        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]

        ref_tag = os.path.basename(ref_path)
        ref_tag = os.path.splitext(ref_tag)[0]    
        #meant to be called on the split-file version.  PP script will not merge gene maps.
        
        sample_tag = sample_root_name + "_" + ref_tag
        sam_name =  sample_tag + ".sam"
        unscanned_name = sample_tag + ".fasta"
        gene_map_name = sample_tag + "_gene_map.tsv"
        genes_name = sample_tag + "_mapped_genes.fna"
        self.make_folder(self.path_obj.GA_BWA_pp_path)
        self.make_folder(self.path_obj.GA_BWA_u_path)
        
        reads_in    = query_file
        bwa_in      = os.path.join(self.path_obj.GA_BWA_run_path, sam_name)
        reads_out = ""
        if(self.path_obj.GA_DB_mode == "multi"):
            print(dt.today(), "BWA_pp running in split-mode")
            reads_out   = os.path.join(self.path_obj.GA_BWA_u_path, unscanned_name)
        else:
            print(dt.today(), "BWA_pp running in single-mode")
            reads_out = os.path.join(self.path_obj.GA_BWA_final_path, unscanned_name)
        

        map_read_bwa = ">&2 echo " + str(dt.today()) + " GA BWA PP generic: " + sample_root_name + " | "
        map_read_bwa += self.path_obj.Python + " "
        map_read_bwa += self.path_obj.Map_reads_gene_BWA + " "
        map_read_bwa += str(self.path_obj.BWA_cigar_cutoff) + " "
        map_read_bwa += ref_path + " "
        if(self.sequence_contigs == "None"):
            map_read_bwa += "None" + " "
        else:        
            map_read_bwa += os.path.join(self.path_obj.contigs_final_path, "contig_map.tsv") + " "  # IN
        map_read_bwa += os.path.join(self.path_obj.GA_BWA_final_path, gene_map_name) + " "  # OUT
        map_read_bwa += os.path.join(self.path_obj.GA_BWA_final_path, genes_name ) + " " #OUT
        map_read_bwa += reads_in + " "
        map_read_bwa += bwa_in + " "
        map_read_bwa += reads_out



        

        make_marker = ">&2 echo bwa pp complete: " + sample_tag + " | " 
        make_marker += "touch" + " " + marker_path

        COMMANDS_Annotate_BWA = [
            map_read_bwa + " && " + make_marker
        ]

        return COMMANDS_Annotate_BWA


 

    def create_BWA_copy_contig_map_command(self, marker_path):
        
    
        copy_contig_map = ">&2 echo " + str(dt.today()) + " copy contig map | "
        copy_contig_map += "cp " + os.path.join(self.path_obj.contigs_final_path, "contig_map.tsv") + " " + os.path.join(self.path_obj.GA_BWA_final_path, "contig_map.tsv")
        
        make_marker = "touch" + " " + marker_path
        
        return [copy_contig_map + " && " + make_marker]

    def create_merge_BWA_fasta_command(self, query_file, marker_path):
        #merges all FASTAs exported by BWA from the split-DB runs.
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]

        merge_bwa_fastas = ">&2 echo " + str(dt.today()) + " GA BWA merge leftover reads " + sample_root_name + " | "
        merge_bwa_fastas += self.path_obj.Python + " "
        merge_bwa_fastas += self.path_obj.GA_merge_fasta + " "
        merge_bwa_fastas += self.path_obj.GA_BWA_u_path + " " 
        merge_bwa_fastas += sample_root_name + " " 
        merge_bwa_fastas += self.path_obj.GA_BWA_final_path

        make_marker = "touch" + " " + marker_path
        return [merge_bwa_fastas + " && " + make_marker]

    def create_BLAT_annotate_command_v2(self, stage_name, query_file, db_path, fasta_db, marker_file):
        
        #takes in a sample query file (expecting a segment of the whole GA data, after BWA
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        
        subfolder   = os.path.join(self.Output_Path, stage_name)
        data_folder = os.path.join(subfolder, "data")
        blat_folder = os.path.join(data_folder, "0_blat")
        jobs_folder = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)
        self.make_folder(jobs_folder)

        blat_command = ">&2 echo " + str(dt.today()) + " BLAT annotation for " + sample_root_name + " " + fasta_db + " | "
        blat_command += self.path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_command += os.path.join(db_path, fasta_db) + " "
        blat_command += query_file
        blat_command += " -fine -q=rna -t=dna -out=blast8 -threads=40" + " "
        blat_command += os.path.join(blat_folder, sample_root_name + "_" + fasta_db + ".blatout")
         
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        if(os.path.getsize(query_file) > 0):
            return [blat_command + " && " + make_marker]
        else:
            dummy_blat_command = ">&2 echo " + str(dt.today()) + " Not running BLAT command on empty file: " + query_file
            return [dummy_blat_command]
        
        
    def create_BLAT_cat_command_v2(self, stage_name, query_file, marker_file):
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        # This merges each blatout file based on the sample's name
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        blat_folder         = os.path.join(data_folder, "0_blat")
        #blat_merge_folder   = os.path.join(data_folder, "1_blat_merge")
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        #self.make_folder(blat_merge_folder)
        self.make_folder(jobs_folder)

        cat_command = ">&2 echo " + str(dt.today()) + " combining and deleting BLATout | "
        cat_command += "for f in " + os.path.join(blat_folder, sample_root_name + "_*.blatout") + ";  do cat $f >> " + os.path.join(blat_merge_folder, sample_root_name + ".blatout") + " && rm $f; done"
        
        make_marker = ">&2 echo completed BLAT cat: " + marker_file + " | " 
        make_marker += "touch" + " "
        make_marker += (os.path.join(jobs_folder, marker_file))
        
        return [
            cat_command + " && " + make_marker
            #cleanup_command
        ]
        

    def create_BLAT_pp_command_v2(self, stage_name, query_file, dependency_stage_name, ref_file, marker_file):
        # this call is meant to be run after the BLAT calls have been completed.
        #aug 16, 2021: modded to consider the split-chocophlan 
        #Dec 14 2021: this one remains the old variant for back compatibility
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        blat_folder         = os.path.join(data_folder, "0_blat")
        final_folder        = os.path.join(subfolder, "final_results")
        dep_loc             = os.path.join(self.Output_Path, dependency_stage_name, "final_results")  # implied to be BWA
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        blat_pp = ">&2 echo " + str(dt.today()) + " BLAT post-processing " + sample_root_name + " | "
        blat_pp += self.path_obj.Python + " "
        blat_pp += self.path_obj.Map_reads_gene_BLAT + " "
        blat_pp += str(self.path_obj.BLAT_identity_cutoff) + " "
        blat_pp += str(self.path_obj.BLAT_length_cutoff) + " "
        blat_pp += str(self.path_obj.BLAT_score_cutoff) + " "
        blat_pp += ref_file + " " 
        
        if(self.sequence_contigs == "None"):
            blat_pp += "None" + " "
        else:
            blat_pp += os.path.join(dep_loc, "contig_map.tsv") + " "
        blat_pp += os.path.join(final_folder, sample_root_name + "_mapped_genes.fna") + " "
        blat_pp += os.path.join(final_folder, sample_root_name + "_gene_map.tsv") + " "
        blat_pp += query_file + " "
        blat_pp += os.path.join(blat_folder, sample_root_name + ".blatout") + " "
        blat_pp += os.path.join(final_folder, sample_root_name + ".fasta") + " "
        
        make_marker = ">&2 echo BLAT pp complete: " + marker_file + " | "
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_Annotate_BLAT_Post = [blat_pp + " && " + make_marker]

        return COMMANDS_Annotate_BLAT_Post
        
    def create_BLAT_pp_command_v3(self, stage_name, reads_in, dependency_stage_name, ref_file, marker_file):
        # this call is meant to be run after the BLAT calls have been completed.
        #aug 16, 2021: modded to consider the split-chocophlan
        #oct 22, 2021: modded to consider that we now use a compact form of splitting to cut down on file numbers
        sample_root_name = os.path.basename(reads_in)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        sample_ref_root_name = os.path.basename(ref_file)
        no_ext_ref_root_name = sample_ref_root_name.strip(".fasta")
        
        
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        blat_folder         = os.path.join(data_folder, "0_blat")
        pp_folder           = os.path.join(data_folder, "1_pp")
        final_folder        = os.path.join(subfolder, "final_results")
        dep_loc             = os.path.join(self.Output_Path, dependency_stage_name, "final_results")  # implied to be BWA
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)
        self.make_folder(pp_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        blat_pp = ">&2 echo " + str(dt.today()) + " BLAT post-processing " + sample_root_name + " | "
        blat_pp += self.path_obj.Python + " "
        blat_pp += self.path_obj.Map_reads_gene_BLAT + " "
        blat_pp += str(self.path_obj.BLAT_identity_cutoff) + " "
        blat_pp += str(self.path_obj.BLAT_length_cutoff) + " "
        blat_pp += str(self.path_obj.BLAT_score_cutoff) + " "
        blat_pp += ref_file + " " 
        
        if(self.sequence_contigs == "None"):
            blat_pp += "None" + " "
        else:
            blat_pp += os.path.join(dep_loc, "contig_map.tsv") + " "
        blat_pp += os.path.join(final_folder, sample_root_name + "_" + no_ext_ref_root_name + "_mapped_genes.fna") + " "
        blat_pp += os.path.join(final_folder, sample_root_name + "_" + no_ext_ref_root_name + "_gene_map.tsv") + " "
        blat_pp += reads_in + " "
        blat_pp += os.path.join(blat_folder, sample_root_name + "_" + sample_ref_root_name+ ".blatout") + " "
        blat_pp += os.path.join(pp_folder, sample_root_name + "_" + no_ext_ref_root_name + ".fasta") + " "

        make_marker = ">&2 echo BLAT pp complete: " + marker_file + " | "
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_Annotate_BLAT_Post = [blat_pp + " && " + make_marker]

        return COMMANDS_Annotate_BLAT_Post        

    def create_BLAT_copy_contig_map_command(self, stage_name, dependency_stage_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        final_folder    = os.path.join(subfolder, "final_results")
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(jobs_folder)
        self.make_folder(final_folder)
    
        copy_contig_map = ">&2 echo " + str(dt.today()) + " copy contig map | "
        copy_contig_map += "cp " + os.path.join(dep_loc, "contig_map.tsv") + " " + os.path.join(final_folder, "contig_map.tsv")
        
        make_marker = ">&2 echo copy contig map done | "
        make_marker += "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [copy_contig_map + " && " + make_marker]

    def create_BLAT_merge_fasta_command(self, stage_name, sample_root_name, marker_file):
        
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        blat_folder         = os.path.join(data_folder, "0_blat")
        pp_folder           = os.path.join(data_folder, "1_pp")
        final_folder        = os.path.join(subfolder, "final_results")
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)
        self.make_folder(pp_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        merge_blat_fastas = ">&2 echo " + str(dt.today()) + " GA BLAT merge leftover reads " + sample_root_name + " | "
        merge_blat_fastas += self.path_obj.Python + " "
        merge_blat_fastas += self.path_obj.GA_merge_fasta + " "
        merge_blat_fastas += pp_folder + " " 
        merge_blat_fastas += sample_root_name + " " 
        merge_blat_fastas += final_folder

        make_marker = ">&2 echo merge BLAT leftover fastas: " + marker_file + " | " 
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)

        return [merge_blat_fastas + " && " + make_marker]

        
    def create_DIAMOND_annotate_command_v2(self, stage_name, query_file, marker_file):
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
    
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        #dep_loc             = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        diamond_folder      = os.path.join(data_folder, "0_diamond")
        main_temp_folder    = os.path.join(data_folder, sample_root_name + "_diamond_temp")
        temp_folder         = os.path.join(main_temp_folder, "temp")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(diamond_folder)
        self.make_folder(main_temp_folder)
        self.make_folder(temp_folder)
        self.make_folder(jobs_folder)
        
        diamond_annotate = ">&2 echo " + str(dt.today()) + " GA DIAMOND " + sample_root_name + " | "
        diamond_annotate += self.path_obj.DIAMOND
        diamond_annotate += " blastx -p " + self.threads_str
        diamond_annotate += " -d " + self.path_obj.Prot_DB
        diamond_annotate += " -q " + query_file 
        diamond_annotate += " -o " + os.path.join(diamond_folder, sample_root_name + ".dmdout")
        diamond_annotate += " -f 6 -t " + temp_folder #section_temp_folder
        diamond_annotate += " -k 10 --id 85 --query-cover 65 --min-score 60 --unal 1"

        #make_marker = ">&2 echo marking DIAMOND complete: " + sample_root_name + " | "
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        return [diamond_annotate + " && " + make_marker]


   
    def create_DIAMOND_pp_command_v2(self, stage_name, dependency_stage_name, query_file, marker_file):
    
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        # the command just calls the merger program
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")  # implied to be blat pp
        diamond_folder  = os.path.join(data_folder, "0_diamond/")
        final_folder    = os.path.join(subfolder, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        diamond_pp = ">&2 echo " + str(dt.today()) + " DIAMOND post process " + sample_root_name + " | "
        diamond_pp += self.path_obj.Python + " "
        diamond_pp += self.path_obj.Map_reads_prot_DMND + " "
        diamond_pp += str(self.path_obj.DIAMOND_identity_cutoff) + " "
        diamond_pp += str(self.path_obj.DIAMOND_length_cutoff) + " "
        diamond_pp += str(self.path_obj.DIAMOND_score_cutoff) + " "
        diamond_pp += self.path_obj.Prot_DB_reads + " "                # IN
        if(self.sequence_contigs == "None"):
            diamond_pp += "None" + " "
        else:
            diamond_pp += os.path.join(dep_loc, "contig_map.tsv") + " "         # IN
        diamond_pp += os.path.join(final_folder, sample_root_name + "_diamond_gene_map.tsv") + " "      # OUT
        diamond_pp += os.path.join(final_folder, sample_root_name + "_diamond_proteins.faa") + " "      # OUT
        
        diamond_pp += query_file + " "                                                  # IN
        diamond_pp += os.path.join(diamond_folder, sample_root_name + ".dmdout") + " "  # IN
        diamond_pp += os.path.join(final_folder, sample_root_name + ".fasta") + " "     # OUT
        
        make_marker = ">&2 echo diamond pp complete: " + sample_root_name + " | "
        make_marker += "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_Annotate_Diamond_Post = [
            diamond_pp + " && " + make_marker
        ]

        return COMMANDS_Annotate_Diamond_Post 



    def create_GA_final_merge_command(self, current_stage_name, dep_0_name, dep_1_name, dep_2_name, dep_3_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, current_stage_name)
        data_folder     = os.path.join(subfolder, "data")
        final_folder    = os.path.join(subfolder, "final_results")
        dep_0_path      = os.path.join(self.Output_Path, dep_0_name, "final_results")   #assemble-contigs
        dep_1_path      = os.path.join(self.Output_Path, dep_1_name, "final_results")   #bwa
        dep_2_path      = os.path.join(self.Output_Path, dep_2_name, "final_results")   #blat
        dep_3_path      = os.path.join(self.Output_Path, dep_3_name, "final_results")   #dmd
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        
        final_merge_fastq = self.path_obj.Python + " "
        final_merge_fastq += self.path_obj.GA_final_merge_fasta + " "
        final_merge_fastq += dep_0_path + " "
        final_merge_fastq += dep_3_path + " "
        final_merge_fastq += self.read_mode + " "
        final_merge_fastq += final_folder
        
        final_merge_proteins = self.path_obj.Python + " "
        final_merge_proteins += self.path_obj.GA_final_merge_proteins + " "
        final_merge_proteins += dep_1_path + " "
        final_merge_proteins += dep_2_path + " "
        final_merge_proteins += dep_3_path + " "
        final_merge_proteins += final_folder
        
        final_merge_maps = self.path_obj.Python + " "
        final_merge_maps += self.path_obj.GA_final_merge_maps + " "
        final_merge_maps += dep_1_path + " "
        final_merge_maps += dep_2_path + " "
        final_merge_maps += dep_3_path + " "
        final_merge_maps += final_folder
        
        
        
        make_marker_p = ">&2 echo " + str(dt.today()) + " GA final merge | "
        make_marker_p += "touch" + " "
        make_marker_p += os.path.join(jobs_folder, marker_file + "_proteins")
        
        make_marker_f = ">&2 echo " + str(dt.today()) + " GA final merge | "
        make_marker_f += "touch" + " "
        make_marker_f += os.path.join(jobs_folder, marker_file + "_fastq")
        
        make_marker_m = ">&2 echo " + str(dt.today()) + " GA final merge | "
        make_marker_m += "touch" + " "
        make_marker_m += os.path.join(jobs_folder, marker_file + "_maps")
        
        
        COMMANDS_ga_final_merge = [
            final_merge_maps + " && " + make_marker_m,
            final_merge_fastq + " && " + make_marker_f,
            final_merge_proteins + " && " + make_marker_p
        ]
        
        return COMMANDS_ga_final_merge
    
    def create_GA_pre_scan_kraken2_command(self, operating_mode, marker_path):
        
        self.make_folder(self.path_obj.GA_pre_scan_top_path)
        self.make_folder(self.path_obj.GA_pre_scan_data_path)
        self.make_folder(self.path_obj.GA_pre_scan_jobs_path)
        self.make_folder(self.path_obj.GA_pre_scan_kraken_path)
        

        if(operating_mode == "contigs"):
            kraken2_c = ">&2 echo Kraken2 on contigs | "
            kraken2_c += self.path_obj.kraken2 + " "
            kraken2_c += "--db " + self.path_obj.kraken2_db + " "
            kraken2_c += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_c += os.path.join(self.path_obj.contigs_final_path, "contigs.fasta") + " "
            kraken2_c += "--output " + os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_c_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_c + " && " + make_marker]
            
        elif(operating_mode == "singletons"):
            kraken2_s = ">&2 echo Kraken2 on singletons | "
            kraken2_s += self.path_obj.kraken2 + " "
            kraken2_s += "--db " + self.path_obj.kraken2_db + " "
            kraken2_s += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_s += os.path.join(self.path_obj.contigs_final_path, "singletons.fastq") + " " 
            kraken2_s += "--output " + os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_s_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_s + " && " + make_marker]
            
        elif(operating_mode == "paired"):
            kraken2_p = ">&2 echo Kraken2 on paired | " 
            kraken2_p += self.path_obj.kraken2 + " "
            kraken2_p += "--db " + self.path_obj.kraken2_db +  " "
            kraken2_p += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_p += "--paired " + os.path.join(self.path_obj.contigs_final_path, "pair_1.fastq") + " " + os.path.join(self.path_obj.contigs_final_path, "pair_2.fastq") + " "
            kraken2_p += "--output " + os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_p_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_p + " && " + make_marker]

    def create_TA_kraken2_command(self, operating_mode, marker_path):
        
        self.make_folder(self.path_obj.TA_top_path)
        self.make_folder(self.path_obj.TA_data_path)
        self.make_folder(self.path_obj.TA_jobs_path)
        self.make_folder(self.path_obj.TA_kraken_path)
        

        if(operating_mode == "contigs"):
            kraken2_c = ">&2 echo Kraken2 on contigs | "
            kraken2_c += self.path_obj.kraken2 + " "
            kraken2_c += "--db " + self.path_obj.kraken2_db + " "
            kraken2_c += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_c += os.path.join(self.path_obj.contigs_final_path, "contigs.fasta") + " "
            kraken2_c += "--output " + os.path.join(self.path_obj.TA_kraken_path, "kraken2_c_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_c + " && " + make_marker]
            
        elif(operating_mode == "singletons"):
            kraken2_s = ">&2 echo Kraken2 on singletons | "
            kraken2_s += self.path_obj.kraken2 + " "
            kraken2_s += "--db " + self.path_obj.kraken2_db + " "
            kraken2_s += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_s += os.path.join(self.path_obj.contigs_final_path, "singletons.fastq") + " " 
            kraken2_s += "--output " + os.path.join(self.path_obj.TA_kraken_path, "kraken2_s_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_s + " && " + make_marker]
            
        elif(operating_mode == "paired"):
            kraken2_p = ">&2 echo Kraken2 on paired | " 
            kraken2_p += self.path_obj.kraken2 + " "
            kraken2_p += "--db " + self.path_obj.kraken2_db +  " "
            kraken2_p += "--threads " + str(self.path_obj.num_threads) + " "
            kraken2_p += "--paired " + os.path.join(self.path_obj.contigs_final_path, "pair_1.fastq") + " " + os.path.join(self.path_obj.contigs_final_path, "pair_2.fastq") + " "
            kraken2_p += "--output " + os.path.join(self.path_obj.TA_kraken_path, "kraken2_p_report.txt")
            
            make_marker = "touch " + marker_path
            
            return [kraken2_p + " && " + make_marker]
        
    def create_GA_pre_scan_kraken2_pp_command(self, marker_path):

        cat_kraken2 = ">&2 echo merging kraken2 reports | "
        cat_kraken2 += "cat "
        cat_kraken2 += os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_s_report.txt") + " "
        if(self.sequence_contigs != "None"):
            cat_kraken2 += os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_c_report.txt") + " "
        if(self.read_mode == "paired"):
            cat_kraken2 += os.path.join(self.path_obj.GA_pre_scan_kraken_path, "kraken2_p_report.txt") + " "
        cat_kraken2 += "> " + os.path.join(self.path_obj.GA_pre_scan_kraken_path, "merged_kraken2.txt")
        
        make_marker = "touch" + " " + marker_path
        
        return [cat_kraken2 + " && " + make_marker]
                    

    def create_TA_kraken2_pp_command(self, current_stage_name, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        kraken2_folder            = os.path.join(data_folder, "1_kraken2")
        jobs_folder             = os.path.join(data_folder, "jobs")

        cat_kraken2 = ">&2 echo merging kraken2 reports | "
        cat_kraken2 += "cat "
        cat_kraken2 += os.path.join(kraken2_folder, "kraken2_s_report.txt") + " "
        if(self.sequence_contigs != "None"):
            cat_kraken2 += os.path.join(kraken2_folder, "kraken2_c_report.txt") + " "
        if(self.read_mode == "paired"):
            cat_kraken2 += os.path.join(kraken2_folder, "kraken2_p_report.txt") + " "
        cat_kraken2 += "> " + os.path.join(kraken2_folder, "merged_kraken2.txt")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [cat_kraken2 + " && " + make_marker]
        
    def create_GA_pre_scan_centrifuge_command(self, operating_mode, marker_path):
        
        self.make_folder(self.path_obj.GA_pre_scan_centr_path)
        
        singletons_extension = os.path.splitext(self.sequence_single)[1]
        
        if(operating_mode == "contigs"):
            patch_contig_name = self.path_obj.Python + " "
            patch_contig_name += self.path_obj.ta_contig_name_convert + " "
            if(self.tutorial_keyword == "TA"):
                patch_contig_name += self.sequence_contigs + " "
            else:
                patch_contig_name += os.path.join(self.path_obj.contigs_final_path, "contigs.fasta") + " "
            patch_contig_name += os.path.join(self.path_obj.GA_pre_scan_centr_path, "contigs_renamed.fasta")
        
            centrifuge_on_contigs = ">&2 echo centrifuge on contigs | "
            centrifuge_on_contigs += self.path_obj.Centrifuge
            centrifuge_on_contigs += " -f -x " + self.path_obj.Centrifuge_db
            centrifuge_on_contigs += " -U " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "contigs_renamed.fasta")
            centrifuge_on_contigs += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_contigs += " --phred" + self.Qual_str
            centrifuge_on_contigs += " -p 6"
            centrifuge_on_contigs += " -S " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "raw_contigs.tsv")
            centrifuge_on_contigs += " --report-file " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "raw_contigs.txt")
            
            back_convert_report = self.path_obj.Python + " "
            back_convert_report += self.path_obj.ta_contig_name_convert + " "
            back_convert_report += os.path.join(self.path_obj.GA_pre_scan_centr_path, "raw_contigs.tsv") + " "
            back_convert_report += os.path.join(self.path_obj.GA_pre_scan_centr_path, "contigs.tsv")
            
            make_marker = "touch" + " " + marker_path
            return [patch_contig_name + " && " + centrifuge_on_contigs + " && " + back_convert_report + " && " +  make_marker]

            
        elif(operating_mode == "reads"):
            centrifuge_on_reads = ">&2 echo centrifuge on reads | "
            centrifuge_on_reads += self.path_obj.Centrifuge
            centrifuge_on_reads += " -x " + self.path_obj.Centrifuge_db
            
            if(self.tutorial_keyword == "TA"):
                if(singletons_extension == ".fa" or singletons_extension == ".fasta"):
                    centrifuge_on_reads += " -f -U " + self.sequence_single
                else:
                    centrifuge_on_reads += " -U " + self.sequence_single
                if self.read_mode == "paired":
                    centrifuge_on_reads += " -1 " + self.sequence_path_1
                    centrifuge_on_reads += " -2 " + self.sequence_path_2
            else:
                centrifuge_on_reads += " -U " + os.path.join(self.path_obj.contigs_final_path, "singletons.fastq")
                if self.read_mode == "paired":
                    centrifuge_on_reads += " -1 " + os.path.join(self.path_obj.contigs_final_path, "pair_1.fastq")
                    centrifuge_on_reads += " -2 " + os.path.join(self.path_obj.contigs_final_path, "pair_2.fastq")
            centrifuge_on_reads += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_reads += " --phred" + self.Qual_str
            centrifuge_on_reads += " -p 6"
            centrifuge_on_reads += " -S " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "reads.tsv")
            centrifuge_on_reads += " --report-file " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "reads.txt")

            make_marker = "touch" + " " + marker_path
            return [centrifuge_on_reads + " && " + make_marker]
            
    
        
    def create_TA_centrifuge_command(marker_path):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        rRNA_folder             = os.path.join(self.Output_Path, rRNA_stage, "final_results", "other")
        assemble_contigs_folder = os.path.join(self.Output_Path, assemble_contigs_stage, "final_results")
        centrifuge_folder       = os.path.join(data_folder, "2_centrifuge")
        jobs_folder             = os.path.join(data_folder, "jobs")
        final_folder            = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(centrifuge_folder)
        self.make_folder(jobs_folder)
        self.make_folder(final_folder)
        
        singletons_extension = os.path.splitext(self.sequence_single)[1]
        
        if(operating_mode == "contigs"):
            patch_contig_name = self.path_obj.Python + " "
            patch_contig_name += self.path_obj.ta_contig_name_convert + " "
            if(self.tutorial_keyword == "TA"):
                patch_contig_name += self.sequence_contigs + " "
            else:
                patch_contig_name += os.path.join(assemble_contigs_folder, "contigs.fasta") + " "
            patch_contig_name += os.path.join(centrifuge_folder, "contigs_renamed.fasta")
        
            centrifuge_on_contigs = ">&2 echo centrifuge on contigs | "
            centrifuge_on_contigs += self.path_obj.Centrifuge
            centrifuge_on_contigs += " -f -x " + self.path_obj.Centrifuge_db
            centrifuge_on_contigs += " -U " + os.path.join(centrifuge_folder, "contigs_renamed.fasta")
            centrifuge_on_contigs += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_contigs += " --phred" + self.Qual_str
            centrifuge_on_contigs += " -p 6"
            centrifuge_on_contigs += " -S " + os.path.join(centrifuge_folder, "raw_contigs.tsv")
            centrifuge_on_contigs += " --report-file " + os.path.join(centrifuge_folder, "raw_contigs.txt")
            
            back_convert_report = self.path_obj.Python + " "
            back_convert_report += self.path_obj.ta_contig_name_convert + " "
            back_convert_report += os.path.join(centrifuge_folder, "raw_contigs.tsv") + " "
            back_convert_report += os.path.join(centrifuge_folder, "contigs.tsv")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [patch_contig_name + " && " + centrifuge_on_contigs + " && " + back_convert_report + " && " +  make_marker]

            
        elif(operating_mode == "reads"):
            centrifuge_on_reads = ">&2 echo centrifuge on reads | "
            centrifuge_on_reads += self.path_obj.Centrifuge
            centrifuge_on_reads += " -x " + self.path_obj.Centrifuge_db
            
            if(self.tutorial_keyword == "TA"):
                if(singletons_extension == ".fa" or singletons_extension == ".fasta"):
                    centrifuge_on_reads += " -f -U " + self.sequence_single
                else:
                    centrifuge_on_reads += " -U " + self.sequence_single
                if self.read_mode == "paired":
                    centrifuge_on_reads += " -1 " + self.sequence_path_1
                    centrifuge_on_reads += " -2 " + self.sequence_path_2
            else:
                centrifuge_on_reads += " -U " + os.path.join(assemble_contigs_folder, "singletons.fastq")
                if self.read_mode == "paired":
                    centrifuge_on_reads += " -1 " + os.path.join(assemble_contigs_folder, "pair_1.fastq")
                    centrifuge_on_reads += " -2 " + os.path.join(assemble_contigs_folder, "pair_2.fastq")
            centrifuge_on_reads += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_reads += " --phred" + self.Qual_str
            centrifuge_on_reads += " -p 6"
            centrifuge_on_reads += " -S " + os.path.join(centrifuge_folder, "reads.tsv")
            centrifuge_on_reads += " --report-file " + os.path.join(centrifuge_folder, "reads.txt")

            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
        
            return [centrifuge_on_reads + " && " + make_marker]
            
        elif(operating_mode == "rRNA"):
        
            centrifuge_on_rRNA = ">&2 echo centrifuge on rRNA | "
            centrifuge_on_rRNA += self.path_obj.Centrifuge
            centrifuge_on_rRNA += " -x " + self.path_obj.Centrifuge_db
            centrifuge_on_rRNA += " -U " + os.path.join(rRNA_folder, "singletons_other.fastq")
            if self.read_mode == "paired":
                centrifuge_on_rRNA += " -1 " + os.path.join(rRNA_folder, "pair_1_other.fastq")
                centrifuge_on_rRNA += " -2 " + os.path.join(rRNA_folder, "pair_2_other.fastq")
            centrifuge_on_rRNA += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_rRNA += " --phred" + self.Qual_str
            centrifuge_on_rRNA += " -p 6"
            centrifuge_on_rRNA += " -S " + os.path.join(final_folder, "other.tsv")
            centrifuge_on_rRNA += " --report-file " + os.path.join(final_folder, "other.txt")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [centrifuge_on_rRNA + " &&  " + make_marker]
    
    def create_GA_pre_scan_centrifuge_pp_command(self, marker_path):
        cat_centrifuge = ">&2 echo combining all centrifuge results | "
        cat_centrifuge += "cat "
        cat_centrifuge += os.path.join(self.path_obj.GA_pre_scan_centr_path, "reads.tsv") + " "
        if(self.sequence_contigs != "None"):
            cat_centrifuge += os.path.join(self.path_obj.GA_pre_scan_centr_path, "contigs.tsv")
        cat_centrifuge += " > " + os.path.join(self.path_obj.GA_pre_scan_centr_path, "merged_centrifuge.tsv")

        make_marker = "touch" + " " + marker_path
        
        return [cat_centrifuge + " && " + make_marker]


    def create_TA_centrifuge_pp_command(self, current_stage_name, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        centrifuge_folder       = os.path.join(data_folder, "2_centrifuge")
        jobs_folder             = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(centrifuge_folder)
        self.make_folder(jobs_folder)        
    
        cat_centrifuge = ">&2 echo combining all centrifuge results | "
        cat_centrifuge += "cat "
        cat_centrifuge += os.path.join(centrifuge_folder, "reads.tsv") + " "
        if(self.sequence_contigs != "None"):
            cat_centrifuge += os.path.join(centrifuge_folder, "contigs.tsv")
        cat_centrifuge += " > " + os.path.join(centrifuge_folder, "merged_centrifuge.tsv")

        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [cat_centrifuge + " && " + make_marker]
    
    def create_TA_taxon_pull_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        final_merge_folder      = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        ga_taxa_folder          = os.path.join(data_folder, "0_gene_taxa")
        jobs_folder             = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(ga_taxa_folder)
        self.make_folder(jobs_folder)
        

        get_taxa_from_gene = ">&2 echo get taxa from gene | "
        get_taxa_from_gene += self.path_obj.Python + " "
        get_taxa_from_gene += self.path_obj.Annotated_taxid + " "  # SLOW STEP
        get_taxa_from_gene += os.path.join(final_merge_folder, "gene_map.tsv") + " "
        get_taxa_from_gene += self.path_obj.accession2taxid + " "
        get_taxa_from_gene += os.path.join(ga_taxa_folder, "ga_taxon.tsv")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [get_taxa_from_gene + " && " + make_marker]
        
    def create_GA_pre_scan_wevote_combine_command(self, marker_path):

        self.make_folder(self.path_obj.GA_pre_scan_wevote_path)

        wevote_combine = ">&2 echo combining classification outputs for wevote | "
        wevote_combine += self.path_obj.Python + " "
        wevote_combine += self.path_obj.Classification_combine + " "
        wevote_combine += os.path.join(self.path_obj.contigs_final_path, "contig_map.tsv")
        wevote_combine += " " + os.path.join(self.path_obj.GA_pre_scan_wevote_path, "wevote_input.csv") + " "
        wevote_combine += "none" + " "
        wevote_combine += "none" + " "
        wevote_combine += "none" + " "
        wevote_combine += os.path.join(self.path_obj.GA_pre_scan_kraken_path, "merged_kraken2.txt") + " "
        wevote_combine += os.path.join(self.path_obj.GA_pre_scan_centr_path, "merged_centrifuge.tsv")  

        wevote_call = ">&2 echo Running WEVOTE | "
        wevote_call += self.path_obj.WEVOTE
        wevote_call += " -i " + os.path.join(self.path_obj.GA_pre_scan_wevote_path, "wevote_input.csv")
        wevote_call += " -d " + self.path_obj.WEVOTEDB
        wevote_call += " -p " + os.path.join(self.path_obj.GA_pre_scan_wevote_path, "wevote")
        wevote_call += " -n " + self.threads_str
        wevote_call += " -k " + "2"
        wevote_call += " -a " + "0"
        wevote_call += " -s " + "0"
        
        wevote_collect = ">&2 echo gathering WEVOTE results | "
        wevote_collect += self.path_obj.Python + " "
        wevote_collect += self.path_obj.Wevote_parser + " "
        wevote_collect += os.path.join(self.path_obj.GA_pre_scan_wevote_path, "wevote_WEVOTE_Details.txt") + " "
        wevote_collect += os.path.join(self.path_obj.GA_pre_scan_wevote_path, "taxonomic_classifications.tsv")
        
        make_marker = "touch" + " " + marker_path
        
        return [wevote_combine + " && " + wevote_call + " && " + wevote_collect +  " && " + make_marker]
        

    def create_TA_wevote_combine_command(self, current_stage_name, assemble_contigs_stage, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        assemble_contigs_folder = os.path.join(self.Output_Path, assemble_contigs_stage, "final_results")
        #kaiju_folder            = os.path.join(data_folder, "1_kaiju")
        kraken2_folder          = os.path.join(data_folder, "1_kraken2")
        centrifuge_folder       = os.path.join(data_folder, "2_centrifuge")
        wevote_folder           = os.path.join(data_folder, "3_wevote")
        final_folder            = os.path.join(subfolder, "final_results")
        jobs_folder             = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(wevote_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        
        wevote_combine = ">&2 echo combining classification outputs for wevote | "
        wevote_combine += self.path_obj.Python + " "
        wevote_combine += self.path_obj.Classification_combine + " "
        wevote_combine += os.path.join(assemble_contigs_folder, "contig_map.tsv")
        wevote_combine += " " + os.path.join(wevote_folder, "wevote_input.csv") + " "
        wevote_combine += "none" + " "
        wevote_combine += "none" + " "
        wevote_combine += "none" + " "
        wevote_combine += os.path.join(kraken2_folder, "merged_kraken2.txt") + " "
        wevote_combine += os.path.join(centrifuge_folder, "merged_centrifuge.tsv")  

        wevote_call = ">&2 echo Running WEVOTE | "
        wevote_call += self.path_obj.WEVOTE
        wevote_call += " -i " + os.path.join(wevote_folder, "wevote_input.csv")
        wevote_call += " -d " + self.path_obj.WEVOTEDB
        wevote_call += " -p " + os.path.join(wevote_folder, "wevote")
        wevote_call += " -n " + self.threads_str
        wevote_call += " -k " + "2"
        wevote_call += " -a " + "0"
        wevote_call += " -s " + "0"
        
        wevote_collect = ">&2 echo gathering WEVOTE results | "
        wevote_collect += self.path_obj.Python + " "
        wevote_collect += self.path_obj.Wevote_parser + " "
        wevote_collect += os.path.join(wevote_folder, "wevote_WEVOTE_Details.txt") + " "
        wevote_collect += os.path.join(wevote_folder, "taxonomic_classifications.tsv")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [wevote_combine + " && " + wevote_call + " && " + wevote_collect +  " && " + make_marker]
        
    
    def create_TA_final_command(self, current_stage_name, assemble_contigs_stage, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        assemble_contigs_folder = os.path.join(self.Output_Path, assemble_contigs_stage, "final_results")
        ga_taxa_folder          = os.path.join(data_folder, "0_gene_taxa")
        kraken2_folder            = os.path.join(data_folder, "1_kraken2")
        centrifuge_folder       = os.path.join(data_folder, "2_centrifuge")
        wevote_folder           = os.path.join(data_folder, "3_wevote")
        final_folder            = os.path.join(subfolder, "final_results")
        jobs_folder             = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(ga_taxa_folder)
        self.make_folder(kraken2_folder)
        self.make_folder(centrifuge_folder)
        self.make_folder(wevote_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        wevote_combine = ">&2 echo combining classification outputs for wevote | "
        wevote_combine += self.path_obj.Python + " "
        wevote_combine += self.path_obj.Classification_combine + " "
        wevote_combine += os.path.join(assemble_contigs_folder, "contig_map.tsv")
        wevote_combine += " " + os.path.join(wevote_folder, "wevote_ensemble.csv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(kraken2_folder, "merged_kraken2.txt") + " "
        wevote_combine += os.path.join(centrifuge_folder, "merged_centrifuge.tsv")        

        wevote_call = ">&2 echo Running WEVOTE | "
        wevote_call += self.path_obj.WEVOTE
        wevote_call += " -i " + os.path.join(wevote_folder, "wevote_ensemble.csv")
        wevote_call += " -d " + self.path_obj.WEVOTEDB
        wevote_call += " -p " + os.path.join(wevote_folder, "wevote")
        wevote_call += " -n " + self.threads_str
        wevote_call += " -k " + "2"
        wevote_call += " -a " + "0"
        wevote_call += " -s " + "0"
        
        wevote_collect = ">&2 echo gathering WEVOTE results | "
        wevote_collect += self.path_obj.Python + " "
        wevote_collect += self.path_obj.Wevote_parser + " "
        wevote_collect += os.path.join(wevote_folder, "wevote_WEVOTE_Details.txt") + " "
        wevote_collect += os.path.join(final_folder, "taxonomic_classifications.tsv")
        
        constrain = ">&2 echo Constraining the Taxonomic Annotation | " 
        constrain += self.path_obj.Python + " " + self.path_obj.Constrain_classification + " "
        constrain += self.path_obj.target_rank + " "
        constrain += os.path.join(final_folder, "taxonomic_classifications.tsv") + " "
        constrain += self.path_obj.nodes + " "
        constrain += self.path_obj.names + " "
        constrain += os.path.join(final_folder, "constrain_classification.tsv")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
   
        return [wevote_combine + " && " + wevote_call + " && " + wevote_collect + " && " + constrain + " && " + make_marker]
        
      

    def create_EC_DETECT_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        detect_folder       = os.path.join(data_folder, "0_detect")
        jobs_folder         = os.path.join(data_folder, "jobs")
    
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(detect_folder)
        self.make_folder(jobs_folder)
        
        
        detect_protein = ">&2 echo running detect on split file | "
        detect_protein += self.path_obj.Python + " "
        detect_protein += self.path_obj.Detect + " "
        detect_protein += os.path.join(final_merge_folder,"all_proteins.faa")
        detect_protein += " --output_file " + os.path.join(detect_folder, "proteins.detect")
        detect_protein += " --fbeta " + os.path.join(detect_folder, "proteins.fbeta")
        detect_protein += " --db " + self.path_obj.DetectDB
        detect_protein += " --blastp " + self.path_obj.Blastp
        detect_protein += " --needle " + self.path_obj.Needle
        detect_protein += " --dump_dir " + detect_folder 
        detect_protein += " --n_count" + " " + str(self.path_obj.DETECT_job_limit)
        detect_protein += " --mem_limit" + " " + str(self.path_obj.DETECT_mem_threshold) 
        detect_protein += " --job_delay" + " " + str(self.path_obj.DETECT_job_delay)
        detect_protein += " >> " + os.path.join(detect_folder, "detect_out.txt") + " 2>&1"

        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_DETECT = [
            detect_protein + " && " + make_marker
        ]

        return COMMANDS_DETECT

    def create_EC_PRIAM_split_command(self, current_stage_name, ga_final_merge_stage, split_folder, marker_file):
        #used to split the proteins file to run PRIAM
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(split_folder)
        self.make_folder(jobs_folder)

        split_command = self.path_obj.Python + " "
        split_command += self.path_obj.File_splitter + " "
        split_command += os.path.join(final_merge_folder, "all_proteins.faa") + " "
        split_command += os.path.join(split_folder, "protein_split") + " "
        split_command += str(self.path_obj.EC_chunksize)
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [split_command + " && " + make_marker]
        
    """
    def create_EC_PRIAM_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        #april 06, 2021: This one's a little tricky.  PRIAM has a user-prompt (and no args) to auto-resume.  
        #We must feed it the bash "Yes" in order to activate it.  So, mind the mess
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        PRIAM_folder        = os.path.join(data_folder, "1_priam")
        jobs_folder    = os.path.join(data_folder, "jobs")

        PRIAM_command = ">&2 echo running PRIAM | "
        
        if(os.path.exists(PRIAM_folder)):
            PRIAM_command += "yes | "
        else:
            self.make_folder(PRIAM_folder)
        self.make_folder(jobs_folder)
        
        
        PRIAM_command += self.path_obj.Java + " "
        PRIAM_command += self.path_obj.Priam
        PRIAM_command += " -n " + "proteins_priam" + " "
        PRIAM_command += " -i " + os.path.join(final_merge_folder, "all_proteins.faa")
        PRIAM_command += " -p " + self.path_obj.PriamDB
        PRIAM_command += " -o " + PRIAM_folder
        PRIAM_command += " --np " + self.threads_str
        PRIAM_command += " --bh --cc --cg --bp --bd "
        PRIAM_command += self.path_obj.BLAST_dir
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_PRIAM = [
            PRIAM_command + " && " + make_marker
        ]

        return COMMANDS_PRIAM
    """

    def create_EC_PRIAM_command_v2(self, current_stage_name, ga_final_merge_stage, priam_out_folder, split_file, id, marker_file):
        #april 06, 2021: This one's a little tricky.  PRIAM has a user-prompt (and no args) to auto-resume.  
        #We must feed it the bash "Yes" in order to activate it.  So, mind the mess
        #dec 05, 2022: now designed to run on split protein files
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        PRIAM_folder        = os.path.join(data_folder, "1_priam")
        jobs_folder    = os.path.join(data_folder, "jobs")

        PRIAM_command = ">&2 echo running PRIAM | "
        
        if(os.path.exists(PRIAM_folder)):
            PRIAM_command += "yes | "
        else:
            self.make_folder(PRIAM_folder)
        self.make_folder(jobs_folder)
        
        
        PRIAM_command += self.path_obj.Java + " "
        PRIAM_command += self.path_obj.Priam
        PRIAM_command += " -n " + "split_" + str(id) 
        PRIAM_command += " -i " + split_file
        PRIAM_command += " -p " + self.path_obj.PriamDB
        PRIAM_command += " -o " + priam_out_folder
        PRIAM_command += " --np " + self.threads_str
        PRIAM_command += " --bh --cc --bp --bd "
        PRIAM_command += self.path_obj.BLAST_dir
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_PRIAM = [
            PRIAM_command + " && " + make_marker
        ]

        return COMMANDS_PRIAM
        
    def create_EC_PRIAM_cat_command(self, current_stage_name, marker_file):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        PRIAM_folder        = os.path.join(data_folder, "1_priam")
        jobs_folder    = os.path.join(data_folder, "jobs")
        
        cat_command = "for i in $(ls " + PRIAM_folder + " | grep split); do cat " + PRIAM_folder + "/$i/PRIAM_$i/ANNOTATION/sequenceECs.txt >> " + PRIAM_folder + "/all_sequenceECs.txt; done"
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [cat_command + " && " + make_marker]
        
        
    def create_EC_DIAMOND_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        diamond_ea_folder   = os.path.join(data_folder, "2_diamond")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        self.make_folder(diamond_ea_folder)
        self.make_folder(jobs_folder)
        
        diamond_ea_command = ">&2 echo running Diamond enzyme annotation | "
        diamond_ea_command += self.path_obj.DIAMOND + " blastp"
        diamond_ea_command += " -p " + self.threads_str
        diamond_ea_command += " --query " + os.path.join(final_merge_folder, "all_proteins.faa")
        diamond_ea_command += " --db " + self.path_obj.SWISS_PROT
        diamond_ea_command += " --outfmt " + "6 qseqid sseqid length qstart qend sstart send evalue bitscore qcovhsp slen pident"
        diamond_ea_command += " --out " + os.path.join(diamond_ea_folder, "proteins.blastout")
        diamond_ea_command += " --evalue 0.0000000001"
        #diamond_ea_command += " --max-target-seqs 1"
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        COMMANDS_DIAMOND_EC = [
            diamond_ea_command + " && " + make_marker
        ]
        
        return COMMANDS_DIAMOND_EC
        
    def create_EC_postprocess_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        detect_folder       = os.path.join(data_folder, "0_detect")
        PRIAM_folder        = os.path.join(data_folder, "1_priam")
        diamond_ea_folder   = os.path.join(data_folder, "2_diamond")
        final_folder        = os.path.join(subfolder, "final_results")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        #combine_detect = "cat " + os.path.join(detect_folder, "protein_*.toppred")
        #combine_detect += " > " + os.path.join(detect_folder, "proteins.toppred")

        postprocess_command = ">&2 echo combining enzyme annotation output | "
        postprocess_command += self.path_obj.Python + " "
        postprocess_command += self.path_obj.EC_Annotation_Post + " "
        postprocess_command += os.path.join(detect_folder, "proteins.fbeta") + " "
        postprocess_command += os.path.join(PRIAM_folder, "all_sequenceECs.txt") + " "
        postprocess_command += os.path.join(diamond_ea_folder, "proteins.blastout") + " "
        postprocess_command += self.path_obj.SWISS_PROT_map + " "
        postprocess_command += os.path.join(final_merge_folder, "gene_map.tsv") + " "
        postprocess_command += self.path_obj.enzyme_db + " "
        postprocess_command += os.path.join(final_folder, "proteins.ECs_All") + " "
        postprocess_command += os.path.join(final_folder, "lq_proteins.ECs_All")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_EC_Postprocess = [
            #combine_detect,
            postprocess_command + " && " + make_marker
        ]

        return COMMANDS_EC_Postprocess

        

        
    def create_output_copy_gene_map_command(self, current_stage_name, ga_final_merge_stage):
        #just copies the gene map over to the output
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        ga_final_merge_folder   = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        final_folder            = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        gene_map_location = os.path.join(ga_final_merge_folder, "gene_map.tsv")
        
        copy_gene_map = ">&2 echo copying gene map | "
        copy_gene_map += "cp " + os.path.join(ga_final_merge_folder, "gene_map.tsv") + " "
        copy_gene_map += os.path.join(final_folder, "gene_map.tsv")
        
        return[copy_gene_map]

        
    def create_output_network_generation_command(self, current_stage_name, ga_final_merge_stage, taxonomic_annotation_stage, enzyme_annotation_stage):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        ga_final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        ta_folder           = os.path.join(self.Output_Path, taxonomic_annotation_stage, "final_results")
        ea_folder           = os.path.join(self.Output_Path, enzyme_annotation_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        #gene_map_location = os.path.join(final_folder, "final_gene_map.tsv")
        gene_map_location = os.path.join(ga_final_merge_folder, "gene_map.tsv")
        
        network_generation = ">&2 echo Generating RPKM and Cytoscape network | "
        network_generation += self.path_obj.Python + " "
        network_generation += self.path_obj.RPKM + " "
        network_generation += str(self.path_obj.RPKM_cutoff) + " "
        network_generation += "None" + " "
        network_generation += self.path_obj.nodes + " "
        network_generation += self.path_obj.names + " "
        network_generation += gene_map_location + " "
        network_generation += os.path.join(ta_folder, "taxonomic_classifications.tsv") + " "
        network_generation += os.path.join(ea_folder, "proteins.ECs_All") + " "
        network_generation += self.path_obj.show_unclassified + " "
        network_generation += os.path.join(final_folder, "RPKM_table.tsv") + " "
        network_generation += os.path.join(final_folder, "Cytoscape_network.tsv") + " "
        
        
        
        flatten_rpkm = ">&2 echo Reformat RPKM for EC heatmap | "
        flatten_rpkm += self.path_obj.Python + " "
        flatten_rpkm += self.path_obj.format_RPKM + " "
        flatten_rpkm += os.path.join(final_folder, "RPKM_table.tsv") + " "
        flatten_rpkm += os.path.join(final_folder, "EC_heatmap_RPKM.tsv")
        
        return [network_generation, flatten_rpkm]
        
    def create_output_unique_hosts_singletons_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path      = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder         = os.path.join(self.Output_Path, host_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        unique_hosts_folder = os.path.join(data_folder, "1_unique_hosts")
        full_hosts_folder   = os.path.join(data_folder, "2_full_hosts")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(data_folder)
        self.make_folder(unique_hosts_folder)
        self.make_folder(full_hosts_folder)
        self.make_folder(final_folder)
        
        
        get_unique_host_reads_singletons = ">&2 echo get singleton host reads for stats | "
        get_unique_host_reads_singletons += self.path_obj.Python + " "
        get_unique_host_reads_singletons += self.path_obj.get_unique_host_reads + " "
        get_unique_host_reads_singletons += os.path.join(host_folder, "singletons.fastq") + " "
        get_unique_host_reads_singletons += os.path.join(self.path_obj.qc_final_path, "singletons.fastq") + " "
        get_unique_host_reads_singletons += os.path.join(unique_hosts_folder, "singletons_hosts.fastq")
        
        
        repop_singletons_hosts = ">&2 echo repopulating singletons hosts | " 
        repop_singletons_hosts += self.path_obj.Python + " "
        repop_singletons_hosts += self.path_obj.duplicate_repopulate + " "
        if(self.read_mode == "single"):
            repop_singletons_hosts += os.path.join(self.path_obj.qc_final_path, "singletons_hq.fastq") + " "
        else:
            repop_singletons_hosts += os.path.join(self.path_obj.qc_final_path, "singletons_with_duplicates.fastq") + " "
        repop_singletons_hosts += os.path.join(unique_hosts_folder, "singletons_hosts.fastq") + " "
        repop_singletons_hosts += os.path.join(self.path_obj.qc_final_path, "singletons_unique.fastq.clstr") + " "
        repop_singletons_hosts += os.path.join(full_hosts_folder, "singletons_full_hosts.fastq")
        
        return [get_unique_host_reads_singletons, repop_singletons_hosts]
        
    def create_output_unique_hosts_pair_1_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path      = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder         = os.path.join(self.Output_Path, host_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        unique_hosts_folder = os.path.join(data_folder, "1_unique_hosts")
        full_hosts_folder   = os.path.join(data_folder, "2_full_hosts")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)

        self.make_folder(data_folder)
        self.make_folder(unique_hosts_folder)
        self.make_folder(full_hosts_folder)
        self.make_folder(final_folder)
        
        get_unique_host_reads_pair_1 = ">&2 echo get pair 1 host reads for stats | " 
        get_unique_host_reads_pair_1 += self.path_obj.Python + " "
        get_unique_host_reads_pair_1 += self.path_obj.get_unique_host_reads + " "
        get_unique_host_reads_pair_1 += os.path.join(host_folder, "pair_1.fastq") + " "
        get_unique_host_reads_pair_1 += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
        get_unique_host_reads_pair_1 += os.path.join(unique_hosts_folder, "pair_1_hosts.fastq")
        
        repop_pair_1_hosts = ">&2 echo repopulating pair 1 hosts | " 
        repop_pair_1_hosts += self.path_obj.Python + " "
        repop_pair_1_hosts += self.path_obj.duplicate_repopulate + " "
        repop_pair_1_hosts += os.path.join(self.path_obj.qc_final_path, "pair_1_match.fastq") + " "
        repop_pair_1_hosts += os.path.join(unique_hosts_folder, "pair_1_hosts.fastq") + " "
        repop_pair_1_hosts += os.path.join(self.path_obj.qc_final_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_hosts += os.path.join(full_hosts_folder, "pair_1_full_hosts.fastq")
        
        return [get_unique_host_reads_pair_1, repop_pair_1_hosts]
        
    def create_output_unique_hosts_pair_2_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path      = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder         = os.path.join(self.Output_Path, host_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        unique_hosts_folder = os.path.join(data_folder, "1_unique_hosts")
        full_hosts_folder   = os.path.join(data_folder, "2_full_hosts")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(data_folder)
        self.make_folder(full_hosts_folder)
        self.make_folder(final_folder)
        
        get_unique_host_reads_pair_2 = ">&2 echo get pair 2 host reads for stats | " 
        get_unique_host_reads_pair_2 += self.path_obj.Python + " "
        get_unique_host_reads_pair_2 += self.path_obj.get_unique_host_reads + " "
        get_unique_host_reads_pair_2 += os.path.join(host_folder, "pair_2.fastq") + " "
        get_unique_host_reads_pair_2 += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
        get_unique_host_reads_pair_2 += os.path.join(unique_hosts_folder, "pair_2_hosts.fastq")
        
        repop_pair_2_hosts = ">&2 echo repopulating pair 2 hosts | " 
        repop_pair_2_hosts += self.path_obj.Python + " "
        repop_pair_2_hosts += self.path_obj.duplicate_repopulate + " "
        repop_pair_2_hosts += os.path.join(self.path_obj.qc_final_path, "pair_2_match.fastq") + " "
        repop_pair_2_hosts += os.path.join(unique_hosts_folder, "pair_2_hosts.fastq") + " "
        repop_pair_2_hosts += os.path.join(self.path_obj.qc_final_path, "pair_1_unique.fastq.clstr") + " " #we do this based on pairs now
        repop_pair_2_hosts += os.path.join(full_hosts_folder, "pair_2_full_hosts.fastq")
        
        return [get_unique_host_reads_pair_2, repop_pair_2_hosts]
#-------------------------------------------------------------------------------------------
    def create_output_unique_vectors_singletons_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path          = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder             = os.path.join(self.Output_Path, host_stage, "final_results")
        vectors_folder          = os.path.join(self.Output_Path, vectors_stage, "final_results")
        data_folder             = os.path.join(subfolder, "data")
        unique_vectors_folder   = os.path.join(data_folder, "3_unique_vectors")
        full_vectors_folder     = os.path.join(data_folder, "4_full_vectors")
        final_folder            = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(unique_vectors_folder)
        self.make_folder(full_vectors_folder)
        self.make_folder(final_folder)
        
        get_unique_vectors_reads_singletons = ">&2 echo get singleton vectors reads for stats | "
        get_unique_vectors_reads_singletons += self.path_obj.Python + " "
        get_unique_vectors_reads_singletons += self.path_obj.get_unique_host_reads + " "
            
        
        if(self.no_host_flag):
            get_unique_vectors_reads_singletons += os.path.join(vectors_folder, "singletons.fastq") + " "
            get_unique_vectors_reads_singletons += os.path.join(self.path_obj.qc_final_path, "singletons.fastq") + " "
            get_unique_vectors_reads_singletons += os.path.join(unique_vectors_folder, "singletons_vectors.fastq")
            
        else:
        
            get_unique_vectors_reads_singletons += os.path.join(vectors_folder, "singletons.fastq") + " "
            get_unique_vectors_reads_singletons += os.path.join(host_folder, "singletons.fastq") + " "
            get_unique_vectors_reads_singletons += os.path.join(unique_vectors_folder, "singletons_vectors.fastq")
            
            
        repop_singletons_vectors = ">&2 echo repopulating singletons vectors | " 
        repop_singletons_vectors += self.path_obj.Python + " "
        repop_singletons_vectors += self.path_obj.duplicate_repopulate + " "
        if(self.read_mode == "single"):
            repop_singletons_vectors += os.path.join(self.path_obj.qc_final_path, "singletons_hq.fastq") + " "
        else:
            repop_singletons_vectors += os.path.join(self.path_obj.qc_final_path, "singletons_with_duplicates.fastq") + " "
        repop_singletons_vectors += os.path.join(unique_vectors_folder, "singletons_vectors.fastq") + " "
        repop_singletons_vectors += os.path.join(self.path_obj.qc_final_path, "singletons_unique.fastq.clstr") + " "
        repop_singletons_vectors += os.path.join(full_vectors_folder, "singletons_full_vectors.fastq")
        
        return [get_unique_vectors_reads_singletons, repop_singletons_vectors]
        
    def create_output_unique_vectors_pair_1_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path          = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder             = os.path.join(self.Output_Path, host_stage, "final_results")
        vectors_folder          = os.path.join(self.Output_Path, vectors_stage, "final_results")
        data_folder             = os.path.join(subfolder, "data")
        unique_vectors_folder   = os.path.join(data_folder, "3_unique_vectors")
        full_vectors_folder     = os.path.join(data_folder, "4_full_vectors")
        final_folder            = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(unique_vectors_folder)
        self.make_folder(full_vectors_folder)
        self.make_folder(final_folder)
        
        get_unique_vectors_reads_pair_1 = ">&2 echo get pair 1 vector reads for stats | " 
        get_unique_vectors_reads_pair_1 += self.path_obj.Python + " "
        get_unique_vectors_reads_pair_1 += self.path_obj.get_unique_host_reads + " "
        
        if(self.no_host_flag):
            get_unique_vectors_reads_pair_1 += os.path.join(vectors_folder, "pair_1.fastq") + " "
            get_unique_vectors_reads_pair_1 += os.path.join(self.path_obj.qc_final_path, "pair_1.fastq") + " "
            get_unique_vectors_reads_pair_1 += os.path.join(unique_vectors_folder, "pair_1_vectors.fastq")

        else:
            get_unique_vectors_reads_pair_1 += os.path.join(vectors_folder, "pair_1.fastq") + " "
            get_unique_vectors_reads_pair_1 += os.path.join(host_folder, "pair_1.fastq") + " "
            get_unique_vectors_reads_pair_1 += os.path.join(unique_vectors_folder, "pair_1_vectors.fastq")
        
        repop_pair_1_vectors = ">&2 echo repopulating pair 1 vectors | " 
        repop_pair_1_vectors += self.path_obj.Python + " "
        repop_pair_1_vectors += self.path_obj.duplicate_repopulate + " "
        repop_pair_1_vectors += os.path.join(self.path_obj.qc_final_path, "pair_1_match.fastq") + " "
        repop_pair_1_vectors += os.path.join(unique_vectors_folder, "pair_1_vectors.fastq") + " "
        repop_pair_1_vectors += os.path.join(self.path_obj.qc_final_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_vectors += os.path.join(full_vectors_folder, "pair_1_full_vectors.fastq")
        
        return [get_unique_vectors_reads_pair_1, repop_pair_1_vectors]
        
    def create_output_unique_vectors_pair_2_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path          = os.path.join(self.Output_Path, quality_stage, "final_results")
        host_folder             = os.path.join(self.Output_Path, host_stage, "final_results")
        vectors_folder          = os.path.join(self.Output_Path, vectors_stage, "final_results")
        data_folder             = os.path.join(subfolder, "data")
        unique_vectors_folder   = os.path.join(data_folder, "3_unique_vectors")
        full_vectors_folder     = os.path.join(data_folder, "4_full_vectors")
        final_folder            = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(unique_vectors_folder)
        self.make_folder(full_vectors_folder)
        self.make_folder(final_folder)
        
        get_unique_vectors_reads_pair_2 = ">&2 echo get pair 2 vector reads for stats | " 
        get_unique_vectors_reads_pair_2 += self.path_obj.Python + " "
        get_unique_vectors_reads_pair_2 += self.path_obj.get_unique_host_reads + " "
        
        if(self.no_host_flag):
            get_unique_vectors_reads_pair_2 += os.path.join(vectors_folder, "pair_2.fastq") + " "
            get_unique_vectors_reads_pair_2 += os.path.join(self.path_obj.qc_final_path, "pair_2.fastq") + " "
            get_unique_vectors_reads_pair_2 += os.path.join(unique_vectors_folder, "pair_2_vectors.fastq")
        else:
            get_unique_vectors_reads_pair_2 += os.path.join(vectors_folder, "pair_2.fastq") + " "
            get_unique_vectors_reads_pair_2 += os.path.join(host_folder, "pair_2.fastq") + " "
            get_unique_vectors_reads_pair_2 += os.path.join(unique_vectors_folder, "pair_2_vectors.fastq")
            
        repop_pair_2_vectors = ">&2 echo repopulating pair 2 vectors | " 
        repop_pair_2_vectors += self.path_obj.Python + " "
        repop_pair_2_vectors += self.path_obj.duplicate_repopulate + " "
        repop_pair_2_vectors += os.path.join(self.path_obj.qc_final_path, "pair_2_match.fastq") + " "
        repop_pair_2_vectors += os.path.join(unique_vectors_folder, "pair_2_vectors.fastq") + " "
        repop_pair_2_vectors += os.path.join(self.path_obj.qc_final_path, "pair_1_unique.fastq.clstr") + " " #we do this based on pairs now
        repop_pair_2_vectors += os.path.join(full_vectors_folder, "pair_2_full_vectors.fastq")
        
        return [get_unique_vectors_reads_pair_2, repop_pair_2_vectors]
        

        
    def create_output_per_read_scores_command(self, current_stage_name, quality_stage):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path      = os.path.join(self.Output_Path, quality_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        per_read_scores = ">&2 echo collecting per-read quality | " 
        per_read_scores += self.path_obj.Python + " "
        per_read_scores += self.path_obj.read_quality_metrics + " "
        if(self.read_mode == "single"):
            per_read_scores += "single" + " "
            per_read_scores += self.sequence_single + " "
            per_read_scores += os.path.join(self.path_obj.qc_final_path, "singletons_hq.fastq") + " "
            per_read_scores += os.path.join(final_folder)
            
        elif(self.read_mode == "paired"):
            per_read_scores += "paired" + " " 
            per_read_scores += self.sequence_path_1 + " "
            per_read_scores += self.sequence_path_2 + " "
            per_read_scores += os.path.join(self.path_obj.qc_final_path, "pair_1_match.fastq") + " "
            per_read_scores += os.path.join(self.path_obj.qc_final_path, "pair_2_match.fastq") + " "
            per_read_scores += os.path.join(self.path_obj.qc_final_path, "singletons_with_duplicates.fastq") + " "
            per_read_scores += os.path.join(final_folder)
            
        return [per_read_scores]
        
    def create_output_copy_taxa_command(self, current_stage_name, taxa_stage):
        subfolder       = os.path.join(self.Output_Path, current_stage_name)
        taxa_folder     = os.path.join(self.Output_Path, taxa_stage, "final_results")
        final_folder    = os.path.join(subfolder, "final_results")
        
        self.make_folder(subfolder)
        self.make_folder(final_folder)
        

        copy_taxa = ">&2 echo " + str(dt.today()) + " copying taxa data | " 
        copy_taxa += "cp" + " "
        copy_taxa += os.path.join(taxa_folder, "constrain_classification.tsv") + " "
        copy_taxa += os.path.join(final_folder, "taxa_classifications.tsv")
        
        return [copy_taxa]
        
        
    def create_output_contig_stats_command(self, current_stage_name, contig_stage):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        contig_folder       = os.path.join(self.Output_Path, contig_stage, "final_results")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        contig_stats = ">&2 echo " + str(dt.today()) + " collecting contig stats | " 
        contig_stats += self.path_obj.Python + " "
        contig_stats += self.path_obj.contig_stats + " "
        contig_stats += os.path.join(contig_folder, "contigs.fasta") + " "
        contig_stats += os.path.join(final_folder, "contig_stats.txt")
        
        return [contig_stats]
        
    def create_output_EC_heatmap_command(self, current_stage_name):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        EC_heatmap = ">&2 echo " + str(dt.today()) + " forming EC heatmap | "
        EC_heatmap += self.path_obj.Python + " "
        EC_heatmap += self.path_obj.ec_heatmap + " "
        EC_heatmap += self.path_obj.EC_pathway + " "
        EC_heatmap += os.path.join(final_folder, "EC_heatmap_RPKM.tsv") + " "
        EC_heatmap += self.path_obj.path_to_superpath + " "
        EC_heatmap += final_folder
        
        return [EC_heatmap]
        
        
        
    def create_output_read_count_command(self, current_stage_name, quality_stage, repopulation_stage, ga_final_merge_stage, enzyme_annotation_stage):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        self.path_obj.qc_final_path      = os.path.join(self.Output_Path, quality_stage, "final_results")
        repopulation_folder = os.path.join(self.Output_Path, repopulation_stage, "final_results")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        ea_folder           = os.path.join(self.Output_Path, enzyme_annotation_stage, "final_results")
        full_hosts_folder   = os.path.join(data_folder, "2_full_hosts")
        full_vectors_folder = os.path.join(data_folder, "4_full_vectors")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(full_hosts_folder)
        self.make_folder(full_vectors_folder)
        self.make_folder(final_folder)
        gene_map_location = os.path.join(final_folder, "gene_map.tsv")
        
        read_counts = ">&2 echo " + str(dt.today()) + " generating read count table | "
        read_counts += self.path_obj.Python + " "
        read_counts += self.path_obj.read_count + " "
        if self.read_mode == "single":
            read_counts += self.sequence_single + " "
            
        elif self.read_mode == "paired":
            read_counts += self.sequence_path_1 + " "
        read_counts += self.path_obj.qc_final_path + " "
        read_counts += full_hosts_folder + " "
        read_counts += full_vectors_folder + " "
        read_counts += repopulation_folder + " "
        read_counts += final_merge_folder + " "
        read_counts += ea_folder + " "
        read_counts += os.path.join(final_folder, "read_count.tsv") + " "
        read_counts += self.read_mode
        
        return [read_counts]
        
        
    def create_output_taxa_groupby_command(self, current_stage_name):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_folder        = os.path.join(subfolder, "final_results")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        taxa_groupby = ">&2 echo making Taxa summary | " 
        taxa_groupby += self.path_obj.Python + " "
        taxa_groupby += self.path_obj.taxa_table + " "
        taxa_groupby += os.path.join(final_folder, "taxa_classifications.tsv") + " "
        taxa_groupby += os.path.join(final_folder, "taxa_summary.tsv")
        
        return [taxa_groupby]
