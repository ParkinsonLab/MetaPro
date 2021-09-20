# The functions here generate the pipeline commands.
# The functions here generate the pipeline commands.
# Each command module is made up of sub stages that are used to get the final result.

import os
import os.path
import MetaPro_paths as mpp
import subprocess as sp
import time
import sys
from datetime import datetime as dt

class mt_pipe_commands:
    # --------------------------------------------------------------------
    # constructor:
    # there should only be one of these objects used for an entire pipeline.
    def __init__(self, no_host, Config_path, Quality_score=33, Thread_count=8, tutorial_keyword = None, sequence_path_1=None, sequence_path_2=None, sequence_single=None, sequence_contigs = None):

        self.tool_path_obj = mpp.tool_path_obj(Config_path)
        self.no_host_flag = no_host
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
        self.Threads_str = str(Thread_count)
        
        
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
                
                
    def create_quality_control_command(self, stage_name):
        subfolder                   = os.path.join(self.Output_Path, stage_name)
        data_folder                 = os.path.join(subfolder, "data")
        sorted_read_folder          = os.path.join(data_folder, "0_sorted_raw_input")
        adapter_folder              = os.path.join(data_folder, "1_adapter_removal")
        tag_remove_folder           = os.path.join(data_folder, "2_tag_remove")
        vsearch_merge_folder        = os.path.join(data_folder, "3_vsearch_pair_merge")
        vsearch_filter_folder       = os.path.join(data_folder, "4_quality_filter")
        orphan_read_filter_folder   = os.path.join(data_folder, "5_orphan_read_filter")
        cdhit_folder                = os.path.join(data_folder, "6_remove_duplicates")
        final_folder                = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(sorted_read_folder)
        self.make_folder(adapter_folder)
        self.make_folder(tag_remove_folder)
        self.make_folder(vsearch_merge_folder)
        self.make_folder(vsearch_filter_folder)
        self.make_folder(orphan_read_filter_folder)
        self.make_folder(cdhit_folder)
        self.make_folder(final_folder)
        
        sort_pair_1 = ">&2 echo Sorting pair 1 | "
        sort_pair_1 += self.tool_path_obj.Python + " "
        sort_pair_1 += self.tool_path_obj.sort_reads + " "
        sort_pair_1 += self.sequence_path_1 + " "
        sort_pair_1 += os.path.join(sorted_read_folder, "pair_1_sorted.fastq") + " "
        sort_pair_1 += "forward"

        sort_pair_2 = ">&2 echo Sorting pair 2 | "
        sort_pair_2 += self.tool_path_obj.Python + " "
        sort_pair_2 += self.tool_path_obj.sort_reads + " "
        sort_pair_2 += self.sequence_path_2 + " "
        sort_pair_2 += os.path.join(sorted_read_folder, "pair_2_sorted.fastq") + " "
        sort_pair_2 += "reverse"

        adapter_removal_line = ">&2 echo Removing adapters | "
        adapter_removal_line += self.tool_path_obj.AdapterRemoval
        if self.read_mode == "single":
            adapter_removal_line += " --file1 " + self.sequence_single
        elif self.read_mode == "paired":
            adapter_removal_line += " --file1 " + os.path.join(sorted_read_folder, "pair_1_sorted.fastq")
            adapter_removal_line += " --file2 " + os.path.join(sorted_read_folder, "pair_2_sorted.fastq")
        adapter_removal_line += " --qualitybase " + str(self.Qual_str)
        if(self.Qual_str == "33"):
            adapter_removal_line += " --qualitymax 75"
        adapter_removal_line += " --threads " + self.Threads_str
        adapter_removal_line += " --minlength " + str(self.tool_path_obj.adapterremoval_minlength)
        adapter_removal_line += " --basename " + adapter_folder
        adapter_removal_line += "_AdapterRemoval"
        adapter_removal_line += " --trimqualities "
        if self.read_mode == "single":
            adapter_removal_line += " --output1 " + os.path.join(adapter_folder, "singletons_adptr_rem.fastq")
        elif self.read_mode == "paired":
            adapter_removal_line += " --output1 " + os.path.join(adapter_folder, "pair_1_adptr_rem.fastq")
            adapter_removal_line += " --output2 " + os.path.join(adapter_folder, "pair_2_adptr_rem.fastq")
            adapter_removal_line += " --singleton " + os.path.join(adapter_folder, "singletons_adptr_rem.fastq")

        #Sort-reads introduces tags at the read-level of the 
        tag_remove_pair_1 = ">&2 echo Remove tags pair 1 | "
        tag_remove_pair_1 += self.tool_path_obj.Python + " "
        tag_remove_pair_1 += self.tool_path_obj.remove_tag + " "
        tag_remove_pair_1 += os.path.join(adapter_folder, "pair_1_adptr_rem.fastq") + " "
        tag_remove_pair_1 += os.path.join(tag_remove_folder, "pair_1_no_tags.fastq")
        
        tag_remove_pair_2 = ">&2 echo Remove tags pair 2 | "
        tag_remove_pair_2 += self.tool_path_obj.Python + " "
        tag_remove_pair_2 += self.tool_path_obj.remove_tag + " "
        tag_remove_pair_2 += os.path.join(adapter_folder, "pair_2_adptr_rem.fastq") + " "
        tag_remove_pair_2 += os.path.join(tag_remove_folder, "pair_2_no_tags.fastq")

        tag_remove_singletons =  ">&2 echo Remove tags singletons | " 
        tag_remove_singletons += self.tool_path_obj.Python + " "
        tag_remove_singletons += self.tool_path_obj.remove_tag + " "
        tag_remove_singletons += os.path.join(adapter_folder, "singletons_adptr_rem.fastq") + " "
        tag_remove_singletons += os.path.join(tag_remove_folder, "singletons_no_tags.fastq")
        # tries to merge the cleaned pairs
        # rejects get sent out
        vsearch_merge = ">&2 echo " + "Vsearch Merge pairs | "
        vsearch_merge += self.tool_path_obj.vsearch
        vsearch_merge += " --fastq_mergepairs " + os.path.join(tag_remove_folder, "pair_1_no_tags.fastq")
        vsearch_merge += " --reverse " + os.path.join(tag_remove_folder, "pair_2_no_tags.fastq")
        vsearch_merge += " --fastq_ascii " + str(self.Qual_str)
        vsearch_merge += " --fastqout " + os.path.join(vsearch_merge_folder, "merge_success.fastq")
        vsearch_merge += " --fastqout_notmerged_fwd " + os.path.join(vsearch_merge_folder, "pair_1_merge_reject.fastq")
        vsearch_merge += " --fastqout_notmerged_rev " + os.path.join(vsearch_merge_folder, "pair_2_merge_reject.fastq")

        # concatenate the merge overlaps with the singletons
        cat_glue = ">&2 echo concatenating singletons | "
        cat_glue += "cat "
        cat_glue += os.path.join(vsearch_merge_folder, "merge_success.fastq") + " "
        cat_glue += os.path.join(tag_remove_folder, "singletons_no_tags.fastq")
        cat_glue += " > " + os.path.join(vsearch_merge_folder, "singletons.fastq")

        # Filter out low-quality reads
        # start with the singles / merged sections
        
        vsearch_filter_0 = ">&2 echo low-quality filter on singletons | "
        vsearch_filter_0 += self.tool_path_obj.vsearch
        if self.read_mode == "single":
            vsearch_filter_0 += " --fastq_filter " + os.path.join(adapter_folder, "singletons_adptr_rem.fastq")
        elif self.read_mode == "paired":
            vsearch_filter_0 += " --fastq_filter " + os.path.join(vsearch_merge_folder, "singletons.fastq")
        vsearch_filter_0 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_0 += " --fastq_maxee " + "2.0"
        vsearch_filter_0 += " --fastqout " + os.path.join(vsearch_filter_folder, "singletons_hq.fastq")

        # then move onto the standalones in pair 1
        vsearch_filter_1 = ">&2 echo low-quality filter on pair 1 | "
        vsearch_filter_1 += self.tool_path_obj.vsearch
        vsearch_filter_1 += " --fastq_filter " + os.path.join(vsearch_merge_folder, "pair_1_merge_reject.fastq")
        vsearch_filter_1 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_1 += " --fastq_maxee " + "2.0"
        vsearch_filter_1 += " --fastqout " + os.path.join(vsearch_filter_folder, "pair_1_hq.fastq")

        vsearch_filter_2 = ">&2 echo low-quality filter on pair 2 | "
        vsearch_filter_2 += self.tool_path_obj.vsearch
        vsearch_filter_2 += " --fastq_filter " + os.path.join(vsearch_merge_folder, "pair_2_merge_reject.fastq")
        vsearch_filter_2 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_2 += " --fastq_maxee " + "2.0"
        vsearch_filter_2 += " --fastqout " + os.path.join(vsearch_filter_folder, "pair_2_hq.fastq")

        # redistribute data into singletons, or paired-reads
        orphan_read_filter = ">&2 echo moving newly orphaned reads | "
        orphan_read_filter += self.tool_path_obj.Python + " "
        orphan_read_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_read_filter += os.path.join(vsearch_filter_folder, "pair_1_hq.fastq") + " "
        orphan_read_filter += os.path.join(vsearch_filter_folder, "pair_2_hq.fastq") + " "
        orphan_read_filter += os.path.join(vsearch_filter_folder, "singletons_hq.fastq") + " "
        orphan_read_filter += os.path.join(orphan_read_filter_folder, "pair_1_match.fastq") + " "
        orphan_read_filter += os.path.join(orphan_read_filter_folder, "pair_2_match.fastq") + " "
        orphan_read_filter += os.path.join(orphan_read_filter_folder, "singletons_with_duplicates.fastq")

        # remove duplicates (to shrink the data size)
        cdhit_singletons = ">&2 echo removing singleton duplicates | "
        cdhit_singletons += self.tool_path_obj.cdhit_dup + " -i "
        if self.read_mode == "single":
            cdhit_singletons += os.path.join(vsearch_filter_folder, "singletons_hq.fastq")
        elif self.read_mode == "paired":
            cdhit_singletons += os.path.join(orphan_read_filter_folder, "singletons_with_duplicates.fastq")
        cdhit_singletons += " -o " + os.path.join(cdhit_folder, "singletons_unique.fastq")

        # remove duplicates in the pairs
        cdhit_paired = ">&2 echo remove duplicates from paired | "
        cdhit_paired += self.tool_path_obj.cdhit_dup + " "
        cdhit_paired += "-i"    + " " + os.path.join(orphan_read_filter_folder, "pair_1_match.fastq") + " "
        cdhit_paired += "-i2"   + " " + os.path.join(orphan_read_filter_folder, "pair_2_match.fastq") + " "
        cdhit_paired += "-o"    + " " + os.path.join(cdhit_folder, "pair_1_unique.fastq") + " "
        cdhit_paired += "-o2"   + " " + os.path.join(cdhit_folder, "pair_2_unique.fastq")

        #move data to appropriate places
        copy_singletons = "cp " + os.path.join(cdhit_folder, "singletons_unique.fastq") + " "
        copy_singletons += os.path.join(final_folder, "singletons.fastq")

        copy_pair_1 = "cp " + os.path.join(cdhit_folder, "pair_1_unique.fastq") + " "
        copy_pair_1 += os.path.join(final_folder, "pair_1.fastq")

        copy_pair_2 = "cp " + os.path.join(cdhit_folder, "pair_2_unique.fastq") + " "
        copy_pair_2 += os.path.join(final_folder, "pair_2.fastq")
        
        # move these particular files to final_folder because they'll be needed by another stage.
        copy_duplicate_singletons = "cp "
        if(self.read_mode == "single"):
            copy_duplicate_singletons += os.path.join(vsearch_filter_folder, "singletons_hq.fastq") + " "
            copy_duplicate_singletons += os.path.join(final_folder, "singletons_hq.fastq")
        else:
            copy_duplicate_singletons += os.path.join(orphan_read_filter_folder, "singletons_with_duplicates.fastq") + " "
            copy_duplicate_singletons += os.path.join(final_folder, "singletons_with_duplicates.fastq")

        copy_pair_1_match = "cp " + os.path.join(orphan_read_filter_folder, "pair_1_match.fastq") + " "
        copy_pair_1_match += os.path.join(final_folder, "pair_1_match.fastq")

        copy_pair_2_match = "cp " + os.path.join(orphan_read_filter_folder, "pair_2_match.fastq") + " "
        copy_pair_2_match += os.path.join(final_folder, "pair_2_match.fastq")

        copy_singletons_cluster = "cp " + os.path.join(cdhit_folder, "singletons_unique.fastq.clstr") + " "
        copy_singletons_cluster += os.path.join(final_folder, "singletons_unique.fastq.clstr")

        copy_paired_cluster = "cp " + os.path.join(cdhit_folder, "pair_1_unique.fastq.clstr") + " "
        copy_paired_cluster += os.path.join(final_folder, "pair_1_unique.fastq.clstr")

        if self.read_mode == "single":
            COMMANDS_qual = [
                adapter_removal_line,
                vsearch_filter_0,
                cdhit_singletons,
                copy_singletons,
                copy_duplicate_singletons,
                copy_singletons_cluster
            ]
        elif self.read_mode == "paired":
            COMMANDS_qual = [
                sort_pair_1,
                sort_pair_2,
                adapter_removal_line,
                tag_remove_pair_1,
                tag_remove_pair_2,
                tag_remove_singletons,
                vsearch_merge,
                cat_glue,
                vsearch_filter_0,
                vsearch_filter_1,
                vsearch_filter_2,
                orphan_read_filter,
                cdhit_singletons,
                cdhit_paired,
                copy_singletons,
                copy_pair_1,
                copy_pair_2,
                copy_duplicate_singletons,
                copy_singletons_cluster,
                copy_pair_1_match,
                copy_paired_cluster,
                copy_pair_2_match
            ]

        return COMMANDS_qual

    def create_host_filter_command(self, stage_name, dependency_name):
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, dependency_name, "final_results")
        host_removal_folder = os.path.join(data_folder, "0_remove_host")
        blat_hr_folder      = os.path.join(data_folder, "1_blat_host")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(host_removal_folder)
        self.make_folder(blat_hr_folder)
        self.make_folder(final_folder)

        Host_Contaminants = os.path.join(host_removal_folder, "host_contaminents_seq.fasta")
        copy_host = ">&2 echo Copy the host file over | "
        copy_host += "cp " + self.tool_path_obj.Host + " " + Host_Contaminants

        # craft a BWA index for the host sequences
        bwa_hr_prep = ">&2 echo make host contaminants index for BWA | "
        bwa_hr_prep += self.tool_path_obj.BWA + " index -a bwtsw " + Host_Contaminants

        samtools_hr_prep = ">&2 echo SAMTOOLS host contaminant prep | "
        samtools_hr_prep += self.tool_path_obj.SAMTOOLS + " faidx " + Host_Contaminants

        # host removal on unique singletons
        bwa_hr_singletons = ">&2 echo BWA host remove on singletons | "
        bwa_hr_singletons += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_singletons += self.Threads_str + " "
        bwa_hr_singletons += Host_Contaminants + " "
        bwa_hr_singletons += os.path.join(quality_folder, "singletons.fastq") + " " 
        bwa_hr_singletons += ">" + " "
        bwa_hr_singletons += os.path.join(host_removal_folder, "singletons_no_host.sam")
        
        #Tutorial-use only.  
        bwa_hr_tut_singletons = ">&2 echo BWA host remove on singletons | "
        bwa_hr_tut_singletons += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_tut_singletons += self.Threads_str + " "
        bwa_hr_tut_singletons += Host_Contaminants + " "
        bwa_hr_tut_singletons += self.sequence_single 
        bwa_hr_tut_singletons += " > " + os.path.join(host_removal_folder, "singletons_no_host.sam")

        # annoying type conversion pt 1
        samtools_hr_singletons_sam_to_bam = ">&2 echo convert singletons host reads | "
        samtools_hr_singletons_sam_to_bam += self.tool_path_obj.SAMTOOLS
        samtools_hr_singletons_sam_to_bam += " view -bS " + os.path.join(host_removal_folder, "singletons_no_host.sam")
        samtools_hr_singletons_sam_to_bam += " > " + os.path.join(host_removal_folder, "singletons_no_host.bam")
        # annoying type conversion pt 2
        samtools_no_host_singletons_bam_to_fastq = self.tool_path_obj.SAMTOOLS + " fastq -n -f 4" + " -0 "
        samtools_no_host_singletons_bam_to_fastq += os.path.join(host_removal_folder, "singletons_no_host.fastq") + " "
        samtools_no_host_singletons_bam_to_fastq += os.path.join(host_removal_folder, "singletons_no_host.bam")

        # apparently, we're to keep the host separation
        samtools_host_singletons_bam_to_fastq = self.tool_path_obj.SAMTOOLS + " fastq -n -F 4" + " -0 "
        samtools_host_singletons_bam_to_fastq += os.path.join(host_removal_folder, "singletons_host_only.fastq") + " "
        samtools_host_singletons_bam_to_fastq += os.path.join(host_removal_folder, "singletons_no_host.bam")


        
        bwa_hr_paired = ">&2 echo bwa host-removal on paired | " 
        bwa_hr_paired += self.tool_path_obj.BWA + " "
        bwa_hr_paired += "mem" + " "  + "-t" + " " + self.Threads_str + " "
        bwa_hr_paired += Host_Contaminants + " "
        bwa_hr_paired += os.path.join(quality_folder, "pair_1.fastq") + " "
        bwa_hr_paired += os.path.join(quality_folder, "pair_2.fastq") + " "
        bwa_hr_paired += ">" + " "
        bwa_hr_paired += os.path.join(host_removal_folder, "paired_on_host.sam")
        
        #Tutorial-use only
        bwa_hr_tut_paired = ">&2 echo bwa host-removal on paired | " 
        bwa_hr_tut_paired += self.tool_path_obj.BWA + " "
        bwa_hr_tut_paired += "mem" + " "  + "-t" + " " + self.Threads_str + " "
        bwa_hr_tut_paired += Host_Contaminants + " "
        bwa_hr_tut_paired += self.sequence_path_1 + " "
        bwa_hr_tut_paired += self.sequence_path_2 + " "
        bwa_hr_tut_paired += ">" + " "
        bwa_hr_tut_paired += os.path.join(host_removal_folder, "paired_on_host.sam")
        
        
        bwa_hr_filter_paired = ">&2 echo BWA host-removal PP on paired | "
        bwa_hr_filter_paired += self.tool_path_obj.Python + " "
        bwa_hr_filter_paired += self.tool_path_obj.bwa_read_sorter + " "
        bwa_hr_filter_paired += "paired" + " "
        bwa_hr_filter_paired += self.tool_path_obj.filter_stringency + " "
        bwa_hr_filter_paired += os.path.join(host_removal_folder, "paired_on_host.sam") + " "
        bwa_hr_filter_paired += os.path.join(quality_folder, "pair_1.fastq") + " "
        bwa_hr_filter_paired += os.path.join(quality_folder, "pair_2.fastq") + " "
        bwa_hr_filter_paired += os.path.join(host_removal_folder, "pair_1_no_host.fastq") + " "
        bwa_hr_filter_paired += os.path.join(host_removal_folder, "pair_2_no_host.fastq") + " "
        bwa_hr_filter_paired += os.path.join(host_removal_folder, "pair_1_host_only.fastq") + " "
        bwa_hr_filter_paired += os.path.join(host_removal_folder, "pair_2_host_only.fastq")

        # blat prep
        make_blast_db_host = ">&2 echo Make BLAST db for host contaminants | "
        make_blast_db_host += self.tool_path_obj.Makeblastdb + " -in " + Host_Contaminants + " -dbtype nucl"

        vsearch_filter_3 = ">&2 echo Convert singletons for BLAT | "
        vsearch_filter_3 += self.tool_path_obj.vsearch
        vsearch_filter_3 += " --fastq_filter " + os.path.join(host_removal_folder, "singletons_no_host.fastq")
        vsearch_filter_3 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_3 += " --fastaout " + os.path.join(host_removal_folder, "singletons_no_host.fasta")

        vsearch_filter_4 = ">&2 echo Convert pair 1 for BLAT | "
        vsearch_filter_4 += self.tool_path_obj.vsearch
        vsearch_filter_4 += " --fastq_filter " + os.path.join(host_removal_folder, "pair_1_no_host.fastq")
        vsearch_filter_4 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_4 += " --fastaout " + os.path.join(host_removal_folder, "pair_1_no_host.fasta")

        vsearch_filter_5 = ">&2 echo Convert pair 2 for BLAT | "
        vsearch_filter_5 += self.tool_path_obj.vsearch
        vsearch_filter_5 += " --fastq_filter " + os.path.join(host_removal_folder, "pair_2_no_host.fastq")
        vsearch_filter_5 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_5 += " --fastaout " + os.path.join(host_removal_folder, "pair_2_no_host.fasta")

        blat_hr_singletons = ">&2 echo BLAT host singletons | "
        blat_hr_singletons += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_hr_singletons += Host_Contaminants + " "
        blat_hr_singletons += os.path.join(host_removal_folder, "singletons_no_host.fasta")
        blat_hr_singletons += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_singletons += " " + os.path.join(host_removal_folder, "singletons_no_host.blatout")

        blat_hr_pair_1 = ">&2 echo BLAT host pair 1 | "
        blat_hr_pair_1 += self.tool_path_obj.BLAT
        blat_hr_pair_1 += " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " "
        blat_hr_pair_1 += os.path.join(host_removal_folder, "pair_1_no_host.fasta")
        blat_hr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_1 += " " + os.path.join(host_removal_folder, "pair_1_no_host.blatout")

        blat_hr_pair_2 = ">&2 echo BLAT host pair 2 | "
        blat_hr_pair_2 += self.tool_path_obj.BLAT
        blat_hr_pair_2 += " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " "
        blat_hr_pair_2 += os.path.join(host_removal_folder, "pair_2_no_host.fasta")
        blat_hr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_2 += " " + os.path.join(host_removal_folder, "pair_2_no_host.blatout")

        # HR BLAT
        hr_singletons = ">&2 echo BLAT contaminant singletons | "
        hr_singletons += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        hr_singletons += "single" +  " "
        hr_singletons += self.tool_path_obj.filter_stringency + " "
        hr_singletons += os.path.join(host_removal_folder, "singletons_no_host.fastq") + " "  # in
        hr_singletons += os.path.join(host_removal_folder, "singletons_no_host.blatout") + " "  # in
        hr_singletons += os.path.join(blat_hr_folder, "singletons_no_host.fastq") + " "  # out
        hr_singletons += os.path.join(blat_hr_folder, "singletons_host_only.fastq")  # out

        hr_paired = ">&2 echo BLAT contaminant paired | "
        hr_paired += self.tool_path_obj.Python + " "
        hr_paired += self.tool_path_obj.BLAT_Contaminant_Filter + " "
        hr_paired += "paired" + " "
        hr_paired += self.tool_path_obj.filter_stringency + " "
        hr_paired += os.path.join(host_removal_folder, "pair_1_no_host.fastq") + " "
        hr_paired += os.path.join(host_removal_folder, "pair_2_no_host.fastq") + " "
        hr_paired += os.path.join(host_removal_folder, "pair_1_no_host.blatout") + " "
        hr_paired += os.path.join(host_removal_folder, "pair_2_no_host.blatout") + " "
        hr_paired += os.path.join(blat_hr_folder, "pair_1_no_host.fastq") + " "
        hr_paired += os.path.join(blat_hr_folder, "pair_2_no_host.fastq") + " "
        hr_paired += os.path.join(blat_hr_folder, "pair_1_host_only.fastq") + " "
        hr_paired += os.path.join(blat_hr_folder, "pair_2_host_only.fastq")
        

        
        #-----------------------------
        #orphan correction
        
        

        copy_singletons = "cp " + os.path.join(blat_hr_folder, "singletons_no_host.fastq") + " "
        copy_singletons += os.path.join(final_folder, "singletons.fastq")

        copy_pair_1 = "cp " + os.path.join(blat_hr_folder, "pair_1_no_host.fastq") + " "
        copy_pair_1 += os.path.join(final_folder, "pair_1.fastq")

        copy_pair_2 = "cp " + os.path.join(blat_hr_folder, "pair_2_no_host.fastq") + " "
        copy_pair_2 += os.path.join(final_folder, "pair_2.fastq")
        
        
        
        if(self.tutorial_keyword is None):
            if self.read_mode == "single":
                COMMANDS_host = [
                    copy_host,
                    bwa_hr_prep,
                    samtools_hr_prep,
                    bwa_hr_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    make_blast_db_host,
                    vsearch_filter_3,
                    blat_hr_singletons,
                    hr_singletons,
                    copy_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_host = [
                    copy_host,
                    bwa_hr_prep,
                    samtools_hr_prep,
                    bwa_hr_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    bwa_hr_paired,
                    bwa_hr_filter_paired,
                    make_blast_db_host,
                    vsearch_filter_3,
                    vsearch_filter_4,
                    vsearch_filter_5,
                    blat_hr_singletons,
                    blat_hr_pair_1,
                    blat_hr_pair_2,
                    hr_singletons,
                    hr_paired,
                    copy_singletons,
                    copy_pair_1,
                    copy_pair_2
                ]
        else:
            print(dt.today(), "Host filter operating in tutorial-mode")
            if self.read_mode == "single":
                COMMANDS_host = [
                    copy_host,
                    bwa_hr_prep,
                    samtools_hr_prep,
                    bwa_hr_tut_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    make_blast_db_host,
                    vsearch_filter_3,
                    blat_hr_singletons,
                    hr_singletons,
                    copy_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_host = [
                    copy_host,
                    bwa_hr_prep,
                    samtools_hr_prep,
                    bwa_hr_tut_singletons,
                    samtools_hr_singletons_sam_to_bam,
                    samtools_no_host_singletons_bam_to_fastq,
                    samtools_host_singletons_bam_to_fastq,
                    bwa_hr_tut_paired,
                    bwa_hr_filter_paired,
                    make_blast_db_host,
                    vsearch_filter_3,
                    vsearch_filter_4,
                    vsearch_filter_5,
                    blat_hr_singletons,
                    blat_hr_pair_1,
                    blat_hr_pair_2,
                    hr_singletons,
                    hr_paired,
                    copy_singletons,
                    copy_pair_1,
                    copy_pair_2
                ]

                
        return COMMANDS_host

    def create_vector_filter_command(self, stage_name, dependency_name):
        # why do we leave all the interim files intact?
        # because science needs repeatable data, and the process needs to be able to start at any point
        subfolder                       = os.path.join(self.Output_Path, stage_name)
        data_folder                     = os.path.join(subfolder, "data")
        dependency_folder               = os.path.join(self.Output_Path, dependency_name, "final_results")
        vector_removal_folder           = os.path.join(data_folder, "0_vector_removal")
        blat_containment_vector_folder  = os.path.join(data_folder, "1_blat_containment_vr")
        final_folder                    = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(vector_removal_folder)
        self.make_folder(blat_containment_vector_folder)
        self.make_folder(final_folder)

        Vector_Contaminants = os.path.join(vector_removal_folder, "vector_contaminants_seq.fasta")

        copy_vector = ">&2 echo copy vector prep | "
        copy_vector += "cp " + self.tool_path_obj.UniVec_Core + " " + Vector_Contaminants

        bwa_vr_prep = ">&2 echo BWA vector prep | "
        bwa_vr_prep += self.tool_path_obj.BWA + " index -a bwtsw " + Vector_Contaminants

        samtools_vr_prep = ">&2 echo samtools vector prep | "
        samtools_vr_prep += self.tool_path_obj.SAMTOOLS + " faidx " + Vector_Contaminants

        bwa_vr_singletons = ">&2 echo BWA vector oprhans | "
        bwa_vr_singletons += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_singletons += Vector_Contaminants + " "
        bwa_vr_singletons += os.path.join(dependency_folder, "singletons.fastq")
        bwa_vr_singletons += " > " + os.path.join(vector_removal_folder, "singletons_no_vectors.sam")
        
        
        bwa_vr_tut_singletons = ">&2 echo BWA vector oprhans TUTORIAL MODE | "
        bwa_vr_tut_singletons += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_tut_singletons += Vector_Contaminants + " "
        bwa_vr_tut_singletons += self.sequence_single
        bwa_vr_tut_singletons += " > " + os.path.join(vector_removal_folder, "singletons_no_vectors.sam")

        samtools_no_vector_singletons_convert = ">&2 echo samtools vector oprhans pt 1 | "
        samtools_no_vector_singletons_convert += self.tool_path_obj.SAMTOOLS + " view -bS "
        samtools_no_vector_singletons_convert += os.path.join(vector_removal_folder, "singletons_no_vectors.sam")
        samtools_no_vector_singletons_convert += " > " + os.path.join(vector_removal_folder, "singletons_no_vectors.bam")

        samtools_no_vector_singletons_export = ">&2 echo samtools vector singletons pt 2 | "
        samtools_no_vector_singletons_export += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_vector_singletons_export += " -0 " + os.path.join(vector_removal_folder, "singletons_no_vectors.fastq") + " "
        samtools_no_vector_singletons_export += os.path.join(vector_removal_folder, "singletons_no_vectors.bam")

        samtools_vector_singletons_export = ">&2 echo samtools vector singletons pt 3 | "
        samtools_vector_singletons_export += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_vector_singletons_export += " -0 " + os.path.join(vector_removal_folder, "singletons_vectors_only.fastq") + " "
        samtools_vector_singletons_export += os.path.join(vector_removal_folder, "singletons_no_vectors.bam")

        bwa_vr_paired = ">&2 echo bwa vector paired | "
        bwa_vr_paired += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_paired += Vector_Contaminants + " "
        bwa_vr_paired += os.path.join(dependency_folder, "pair_1.fastq") + " "
        bwa_vr_paired += os.path.join(dependency_folder, "pair_2.fastq") + " "
        bwa_vr_paired += " > " + os.path.join(vector_removal_folder, "paired_on_vectors.sam")

        bwa_vr_tut_paired = ">&2 echo bwa vector paired TUTORIAL MODE | "
        bwa_vr_tut_paired += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_tut_paired += Vector_Contaminants + " "
        bwa_vr_tut_paired += self.sequence_path_1 + " "
        bwa_vr_tut_paired += self.sequence_path_2 + " "
        bwa_vr_tut_paired += " > " + os.path.join(vector_removal_folder, "paired_on_vectors.sam")
        
        bwa_vr_filter_paired = ">&2 echo BWA vector filter on paired | "
        bwa_vr_filter_paired += self.tool_path_obj.Python + " "
        bwa_vr_filter_paired += self.tool_path_obj.bwa_read_sorter + " "
        bwa_vr_filter_paired += "paired" + " "
        bwa_vr_filter_paired += self.tool_path_obj.filter_stringency + " "
        bwa_vr_filter_paired += os.path.join(vector_removal_folder, "paired_on_vectors.sam") + " "
        bwa_vr_filter_paired += os.path.join(dependency_folder, "pair_1.fastq") + " "
        bwa_vr_filter_paired += os.path.join(dependency_folder, "pair_2.fastq") + " "
        bwa_vr_filter_paired += os.path.join(vector_removal_folder, "pair_1_no_vectors.fastq") + " "
        bwa_vr_filter_paired += os.path.join(vector_removal_folder, "pair_2_no_vectors.fastq") + " "
        bwa_vr_filter_paired += os.path.join(vector_removal_folder, "pair_1_vectors_only.fastq") + " "
        bwa_vr_filter_paired += os.path.join(vector_removal_folder, "pair_2_vectors_only.fastq")


        make_blast_db_vector = ">&2 echo BLAST make db vectors | "
        make_blast_db_vector += self.tool_path_obj.Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl"

        vsearch_filter_6 = ">&2 echo convert vector singletons for BLAT | "
        vsearch_filter_6 += self.tool_path_obj.vsearch
        vsearch_filter_6 += " --fastq_filter " + os.path.join(vector_removal_folder, "singletons_no_vectors.fastq")
        vsearch_filter_6 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_6 += " --fastaout " + os.path.join(vector_removal_folder, "singletons_no_vectors.fasta")

        vsearch_filter_7 = ">&2 echo convert vector pair 1 for BLAT | "
        vsearch_filter_7 += self.tool_path_obj.vsearch
        vsearch_filter_7 += " --fastq_filter " + os.path.join(vector_removal_folder, "pair_1_no_vectors.fastq")
        vsearch_filter_7 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_7 += " --fastaout " + os.path.join(vector_removal_folder, "pair_1_no_vectors.fasta")

        vsearch_filter_8 = ">&2 echo convert vector pair 2 for BLAT | "
        vsearch_filter_8 += self.tool_path_obj.vsearch
        vsearch_filter_8 += " --fastq_filter " + os.path.join(vector_removal_folder, "pair_2_no_vectors.fastq")
        vsearch_filter_8 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_8 += " --fastaout " + os.path.join(vector_removal_folder, "pair_2_no_vectors.fasta")

        blat_vr_singletons = ">&2 echo BLAT vector singletons | "
        blat_vr_singletons += self.tool_path_obj.BLAT
        blat_vr_singletons += " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_singletons += Vector_Contaminants + " "
        blat_vr_singletons += os.path.join(vector_removal_folder, "singletons_no_vectors.fasta")
        blat_vr_singletons += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_singletons += os.path.join(vector_removal_folder, "singletons_no_vectors.blatout")

        blat_vr_pair_1 = ">&2 echo BLAT vector pair 1 | "
        blat_vr_pair_1 += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_1 += Vector_Contaminants + " "
        blat_vr_pair_1 += os.path.join(vector_removal_folder, "pair_1_no_vectors.fasta")
        blat_vr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_pair_1 += os.path.join(vector_removal_folder, "pair_1_no_vectors.blatout")

        blat_vr_pair_2 = ">&2 echo BLAT vector pair 2 | "
        blat_vr_pair_2 += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_2 += Vector_Contaminants + " "
        blat_vr_pair_2 += os.path.join(vector_removal_folder, "pair_2_no_vectors.fasta")
        blat_vr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_pair_2 += os.path.join(vector_removal_folder, "pair_2_no_vectors.blatout")

        blat_filter_vector_singletons = ">&2 echo BLAT contaminant singletons | "
        blat_filter_vector_singletons += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        blat_filter_vector_singletons += "single" + " "
        blat_filter_vector_singletons += self.tool_path_obj.filter_stringency + " "
        blat_filter_vector_singletons += os.path.join(vector_removal_folder, "singletons_no_vectors.fastq") + " "  # in
        blat_filter_vector_singletons += os.path.join(vector_removal_folder, "singletons_no_vectors.blatout") + " "  # in
        blat_filter_vector_singletons += os.path.join(blat_containment_vector_folder, "singletons_no_vectors.fastq") + " "  # out
        blat_filter_vector_singletons += os.path.join(blat_containment_vector_folder, "singletons_vectors_only.fastq")  # out

        blat_filter_vector_paired = ">&2 echo BLAT contaminant pair 1 | "
        blat_filter_vector_paired += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        blat_filter_vector_paired += "paired" + " "
        blat_filter_vector_paired += self.tool_path_obj.filter_stringency + " "
        blat_filter_vector_paired += os.path.join(vector_removal_folder, "pair_1_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(vector_removal_folder, "pair_2_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(vector_removal_folder, "pair_1_no_vectors.blatout") + " "
        blat_filter_vector_paired += os.path.join(vector_removal_folder, "pair_2_no_vectors.blatout") + " "
        blat_filter_vector_paired += os.path.join(blat_containment_vector_folder, "pair_1_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(blat_containment_vector_folder, "pair_2_no_vectors.fastq") + " "
        blat_filter_vector_paired += os.path.join(blat_containment_vector_folder, "pair_1_vectors_only.fastq") + " "
        blat_filter_vector_paired += os.path.join(blat_containment_vector_folder, "pair_2_vectors_only.fastq")

        copy_singletons = "cp " + os.path.join(blat_containment_vector_folder, "singletons_no_vectors.fastq") + " "
        copy_singletons += os.path.join(final_folder, "singletons.fastq")

        copy_pair_1 = "cp " + os.path.join(blat_containment_vector_folder, "pair_1_no_vectors.fastq") + " "
        copy_pair_1 += os.path.join(final_folder, "pair_1.fastq")

        copy_pair_2 = "cp " + os.path.join(blat_containment_vector_folder, "pair_2_no_vectors.fastq") + " "
        copy_pair_2 += os.path.join(final_folder, "pair_2.fastq")
        
        if(self.tutorial_keyword == "vectors"):
            if self.read_mode == "single":
                COMMANDS_vector = [
                    copy_vector,
                    bwa_vr_prep,
                    samtools_vr_prep,
                    bwa_vr_tut_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    make_blast_db_vector,
                    vsearch_filter_6,
                    blat_vr_singletons,
                    blat_filter_vector_singletons,
                    copy_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    copy_vector,
                    bwa_vr_prep,
                    samtools_vr_prep,
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
                    copy_singletons,
                    copy_pair_1,
                    copy_pair_2
                ]
        else:    
            if self.read_mode == "single":
                COMMANDS_vector = [
                    copy_vector,
                    bwa_vr_prep,
                    samtools_vr_prep,
                    bwa_vr_singletons,
                    samtools_no_vector_singletons_convert,
                    samtools_no_vector_singletons_export,
                    samtools_vector_singletons_export,
                    make_blast_db_vector,
                    vsearch_filter_6,
                    blat_vr_singletons,
                    blat_filter_vector_singletons,
                    copy_singletons
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    copy_vector,
                    bwa_vr_prep,
                    samtools_vr_prep,
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
                    copy_singletons,
                    copy_pair_1,
                    copy_pair_2
                ]    

        return COMMANDS_vector

   

        
        
    def create_rRNA_filter_prep_command_v3(self, stage_name, category, dependency_name, marker_file):
        #split the data into tiny shards.  called once
        dep_loc                 = os.path.join(self.Output_Path, dependency_name, "final_results")
        subfolder               = os.path.join(self.Output_Path, stage_name)
        data_folder             = os.path.join(subfolder, "data")
        
        jobs_folder             = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(jobs_folder)
        
        
        if(self.tutorial_keyword == "rRNA"):
            category = "singletons"
            split_folder            = os.path.join(data_folder, category + "_fastq")
            self.make_folder(split_folder)
            split_tut_single_fastq = ">&2 echo splitting fastq for " + category + " | " 
            split_tut_single_fastq += self.tool_path_obj.Python + " "
            split_tut_single_fastq += self.tool_path_obj.File_splitter + " "
            split_tut_single_fastq += self.sequence_single + " " #os.path.join(dep_loc, category + ".fastq") + " "
            split_tut_single_fastq += os.path.join(split_folder, category) + " "
            split_tut_single_fastq += str(self.tool_path_obj.rRNA_chunksize)

            if(self.read_mode == "paired"):
                category = "pair_1"
                split_folder            = os.path.join(data_folder, category + "_fastq")
                self.make_folder(split_folder)
                split_tut_pair_1_fastq = ">&2 echo splitting fastq for " + category + " | " 
                split_tut_pair_1_fastq += self.tool_path_obj.Python + " "
                split_tut_pair_1_fastq += self.tool_path_obj.File_splitter + " "
                split_tut_pair_1_fastq += self.sequence_path_1 + " "#os.path.join(dep_loc, category + ".fastq") + " "
                split_tut_pair_1_fastq += os.path.join(split_folder, category) + " "
                split_tut_pair_1_fastq += str(self.tool_path_obj.rRNA_chunksize)
                
                category = "pair_2"
                split_folder            = os.path.join(data_folder, category + "_fastq")
                self.make_folder(split_folder)
                split_tut_pair_2_fastq = ">&2 echo splitting fastq for " + category + " | " 
                split_tut_pair_2_fastq += self.tool_path_obj.Python + " "
                split_tut_pair_2_fastq += self.tool_path_obj.File_splitter + " "
                split_tut_pair_2_fastq += self.sequence_path_2 + " "#os.path.join(dep_loc, category + ".fastq") + " "
                split_tut_pair_2_fastq += os.path.join(split_folder, category) + " "
                split_tut_pair_2_fastq += str(self.tool_path_obj.rRNA_chunksize)
            
                return [split_tut_single_fastq + " && " + split_tut_pair_1_fastq + " && " + split_tut_pair_2_fastq]
            else:
                return [split_tut_single_fastq]
        
        else:
            split_folder            = os.path.join(data_folder, category + "_fastq")
            self.make_folder(split_folder)
            
            
            split_fastq = ">&2 echo splitting fastq for " + category + " | " 
            split_fastq += self.tool_path_obj.Python + " "
            split_fastq += self.tool_path_obj.File_splitter + " "
            split_fastq += os.path.join(dep_loc, category + ".fastq") + " "
            split_fastq += os.path.join(split_folder, category) + " "
            split_fastq += str(self.tool_path_obj.rRNA_chunksize)
            
        
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [split_fastq + " && " + make_marker]

    def create_rRNA_filter_convert_fastq_command(self, stage_name, category, fastq_name, marker_file):
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        file_name           = fastq_name.split(".")[0]
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        fastq_seqs          = os.path.join(fastq_folder, fastq_name)
        
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")

        self.make_folder(jobs_folder)
        
        #if(self.tutorial_keyword == "rRNA"):
        #    self.make_folder(tut_fasta_folder)
        #else:
        self.make_folder(fasta_folder)
            
        
        convert_fastq_to_fasta = ">&2 echo " + " converting " + file_name + " file to fasta | "
        convert_fastq_to_fasta += self.tool_path_obj.vsearch
        convert_fastq_to_fasta += " --fastq_filter " + fastq_seqs
        convert_fastq_to_fasta += " --fastq_ascii " + self.Qual_str
        convert_fastq_to_fasta += " --fastaout " + fasta_seqs
        
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
 
        return [convert_fastq_to_fasta + " && " + make_marker]
    
    def create_rRNA_filter_barrnap_arc_command(self, stage_name, category, fastq_name, marker_file):
        # called by each split file
        # category -> singletons, pair 1, pair 2
        # file name -> the specific split section of the category (the fastq segments)
        # stage_name -> "rRNA_Filter"
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        file_name           = fastq_name.split(".")[0]
        Barrnap_arc_out     = os.path.join(Barrnap_out_folder, file_name + "_arc.barrnap_out")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        
        
        #if(self.tutorial_keyword == "rRNA"):
        #    self.make_folder(tut_fasta_folder)
        #    self.make_folder(tut_Barrnap_out_folder)
        #else:
        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)
    
        
        Barrnap_archaea = ">&2 echo running Barrnap on " + file_name + " file: arc | "
        Barrnap_archaea += self.tool_path_obj.Barrnap
        Barrnap_archaea += " --quiet --reject 0.01 --kingdom " + "arc"
        Barrnap_archaea += " --threads " + self.Threads_str
        Barrnap_archaea += " " + fasta_seqs
        Barrnap_archaea += " >> " + Barrnap_arc_out
        

    
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
  

        return [Barrnap_archaea + " && " + make_marker]
         
        
    def create_rRNA_filter_barrnap_bac_command(self, stage_name, category, fastq_name, marker_file):

        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        file_name           = fastq_name.split(".")[0]
        Barrnap_bac_out     = os.path.join(Barrnap_out_folder, file_name + "_bac.barrnap_out")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        
        #if(self.tutorial_keyword == "rRNA"):
        #    self.make_folder(tut_fasta_folder)
        #    self.make_folder(tut_Barrnap_out_folder)
        #else:
        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)
    
        

        Barrnap_bacteria = ">&2 echo Running Barrnap on " + file_name + " file:  bac | "
        Barrnap_bacteria += self.tool_path_obj.Barrnap
        Barrnap_bacteria += " --quiet --reject 0.01 --kingdom " + "bac"
        Barrnap_bacteria += " --threads " + self.Threads_str
        Barrnap_bacteria += " " + fasta_seqs
        Barrnap_bacteria += " >> " + Barrnap_bac_out
  
        
        make_marker = "touch"  + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        # if(self.tutorial_keyword == "rRNA"):
            # if(self.read_mode == "single"):
                # return [Barrnap_tut_singletons_bacteria]
            # elif(self.read_mode == "paired"):
                # return [Barrnap_tut_singletons_bacteria + " && " + Barrnap_tut_pair_1_bacteria + " && " + Barrnap_tut_pair_2_bacteria]
        # else:
        return [Barrnap_bacteria + " && " + make_marker]
        
    def create_rRNA_filter_barrnap_euk_command(self, stage_name, category, fastq_name, marker_file):

        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        file_name           = fastq_name.split(".")[0]
        Barrnap_euk_out     = os.path.join(Barrnap_out_folder, file_name + "_euk.barrnap_out")
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        

        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)

        Barrnap_eukaryote = ">&2 echo Running Barrnap on " + file_name + " file: euk | "
        Barrnap_eukaryote += self.tool_path_obj.Barrnap
        Barrnap_eukaryote += " --quiet --reject 0.01 --kingdom " + "euk"
        Barrnap_eukaryote += " --threads " + self.Threads_str
        Barrnap_eukaryote += " " + fasta_seqs
        Barrnap_eukaryote += " >> " + Barrnap_euk_out
        

        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
    
    
        return [Barrnap_eukaryote + " && " + make_marker]
        
    def create_rRNA_filter_barrnap_mit_command(self, stage_name, category, fastq_name, marker_file):
        #designed to run on a single split sample.
        #expected to be merged later with all the other runs of the same fastq name
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        file_name           = fastq_name.split(".")[0]
        Barrnap_mit_out     = os.path.join(Barrnap_out_folder, file_name + "_mit.barrnap_out")
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        

        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)

        Barrnap_mitochondria = ">&2 echo Running Barrnap on " + file_name + " file: mito | " 
        Barrnap_mitochondria += self.tool_path_obj.Barrnap
        Barrnap_mitochondria += " --quiet --reject 0.01 --kingdom " + "mito"
        Barrnap_mitochondria += " --threads " + self.Threads_str
        Barrnap_mitochondria += " " + fasta_seqs
        Barrnap_mitochondria += " >> " + Barrnap_mit_out
        


        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [Barrnap_mitochondria + " && " + make_marker]
        
        
        
        
        
        
    def create_rRNA_filter_barrnap_cat_command(self, stage_name, category, fastq_name, marker_file):
        #this is expected to run on each sample split
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        infernal_out_folder = os.path.join(data_folder, category + "_infernal")
        file_name           = fastq_name.split(".")[0]
        Barrnap_arc_out     = os.path.join(Barrnap_out_folder, file_name + "_arc.barrnap_out")
        Barrnap_bac_out     = os.path.join(Barrnap_out_folder, file_name + "_bac.barrnap_out")
        Barrnap_euk_out     = os.path.join(Barrnap_out_folder, file_name + "_euk.barrnap_out")
        Barrnap_mit_out     = os.path.join(Barrnap_out_folder, file_name + "_mit.barrnap_out")
        infernal_out        = os.path.join(infernal_out_folder, file_name + ".infernal_out")
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        Barrnap_out         = os.path.join(Barrnap_out_folder, file_name + ".barrnap_out")

        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        

        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(infernal_out_folder)
        self.make_folder(jobs_folder)
   
        

        
        #combine the arc, bac, euk, mit files into 1
        cat_command = ">&2 echo Combining files for:" + file_name + " | "
        cat_command += "cat" + " "
        cat_command += Barrnap_arc_out + " " + Barrnap_bac_out + " " + Barrnap_euk_out + " " + Barrnap_mit_out + " "
        cat_command += ">>" + " " + Barrnap_out
        
        rm_arc = ">&2 echo delete arc: " + file_name + " | "
        rm_arc += "rm" + " "
        rm_arc += Barrnap_arc_out
        
        rm_bac = ">&2 echo delete bac: " + file_name + " | "
        rm_bac += "rm" + " "
        rm_bac += Barrnap_bac_out
        
        rm_euk = ">&2 echo delete euk: " + file_name + " | "
        rm_euk += "rm" + " "
        rm_euk += Barrnap_euk_out
        
        rm_mit = ">&2 echo delete mit: " + file_name + " | "
        rm_mit += "rm" + " "
        rm_mit += Barrnap_mit_out
        
        make_marker = "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [cat_command + " && " + make_marker + " && " + rm_arc  + " && " + rm_bac  + " && " +  rm_euk  + " && " +  rm_mit]
        
        
    def create_rRNA_filter_barrnap_pp_command(self, stage_name, category, fastq_name, marker_file):
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap")
        infernal_out_folder = os.path.join(data_folder, category + "_infernal")
        mRNA_folder         = os.path.join(data_folder, category + "_barrnap_mRNA")
        rRNA_folder         = os.path.join(data_folder, category + "_barrnap_other")
        file_name           = fastq_name.split(".")[0]
        Barrnap_out         = os.path.join(Barrnap_out_folder, file_name + ".barrnap_out")
        fastq_seqs          = os.path.join(fastq_folder, fastq_name)
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        tut_mRNA_folder     = os.path.join(data_folder, "tutorial_barrnap_mRNA")
        tut_rRNA_folder     = os.path.join(data_folder, "tutorial_barrnap_rRNA")
        tut_fasta_folder    = os.path.join(data_folder, "tutorial_fasta")
        tut_Barrnap_out_folder = os.path.join(data_folder, "tutorial_barrnap")
        # if(self.tutorial_keyword == "rRNA"):
            # self.make_folder(tut_mRNA_folder)
            # self.make_folder(tut_rRNA_folder)
        # else:
        self.make_folder(fasta_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(infernal_out_folder)
        self.make_folder(mRNA_folder)
        self.make_folder(rRNA_folder)
        self.make_folder(jobs_folder)
        
        
        Barrnap_pp = ">&2 echo Running Barrnap pp scripts | "
        Barrnap_pp += self.tool_path_obj.Python + " "
        Barrnap_pp += self.tool_path_obj.barrnap_post + " "
        Barrnap_pp += Barrnap_out + " "
        Barrnap_pp += fastq_seqs + " "
        Barrnap_pp += mRNA_folder + " "
        Barrnap_pp += rRNA_folder + " "
        Barrnap_pp += file_name + "_barrnap"
        

        
        
        #make_marker = ">&2 echo " + file_name + "_barrnap Marking job completed | " 
        make_marker = "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [Barrnap_pp + " && " + make_marker]
        
        
    def create_rRNA_filter_infernal_prep_command(self, stage_name, category, fastq_name, root_name, marker_file):
        #expecting full file name in fastq_name
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap_mRNA_fasta")
        infernal_out_folder = os.path.join(data_folder, category + "_infernal")
        mRNA_folder         = os.path.join(data_folder, category + "_barrnap_mRNA")
        file_name           = fastq_name.split(".")[0]
        Barrnap_out         = os.path.join(Barrnap_out_folder, file_name + ".barrnap_out")
        infernal_out        = os.path.join(infernal_out_folder, file_name + ".infernal_out")
        jobs_folder         = os.path.join(data_folder, "jobs")
        fastq_seqs          = os.path.join(fastq_folder, fastq_name)
        
        fasta_seqs          = os.path.join(fasta_folder, file_name + ".fasta")
        
        tut_Barrnap_mRNA_folder     = os.path.join(data_folder, "tutorial_barrnap_mRNA")
        tut_infernal_input_folder   = os.path.join(data_folder, "tutorial_infernal_input")
        # if(self.tutorial_keyword == "rRNA"):
            # self.make_folder(tut_Barrnap_mRNA_folder)
            # self.make_folder(tut_infernal_input_folder)
        # else:
        self.make_folder(infernal_out_folder)
        self.make_folder(mRNA_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)
        
        convert_fastq_to_fasta_barrnap = ">&2 echo converting barrnap fastq to fasta:" + file_name + " | "
        convert_fastq_to_fasta_barrnap += self.tool_path_obj.vsearch
        convert_fastq_to_fasta_barrnap += " --fastq_filter " + os.path.join(mRNA_folder, fastq_name)
        convert_fastq_to_fasta_barrnap += " --fastq_ascii " + self.Qual_str
        convert_fastq_to_fasta_barrnap += " --fastaout " + os.path.join(Barrnap_out_folder, file_name + ".fasta")
        

        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [convert_fastq_to_fasta_barrnap + " && " + make_marker]

    def create_rRNA_filter_infernal_command(self, stage_name, category, file_name, marker_file):
    
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        Barrnap_out_folder  = os.path.join(data_folder, category + "_barrnap_mRNA_fasta")
        infernal_out_folder = os.path.join(data_folder, category + "_infernal")
        mRNA_folder         = os.path.join(data_folder, category + "_barrnap_mRNA")
        Barrnap_out         = os.path.join(Barrnap_out_folder, file_name + ".barrnap_out")
        infernal_out        = os.path.join(infernal_out_folder, file_name + ".infernal_out")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        tut_Barrnap_mRNA_folder     = os.path.join(data_folder, "tutorial_barrnap_mRNA")
        tut_infernal_input_folder   = os.path.join(data_folder, "tutorial_infernal_input")
        # if(self.tutorial_keyword == "rRNA"):
            # self.make_folder(tut_Barrnap_mRNA_folder)
            # self.make_folder(tut_infernal_input_folder)
            # self.make_folder(infernal_out_folder)
        # else:

        self.make_folder(infernal_out_folder)
        self.make_folder(mRNA_folder)
        self.make_folder(Barrnap_out_folder)
        self.make_folder(jobs_folder)
        

        infernal_command = ">&2 echo " + str(dt.today()) + " running infernal on " + file_name + " file | "
        infernal_command += self.tool_path_obj.Infernal
        infernal_command += " -o /dev/null --tblout "
        infernal_command += infernal_out
        #infernal_command += " --cpu " + self.Threads_str -> lined nerf'd because infernal's parallelism is not good
        infernal_command += " --cpu 1"
        infernal_command += " --anytrunc --rfam -E 0.001 "
        infernal_command += self.tool_path_obj.Rfam + " "
        infernal_command += os.path.join(Barrnap_out_folder, file_name + "_barrnap_mRNA.fasta")
  
        
        #make_marker = ">&2 echo " + file_name + "_infernal Marking job completed | " 
        make_marker = "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)
        

        return [infernal_command + " && " + make_marker]
        
    
    def create_rRNA_filter_splitter_command(self, stage_name, category, file_name, marker_file):
    #file name expected to have no extensions.  eg: pair_1_0
    #expected to be called for each category (pair1, singletons).  not pair 2.  paired data is handled in combination
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        fasta_folder        = os.path.join(data_folder, category + "_fasta")
        fastq_folder        = os.path.join(data_folder, category + "_fastq")
        infernal_out_folder = os.path.join(data_folder, category + "_infernal")
        mRNA_barrnap_folder = os.path.join(data_folder, category + "_mRNA")
        mRNA_infernal_folder= os.path.join(data_folder, category + "_infernal_mRNA")
        rRNA_folder         = os.path.join(data_folder, category + "_infernal_rRNA")
        infernal_out        = os.path.join(infernal_out_folder, file_name + ".infernal_out")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        file_name_code = file_name.split("_")[-1]
        self.make_folder(mRNA_infernal_folder)
        self.make_folder(jobs_folder)
        
        
        if(category == "pair_1"):            
            infernal_pair_1_out_folder = os.path.join(data_folder, "pair_1_infernal")
            infernal_pair_2_out_folder = os.path.join(data_folder, "pair_2_infernal")
            
            infernal_mRNA_pair_1_folder = os.path.join(data_folder, "pair_1_infernal_mRNA")
            infernal_mRNA_pair_2_folder = os.path.join(data_folder, "pair_2_infernal_mRNA")
            
            infernal_rRNA_pair_1_folder = os.path.join(data_folder, "pair_1_infernal_other")
            infernal_rRNA_pair_2_folder = os.path.join(data_folder, "pair_2_infernal_other")
            
            Barrnap_pair_1_out_folder = os.path.join(data_folder, "pair_1_barrnap")
            Barrnap_pair_2_out_folder = os.path.join(data_folder, "pair_2_barrnap")
            
            self.make_folder(infernal_mRNA_pair_1_folder)
            self.make_folder(infernal_mRNA_pair_2_folder)
            self.make_folder(infernal_rRNA_pair_1_folder)
            self.make_folder(infernal_rRNA_pair_2_folder)
            self.make_folder(Barrnap_pair_1_out_folder)
            self.make_folder(Barrnap_pair_2_out_folder)
            
            rRNA_filtration = ">&2 echo extracting mRNA with infernal report: " + file_name + " | "
            rRNA_filtration += self.tool_path_obj.Python + " "
            rRNA_filtration += self.tool_path_obj.rRNA_filter + " "
            rRNA_filtration += self.tool_path_obj.filter_stringency + " "
            rRNA_filtration += "paired" + " "
            rRNA_filtration += os.path.join(infernal_pair_1_out_folder, "pair_1_" + file_name_code + ".infernal_out") + " "
            rRNA_filtration += os.path.join(infernal_pair_2_out_folder, "pair_2_" + file_name_code + ".infernal_out") + " "
            rRNA_filtration += os.path.join(Barrnap_pair_1_out_folder, "pair_1_" + file_name_code + ".barrnap_out") + " "
            rRNA_filtration += os.path.join(Barrnap_pair_2_out_folder, "pair_2_" + file_name_code + ".barrnap_out") + " "
            rRNA_filtration += os.path.join(data_folder, "pair_1_fastq", "pair_1_" + file_name_code + ".fastq") + " "
            rRNA_filtration += os.path.join(data_folder, "pair_2_fastq", "pair_2_" + file_name_code + ".fastq") + " "
            rRNA_filtration += os.path.join(infernal_mRNA_pair_1_folder, "pair_1_" + file_name_code + "_infernal_mRNA.fastq") + " "
            rRNA_filtration += os.path.join(infernal_mRNA_pair_2_folder, "pair_2_" + file_name_code + "_infernal_mRNA.fastq") + " "
            rRNA_filtration += os.path.join(infernal_rRNA_pair_1_folder, "pair_1_" + file_name_code + "_infernal_other.fastq") + " "
            rRNA_filtration += os.path.join(infernal_rRNA_pair_2_folder, "pair_2_" + file_name_code + "_infernal_other.fastq")
            
            
            
        elif(category == "singletons"):    
            infernal_mRNA_singletons_folder = os.path.join(data_folder, "singletons_infernal_mRNA")
            infernal_rRNA_singletons_folder = os.path.join(data_folder, "singletons_infernal_other")
            Barrnap_singletons_out_folder   = os.path.join(data_folder, "singletons_barrnap")
            infernal_singletons_out_folder  = os.path.join(data_folder, "singletons_infernal")
            
            self.make_folder(infernal_mRNA_singletons_folder)
            self.make_folder(infernal_rRNA_singletons_folder)
            
            rRNA_filtration = ">&2 echo extracting mRNA with infernal report: " + file_name + " | "
            rRNA_filtration += self.tool_path_obj.Python + " "
            rRNA_filtration += self.tool_path_obj.rRNA_filter + " "
            rRNA_filtration += self.tool_path_obj.filter_stringency + " "
            rRNA_filtration += "single" + " "
            rRNA_filtration += os.path.join(infernal_out_folder, "singletons_" + file_name_code + ".infernal_out") + " "
            rRNA_filtration += os.path.join(Barrnap_singletons_out_folder, "singletons_" + file_name_code + ".barrnap_out") + " "
            rRNA_filtration += os.path.join(data_folder, "singletons_fastq", "singletons_" + file_name_code + ".fastq") + " "
            rRNA_filtration += os.path.join(infernal_mRNA_singletons_folder, "singletons_" + file_name_code + "_infernal_mRNA.fastq") + " "
            rRNA_filtration += os.path.join(infernal_rRNA_singletons_folder, "singletons_" + file_name_code + "_infernal_rRNA.fastq")
            
        else:
            rRNA_filtration = ">&2 echo rRNA filtration ignored for pair 2 data: " + file_name
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
            
        return [rRNA_filtration + " && " + make_marker]

    
    
    def create_rRNA_filter_final_cat_command(self, stage_name, category, marker_file):
        # 
        # Cat then final filter.  
        # operates in sections
        subfolder               = os.path.join(self.Output_Path, stage_name)
        data_folder             = os.path.join(subfolder, "data")
        infernal_mRNA_folder    = os.path.join(data_folder, category + "_infernal_mRNA")
        infernal_rRNA_folder    = os.path.join(data_folder, category + "_infernal_other")
        final_folder            = os.path.join(subfolder, "final_results")
        final_rRNA_folder       = os.path.join(final_folder, "other")
        final_mRNA_folder       = os.path.join(final_folder, "mRNA")
        jobs_folder             = os.path.join(data_folder, "jobs")
        
        #self.make_folder(merged_infernal_mRNA_folder)
        self.make_folder(final_folder)
        self.make_folder(final_rRNA_folder)
        self.make_folder(final_mRNA_folder)
        
        cat_infernal_mRNA = ">&2 echo merging mRNA | "
        cat_infernal_mRNA += "for f in" + " "
        cat_infernal_mRNA += infernal_mRNA_folder
        cat_infernal_mRNA += "/*; do cat \"$f\" >>" + " "
        cat_infernal_mRNA += os.path.join(final_mRNA_folder, category + ".fastq")
        cat_infernal_mRNA += "; done"
        
        cat_infernal_rRNA = ">&2 echo merging infernal rRNA | " 
        cat_infernal_rRNA += "for f in" + " "
        cat_infernal_rRNA += infernal_rRNA_folder
        cat_infernal_rRNA += "/*; do cat \"$f\" >>" + " " 
        cat_infernal_rRNA += os.path.join(final_rRNA_folder, category + "_other.fastq")
        cat_infernal_rRNA += "; done"

        make_marker = "touch" +  " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        
        COMMANDS_rRNA_post = [cat_infernal_mRNA + " && " + cat_infernal_rRNA + " && " + make_marker]

        return COMMANDS_rRNA_post

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
        repop_singletons += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        #the reference data to be drawn from 
        if self.read_mode == "single":
            repop_singletons += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        
        repop_singletons += os.path.join(dep_loc, "mRNA", "singletons.fastq") + " "  # in -> rRNA filtration output
        repop_singletons += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output

        
        if self.read_mode == "single":
            repop_singletons += os.path.join(final_folder, "singletons.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(repop_folder, "singletons.fastq")  # out
            
            

        repop_singletons_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulations singletons rRNA | "
        repop_singletons_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
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
        repop_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1 += os.path.join(hq_path, "pair_1_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_1 += self.sequence_path_1 + " "
        else:
            repop_pair_1 += os.path.join(dep_loc, "mRNA", "pair_1.fastq") + " "
        repop_pair_1 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1 += os.path.join(repop_folder, "pair_1.fastq")

        repop_pair_1_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 rRNA | "
        repop_pair_1_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1_rRNA += os.path.join(hq_path, "pair_1_match.fastq") + " "
        repop_pair_1_rRNA += os.path.join(dep_loc, "other", "pair_1_other.fastq") + " "
        repop_pair_1_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq")

        repop_pair_2 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2 += os.path.join(hq_path, "pair_2_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_2 += self.sequence_path_2 + " "
        else:
            repop_pair_2 += os.path.join(dep_loc, "mRNA", "pair_2.fastq") + " "
        repop_pair_2 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2 += os.path.join(repop_folder, "pair_2.fastq")

        repop_pair_2_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2_rRNA += os.path.join(hq_path, "pair_2_match.fastq") + " "
        repop_pair_2_rRNA += os.path.join(dep_loc, "other", "pair_2_other.fastq") + " "
        repop_pair_2_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq")

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.tool_path_obj.Python + " "
        singleton_repop_filter += self.tool_path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.tool_path_obj.Python + " "
        singleton_repop_filter_rRNA += self.tool_path_obj.orphaned_read_filter + " "
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
        repop_singletons += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        #the reference data to be drawn from 
        if self.read_mode == "single":
            repop_singletons += os.path.join(singleton_path, "singletons_hq.fastq") + " "
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(hq_path, "singletons_with_duplicates.fastq") + " "
        
        repop_singletons += os.path.join(dep_loc, "mRNA", "singletons.fastq") + " "  # in -> rRNA filtration output
        repop_singletons += os.path.join(cluster_path, "singletons_unique.fastq.clstr") + " "  # in -> duplicates filter output

        
        if self.read_mode == "single":
            repop_singletons += os.path.join(final_folder, "singletons.fastq")  # out
        elif self.read_mode == "paired":
            repop_singletons += os.path.join(repop_folder, "singletons.fastq")  # out
            
            

        repop_singletons_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulations singletons rRNA | "
        repop_singletons_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
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
        repop_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1 += os.path.join(hq_path, "pair_1_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_1 += self.sequence_path_1 + " "
        else:
            repop_pair_1 += os.path.join(dep_loc, "mRNA", "pair_1.fastq") + " "
        repop_pair_1 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1 += os.path.join(repop_folder, "pair_1.fastq")

        repop_pair_1_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 rRNA | "
        repop_pair_1_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1_rRNA += os.path.join(hq_path, "pair_1_match.fastq") + " "
        repop_pair_1_rRNA += os.path.join(dep_loc, "other", "pair_1_other.fastq") + " "
        repop_pair_1_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_rRNA += os.path.join(repop_folder, "pair_1_rRNA.fastq")

        repop_pair_2 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2 += os.path.join(hq_path, "pair_2_match.fastq") + " "
        if(self.tutorial_keyword == tut_keyword):
            repop_pair_2 += self.sequence_path_2 + " "
        else:
            repop_pair_2 += os.path.join(dep_loc, "mRNA", "pair_2.fastq") + " "
        repop_pair_2 += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2 += os.path.join(repop_folder, "pair_2.fastq")

        repop_pair_2_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2_rRNA += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2_rRNA += os.path.join(hq_path, "pair_2_match.fastq") + " "
        repop_pair_2_rRNA += os.path.join(dep_loc, "other", "pair_2_other.fastq") + " "
        repop_pair_2_rRNA += os.path.join(cluster_path, "pair_1_unique.fastq.clstr") + " "
        repop_pair_2_rRNA += os.path.join(repop_folder, "pair_2_rRNA.fastq")

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.tool_path_obj.Python + " "
        singleton_repop_filter += self.tool_path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.tool_path_obj.Python + " "
        singleton_repop_filter_rRNA += self.tool_path_obj.orphaned_read_filter + " "
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
        singleton_repop_filter += self.tool_path_obj.Python + " "
        singleton_repop_filter += self.tool_path_obj.orphaned_read_filter + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(repop_folder, "singletons.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_1.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "pair_2.fastq") + " "
        singleton_repop_filter += os.path.join(final_folder, "singletons.fastq")
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.tool_path_obj.Python + " "
        singleton_repop_filter_rRNA += self.tool_path_obj.orphaned_read_filter + " "
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
        spades += self.tool_path_obj.Python + " "
        spades += self.tool_path_obj.Spades + " --rna"
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
        disassemble_contigs += self.tool_path_obj.MetaGeneMark + " -o " + mgm_report + " "
        disassemble_contigs += "-D " + post_mgm_contig + " "
        disassemble_contigs += "-m " + self.tool_path_obj.mgm_model + " "
        disassemble_contigs += os.path.join(spades_folder, "contigs.fasta")
        
        remove_whitespace = ">&2 echo Removing whitespace from fasta | " 
        remove_whitespace += self.tool_path_obj.Python + " " + self.tool_path_obj.remove_gaps_in_fasta + " "
        remove_whitespace += post_mgm_contig + " "
        remove_whitespace += final_contigs
        
        #BWA-ing against the final contigs gives us a proper contig-segment -> read map. 
        bwa_index = self.tool_path_obj.BWA + " index -a bwtsw " + final_contigs
        
        
        # Build a report of what was consumed by contig transmutation (assemble/disassemble)
        bwa_paired_contigs = ">&2 echo BWA pair contigs | "
        bwa_paired_contigs += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_paired_contigs += final_contigs + " "
        bwa_paired_contigs += os.path.join(dep_loc, "pair_1.fastq") + " "
        bwa_paired_contigs += os.path.join(dep_loc, "pair_2.fastq") + " "
        bwa_paired_contigs += ">" + " " 
        bwa_paired_contigs += os.path.join(bwa_folder, "paired_on_contigs.sam")

        bwa_singletons_contigs = ">&2 echo BWA singleton contigs | "
        bwa_singletons_contigs += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_singletons_contigs += final_contigs + " "
        bwa_singletons_contigs += os.path.join(dep_loc, "singletons.fastq")
        bwa_singletons_contigs += " > " + os.path.join(bwa_folder, "singletons_on_contigs.sam")
        
        make_contig_map = ">&2 echo Making contig map | " 
        make_contig_map += self.tool_path_obj.Python + " "
        make_contig_map += self.tool_path_obj.Map_contig + " "
        make_contig_map += self.read_mode + " "
        make_contig_map += dep_loc + " "
        make_contig_map += final_folder + " "
        make_contig_map += os.path.join(bwa_folder, "singletons_on_contigs.sam") + " "
        if(self.read_mode == "paired"):
            make_contig_map += os.path.join(bwa_folder, "paired_on_contigs.sam")
            
            
        flush_bad_contigs = ">&2 echo flush bad contigs | " 
        flush_bad_contigs += self.tool_path_obj.Python + " "
        flush_bad_contigs += self.tool_path_obj.flush_bad_contigs + " "
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

   
    def create_split_ga_fastq_data_command(self, stage_name, dependency_stage_name, category, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        split_folder    = os.path.join(data_folder, "0_read_split", category)
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(split_folder)
        self.make_folder(jobs_folder)
        
        if(self.tutorial_keyword == "GA"):
            if(category == "pair_1"):
            
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.tool_path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_path_1 + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(split_folder, category + "_")
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
                
            elif(category == "pair_2"):
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.tool_path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_path_2 + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(split_folder, category + "_")
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
            elif(category == "singletons"):
                split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
                split_fastq += "split -l " + str(int(self.tool_path_obj.GA_chunksize) * 4) + " "        
                split_fastq += self.sequence_single + " "
                split_fastq += "--additional-suffix .fastq" + " "
                split_fastq += "-d" + " "
                split_fastq += os.path.join(split_folder, category + "_")
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fastq = [
                    split_fastq + " && " + make_marker
                ]
     
        else:
            split_fastq = ">&2 echo splitting fastq for " + category + " GA | "
            split_fastq += "split -l " + str(int(self.tool_path_obj.GA_chunksize) * 4) + " "        
            split_fastq += os.path.join(dep_loc, category + ".fastq") + " "
            split_fastq += "--additional-suffix .fastq" + " "
            split_fastq += "-d" + " "
            split_fastq += os.path.join(split_folder, category + "_")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            COMMANDS_GA_prep_fastq = [
                split_fastq + " && " + make_marker
            ]
            
        return COMMANDS_GA_prep_fastq

    def create_split_ga_fasta_data_command(self, stage_name, dependency_stage_name, category, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        split_folder    = os.path.join(data_folder, "0_read_split", category)
        dep_folder      = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(split_folder)
        self.make_folder(jobs_folder)
        
        
        if(self.tutorial_keyword == "GA"):
            if(category == "singletons"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.tool_path_obj.Python + " "    
                split_fasta += self.tool_path_obj.File_splitter + " "
                split_fasta += self.sequence_single + " "
                split_fasta += os.path.join(split_folder, category) + " "
                split_fasta += str(self.tool_path_obj.GA_chunksize)
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
                
            elif(category == "contigs"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.tool_path_obj.Python + " "    
                split_fasta += self.tool_path_obj.File_splitter + " "
                split_fasta += self.sequence_contigs + " "
                split_fasta += os.path.join(split_folder, category) + " "
                split_fasta += str(self.tool_path_obj.GA_chunksize)
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
            elif(category == "pair_1"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.tool_path_obj.Python + " "    
                split_fasta += self.tool_path_obj.File_splitter + " "
                split_fasta += self.sequence_path_1 + " "
                split_fasta += os.path.join(split_folder, category) + " "
                split_fasta += str(self.tool_path_obj.GA_chunksize)
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
            elif(category == "pair_2"):
                split_fasta = ">&2 echo splitting fasta for " + category + " | "
                split_fasta += self.tool_path_obj.Python + " "    
                split_fasta += self.tool_path_obj.File_splitter + " "
                split_fasta += self.sequence_path_2 + " "
                split_fasta += os.path.join(split_folder, category) + " "
                split_fasta += str(self.tool_path_obj.GA_chunksize)
                
                make_marker = "touch" + " "
                make_marker += os.path.join(jobs_folder, marker_file)
                
                COMMANDS_GA_prep_fasta = [
                    split_fasta + " && " + make_marker
                ]
                
        else:
            split_fasta = ">&2 echo splitting fasta for " + category + " | "
            split_fasta += self.tool_path_obj.Python + " "    
            split_fasta += self.tool_path_obj.File_splitter + " "
            split_fasta += os.path.join(dep_folder, category +".fasta") + " "
            split_fasta += os.path.join(split_folder, category) + " "
            split_fasta += str(self.tool_path_obj.GA_chunksize)
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            COMMANDS_GA_prep_fasta = [
                split_fasta + " && " + make_marker
            ]
        
        return COMMANDS_GA_prep_fasta


    def create_BWA_annotate_command_v2(self, stage_name, ref_path, ref_tag, query_file, marker_file):
        # meant to be called multiple times: query file is a split file
        # aug 10, 2021: changed ref path to accomodate new split-chocophlan
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        bwa_folder      = os.path.join(data_folder, "1_bwa")
        jobs_folder     = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(jobs_folder)

        file_tag = os.path.basename(query_file)
        file_tag = os.path.splitext(file_tag)[0]
        
        bwa_job = ">&2 echo " + str(dt.today()) + " BWA on " + file_tag + " | "
        bwa_job += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_job += ref_path + " " #self.tool_path_obj.DNA_DB + " "
        #bwa_job += os.path.join(dep_loc, section_file) + " | "
        bwa_job += query_file + " | "
        bwa_job += self.tool_path_obj.SAMTOOLS + " view "
        bwa_job += "> " + os.path.join(bwa_folder, file_tag +"_" + ref_tag + ".sam")
        
        #make_marker = ">&2 echo marking BWA job complete: " + file_tag + " | "
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_BWA = [
            bwa_job + " && " + make_marker
        ]

        return COMMANDS_BWA
        
        
        
        
    def create_BWA_pp_command_v2(self, stage_name, dependency_stage_name, ref_tag, ref_path, query_file, marker_file):
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
            
        
        #meant to be called on the split-file version.  PP script will not merge gene maps.
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        bwa_folder      = os.path.join(data_folder, "1_bwa")
        split_folder    = os.path.join(data_folder, "0_read_split")
        final_folder    = os.path.join(subfolder, "final_results")
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        
        reads_in    = query_file
        bwa_in      = os.path.join(bwa_folder, sample_root_name + "_" + ref_tag + ".sam")
        reads_out = os.path.join(final_folder, sample_root_name + "_" + ref_tag + ".fasta")
            
        map_read_bwa = ">&2 echo " + str(dt.today()) + " GA BWA PP generic: " + sample_root_name + " | "
        map_read_bwa += self.tool_path_obj.Python + " "
        map_read_bwa += self.tool_path_obj.Map_reads_gene_BWA + " "
        map_read_bwa += str(self.tool_path_obj.BWA_cigar_cutoff) + " "
        map_read_bwa += ref_path + " " #self.tool_path_obj.DNA_DB + " "  # IN
        if(self.sequence_contigs == "None"):
            map_read_bwa += "None" + " "
        else:        
            map_read_bwa += os.path.join(dep_loc, "contig_map.tsv") + " "  # IN
        map_read_bwa += os.path.join(final_folder, sample_root_name + "_" + ref_tag + "_gene_map.tsv") + " "  # OUT
        map_read_bwa += os.path.join(final_folder, sample_root_name + "_" + ref_tag + "_mapped_genes.fna") + " " #OUT
        map_read_bwa += reads_in + " "
        map_read_bwa += bwa_in + " "
        map_read_bwa += reads_out
        
        make_marker = ">&2 echo bwa pp complete: " + marker_file + " | " 
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_Annotate_BWA = [
            map_read_bwa + " && " + make_marker
            #copy_contig_map
        ]
        return COMMANDS_Annotate_BWA

    def create_BWA_copy_contig_map_command(self, stage_name, dependency_stage_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        bwa_folder      = os.path.join(data_folder, "1_bwa")
        split_folder    = os.path.join(data_folder, "0_read_split")
        final_folder    = os.path.join(subfolder, "final_results")
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
    
        copy_contig_map = ">&2 echo " + str(dt.today()) + " copy contig map | "
        copy_contig_map += "cp " + os.path.join(dep_loc, "contig_map.tsv") + " " + os.path.join(final_folder, "contig_map.tsv")
        
        make_marker = ">&2 echo bwa copy contig map complete: " + marker_file + " | " 
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [copy_contig_map + " && " + make_marker]

    def create_BLAT_annotate_command_v2(self, stage_name, query_file, fasta_db, marker_file):
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
        blat_command += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_command += self.tool_path_obj.DNA_DB_Split + fasta_db + " "
        blat_command += query_file
        blat_command += " -fine -q=rna -t=dna -out=blast8 -threads=2" + " "
        blat_command += os.path.join(blat_folder, sample_root_name + "_" + fasta_db + ".blatout")
        
        #make_marker = ">&2 echo marking BLAT job complete: " + sample_root_name + "_" + fasta_db + " | " 
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
        blat_merge_folder   = os.path.join(data_folder, "1_blat_merge")
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_merge_folder)
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
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]
        
        subfolder           = os.path.join(self.Output_Path, stage_name)
        data_folder         = os.path.join(subfolder, "data")
        blat_folder         = os.path.join(data_folder, "1_blat_merge")
        final_folder        = os.path.join(subfolder, "final_results")
        dep_loc             = os.path.join(self.Output_Path, dependency_stage_name, "final_results")  # implied to be BWA
        jobs_folder         = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        blat_pp = ">&2 echo " + str(dt.today()) + " BLAT post-processing " + sample_root_name + " | "
        blat_pp += self.tool_path_obj.Python + " "
        blat_pp += self.tool_path_obj.Map_reads_gene_BLAT + " "
        blat_pp += str(self.tool_path_obj.BLAT_identity_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_length_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_score_cutoff) + " "
        blat_pp += ref_file + " " #self.tool_path_obj.DNA_DB + " "
        
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
        diamond_annotate += self.tool_path_obj.DIAMOND
        diamond_annotate += " blastx -p " + self.Threads_str
        diamond_annotate += " -d " + self.tool_path_obj.Prot_DB
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
        diamond_pp += self.tool_path_obj.Python + " "
        diamond_pp += self.tool_path_obj.Map_reads_prot_DMND + " "
        diamond_pp += str(self.tool_path_obj.DIAMOND_identity_cutoff) + " "
        diamond_pp += str(self.tool_path_obj.DIAMOND_length_cutoff) + " "
        diamond_pp += str(self.tool_path_obj.DIAMOND_score_cutoff) + " "
        diamond_pp += self.tool_path_obj.Prot_DB_reads + " "                # IN
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
        dep_0_path      = os.path.join(self.Output_Path, dep_0_name, "final_results")
        dep_1_path      = os.path.join(self.Output_Path, dep_1_name, "final_results")
        dep_2_path      = os.path.join(self.Output_Path, dep_2_name, "final_results")
        dep_3_path      = os.path.join(self.Output_Path, dep_3_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        
        final_merge = self.tool_path_obj.Python + " "
        final_merge += self.tool_path_obj.GA_final_merge + " "
        final_merge += dep_3_path + " "
        final_merge += dep_0_path + " "
        final_merge += dep_1_path + " "
        final_merge += dep_2_path + " "
        final_merge += data_folder + " "
        final_merge += final_folder + " "
        final_merge += self.read_mode + " "
        final_merge += jobs_folder
        
        make_marker = ">&2 echo " + str(dt.today()) + " GA final merge | "
        make_marker += "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        
        COMMANDS_ga_final_merge = [
            final_merge + " && " + make_marker
        ]
        
        return COMMANDS_ga_final_merge
        

    def create_TA_kaiju_command(self, current_stage_name, assemble_contigs_stage, operating_mode, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        assemble_contigs_folder = os.path.join(self.Output_Path, assemble_contigs_stage, "final_results")
        kaiju_folder            = os.path.join(data_folder, "1_kaiju")
        jobs_folder             = os.path.join(data_folder, "jobs")
        
        self.make_folder(kaiju_folder)
        self.make_folder(jobs_folder)
        
        if(operating_mode == "contigs"):
            kaiju_on_contigs = ">&2 echo kaiju on contigs | "
            kaiju_on_contigs += self.tool_path_obj.Kaiju
            kaiju_on_contigs += " -t " + self.tool_path_obj.nodes
            kaiju_on_contigs += " -f " + self.tool_path_obj.Kaiju_db
            if(self.tutorial_keyword == "TA"):
                kaiju_on_contigs += " -i " + self.sequence_contigs
            else:
                kaiju_on_contigs += " -i " + os.path.join(assemble_contigs_folder, "contigs.fasta")
            kaiju_on_contigs += " -z " + self.Threads_str
            kaiju_on_contigs += " -o " + os.path.join(kaiju_folder, "contigs.tsv")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [kaiju_on_contigs + " && " + make_marker]

        elif(operating_mode == "singletons"):
            kaiju_on_singletons = ">&2 echo kaiju on singletons | "
            kaiju_on_singletons += self.tool_path_obj.Kaiju
            kaiju_on_singletons += " -t " + self.tool_path_obj.nodes
            kaiju_on_singletons += " -f " + self.tool_path_obj.Kaiju_db
            if(self.tutorial_keyword == "TA"):
                kaiju_on_singletons += " -i " + self.sequence_single
            else:
                kaiju_on_singletons += " -i " + os.path.join(assemble_contigs_folder, "singletons.fastq")
            kaiju_on_singletons += " -z " + self.Threads_str
            kaiju_on_singletons += " -o " + os.path.join(kaiju_folder, "singletons.tsv")

            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [kaiju_on_singletons + " && " + make_marker]
            
        elif(operating_mode == "paired"):
            kaiju_on_paired = ">&2 echo kaiju on pairs | "
            kaiju_on_paired += self.tool_path_obj.Kaiju
            kaiju_on_paired += " -t " + self.tool_path_obj.nodes
            kaiju_on_paired += " -f " + self.tool_path_obj.Kaiju_db

            if(self.tutorial_keyword == "TA"):
                kaiju_on_paired += " -i " + self.sequence_path_1
                kaiju_on_paired += " -j " + self.sequence_path_2

            else:
                kaiju_on_paired += " -i " + os.path.join(assemble_contigs_folder, "pair_1.fastq")
                kaiju_on_paired += " -j " + os.path.join(assemble_contigs_folder, "pair_2.fastq")

            kaiju_on_paired += " -z " + self.Threads_str
            kaiju_on_paired += " -o " + os.path.join(kaiju_folder, "pairs.tsv")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [kaiju_on_paired + " && " + make_marker]
            
            
    def create_TA_kaiju_pp_command(self, current_stage_name, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        kaiju_folder            = os.path.join(data_folder, "1_kaiju")
        jobs_folder             = os.path.join(data_folder, "jobs")

        self.make_folder(jobs_folder)
        
        cat_kaiju = ">&2 echo merging all kaiju results | "
        cat_kaiju += "cat "
        if(self.sequence_contigs != "None"):
            cat_kaiju += os.path.join(kaiju_folder, "contigs.tsv") + " "
        cat_kaiju += os.path.join(kaiju_folder, "singletons.tsv")
        if self.read_mode == "paired":
            cat_kaiju += " " + os.path.join(kaiju_folder, "pairs.tsv")
        cat_kaiju += " > " + os.path.join(kaiju_folder, "merged_kaiju.tsv")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [cat_kaiju + " && " + make_marker]
        
        
    def create_TA_centrifuge_command(self, current_stage_name, rRNA_stage, assemble_contigs_stage, operating_mode, marker_file):
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
            patch_contig_name = self.tool_path_obj.Python + " "
            patch_contig_name += self.tool_path_obj.ta_contig_name_convert + " "
            if(self.tutorial_keyword == "TA"):
                patch_contig_name += self.sequence_contigs + " "
            else:
                patch_contig_name += os.path.join(assemble_contigs_folder, "contigs.fasta") + " "
            patch_contig_name += os.path.join(centrifuge_folder, "contigs_renamed.fasta")
        
            centrifuge_on_contigs = ">&2 echo centrifuge on contigs | "
            centrifuge_on_contigs += self.tool_path_obj.Centrifuge
            centrifuge_on_contigs += " -f -x " + self.tool_path_obj.Centrifuge_db
            centrifuge_on_contigs += " -U " + os.path.join(centrifuge_folder, "contigs_renamed.fasta")
            centrifuge_on_contigs += " --exclude-taxids 2759 -k 1 --tab-fmt-cols " + "score,readID,taxID"
            centrifuge_on_contigs += " --phred" + self.Qual_str
            centrifuge_on_contigs += " -p 6"
            centrifuge_on_contigs += " -S " + os.path.join(centrifuge_folder, "raw_contigs.tsv")
            centrifuge_on_contigs += " --report-file " + os.path.join(centrifuge_folder, "raw_contigs.txt")
            
            back_convert_report = self.tool_path_obj.Python + " "
            back_convert_report += self.tool_path_obj.ta_contig_name_convert + " "
            back_convert_report += os.path.join(centrifuge_folder, "raw_contigs.tsv") + " "
            back_convert_report += os.path.join(centrifuge_folder, "contigs.tsv")
            
            make_marker = "touch" + " "
            make_marker += os.path.join(jobs_folder, marker_file)
            
            return [patch_contig_name + " && " + centrifuge_on_contigs + " && " + back_convert_report + " && " +  make_marker]

            
        elif(operating_mode == "reads"):
            centrifuge_on_reads = ">&2 echo centrifuge on reads | "
            centrifuge_on_reads += self.tool_path_obj.Centrifuge
            centrifuge_on_reads += " -x " + self.tool_path_obj.Centrifuge_db
            
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
            centrifuge_on_rRNA += self.tool_path_obj.Centrifuge
            centrifuge_on_rRNA += " -x " + self.tool_path_obj.Centrifuge_db
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
        get_taxa_from_gene += self.tool_path_obj.Python + " "
        get_taxa_from_gene += self.tool_path_obj.Annotated_taxid + " "  # SLOW STEP
        get_taxa_from_gene += os.path.join(final_merge_folder, "gene_map.tsv") + " "
        get_taxa_from_gene += self.tool_path_obj.accession2taxid + " "
        get_taxa_from_gene += os.path.join(ga_taxa_folder, "ga_taxon.tsv")
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        return [get_taxa_from_gene + " && " + make_marker]

    
    def create_TA_final_command(self, current_stage_name, assemble_contigs_stage, marker_file):
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        assemble_contigs_folder = os.path.join(self.Output_Path, assemble_contigs_stage, "final_results")
        ga_taxa_folder          = os.path.join(data_folder, "0_gene_taxa")
        kaiju_folder            = os.path.join(data_folder, "1_kaiju")
        centrifuge_folder       = os.path.join(data_folder, "2_centrifuge")
        wevote_folder           = os.path.join(data_folder, "3_wevote")
        final_folder            = os.path.join(subfolder, "final_results")
        jobs_folder             = os.path.join(data_folder, "jobs")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(ga_taxa_folder)
        self.make_folder(kaiju_folder)
        self.make_folder(centrifuge_folder)
        self.make_folder(wevote_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)

        wevote_combine = ">&2 echo combining classification outputs for wevote | "
        wevote_combine += self.tool_path_obj.Python + " "
        wevote_combine += self.tool_path_obj.Classification_combine + " "
        wevote_combine += os.path.join(assemble_contigs_folder, "contig_map.tsv")
        wevote_combine += " " + os.path.join(wevote_folder, "wevote_ensemble.csv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(ga_taxa_folder, "ga_taxon.tsv") + " "
        wevote_combine += os.path.join(kaiju_folder, "merged_kaiju.tsv") + " "
        wevote_combine += os.path.join(centrifuge_folder, "merged_centrifuge.tsv")        

        wevote_call = ">&2 echo Running WEVOTE | "
        wevote_call += self.tool_path_obj.WEVOTE
        wevote_call += " -i " + os.path.join(wevote_folder, "wevote_ensemble.csv")
        wevote_call += " -d " + self.tool_path_obj.WEVOTEDB
        wevote_call += " -p " + os.path.join(wevote_folder, "wevote")
        wevote_call += " -n " + self.Threads_str
        wevote_call += " -k " + "2"
        wevote_call += " -a " + "0"
        wevote_call += " -s " + "0"
        
        wevote_collect = ">&2 echo gathering WEVOTE results | "
        wevote_collect += self.tool_path_obj.Python + " "
        wevote_collect += self.tool_path_obj.Wevote_parser + " "
        wevote_collect += os.path.join(wevote_folder, "wevote_WEVOTE_Details.txt") + " "
        wevote_collect += os.path.join(final_folder, "taxonomic_classifications.tsv")
        
        constrain = ">&2 echo Constraining the Taxonomic Annotation | " 
        constrain += self.tool_path_obj.Python + " " + self.tool_path_obj.Constrain_classification + " "
        constrain += self.tool_path_obj.target_rank + " "
        constrain += os.path.join(final_folder, "taxonomic_classifications.tsv") + " "
        constrain += self.tool_path_obj.nodes + " "
        constrain += self.tool_path_obj.names + " "
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
        detect_protein += self.tool_path_obj.Python + " "
        detect_protein += self.tool_path_obj.Detect + " "
        detect_protein += os.path.join(final_merge_folder,"all_proteins.faa")
        detect_protein += " --output_file " + os.path.join(detect_folder, "proteins.detect")
        detect_protein += " --fbeta " + os.path.join(detect_folder, "proteins.fbeta")
        detect_protein += " --db " + self.tool_path_obj.DetectDB
        detect_protein += " --blastp " + self.tool_path_obj.Blastp
        detect_protein += " --needle " + self.tool_path_obj.Needle
        detect_protein += " --dump_dir " + detect_folder 
        detect_protein += " --n_count" + " " + str(self.tool_path_obj.DETECT_job_limit)
        detect_protein += " --mem_limit" + " " + str(self.tool_path_obj.DETECT_mem_threshold) 
        detect_protein += " --job_delay" + " " + str(self.tool_path_obj.DETECT_job_delay)
        detect_protein += " >> " + os.path.join(detect_folder, "detect_out.txt") + " 2>&1"

        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_DETECT = [
            detect_protein + " && " + make_marker
        ]

        return COMMANDS_DETECT

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
        
        
        PRIAM_command += self.tool_path_obj.Java + " "
        PRIAM_command += self.tool_path_obj.Priam
        PRIAM_command += " -n " + "proteins_priam" + " "
        PRIAM_command += " -i " + os.path.join(final_merge_folder, "all_proteins.faa")
        PRIAM_command += " -p " + self.tool_path_obj.PriamDB
        PRIAM_command += " -o " + PRIAM_folder
        PRIAM_command += " --np " + self.Threads_str
        PRIAM_command += " --bh --cc --cg --bp --bd "
        PRIAM_command += self.tool_path_obj.BLAST_dir
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)

        COMMANDS_PRIAM = [
            PRIAM_command + " && " + make_marker
        ]

        return COMMANDS_PRIAM
        
        
    def create_EC_DIAMOND_command(self, current_stage_name, ga_final_merge_stage, marker_file):
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        final_merge_folder  = os.path.join(self.Output_Path, ga_final_merge_stage, "final_results")
        diamond_ea_folder   = os.path.join(data_folder, "2_diamond")
        jobs_folder         = os.path.join(data_folder, "jobs")
        
        self.make_folder(diamond_ea_folder)
        self.make_folder(jobs_folder)
        
        diamond_ea_command = ">&2 echo running Diamond enzyme annotation | "
        diamond_ea_command += self.tool_path_obj.DIAMOND + " blastp"
        diamond_ea_command += " -p " + self.Threads_str
        diamond_ea_command += " --query " + os.path.join(final_merge_folder, "all_proteins.faa")
        diamond_ea_command += " --db " + self.tool_path_obj.SWISS_PROT
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
        postprocess_command += self.tool_path_obj.Python + " "
        postprocess_command += self.tool_path_obj.EC_Annotation_Post + " "
        postprocess_command += os.path.join(detect_folder, "proteins.fbeta") + " "
        postprocess_command += os.path.join(PRIAM_folder, "PRIAM_proteins_priam", "ANNOTATION", "sequenceECs.txt") + " "
        postprocess_command += os.path.join(diamond_ea_folder, "proteins.blastout") + " "
        postprocess_command += self.tool_path_obj.SWISS_PROT_map + " "
        postprocess_command += os.path.join(final_merge_folder, "gene_map.tsv") + " "
        postprocess_command += self.tool_path_obj.enzyme_db + " "
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
        network_generation += self.tool_path_obj.Python + " "
        network_generation += self.tool_path_obj.RPKM + " "
        network_generation += str(self.tool_path_obj.RPKM_cutoff) + " "
        network_generation += "None" + " "
        network_generation += self.tool_path_obj.nodes + " "
        network_generation += self.tool_path_obj.names + " "
        network_generation += gene_map_location + " "
        network_generation += os.path.join(ta_folder, "taxonomic_classifications.tsv") + " "
        network_generation += os.path.join(ea_folder, "proteins.ECs_All") + " "
        network_generation += self.tool_path_obj.show_unclassified + " "
        network_generation += os.path.join(final_folder, "RPKM_table.tsv") + " "
        network_generation += os.path.join(final_folder, "Cytoscape_network.tsv") + " "
        
        
        
        flatten_rpkm = ">&2 echo Reformat RPKM for EC heatmap | "
        flatten_rpkm += self.tool_path_obj.Python + " "
        flatten_rpkm += self.tool_path_obj.format_RPKM + " "
        flatten_rpkm += os.path.join(final_folder, "RPKM_table.tsv") + " "
        flatten_rpkm += os.path.join(final_folder, "EC_heatmap_RPKM.tsv")
        
        return [network_generation, flatten_rpkm]
        
    def create_output_unique_hosts_singletons_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_host_reads_singletons += self.tool_path_obj.Python + " "
        get_unique_host_reads_singletons += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_host_reads_singletons += os.path.join(host_folder, "singletons.fastq") + " "
        get_unique_host_reads_singletons += os.path.join(quality_folder, "singletons.fastq") + " "
        get_unique_host_reads_singletons += os.path.join(unique_hosts_folder, "singletons_hosts.fastq")
        
        
        repop_singletons_hosts = ">&2 echo repopulating singletons hosts | " 
        repop_singletons_hosts += self.tool_path_obj.Python + " "
        repop_singletons_hosts += self.tool_path_obj.duplicate_repopulate + " "
        if(self.read_mode == "single"):
            repop_singletons_hosts += os.path.join(quality_folder, "singletons_hq.fastq") + " "
        else:
            repop_singletons_hosts += os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
        repop_singletons_hosts += os.path.join(unique_hosts_folder, "singletons_hosts.fastq") + " "
        repop_singletons_hosts += os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
        repop_singletons_hosts += os.path.join(full_hosts_folder, "singletons_full_hosts.fastq")
        
        return [get_unique_host_reads_singletons, repop_singletons_hosts]
        
    def create_output_unique_hosts_pair_1_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_host_reads_pair_1 += self.tool_path_obj.Python + " "
        get_unique_host_reads_pair_1 += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_host_reads_pair_1 += os.path.join(host_folder, "pair_1.fastq") + " "
        get_unique_host_reads_pair_1 += os.path.join(quality_folder, "pair_1.fastq") + " "
        get_unique_host_reads_pair_1 += os.path.join(unique_hosts_folder, "pair_1_hosts.fastq")
        
        repop_pair_1_hosts = ">&2 echo repopulating pair 1 hosts | " 
        repop_pair_1_hosts += self.tool_path_obj.Python + " "
        repop_pair_1_hosts += self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1_hosts += os.path.join(quality_folder, "pair_1_match.fastq") + " "
        repop_pair_1_hosts += os.path.join(unique_hosts_folder, "pair_1_hosts.fastq") + " "
        repop_pair_1_hosts += os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_hosts += os.path.join(full_hosts_folder, "pair_1_full_hosts.fastq")
        
        return [get_unique_host_reads_pair_1, repop_pair_1_hosts]
        
    def create_output_unique_hosts_pair_2_command(self, current_stage_name, quality_stage, host_stage):
        #only call if we had hosts to filter
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_host_reads_pair_2 += self.tool_path_obj.Python + " "
        get_unique_host_reads_pair_2 += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_host_reads_pair_2 += os.path.join(host_folder, "pair_2.fastq") + " "
        get_unique_host_reads_pair_2 += os.path.join(quality_folder, "pair_2.fastq") + " "
        get_unique_host_reads_pair_2 += os.path.join(unique_hosts_folder, "pair_2_hosts.fastq")
        
        repop_pair_2_hosts = ">&2 echo repopulating pair 2 hosts | " 
        repop_pair_2_hosts += self.tool_path_obj.Python + " "
        repop_pair_2_hosts += self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2_hosts += os.path.join(quality_folder, "pair_2_match.fastq") + " "
        repop_pair_2_hosts += os.path.join(unique_hosts_folder, "pair_2_hosts.fastq") + " "
        repop_pair_2_hosts += os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " " #we do this based on pairs now
        repop_pair_2_hosts += os.path.join(full_hosts_folder, "pair_2_full_hosts.fastq")
        
        return [get_unique_host_reads_pair_2, repop_pair_2_hosts]
#-------------------------------------------------------------------------------------------
    def create_output_unique_vectors_singletons_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        quality_folder          = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_vectors_reads_singletons += self.tool_path_obj.Python + " "
        get_unique_vectors_reads_singletons += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_vectors_reads_singletons += os.path.join(vectors_folder, "singletons.fastq") + " "
        get_unique_vectors_reads_singletons += os.path.join(host_folder, "singletons.fastq") + " "
        get_unique_vectors_reads_singletons += os.path.join(unique_vectors_folder, "singletons_vectors.fastq")
        
        
        repop_singletons_vectors = ">&2 echo repopulating singletons vectors | " 
        repop_singletons_vectors += self.tool_path_obj.Python + " "
        repop_singletons_vectors += self.tool_path_obj.duplicate_repopulate + " "
        if(self.read_mode == "single"):
            repop_singletons_vectors += os.path.join(quality_folder, "singletons_hq.fastq") + " "
        else:
            repop_singletons_vectors += os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
        repop_singletons_vectors += os.path.join(unique_vectors_folder, "singletons_vectors.fastq") + " "
        repop_singletons_vectors += os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
        repop_singletons_vectors += os.path.join(full_vectors_folder, "singletons_full_vectors.fastq")
        
        return [get_unique_vectors_reads_singletons, repop_singletons_vectors]
        
    def create_output_unique_vectors_pair_1_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        quality_folder          = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_vectors_reads_pair_1 += self.tool_path_obj.Python + " "
        get_unique_vectors_reads_pair_1 += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_vectors_reads_pair_1 += os.path.join(vectors_folder, "pair_1.fastq") + " "
        get_unique_vectors_reads_pair_1 += os.path.join(host_folder, "pair_1.fastq") + " "
        get_unique_vectors_reads_pair_1 += os.path.join(unique_vectors_folder, "pair_1_vectors.fastq")
        
        repop_pair_1_vectors = ">&2 echo repopulating pair 1 vectors | " 
        repop_pair_1_vectors += self.tool_path_obj.Python + " "
        repop_pair_1_vectors += self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1_vectors += os.path.join(quality_folder, "pair_1_match.fastq") + " "
        repop_pair_1_vectors += os.path.join(unique_vectors_folder, "pair_1_vectors.fastq") + " "
        repop_pair_1_vectors += os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " "
        repop_pair_1_vectors += os.path.join(full_vectors_folder, "pair_1_full_vectors.fastq")
        
        return [get_unique_vectors_reads_pair_1, repop_pair_1_vectors]
        
    def create_output_unique_vectors_pair_2_command(self, current_stage_name, quality_stage, host_stage, vectors_stage):
        #only call if we had hosts to filter
        subfolder               = os.path.join(self.Output_Path, current_stage_name)
        data_folder             = os.path.join(subfolder, "data")
        quality_folder          = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        get_unique_vectors_reads_pair_2 += self.tool_path_obj.Python + " "
        get_unique_vectors_reads_pair_2 += self.tool_path_obj.get_unique_host_reads + " "
        get_unique_vectors_reads_pair_2 += os.path.join(vectors_folder, "pair_2.fastq") + " "
        get_unique_vectors_reads_pair_2 += os.path.join(host_folder, "pair_2.fastq") + " "
        get_unique_vectors_reads_pair_2 += os.path.join(unique_vectors_folder, "pair_2_vectors.fastq")
        
        repop_pair_2_vectors = ">&2 echo repopulating pair 2 vectors | " 
        repop_pair_2_vectors += self.tool_path_obj.Python + " "
        repop_pair_2_vectors += self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2_vectors += os.path.join(quality_folder, "pair_2_match.fastq") + " "
        repop_pair_2_vectors += os.path.join(unique_vectors_folder, "pair_2_vectors.fastq") + " "
        repop_pair_2_vectors += os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " " #we do this based on pairs now
        repop_pair_2_vectors += os.path.join(full_vectors_folder, "pair_2_full_vectors.fastq")
        
        return [get_unique_vectors_reads_pair_2, repop_pair_2_vectors]
        

        
    def create_output_per_read_scores_command(self, current_stage_name, quality_stage):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, quality_stage, "final_results")
        data_folder         = os.path.join(subfolder, "data")
        final_folder        = os.path.join(subfolder, "final_results")

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        per_read_scores = ">&2 echo collecting per-read quality | " 
        per_read_scores += self.tool_path_obj.Python + " "
        per_read_scores += self.tool_path_obj.read_quality_metrics + " "
        if(self.read_mode == "single"):
            per_read_scores += "single" + " "
            per_read_scores += self.sequence_single + " "
            per_read_scores += os.path.join(quality_folder, "singletons_hq.fastq") + " "
            per_read_scores += os.path.join(final_folder)
            
        elif(self.read_mode == "paired"):
            per_read_scores += "paired" + " " 
            per_read_scores += self.sequence_path_1 + " "
            per_read_scores += self.sequence_path_2 + " "
            per_read_scores += os.path.join(quality_folder, "pair_1_match.fastq") + " "
            per_read_scores += os.path.join(quality_folder, "pair_2_match.fastq") + " "
            per_read_scores += os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
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
        contig_stats += self.tool_path_obj.Python + " "
        contig_stats += self.tool_path_obj.contig_stats + " "
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
        EC_heatmap += self.tool_path_obj.Python + " "
        EC_heatmap += self.tool_path_obj.ec_heatmap + " "
        EC_heatmap += self.tool_path_obj.EC_pathway + " "
        EC_heatmap += os.path.join(final_folder, "EC_heatmap_RPKM.tsv") + " "
        EC_heatmap += self.tool_path_obj.path_to_superpath + " "
        EC_heatmap += final_folder
        
        return [EC_heatmap]
        
        
        
    def create_output_read_count_command(self, current_stage_name, quality_stage, repopulation_stage, ga_final_merge_stage, enzyme_annotation_stage):
        #only call if we had hosts to filter, and run it after the host regen is complete.
        subfolder           = os.path.join(self.Output_Path, current_stage_name)
        data_folder         = os.path.join(subfolder, "data")
        quality_folder      = os.path.join(self.Output_Path, quality_stage, "final_results")
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
        read_counts += self.tool_path_obj.Python + " "
        read_counts += self.tool_path_obj.read_count + " "
        if self.read_mode == "single":
            read_counts += self.sequence_single + " "
            
        elif self.read_mode == "paired":
            read_counts += self.sequence_path_1 + " "
        read_counts += quality_folder + " "
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
        taxa_groupby += self.tool_path_obj.Python + " "
        taxa_groupby += self.tool_path_obj.taxa_table + " "
        taxa_groupby += os.path.join(final_folder, "taxa_classifications.tsv") + " "
        taxa_groupby += os.path.join(final_folder, "taxa_summary.tsv")
        
        return [taxa_groupby]