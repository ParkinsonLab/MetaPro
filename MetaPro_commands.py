# The functions here generate the pipeline commands.
# The functions here generate the pipeline commands.
# Each command module is made up of sub stages that are used to get the final result.

import os
from datetime import datetime as dt
import time

class mt_pipe_commands:
    # --------------------------------------------------------------------
    # constructor:
    # there should only be one of these objects used for an entire pipeline.
    def __init__(self, config_dict, dir_dict, file_dict, marker_dict): #no_host, config_obj, Quality_score=33, tutorial_keyword = None, self.config_dict["pair_1"]=None, self.config_dict["pair_2"]=None, self.config_dict["single"]=None, sequence_contigs = None):
        
        self.config_dict = config_dict
        self.file_dict = file_dict
        self.dir_dict = dir_dict
        self.marker_dict = marker_dict
        
        #self.config_dict = config_dict
        #self.tool_path_obj = config_obj #mpp.tool_path_obj(Config_path)
        self.no_host_flag = self.config_dict["no_host"]
        # path to the genome sequence file
        
        self.tutorial_keyword = self.config_dict["tutorial_keyword"]

        print("tutorial keyword:", self.tutorial_keyword)
        time.sleep(4)
        self.read_mode = self.config_dict["read_mode"]        

                
        self.q_enc = self.config_dict["q_enc"]
        self.output_path = dir_dict["out"]
        self.threads_str = str(self.config_dict["num_threads"])
        self.thread_count = self.config_dict["num_threads"]
        self.DNA_DB = self.config_dict["DNA_DB"]
        
        
        print("Output filepath:", self.output_path)

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
                
    def create_quality_control_command(self, marker, q_enc):
        sort_pair_1 = ">&2 echo Sorting pair 1 | "
        sort_pair_1 += self.config_dict["Python"] + " "
        sort_pair_1 += self.config_dict["sort_reads"] + " "
        sort_pair_1 += self.file_dict["raw_p1"] + " "
        sort_pair_1 += self.file_dict["qf_sort_p1"] 
        sort_pair_2 = ">&2 echo Sorting pair 2 | "
        sort_pair_2 += self.config_dict["Python"] + " "
        sort_pair_2 += self.config_dict["sort_reads"] + " "
        sort_pair_2 += self.file_dict["raw_p2"] + " "
        sort_pair_2 += self.file_dict["qf_sort_p2"]


        trimgalore_trim = ">&2 echo Removing adapters | "
        trimgalore_trim += self.config_dict["trimgalore"]
        if self.read_mode == "single":
            trimgalore_trim += " " + self.file_dict["raw_s"]
        elif self.read_mode == "paired":
            trimgalore_trim += " --paired"
            trimgalore_trim += " " + self.file_dict["qf_sort_p1"]
            trimgalore_trim += " " + self.file_dict["qf_sort_p2"]
            trimgalore_trim += " --retain_unpaired"
        trimgalore_trim += " --quality " + self.config_dict["adapterremoval_minlength"]
        trimgalore_trim += " --cores " + self.threads_str
        trimgalore_trim += " --length " + self.config_dict["adapterremoval_minlength"]
        trimgalore_trim += " --output_dir " + self.dir_dict["qf_adapt"]
        trimgalore_trim += " --poly-g"
        trimgalore_trim += " --poly-a"
        trimgalore_trim += " --trim_n"
        trimgalore_trim += " --dont_gzip"

        if self.read_mode == "paired":
            # WARNING: filenames below assume TrimGalore's documented naming
            # convention (INPUT_unpaired_1.fq / INPUT_unpaired_2.fq, stem taken
            # from the input file minus .fastq/.fq/.gz) holds for this v2 Rust
            # build exactly as it does for v0.6.x. Not independently confirmed
            # against a real run of --retain_unpaired on this pipeline yet.
            trim_unpaired_1 = os.path.join(self.dir_dict["qf_adapt"], os.path.splitext(os.path.basename(self.file_dict["qf_sort_p1"]))[0] + "_unpaired_1.fq")
            trim_unpaired_2 = os.path.join(self.dir_dict["qf_adapt"], os.path.splitext(os.path.basename(self.file_dict["qf_sort_p2"]))[0] + "_unpaired_2.fq")
            trimgalore_trim += " && cat " + trim_unpaired_1 + " " + trim_unpaired_2 + " > " + self.file_dict["qf_adapt_s"]


        #Sort-reads introduces tags at the read-level of the 
        tag_remove_pair_1 = ">&2 echo Remove tags pair 1 | "
        tag_remove_pair_1 += self.config_dict["Python"] + " "
        tag_remove_pair_1 += self.config_dict["remove_tag"] + " "
        tag_remove_pair_1 += self.file_dict["qf_adapt_p1"] + " "
        tag_remove_pair_1 += self.file_dict["qf_tags_p1"]
        
        tag_remove_pair_2 = ">&2 echo Remove tags pair 2 | "
        tag_remove_pair_2 += self.config_dict["Python"] + " "
        tag_remove_pair_2 += self.config_dict["remove_tag"] + " "
        tag_remove_pair_2 += self.file_dict["qf_adapt_p2"] + " "
        tag_remove_pair_2 += self.file_dict["qf_tags_p2"]

        tag_remove_singletons =  ">&2 echo Remove tags singletons | " 
        tag_remove_singletons += self.config_dict["Python"] + " "
        tag_remove_singletons += self.config_dict["remove_tag"] + " "
        tag_remove_singletons += self.file_dict["qf_adapt_s"] + " "
        tag_remove_singletons += self.file_dict["qf_tags_s"]
        # tries to merge the cleaned pairs
        # rejects get sent out
        vsearch_merge = ">&2 echo " + "Vsearch Merge pairs | "
        vsearch_merge += self.config_dict["vsearch"]
        vsearch_merge += " --fastq_mergepairs " + self.file_dict["qf_tags_p1"]
        vsearch_merge += " --reverse " + self.file_dict["qf_tags_p2"]
        vsearch_merge += " --fastq_ascii " + self.config_dict["q_enc"]
        vsearch_merge += " --fastqout " + self.file_dict["qf_merge_s"]
        vsearch_merge += " --fastqout_notmerged_fwd " + self.file_dict["qf_merge_p1"]
        vsearch_merge += " --fastqout_notmerged_rev " + self.file_dict["qf_merge_p2"]

        # concatenate the merge overlaps with the singletons
        cat_glue = ">&2 echo concatenating singletons | "
        cat_glue += "cat "
        cat_glue += self.file_dict["qf_merge_s"] + " "
        cat_glue += self.file_dict["qf_tags_s"]
        cat_glue += " > " + self.file_dict["qf_merge_s2"]

        # Filter out low-quality reads
        # start with the singles / merged sections
        
        vsearch_filter_0 = ">&2 echo low-quality filter on singletons | "
        vsearch_filter_0 += self.config_dict["vsearch"]
        if self.read_mode == "single":
            vsearch_filter_0 += " --fastq_filter " + self.file_dict["qf_tags_s"]
        elif self.read_mode == "paired":
            vsearch_filter_0 += " --fastq_filter " + self.file_dict["qf_merge_s2"]
        vsearch_filter_0 += " --fastq_ascii " + self.config_dict["q_enc"]
        vsearch_filter_0 += " --fastq_maxee " + "2.0"
        vsearch_filter_0 += " --fastqout " + self.file_dict["qf_hq_s"]

        # then move onto the standalones in pair 1
        vsearch_filter_1 = ">&2 echo low-quality filter on pair 1 | "
        vsearch_filter_1 += self.config_dict["vsearch"]
        vsearch_filter_1 += " --fastq_filter " + self.file_dict["qf_merge_p1"]
        vsearch_filter_1 += " --fastq_ascii " + self.config_dict["q_enc"]
        vsearch_filter_1 += " --fastq_maxee " + "2.0"
        vsearch_filter_1 += " --fastqout " + self.file_dict["qf_hq_p1"]

        vsearch_filter_2 = ">&2 echo low-quality filter on pair 2 | "
        vsearch_filter_2 += self.config_dict["vsearch"]
        vsearch_filter_2 += " --fastq_filter " + self.file_dict["qf_merge_p2"]
        vsearch_filter_2 += " --fastq_ascii " + self.config_dict["q_enc"]
        vsearch_filter_2 += " --fastq_maxee " + "2.0"
        vsearch_filter_2 += " --fastqout " + self.file_dict["qf_hq_p2"]

        # redistribute data into singletons, or paired-reads
        orphan_read_filter = ">&2 echo moving newly orphaned reads | "
        orphan_read_filter += self.config_dict["Python"] + " "
        orphan_read_filter += self.config_dict["orphaned_read_filter"] + " "
        orphan_read_filter += self.file_dict["qf_hq_p1"] + " "
        orphan_read_filter += self.file_dict["qf_hq_p2"] + " "
        orphan_read_filter += self.file_dict["qf_hq_s"] + " "
        orphan_read_filter += self.file_dict["qf_o_p1"] + " "
        orphan_read_filter += self.file_dict["qf_o_p2"] + " "
        orphan_read_filter += self.file_dict["qf_o_s"]

        # remove duplicates (to shrink the data size)
        cdhit_singletons = ">&2 echo removing singleton duplicates | "
        cdhit_singletons += self.config_dict["cdhit_dup"] + " -i "
        if self.read_mode == "single":
            cdhit_singletons += self.file_dict["qf_hq_s"]
        elif self.read_mode == "paired":
            cdhit_singletons += self.file_dict["qf_o_s"]
        cdhit_singletons += " -o " + self.file_dict["qf_u_s"]

        # remove duplicates in the pairs
        cdhit_paired = ">&2 echo remove duplicates from paired | "
        cdhit_paired += self.config_dict["cdhit_dup"] + " "
        cdhit_paired += "-i"    + " " + self.file_dict["qf_o_p1"] + " "
        cdhit_paired += "-i2"   + " " + self.file_dict["qf_o_p2"] + " "
        cdhit_paired += "-o"    + " " + self.file_dict["qf_u_p1"] + " "
        cdhit_paired += "-o2"   + " " + self.file_dict["qf_u_p2"]

        make_marker = "touch " + marker
        
        if self.read_mode == "single":
            COMMANDS_qual = [
                trimgalore_trim,
                tag_remove_singletons,
                vsearch_filter_0,
                cdhit_singletons, 
                make_marker
            ]
        elif self.read_mode == "paired":
            COMMANDS_qual = [
                sort_pair_1 + " & " +
                sort_pair_2,
                trimgalore_trim,
                tag_remove_pair_1 + " & " +
                tag_remove_pair_2,
                tag_remove_singletons,
                ">&2 echo delaying for writes | sleep " + str(self.config_dict["qf_delay"]), 
                vsearch_merge,
                cat_glue,
                vsearch_filter_0 + " & " +
                vsearch_filter_1 + " & " +
                vsearch_filter_2,
                ">&2 echo delaying for writes | sleep " + str(self.config_dict["qf_delay"]), 
                orphan_read_filter,
                ">&2 echo delaying for writes | sleep " + str(self.config_dict["qf_delay"]), 
                cdhit_singletons,
                cdhit_paired + " && "  + make_marker
            ]

        return COMMANDS_qual

    def create_host_filter_command(self, host_id_list, host_count, marker):
        #may 09, 2025: removing python scripts. samtools can just do it all.

        cur_host = host_id_list[host_count]
        print(dt.today(), "COMMANDS| working on host:", cur_host)
        
        host_seq_1 = self.file_dict[cur_host + "_host_p1"]
        host_seq_2 = self.file_dict[cur_host + "_host_p2"]
        host_seq_s = self.file_dict[cur_host + "_host_s"]
        host_seq_o = self.file_dict[cur_host + "_host_o"]
        host_seq_u = self.file_dict[cur_host + "_host_u"]
        out_seq_s = self.file_dict[cur_host + "_no_host_s"]
        out_seq_o = self.file_dict[cur_host + "_no_host_o"]
        out_seq_u = self.file_dict[cur_host + "_no_host_u"]
        out_seq_1 = self.file_dict[cur_host + "_no_host_p1"]
        out_seq_2 = self.file_dict[cur_host + "_no_host_p2"]
        no_host_sam_s = self.file_dict[cur_host + "_no_host_s_sam"]
        
        no_host_sam_p = self.file_dict[cur_host + "_no_host_p_sam"]
        in_seq_1 = "none"
        in_seq_2 = "none"
        in_seq_s = "none"
        if(host_count == 0):
            in_seq_1 = self.file_dict["qf_u_p1"]
            in_seq_2 = self.file_dict["qf_u_p2"]
            in_seq_s = self.file_dict["qf_u_s"]
            
        else:
            prev_host = host_id_list[host_count - 1]
            in_seq_1 = self.file_dict[prev_host + "_no_host_p1"]
            in_seq_2 = self.file_dict[prev_host + "_no_host_p2"]
            in_seq_s = self.file_dict[prev_host + "_no_host_s"]
            
        if(host_count == (len(host_id_list) -1)):
            last_host = host_id_list[-1]
            out_seq_1 = self.file_dict[last_host + "_no_host_p1"]
            out_seq_2 = self.file_dict[last_host + "_no_host_p2"]
            out_seq_s = self.file_dict[last_host + "_no_host_s"]
            
        # host removal on unique singletons
        bt2_hr_s = ">&2 echo bt2 host remove on singletons | "
        bt2_hr_s += self.config_dict["BT2"] + " -p "
        bt2_hr_s += self.threads_str + " "
        bt2_hr_s += "-x " + os.path.join(self.config_dict["Host_db"], cur_host) + " "
        bt2_hr_s += "-U " + in_seq_s + " " 
        bt2_hr_s += "-S " + no_host_sam_s
        
        sam_no_host_s = self.config_dict["samtools"] + " view -f 4 -h "
        sam_no_host_s += no_host_sam_s + " | "
        sam_no_host_s += self.config_dict["samtools"] + " fastq - > "
        sam_no_host_s += out_seq_s
        
        sam_host_s = self.config_dict["samtools"] + " view -F 4 -h "
        sam_host_s += no_host_sam_s + " | " 
        sam_host_s += self.config_dict["samtools"] + " fastq - > "
        sam_host_s += host_seq_s
    

        bt2_hr_paired = ">&2 echo bt2 host-removal on paired | " 
        bt2_hr_paired += self.config_dict["BT2"] + " "
        bt2_hr_paired += "-p " + self.threads_str + " "
        bt2_hr_paired += "-x " + os.path.join(self.config_dict["Host_db"], cur_host) + " "
        bt2_hr_paired += "-1 " + in_seq_1 + " "
        bt2_hr_paired += "-2 " + in_seq_2 + " "
        bt2_hr_paired += "-S " + no_host_sam_p
        
        
        sam_no_host_p = self.config_dict["samtools"] + " view "
        if(self.config_dict["filter_stringency"] == "high"):
            #read and mate unmapped.  12
            sam_no_host_p += "-f 12 -h " 
        else:
            sam_no_host_p += "-f 4 -h " 
        sam_no_host_p += no_host_sam_p + " | "
        sam_no_host_p += self.config_dict["samtools"] + " fastq - "
        sam_no_host_p += "-1 " + out_seq_1 + " "
        sam_no_host_p += "-2 " + out_seq_2 + " "
        sam_no_host_p += "-s " + out_seq_o + " "
        sam_no_host_p += "-0 " + out_seq_u
        
        sam_host_p = self.config_dict["samtools"] + " view "
        if(self.config_dict["filter_stringency"] == "high"):
            sam_host_p += "-F 12 -h " 
        else:
            sam_host_p += "-F 4 -h "
        sam_host_p += no_host_sam_p + " | "
        sam_host_p += self.config_dict["samtools"] + " fastq - "
        sam_host_p += "-1 " + host_seq_1 + " "
        sam_host_p += "-2 " + host_seq_2 + " "
        sam_host_p += "-s " + host_seq_o + " "
        sam_host_p += "-0 " + host_seq_u
        
        # Create empty files for cases where samtools doesn't create some output files
        create_empty_files = "touch " + out_seq_u + " " + out_seq_o + " " + host_seq_o + " " + host_seq_u
        
        # Fixed cat commands using temporary files to avoid overwriting input
        cat_no_host_p = "cat "
        cat_no_host_p += out_seq_u + " "
        cat_no_host_p += out_seq_o + " "
        cat_no_host_p += out_seq_s + " > " + out_seq_s + "_temp"
        cat_no_host_p += " && mv " + out_seq_s + "_temp " + out_seq_s

        cat_host_p = "cat "
        cat_host_p += host_seq_o + " " 
        cat_host_p += host_seq_u + " "
        cat_host_p += host_seq_s + " > " + host_seq_s + "_temp"
        cat_host_p += " && mv " + host_seq_s + "_temp " + host_seq_s
    
        
        make_marker = "touch " + marker

        
        #-----------------------------
        #orphan correction
        
        

        
        if(self.tutorial_keyword is None or self.tutorial_keyword == "None"):
            if self.read_mode == "single":
                COMMANDS_host = [
                    bt2_hr_s,
                    sam_host_s, 
                    sam_no_host_s + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_host = [
                    bt2_hr_s,
                    bt2_hr_paired,
                    sam_no_host_s,
                    sam_host_s,
                    sam_no_host_p,
                    sam_host_p,
                    create_empty_files,
                    cat_host_p,
                    cat_no_host_p + " && " + make_marker

                ]
        elif(self.tutorial_keyword == "host"):
            print(dt.today(), "Host filter operating in tutorial-mode")
            if self.read_mode == "single":
                COMMANDS_host = [
                    bt2_hr_tut_s,
                    samtools_hr_s_sam_to_bam,
                    samtools_no_host_s_bam_to_fastq,
                    samtools_host_s_bam_to_fastq + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_host = [
                    bt2_hr_tut_s,
                    samtools_hr_s_sam_to_bam,
                    samtools_no_host_s_bam_to_fastq,
                    samtools_host_s_bam_to_fastq,
                    bt2_hr_tut_paired,
                    bt2_hr_filter_paired + " && " + make_marker
                ]
        else:
            print("tutorial keyword:", type(self.tutorial_keyword))
            print(self.tutorial_keyword)
                
        return COMMANDS_host

    def create_vector_filter_command(self, marker_file, final_host):
        # why do we leave all the interim files intact?
        # because science needs repeatable data, and the process needs to be able to start at any point
        
        bt2_vr_s = ">&2 echo bt2 vector oprhans | "
        bt2_vr_s += self.config_dict["BT2"] + " -p " + self.threads_str + " "
        bt2_vr_s += "-x " + os.path.join(self.config_dict["vector_db"], self.config_dict["vector_ID"]) + " "
        if(final_host == "none"):
            bt2_vr_s += "-U " + self.file_dict["qf_u_s"] + " "
        else:
            bt2_vr_s += "-U " + self.file_dict[final_host + "_no_host_s"] + " "
        bt2_vr_s += "-S " + self.file_dict["vec_s_sam"]
        
        
        bt2_vr_tut_s = ">&2 echo bt2 vector oprhans TUTORIAL MODE | "
        bt2_vr_tut_s += self.config_dict["BT2"] + " -p " + self.threads_str + " "
        bt2_vr_tut_s += "-x " + os.path.join(self.config_dict["vector_db"], self.config_dict["vector_ID"]) + " "
        bt2_vr_tut_s += "-U " + self.config_dict["single"] + " "
        bt2_vr_tut_s += "-S " + self.file_dict["vec_s_sam"]

        samtools_no_vec_s = ">&2 echo samtools vector oprhans pt 1 | "
        samtools_no_vec_s += self.config_dict["samtools"] + " view -h -f 4 "
        samtools_no_vec_s += self.file_dict["vec_s_sam"] + " | "
        samtools_no_vec_s += self.config_dict["samtools"] + " fastq - > "
        samtools_no_vec_s += self.file_dict["no_vec_s"] + " "
        
        samtools_vec_s = ">&2 echo samtools vector oprhans pt 1 | "
        samtools_vec_s += self.config_dict["samtools"] + " view -h -F 4 "
        samtools_vec_s += self.file_dict["vec_s_sam"] + " | "
        samtools_vec_s += self.config_dict["samtools"] + " fastq - > "
        samtools_vec_s += self.file_dict["vec_s"] + " "
        
        

        bt2_vr_paired = ">&2 echo bt2 vector paired | "
        bt2_vr_paired += self.config_dict["BT2"] + " -p " + self.threads_str + " "
        bt2_vr_paired += " -x " + os.path.join(self.config_dict["vector_db"], self.config_dict["vector_ID"]) + " "
        if(final_host != "none"):
            bt2_vr_paired += "-1 " + self.file_dict[final_host + "_no_host_p1"] + " "
            bt2_vr_paired += "-2 " + self.file_dict[final_host + "_no_host_p2"] + " "
        else:
            bt2_vr_paired += "-1 " + self.file_dict["qf_u_p1"] + " "
            bt2_vr_paired += "-2 " + self.file_dict["qf_u_p2"] + " "
        bt2_vr_paired += "-S " + self.file_dict["vec_p_sam"]

        bt2_vr_tut_paired = ">&2 echo bt2 vector paired TUTORIAL MODE | "
        bt2_vr_tut_paired += self.config_dict["BT2"] + " -p " + self.threads_str + " "
        bt2_vr_tut_paired += "-x " + os.path.join(self.config_dict["vector_db"], self.config_dict["vector_ID"]) + " "
        bt2_vr_tut_paired += "-1 " + self.config_dict["pair_1"] + " "
        bt2_vr_tut_paired += "-2 " + self.config_dict["pair_2"] + " "
        bt2_vr_tut_paired += "-S " + self.file_dict["vec_p_sam"]
        
        bt2_vr_filter_paired = ">&2 echo bt2 vector filter on paired | "
        bt2_vr_filter_paired += self.config_dict["Python"] + " "
        bt2_vr_filter_paired += self.config_dict["bt2_read_sorter"] + " "
        bt2_vr_filter_paired += "paired" + " "
        bt2_vr_filter_paired += self.config_dict["filter_stringency"] + " "
        bt2_vr_filter_paired += self.file_dict["vec_p_sam"] + " "
        if(final_host == "none"):
            bt2_vr_filter_paired += self.file_dict["qf_u_p1"] + " "
            bt2_vr_filter_paired += self.file_dict["qf_u_p2"] + " "

        else:
            bt2_vr_filter_paired += self.file_dict[final_host + "_no_host_p1"] + " "
            bt2_vr_filter_paired += self.file_dict[final_host + "_no_host_p2"] + " "
        bt2_vr_filter_paired += self.file_dict["no_vec_p1"] + " "
        bt2_vr_filter_paired += self.file_dict["no_vec_p2"] + " "
        bt2_vr_filter_paired += self.file_dict["vec_p1"] + " "
        bt2_vr_filter_paired += self.file_dict["vec_p2"]

        make_marker = "touch " + marker_file

        if(self.tutorial_keyword == "vectors" or self.tutorial_keyword == "vector"):
            if self.read_mode == "single":
                COMMANDS_vector = [
                    bt2_vr_tut_s,
                    samtools_no_vec_s,
                    samtools_vec_s + " && " + make_marker
                    
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    bt2_vr_tut_s,
                    samtools_no_vec_s,
                    samtools_vec_s,
                    bt2_vr_tut_paired,
                    bt2_vr_filter_paired + " && " + make_marker
                ]
        else:    
            if self.read_mode == "single":
                COMMANDS_vector = [
                    bt2_vr_s,
                    samtools_no_vec_s,
                    samtools_vec_s + " && " + make_marker
                ]
            elif self.read_mode == "paired":
                COMMANDS_vector = [
                    bt2_vr_s,
                    samtools_no_vec_s,
                    samtools_vec_s,
                    bt2_vr_paired,
                    bt2_vr_filter_paired + " && " + make_marker
                ]    

        return COMMANDS_vector
         
    def create_rRNA_filter_bnap_command(self, fasta_seqs, fastq_in, mRNA_out, rRNA_out, Barrnap_out, marker_file):
        # called by each split file
        # category -> singletons, pair 1, pair 2
        # file name -> the specific split section of the category (the fastq segments)
        # stage_name -> "rRNA_Filter"
        
        
        Barrnap_arc = self.config_dict["Barrnap"]
        Barrnap_arc += " --quiet --reject 0.01 --kingdom " + "arc"
        Barrnap_arc += " --threads " + self.threads_str
        Barrnap_arc += " " + fasta_seqs
        Barrnap_arc += " >> " + Barrnap_out


        Barrnap_bac = self.config_dict["Barrnap"]
        Barrnap_bac += " --quiet --reject 0.01 --kingdom " + "bac"
        Barrnap_bac += " --threads " + self.threads_str
        Barrnap_bac += " " + fasta_seqs
        Barrnap_bac += " >> " + Barrnap_out


        Barrnap_euk = self.config_dict["Barrnap"]
        Barrnap_euk += " --quiet --reject 0.01 --kingdom " + "euk"
        Barrnap_euk += " --threads " + self.threads_str
        Barrnap_euk += " " + fasta_seqs
        Barrnap_euk += " >> " + Barrnap_out

        Barrnap_mito = self.config_dict["Barrnap"]
        Barrnap_mito += " --quiet --reject 0.01 --kingdom " + "mito"
        Barrnap_mito += " --threads " + self.threads_str
        Barrnap_mito += " " + fasta_seqs
        Barrnap_mito += " >> " + Barrnap_out


        Barrnap_pp = ">&2 echo Running Barrnap pp scripts | "
        Barrnap_pp += self.config_dict["Python"] + " "
        Barrnap_pp += self.config_dict["barrnap_post"] + " "
        Barrnap_pp += Barrnap_out + " "
        Barrnap_pp += fastq_in + " "
        Barrnap_pp += mRNA_out + " "
        Barrnap_pp += rRNA_out + " "
        
        make_marker = "touch" + " " + marker_file
  

        return [Barrnap_arc, 
                Barrnap_bac,
                Barrnap_euk,
                Barrnap_mito,
                Barrnap_pp + " && " + make_marker
                ]
              
    def create_rRNA_filter_infernal_command(self, fasta_in, infernal_out, marker_file):
        
        infernal_command = self.config_dict["Infernal"]
        infernal_command += " -o /dev/null --tblout "
        infernal_command += infernal_out
        #infernal_command += " --cpu " + self.threads_str -> lined nerf'd because infernal's parallelism is not good
        infernal_command += " --cpu 1"
        infernal_command += " --anytrunc --rfam -E 0.001 "
        infernal_command += self.config_dict["Rfam"] + " "
        infernal_command += fasta_in

        
        #make_marker = ">&2 echo " + file_name + "_infernal Marking job completed | " 
        make_marker = "touch " + marker_file

        return [infernal_command + " && " + make_marker]
    
    def create_rRNA_inf_pp_pair_command(self, inf_p1, inf_p2, barrnap_p1, barrnap_p2, raw_p1, raw_p2, mRNA_p1, mRNA_p2, other_p1, other_p2, marker_file):
    #file name expected to have no extensions.  eg: pair_1_0
    #expected to be called for each category (pair1, singletons).  not pair 2.  paired data is handled in combination

        inf_pp = self.config_dict["Python"] + " "
        inf_pp += self.config_dict["rRNA_inf_pp"] + " "
        inf_pp += self.config_dict["filter_stringency"] + " "
        inf_pp += "paired" + " "
        inf_pp += inf_p1 + " "
        inf_pp += inf_p2 + " "
        inf_pp += barrnap_p1 + " "
        inf_pp += barrnap_p2 + " "
        inf_pp += raw_p1 + " "
        inf_pp += raw_p2 + " "
        inf_pp += mRNA_p1 + " "
        inf_pp += mRNA_p2 + " "
        inf_pp += other_p1 + " "
        inf_pp += other_p2


        make_marker = "touch " + marker_file
        return [inf_pp  + " && " + make_marker]
    
    def create_rRNA_inf_pp_s_command(self, inf_s, barrnap_s, raw_s, mRNA_s, other_s, marker_file):
       
        inf_pp = self.config_dict["Python"] + " "
        inf_pp += self.config_dict["rRNA_inf_pp"] + " "
        inf_pp += self.config_dict["filter_stringency"] + " "
        inf_pp += "single" + " "
        inf_pp += inf_s + " "
        inf_pp += barrnap_s + " "
        inf_pp += raw_s + " "
        inf_pp += mRNA_s + " "
        inf_pp += other_s
        
        make_marker = "touch" + " " + marker_file
            
        return [inf_pp + " && " + make_marker]
  
    def create_repop_command(self, marker_file):
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


        repop_singletons = ">&2 echo " + str(dt.today()) + " Duplication repopulation singletons mRNA| "
        repop_singletons += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"] + " "
        #the reference data to be drawn from 
        if self.read_mode == "single":
            repop_singletons += self.file_dict["qf_hq_s"] + " "
        elif self.read_mode == "paired":
            repop_singletons += self.file_dict["qf_o_s"] + " "
        repop_singletons += self.file_dict["rRNA_mRNA_s_fq"] + " "  # in -> rRNA filtration output
        repop_singletons += self.file_dict["qf_clstr_s"] + " "  # in -> duplicates filter output
        repop_singletons += self.file_dict["repop_s"] #OUT doesn't matter if it's single or paired.  
        
        #if self.read_mode == "single":
        #    repop_singletons += os.path.join(final_folder, "singletons.fastq")  # out
        #elif self.read_mode == "paired":
        #    repop_singletons += os.path.join(repop_folder, "singletons.fastq")  # out
            
            

        repop_singletons_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulations singletons rRNA | "
        repop_singletons_rRNA += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"] + " "
        if self.read_mode == "single":
            #repop_singletons_rRNA += os.path.join(singleton_path, "singletons_hq.fastq") + " "
            repop_singletons_rRNA += self.file_dict["qf_hq_s"] + " "
        elif self.read_mode == "paired":
            repop_singletons_rRNA += self.file_dict["qf_o_s"] + " "
        repop_singletons_rRNA += self.file_dict["rRNA_other_s_fq"] + " "  # in -> rRNA filtration output
        repop_singletons_rRNA += self.file_dict["qf_clstr_s"] + " "  # in -> duplicates filter output
        repop_singletons_rRNA += self.file_dict["repop_other_s"]  # out


        repop_pair_1 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 mRNA | "
        repop_pair_1 += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"] + " "
        repop_pair_1 += self.file_dict["qf_o_p1"] + " "
        if(self.tutorial_keyword == "repop"):
            repop_pair_1 += self.config_dict["pair_1"] + " "
        else:
            repop_pair_1 += self.file_dict["rRNA_mRNA_p1_fq"] + " "
        repop_pair_1 += self.file_dict["qf_clstr_p1"] + " "
        repop_pair_1 += self.file_dict["repop_p1"]

        repop_pair_1_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 1 rRNA | "
        repop_pair_1_rRNA += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"] + " "
        repop_pair_1_rRNA += self.file_dict["qf_o_p1"] + " "
        repop_pair_1_rRNA += self.file_dict["rRNA_other_p1_fq"] + " "
        repop_pair_1_rRNA += self.file_dict["qf_clstr_p1"] + " "
        repop_pair_1_rRNA += self.file_dict["repop_other_p1"]

        repop_pair_2 = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2 += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"]+ " "
        repop_pair_2 += self.file_dict["qf_o_p2"] + " "
        if(self.tutorial_keyword == "repop"):
            repop_pair_2 += self.config_dict["pair_2"] + " "
        else:
            repop_pair_2 += self.file_dict["rRNA_mRNA_p2_fq"] + " "
        repop_pair_2 += self.file_dict["qf_clstr_p1"] + " "
        repop_pair_2 += self.file_dict["repop_p2"]

        repop_pair_2_rRNA = ">&2 echo " + str(dt.today()) + " Duplication repopulation pair 2 | "
        repop_pair_2_rRNA += self.config_dict["Python"] + " " + self.config_dict["duplicate_repopulate"]+ " "
        repop_pair_2_rRNA += self.file_dict["qf_o_p2"] + " "
        repop_pair_2_rRNA += self.file_dict["rRNA_other_p2_fq"] + " "
        repop_pair_2_rRNA += self.file_dict["qf_clstr_p1"] + " "
        repop_pair_2_rRNA += self.file_dict["repop_other_p2"]

        singleton_repop_filter = ">&2 echo filtering mRNA for new singletons | "
        singleton_repop_filter += self.config_dict["Python"] + " "
        singleton_repop_filter += self.config_dict["orphaned_read_filter"] + " "
        singleton_repop_filter += self.file_dict["repop_p1"] + " "
        singleton_repop_filter += self.file_dict["repop_p2"] + " "
        singleton_repop_filter += self.file_dict["repop_s"]
    
        singleton_repop_filter_rRNA = ">&2 echo filtering rRNA for new singletons | "  
        singleton_repop_filter_rRNA += self.config_dict["Python"] + " "
        singleton_repop_filter_rRNA += self.config_dict["orphaned_read_filter"] + " "
        singleton_repop_filter_rRNA += self.file_dict["repop_other_p1"] + " "
        singleton_repop_filter_rRNA += self.file_dict["repop_other_p2"] + " "
        singleton_repop_filter_rRNA += self.file_dict["repop_other_s"]
        
        make_marker = "touch " + marker_file

        if(self.tutorial_keyword == "repop"):
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons + " && " + make_marker
                ]
            elif (self.read_mode == "paired") or (self.read_mode == "p"):
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_pair_1,
                    repop_pair_2,
                    singleton_repop_filter  + " && " + make_marker
                ]
        
        else:
            if self.read_mode == "single":
                COMMANDS_Repopulate = [
                    repop_singletons,
                    repop_singletons_rRNA  + " && " + make_marker
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
                    singleton_repop_filter_rRNA  + " && " + make_marker
                ]

        return COMMANDS_Repopulate
        
    def create_assemble_contigs_command(self, marker_file):
        
        tut_keyword = "assembly"
        
        # this assembles contigs
        spades = ">&2 echo Spades Contig assembly | "
        spades += self.config_dict["Python"] + " "
        spades += self.config_dict["Spades"] + " --rna"
        if(self.tutorial_keyword == tut_keyword):
            if self.read_mode == "paired":
                spades += " -1 " + self.config_dict["pair_1"]  # in1 (pair 1)
                spades += " -2 " + self.config_dict["pair_2"]  # in2 (pair 2)
            spades += " -s " + self.config_dict["single"]  # in_single (singletons)
        else:
            if self.read_mode == "paired":
                spades += " -1 " + self.file_dict["repop_p1"]  # in1 (pair 1)
                spades += " -2 " + self.file_dict["repop_p2"]  # in2 (pair 2)
            if(os.path.getsize(self.file_dict["repop_s"]) <= 0):
                print(dt.today(), "warning: empty singletons in SPAdes call.  MetaPro will ignore singletons")
            else:
                spades += " -s " + self.file_dict["repop_s"]  # in_single (singletons)
        spades += " -o " + self.dir_dict["contigs_spades"]  # out

        spades += " && touch " + self.file_dict["contigs_spades_done"]

        #if there is no output, bypass contigs. -> But this is a v2 upgrade.  
        #spades_rename = "cp " + self.file_dict["contigs_transcripts"] + " " + self.file_dict["contigs_og_fa"]  # rename output
        
        
        
       
        final_contigs       = self.file_dict["contigs_out_fa"]
        contig_map          = self.file_dict["contigs_map"]
        contigs_idx         = self.file_dict["contigs_idx"]
        
        #-------------------------------------------------------
        #spades does too good of a job sometimes.  Disassemble it into genes.
        disassemble_contigs = ">&2 echo Disassembling contigs | "
        disassemble_contigs += self.config_dict["mgm1"] + " -o " + self.file_dict["contigs_gene_report"] + " "
        disassemble_contigs += "-D " + self.file_dict["contigs_split"] + " "
        disassemble_contigs += "-m " + self.config_dict["mgm_model"] + " "
        disassemble_contigs += self.file_dict["contigs_transcripts"]
        
        
        remove_whitespace = ">&2 echo Removing whitespace from fasta | " 
        remove_whitespace += self.config_dict["Python"] + " " + self.config_dict["remove_gaps_in_fasta"] + " "
        remove_whitespace += self.file_dict["contigs_split"] + " "
        remove_whitespace += final_contigs
        
        #bt2-ing against the final contigs gives us a proper contig-segment -> read map. 
        #bt2_index = self.config_dict["BT2"] + " index -a bwtsw " + final_contigs
        #note: bt2 index MUST be fasta.  not fastq
        bt2_index = self.config_dict["BT2_index"] + " " + final_contigs + " " + contigs_idx 
        
        # Build a report of what was consumed by contig transmutation (assemble/disassemble)
        bt2_paired_contigs = ">&2 echo bt2 pair contigs | "
        #bt2_paired_contigs += self.config_dict["BT2"] + " mem -t " + self.threads_str + " -B 40 -O 60 -E 10 -L 50 "
        # --score-min L,0,-0.2 --mp 40,40 --rdg 10,10 --rfg 10,10 --np 60 --dpad 15 --gbar 4 -L 50 -i S,1,0.75
        # Apr 29 2025: exec order to ditch the scores and go default due to translation conflicts
        bt2_paired_contigs += self.config_dict["BT2"] 
        bt2_paired_contigs += " -x " + contigs_idx + " "
        bt2_paired_contigs += " -1 " + self.file_dict["repop_p1"] + " "
        bt2_paired_contigs += " -2 " + self.file_dict["repop_p2"] + " " 
        bt2_paired_contigs += "| samtools view > " + self.file_dict["contigs_p_sam"]

        bt2_singletons_contigs = ">&2 echo bt2 singleton contigs | "
        bt2_singletons_contigs += self.config_dict["BT2"] 
        #bt2_singletons_contigs += self.config_dict["BT2"] + " mem -t " + self.threads_str + " -B 40 -O 60 -E 10 -L 50 "
        #bt2_singletons_contigs += self.config_dict["BT2"] +  " --score-min L,0,-0.2 --mp 40,40 --rdg 10,10 --rfg 10,10 --np 60 --dpad 15 --gbar 4 -L 50 -i S,1,0.75 "
        bt2_singletons_contigs += " -x " + contigs_idx + " "
        bt2_singletons_contigs += " -U " + self.file_dict["repop_s"]
        bt2_singletons_contigs += " > " + self.file_dict["contigs_s_sam"]
        
        make_contig_map = ">&2 echo Making contig map | " 
        make_contig_map += self.config_dict["Python"] + " "
        make_contig_map += self.config_dict["Map_contig"] + " "
        make_contig_map += self.read_mode + " "
        make_contig_map += self.file_dict["contigs_map"] + " "
        make_contig_map += self.file_dict["repop_p1"] + " "
        make_contig_map += self.file_dict["repop_p2"] + " "
        make_contig_map += self.file_dict["contigs_p1"] + " "
        make_contig_map += self.file_dict["contigs_p2"] + " "
        make_contig_map += self.file_dict["repop_s"] + " "
        make_contig_map += self.file_dict["contigs_s"] + " "
        make_contig_map += self.file_dict["contigs_s_sam"] + " "
        if(self.read_mode == "paired"):
            make_contig_map += self.file_dict["contigs_p_sam"]
            
            
        flush_bad_contigs = ">&2 echo flush bad contigs | " 
        flush_bad_contigs += self.config_dict["Python"] + " "
        flush_bad_contigs += self.config_dict["flush_bad_contigs"] + " "
        flush_bad_contigs += contig_map + " "
        flush_bad_contigs += final_contigs + " "
        flush_bad_contigs += self.file_dict["contigs_out_fa"]
        

        make_marker = "touch " + marker_file

        if self.read_mode == "single":
            COMMANDS_Assemble = [
                spades,
                #spades_rename,
                disassemble_contigs,
                remove_whitespace,
                bt2_index,
                bt2_singletons_contigs,
                make_contig_map,
                flush_bad_contigs, 
                make_marker
                
            ]
        elif self.read_mode == "paired":
            COMMANDS_Assemble = [
                spades,
                #spades_rename,
                disassemble_contigs,
                remove_whitespace,
                bt2_index,
                bt2_paired_contigs,
                bt2_singletons_contigs,
                make_contig_map, 
                flush_bad_contigs, 
                make_marker
            ]

        return COMMANDS_Assemble
    
    def create_GA_pre_scan_taxa_command(self, operating_mode, marker_file):
        #keep this but don't copy the DB files.  refer to the list instead. 
        if(operating_mode == "c"):
            kraken2_c = ">&2 echo Kraken2 on contigs | "
            kraken2_c += self.config_dict["kraken2"] + " "
            kraken2_c += "--db " + self.config_dict["kraken2_db"] + " "
            kraken2_c += "--threads " + str(self.config_dict["num_threads"]) + " "
            kraken2_c += self.file_dict["contigs_out_fa"] + " "
            kraken2_c += "--output " + self.file_dict["ga_ps_k2_report_c"]
            
            make_marker = "touch " + marker_file
            return [kraken2_c + " && " + make_marker]
            
        elif(operating_mode == "s"):
            kraken2_s = ">&2 echo Kraken2 on singletons | "
            kraken2_s += self.config_dict["kraken2"] + " "
            kraken2_s += "--db " + self.config_dict["kraken2_db"] + " "
            kraken2_s += "--threads " + str(self.config_dict["num_threads"]) + " "
            kraken2_s += self.file_dict["contigs_s"] + " " 
            kraken2_s += "--output " + self.file_dict["ga_ps_k2_report_s"]
            
            make_marker = "touch " + marker_file
            
            return [kraken2_s + " && " + make_marker]
            
        elif(operating_mode == "p1"):
            kraken2_p = ">&2 echo Kraken2 on paired | " 
            kraken2_p += self.config_dict["kraken2"] + " "
            kraken2_p += "--db " + self.config_dict["kraken2_db"] +  " "
            kraken2_p += "--threads " + str(self.config_dict["num_threads"]) + " "
            kraken2_p += "--paired " + self.file_dict["contigs_p1"] + " " 
            kraken2_p += self.file_dict["contigs_p2"]+ " "
            kraken2_p += "--output " + self.file_dict["ga_ps_k2_report_p"]
            
            make_marker = "touch " + marker_file
            
            return [kraken2_p + " && " + make_marker]

    def create_GA_pre_scan_command(self, marker_file):

        
        ga_get_lib = ">&2 echo GA pre-scan get libs | "
        ga_get_lib += self.config_dict["Python"] + " "
        ga_get_lib += self.config_dict["GA_pre_scan_get_lib"] + " "
        ga_get_lib += self.file_dict["ga_ps_k2_report_all"] + " "
        ga_get_lib += self.config_dict["taxid_tree"] + " "
        ga_get_lib += self.config_dict["nodes"] + " "
        ga_get_lib += self.file_dict["ga_lib_list"]+ " "
        ga_get_lib += self.file_dict["ga_lib_reject"] + " "
        ga_get_lib += self.config_dict["source_taxa_DB"] + " "
        ga_get_lib += str(self.config_dict["taxa_exist_cutoff"])
        
        make_marker = "touch " + marker_file
        
        return [ga_get_lib + " && " + make_marker]
        
    def create_GA_BT2_command(self, db_path, read_1, read_2, sam_out, marker_file, op_mode):
        # meant to be called multiple times: query file is a split file
        # aug 10, 2021: changed ref path to accomodate new split-chocophlan
        #feb 20, 2025: reiterating call-per-file.
           
        bt2_job = self.config_dict["BT2"] + " -p " + self.threads_str + " "
        bt2_job += " -x " + db_path + " "
        if(op_mode == "p"):
            bt2_job += "-1 " + read_1 + " "
            bt2_job += "-2" + read_2 + " "
        
        elif(read_1.endswith(".fasta")):
            bt2_job += "-f -U " + read_1 + " "
        else:
            bt2_job += "-U " + read_1 + " "
    
        bt2_job += " | samtools view >> " + sam_out

        #make_marker = ">&2 echo marking bt2 job complete: " + file_tag + " | "
        make_marker = "touch" + " " + marker_file
    
        COMMANDS_bt2 = [
            bt2_job + " && " + make_marker
        ]
    
        return COMMANDS_bt2
           
    def create_GA_BT2_pp_command(self, ref_path, gene_map, genes_hit, prot_hit, read_in_1, read_in_2, sam_file, read_out_1, read_out_2, op_mode, marker_file):
        #may 01, 2025: simplified pp call.  
        #run on every lib list section.
        #

        ga_bt2_pp = self.config_dict["Python"] + " "
        ga_bt2_pp += self.config_dict["GA_BT2_pp"] + " "
        ga_bt2_pp += str(self.config_dict["BT2_cigar_cutoff"]) + " "
        ga_bt2_pp += ref_path + " "
        if(not os.path.exists(self.file_dict["contigs_out_fa"])):
            ga_bt2_pp += "None" + " "
        else:        
            ga_bt2_pp += self.file_dict["contigs_map"] + " "  # IN
        ga_bt2_pp += gene_map + " "  # OUT
        ga_bt2_pp += genes_hit + " " #OUT
        ga_bt2_pp += prot_hit + " "
        ga_bt2_pp += read_in_1 + " "
        ga_bt2_pp += read_in_2 + " "

        ga_bt2_pp += sam_file + " "
        ga_bt2_pp += read_out_1 + " "
        ga_bt2_pp += read_out_2 + " "
        ga_bt2_pp += op_mode


        make_marker = "touch " + marker_file


        COMMANDS_bt2_pp = [
            ga_bt2_pp + " && " + make_marker
        ]

        return COMMANDS_bt2_pp

    def create_GA_DMD_command(self, query_file, out_file, block_size, marker_file):
        

        
        diamond_annotate = ">&2 echo " + str(dt.today()) + " GA DIAMOND " + query_file + " | "
        diamond_annotate += self.config_dict["DMD"]
        diamond_annotate += " blastx -p " + self.threads_str
        diamond_annotate += " -d " + self.config_dict["Prot_DB"]
        diamond_annotate += " -q " + query_file 
        diamond_annotate += " -o " + out_file
        diamond_annotate += " --" + str(self.config_dict["DMD_speed"])
        diamond_annotate += " --block-size " + str(block_size)
        diamond_annotate += " -f 6 --tmpdir " + self.dir_dict["GA_DMD_temp"] #section_temp_folder
        diamond_annotate += " -k " + str(self.config_dict["DMD_hit_count"]) 
        diamond_annotate += " --id " + str(self.config_dict["DMD_id_score"]) 
        diamond_annotate += " --query-cover " + str(self.config_dict["DMD_query_cover"])
        diamond_annotate += " --min-score " + str(self.config_dict["DMD_min_score"])
        diamond_annotate += " --unal 1"

        #make_marker = ">&2 echo marking DIAMOND complete: " + sample_root_name + " | "
        make_marker = "touch" + " "
        make_marker += marker_file

        return [diamond_annotate + " && " + make_marker]

    def create_DIAMOND_pp_command_v2(self, reads_in, dmd_in, reads_out, prot_out, prot_map_out, marker_file):
    
        
        diamond_pp = ">&2 echo " + str(dt.today()) + " DMD post process"  + " | "
        diamond_pp += self.config_dict["Python"] + " "
        diamond_pp += self.config_dict["GA_dmd_pp"] + " "
        diamond_pp += str(self.config_dict["DMD_identity_cutoff"]) + " "
        diamond_pp += str(self.config_dict["DMD_length_cutoff"]) + " "
        diamond_pp += str(self.config_dict["DMD_score_cutoff"]) + " "
        diamond_pp += self.config_dict["Prot_DB_reads"] + " "                # IN
        if(not os.path.exists(self.file_dict["contigs_map"]) == "None"):
            diamond_pp += "None" + " "
        else:
            diamond_pp += self.file_dict["contigs_map"] + " "         # IN
        diamond_pp += prot_map_out + " "      # OUT
        diamond_pp += prot_out + " "      # OUT
        
        diamond_pp += reads_in + " "                                                  # IN
        diamond_pp += dmd_in + " "  # IN
        diamond_pp += reads_out + " "     # OUT
        
        make_marker = "touch" + " " + marker_file

        COMMANDS_Annotate_Diamond_Post = [
            diamond_pp + " && " + make_marker
        ]

        return COMMANDS_Annotate_Diamond_Post 

    def create_GA_final_merge_command(self, marker_0, marker_1, marker_2, marker_3):
        
        final_merge_bt2 = self.config_dict["Python"] + " "
        final_merge_bt2 += self.config_dict["GA_final_merge_fasta"] + " "
        final_merge_bt2 += self.dir_dict["GA_BT2_export_genes"] + " "
        final_merge_bt2 += self.file_dict["ga_fm_bt2_genes"]

        final_merge_dmd = self.config_dict["Python"] + " "
        final_merge_dmd += self.config_dict["GA_final_merge_fasta"] + " "
        final_merge_dmd += self.dir_dict["GA_DMD_export_prot"] + " "
        final_merge_dmd += self.file_dict["ga_fm_dmd_prot"]

        final_merge_proteins = self.config_dict["Python"] + " "
        final_merge_proteins += self.config_dict["GA_final_merge_proteins"] + " "
        final_merge_proteins += self.dir_dict["GA_BT2_export_prot"] + " "
        final_merge_proteins += self.dir_dict["GA_DMD_export_prot"] + " "
        final_merge_proteins += self.file_dict["ga_fm_all_prot"]

        final_merge_maps = self.config_dict["Python"] + " "
        final_merge_maps += self.config_dict["GA_final_merge_maps"] + " "
        final_merge_maps += self.dir_dict["GA_BT2_export_maps"] + " "
        final_merge_maps += self.dir_dict["GA_DMD_export_maps"] + " "
        final_merge_maps += self.file_dict["contigs_map"] + " "
        final_merge_maps += self.file_dict["ga_fm_gene_map"]
        
        
        
        
        COMMANDS_ga_final_merge = [
            final_merge_bt2 + " && touch " + marker_0,
            final_merge_dmd + " && touch " + marker_1,
            final_merge_proteins +  " && touch " + marker_2,
            final_merge_maps + " && touch " + marker_3
        
        ]
    
        return COMMANDS_ga_final_merge
    

        subfolder               = os.path.join(self.output_path, current_stage_name)
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
        
    def create_TA_repack_command(self, marker):
        repack = self.config_dict["Python"] + " "
        repack += self.config_dict["TA_apply_names"] + " "
        repack += self.file_dict["ga_ps_k2_report_all"] + " "
        repack += self.config_dict["names"] + " "
        repack += self.config_dict["nodes"] + " "
        repack += self.file_dict["ta_report"]

        make_marker = "touch " + marker 

        COMMANDS_ta_repack = [
            repack + " && " + make_marker
        ]


        return COMMANDS_ta_repack

    def create_deepec_command(self, marker):
        deepec_run = self.config_dict["DeepEC"] + " "
        deepec_run += "-i" + " " + self.file_dict["ga_fm_all_prot"] + " "
        deepec_run += "-o" + " " + self.dir_dict["EC_data"] + " "
        deepec_run += "-b" + " " + str(self.config_dict["deepec_batch_size"]) + " "
        deepec_run += "-cpu" + " " + str(self.config_dict["deepec_cpu_count"])

        make_marker = "touch" + " " + marker
        return [deepec_run + " && " + make_marker]
        
    def create_rpkm_command(self, marker):
       
        move_map = ">&2 echo cp gene map | "
        move_map += "cp" + " " 
        move_map += self.file_dict["ga_fm_gene_map"] + " "
        move_map += self.file_dict["out_gene_map"]

        network_generation = ">&2 echo Generating RPKM and Cytoscape network | "
        network_generation += self.config_dict["Python"] + " "
        network_generation += self.config_dict["RPKM"] + " "
        network_generation += str(self.config_dict["RPKM_cutoff"]) + " "
        network_generation += "None" + " "
        network_generation += self.config_dict["nodes"] + " "
        network_generation += self.config_dict["names"] + " "
        network_generation += self.file_dict["out_gene_map"] + " "
        network_generation += self.file_dict["ga_ps_k2_report_all"]+ " "
        network_generation += self.file_dict["ec_final_report"] + " "
        network_generation += self.config_dict["show_unclassified"] + " "
        network_generation += self.file_dict["out_rpkm"] + " "
        network_generation += self.file_dict["out_cytoscape"]
        
        
        
        flatten_rpkm = ">&2 echo Reformat RPKM for EC heatmap | "
        flatten_rpkm += self.config_dict["Python"] + " "
        flatten_rpkm += self.config_dict["format_RPKM"] + " "
        flatten_rpkm += self.file_dict["out_rpkm"] + " "
        flatten_rpkm += self.file_dict["out_heatmap_rpkm"]

        make_marker = "touch " + marker
        
        return [move_map, network_generation, flatten_rpkm + " && " + make_marker]
        
    def create_output_unique_junk(self):
        command_list = list()
        for item in self.config_dict["Host_IDs"]:
            #no longer needed. we save the host reads.
            #get_unique_host_reads_singletons = ">&2 echo get singleton host reads for stats | "
            #get_unique_host_reads_singletons += self.config_dict["Python"] + " "
            #get_unique_host_reads_singletons += self.config_dict["get_unique_host_reads"] + " "
            #get_unique_host_reads_singletons += self.file_dict[item + "_host_s"] + " "
            #get_unique_host_reads_singletons += self.file_dict["qf_u_s"] + " "
            #get_unique_host_reads_singletons += os.path.join(unique_hosts_folder, "singletons_hosts.fastq")
        
            if not (os.path.exists(self.marker_dict["Out_" + item + "_repop_host_s"])):
                repop_singletons_hosts = ">&2 echo repopulating singletons hosts | " 
                repop_singletons_hosts += self.config_dict["Python"] + " "
                repop_singletons_hosts += self.config_dict["duplicate_repopulate"]+ " "
                if(self.read_mode == "single"):
                    repop_singletons_hosts += self.file_dict["qf_hq_s"] + " " #os.path.join(quality_folder, "singletons_hq.fastq") + " "
                else:
                    repop_singletons_hosts += self.file_dict["qf_o_s"] + " " #os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
                repop_singletons_hosts += self.file_dict[item + "_host_s"] + " " #os.path.join(unique_hosts_folder, "singletons_hosts.fastq") + " "
                repop_singletons_hosts += self.file_dict["qf_clstr_s"] + " " #os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
                repop_singletons_hosts += self.file_dict[item + "_full_host_s"] + " "#os.path.join(full_hosts_folder, "singletons_full_hosts.fastq")
                repop_singletons_hosts += "&& touch " + self.marker_dict["Out_" + item + "_repop_host_s"]

                command_list.append(repop_singletons_hosts)

            if not (os.path.exists(self.marker_dict["Out_" + item + "_repop_host_p1"])):
                repop_pair_1_hosts = ">&2 echo repopulating pair 1 hosts | " 
                repop_pair_1_hosts += self.config_dict["Python"] + " "
                repop_pair_1_hosts += self.config_dict["duplicate_repopulate"]+ " "
                if self.read_mode == "single":
                    repop_pair_1_hosts += self.file_dict["qf_hq_p1"] + " " #os.path.join(quality_folder, "pair_1_match.fastq") + " "
                else:
                    repop_pair_1_hosts += self.file_dict["qf_o_p1"] + " "

                repop_pair_1_hosts += self.file_dict[item + "_host_p1"] + " "#os.path.join(unique_hosts_folder, "pair_1_hosts.fastq") + " "
                repop_pair_1_hosts += self.file_dict["qf_clstr_p1"] + " "#os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " "
                repop_pair_1_hosts += self.file_dict[item + "_full_host_p1"] + " "#os.path.join(full_hosts_folder, "pair_1_full_hosts.fastq")
                repop_pair_1_hosts += "&& touch " + self.marker_dict["Out_" + item + "_repop_host_p1"]
                command_list.append(repop_pair_1_hosts)

            if not (os.path.exists(self.marker_dict["Out_" + item + "_repop_host_p2"])):
                repop_pair_2_hosts = ">&2 echo repopulating pair 2 hosts | " 
                repop_pair_2_hosts += self.config_dict["Python"] + " "
                repop_pair_2_hosts += self.config_dict["duplicate_repopulate"]+ " "
                if self.read_mode == "single":
                    repop_pair_2_hosts += self.file_dict["qf_hq_p2"] + " " #os.path.join(quality_folder, "pair_1_match.fastq") + " "
                else:
                    repop_pair_2_hosts += self.file_dict["qf_o_p2"] + " "

                repop_pair_2_hosts += self.file_dict[item + "_host_p2"] + " "#os.path.join(unique_hosts_folder, "pair_1_hosts.fastq") + " "
                repop_pair_2_hosts += self.file_dict["qf_clstr_p1"] + " "#os.path.join(quality_folder, "pair_1_unique.fastq.clstr") + " "
                repop_pair_2_hosts += self.file_dict[item + "_full_host_p2"] + " "#os.path.join(full_hosts_folder, "pair_1_full_hosts.fastq")
                repop_pair_2_hosts += "&& touch " + self.marker_dict["Out_" + item + "_repop_host_p2"]
                command_list.append(repop_pair_2_hosts)

        if not (os.path.exists(self.marker_dict["Out_vec_repop_s"])):
            repop_vec_s = ">&2 echo repopulating singletons vectors | " 
            repop_vec_s += self.config_dict["Python"] + " "
            repop_vec_s += self.config_dict["duplicate_repopulate"]+ " "
            if(self.read_mode == "single"):
                repop_vec_s += self.file_dict["qf_hq_s"] + " "#os.path.join(quality_folder, "singletons_hq.fastq") + " "
            else:
                repop_vec_s += self.file_dict["qf_o_s"] + " "#os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
            repop_vec_s += self.file_dict["vec_s"] + " "#os.path.join(unique_vectors_folder, "singletons_vectors.fastq") + " "
            repop_vec_s += self.file_dict["qf_clstr_s"] + " "#os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
            repop_vec_s += self.file_dict["out_full_vec_s"] + " " #os.path.join(full_vectors_folder, "singletons_full_vectors.fastq")
            repop_vec_s += "&& touch " + self.marker_dict["Out_vec_repop_s"]
            command_list.append(repop_vec_s)

        if not (os.path.exists(self.marker_dict["Out_vec_repop_p1"])):
            repop_vec_p1 = ">&2 echo repopulating vectors p1 | " 
            repop_vec_p1 += self.config_dict["Python"] + " "
            repop_vec_p1 += self.config_dict["duplicate_repopulate"]+ " "
            if(self.read_mode == "single"):
                repop_vec_p1 += self.file_dict["qf_hq_p1"] + " "#os.path.join(quality_folder, "singletons_hq.fastq") + " "
            else:
                repop_vec_p1 += self.file_dict["qf_o_p1"] + " "#os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
            repop_vec_p1 += self.file_dict["vec_p1"] + " "#os.path.join(unique_vectors_folder, "singletons_vectors.fastq") + " "
            repop_vec_p1 += self.file_dict["qf_clstr_p1"] + " "#os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
            repop_vec_p1 += self.file_dict["out_full_vec_p1"] + " " #os.path.join(full_vectors_folder, "singletons_full_vectors.fastq")
            repop_vec_p1 += "&& touch " + self.marker_dict["Out_vec_repop_p1"]
            command_list.append(repop_vec_p1)

        if not (os.path.exists(self.marker_dict["Out_vec_repop_p2"])):
            repop_vec_p2 = ">&2 echo repopulating vectors p2 | " 
            repop_vec_p2 += self.config_dict["Python"] + " "
            repop_vec_p2 += self.config_dict["duplicate_repopulate"]+ " "
            if(self.read_mode == "single"):
                repop_vec_p2 += self.file_dict["qf_hq_p2"] + " "#os.path.join(quality_folder, "singletons_hq.fastq") + " "
            else:
                repop_vec_p2 += self.file_dict["qf_o_p2"] + " "#os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
            repop_vec_p2 += self.file_dict["vec_p2"] + " "#os.path.join(unique_vectors_folder, "singletons_vectors.fastq") + " "
            repop_vec_p2 += self.file_dict["qf_clstr_p1"] + " "#os.path.join(quality_folder, "singletons_unique.fastq.clstr") + " "
            repop_vec_p2 += self.file_dict["out_full_vec_p2"] + " " #os.path.join(full_vectors_folder, "singletons_full_vectors.fastq")
            repop_vec_p2 += "&& touch " + self.marker_dict["Out_vec_repop_p2"]
            command_list.append(repop_vec_p2)


        return command_list
        

    def create_output_per_read_scores_command(self, marker):
     
        
        per_read_scores = ">&2 echo collecting per-read quality | " 
        per_read_scores += self.config_dict["Python"] + " "
        per_read_scores += self.config_dict["read_quality_metrics"] + " "
        if(self.read_mode == "single"):
            per_read_scores += "single" + " "
            per_read_scores += self.config_dict["single"] + " "
            per_read_scores += self.file_dict["qf_hq_s"] + " "#os.path.join(quality_folder, "singletons_hq.fastq") + " "
            per_read_scores += self.dir_dict["out_export"] #os.path.join(final_folder)
            
        elif(self.read_mode == "paired"):
            per_read_scores += "paired" + " " 
            per_read_scores += self.config_dict["pair_1"] + " "
            per_read_scores += self.config_dict["pair_2"] + " "
            per_read_scores += self.file_dict["qf_o_p1"] + " "#os.path.join(quality_folder, "pair_1_match.fastq") + " "
            per_read_scores += self.file_dict["qf_o_p2"] + " "#os.path.join(quality_folder, "pair_2_match.fastq") + " "
            per_read_scores += self.file_dict["qf_o_s"] + " " #os.path.join(quality_folder, "singletons_with_duplicates.fastq") + " "
            per_read_scores += self.dir_dict["out_export"] #os.path.join(final_folder)
            
        make_marker = "touch " + marker
        return [per_read_scores + " && "+ make_marker]
        
    def create_output_copy_taxa_command(self, marker):
        copy_taxa = ">&2 echo " + str(dt.today()) + " copying taxa data | " 
        copy_taxa += "cp" + " "
        copy_taxa += self.file_dict["ta_report"] + " "
        copy_taxa += self.file_dict["out_taxa_report"]
        #copy_taxa += os.path.join(taxa_folder, "constrain_classification.tsv") + " "
        #copy_taxa += os.path.join(final_folder, "taxa_classifications.tsv")
        
        make_marker = "touch " + marker
        return [copy_taxa + " && " + make_marker]
    
    def create_output_copy_leftovers_command(self, marker):
        copy_rem_s = ">&2 echo " + str(dt.today()) + " copying leftover reads | "
        copy_rem_s += "cp" + " " 
        copy_rem_s += self.file_dict["ga_dmd_rem_s"] + " " + self.file_dict["out_rem_s"]

        copy_rem_c = "cp" + " "
        copy_rem_c += self.file_dict["ga_dmd_rem_c"] + " " + self.file_dict["out_rem_c"]

        copy_rem_p1 = "cp" + " "
        copy_rem_p1 += self.file_dict["ga_dmd_rem_p1"] + " " + self.file_dict["out_rem_p1"]

        copy_rem_p2 = "cp" + " "
        copy_rem_p2 += self.file_dict["ga_dmd_rem_p2"] + " " + self.file_dict["out_rem_p2"]

        return [copy_rem_s, copy_rem_c, copy_rem_p1, copy_rem_p2 + " && touch " + marker]
        
    def create_output_contig_stats_command(self, marker):
        
        contig_stats = ">&2 echo " + str(dt.today()) + " collecting contig stats | " 
        contig_stats += self.config_dict["Python"] + " "
        contig_stats += self.config_dict["contig_stats"] + " "
        contig_stats += self.file_dict["contigs_out_fa"] + " "#os.path.join(contig_folder, "contigs.fasta") + " "
        contig_stats += self.file_dict["out_contig_stats"]#os.path.join(final_folder, "contig_stats.txt")

        make_marker = "touch " + marker
        
        return [contig_stats + " && " + make_marker]
        
    def create_output_EC_heatmap_command(self, marker):
        
        EC_heatmap = ">&2 echo " + str(dt.today()) + " forming EC heatmap | "
        EC_heatmap += self.config_dict["Python"] + " "
        EC_heatmap += self.config_dict["ec_heatmap"] + " "
        EC_heatmap += self.config_dict["EC_pathway"] + " "
        EC_heatmap += self.file_dict["out_heatmap_rpkm"] + " "#os.path.join(final_folder, "EC_heatmap_RPKM.tsv") + " "
        EC_heatmap += self.config_dict["path_to_superpath"] + " "
        EC_heatmap += self.dir_dict["out_export"] #final_folder
        
        make_marker = "touch " + marker
        return [EC_heatmap + " && " + make_marker]
        
        
        
    def create_output_read_count_command(self, marker):
        #new approach. needs to cycle through.

        read_counts = ">&2 echo " + str(dt.today()) + " generating read count table | "
        read_counts += self.config_dict["Python"] + " "
        read_counts += self.config_dict["read_count"] + " "
        if(self.config_dict["read_mode"] == "s"):
            read_counts += self.config_dict["single"] + " "
        else:
            read_counts += self.config_dict["pair_1"] + " "
            
        read_counts += self.file_dict["qf_o_s"] + " "
        read_counts += self.file_dict["qf_o_p1"] + " "
        read_counts += self.file_dict["qf_o_p2"] + " "
        read_counts += self.file_dict["qf_u_s"] + " "
        read_counts += self.file_dict["qf_u_p1"] + " "
        read_counts += self.file_dict["qf_u_p2"] + " "

        read_counts += self.dir_dict["out_hosts"] + " "
        read_counts += self.file_dict["out_full_vec_s"] + " "
        read_counts += self.file_dict["out_full_vec_p1"] + " "
        read_counts += self.file_dict["out_full_vec_p2"] + " "
        read_counts += self.file_dict["repop_other_s"] + " "
        read_counts += self.file_dict["repop_other_p1"] + " "
        read_counts += self.file_dict["repop_other_p2"] + " "
        read_counts += self.file_dict["repop_s"] + " "
        read_counts += self.file_dict["repop_p1"] + " "
        read_counts += self.file_dict["repop_p2"] + " "

        read_counts += self.file_dict["ga_fm_gene_map"] + " "
        read_counts += self.file_dict["ec_final_report"] + " "
        read_counts += self.file_dict["out_read_count"]

        make_marker = "touch " + marker

        return [read_counts + " && " + make_marker]
        
        

    def create_accounting_command(self, stage_name, entries):
        # Builds a subprocess call that appends read-count rows for `stage_name`
        # to the cumulative accounting log (self.file_dict["read_accounting_log"]).
        # entries: list of (role, label, filepath) tuples, role is "input" or "output".
        # No marker is used here -- by design this call re-runs every time its
        # owning stage method executes, including on resumed pipeline runs, so
        # the log can accumulate duplicate rows for a stage across multiple runs.
        accounting_call = ">&2 echo " + str(dt.today()) + " logging read counts for " + stage_name + " | "
        accounting_call += self.config_dict["Python"] + " "
        accounting_call += self.config_dict["read_accounting"] + " "
        accounting_call += stage_name + " "
        accounting_call += self.file_dict["read_accounting_log"]
        for role, label, filepath in entries:
            accounting_call += " " + role + " " + label + " " + filepath

        return [accounting_call]

    def create_output_taxa_groupby_command(self, marker):

        taxa_groupby = ">&2 echo making Taxa summary | " 
        taxa_groupby += self.config_dict["Python"] + " "
        taxa_groupby += self.config_dict["taxa_table"] + " "
        taxa_groupby += self.file_dict["ta_report"] + " "
        taxa_groupby += self.file_dict["out_taxa_groupby"]
        #taxa_groupby += os.path.join(final_folder, "taxa_classifications.tsv") + " "
        #taxa_groupby += os.path.join(final_folder, "taxa_summary.tsv")
        
        make_marker = "touch " + marker

        return [taxa_groupby + " && " + make_marker]