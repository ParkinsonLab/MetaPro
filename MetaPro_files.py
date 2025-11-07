
import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


        

class mpro_file_handler:
    #the file-interconnect. 
    

    def read_lib_list(self, lib_file):
        #reads the lib list and fills in the file dict
        #new format: just lines that are found.
        lib_set = set()
        with open(lib_file, "r") as lib_list:
            for line in lib_list:
                lib_line = line.strip("\n")
                if(os.path.exists(lib_line)):
                    lib_set.add(lib_line)
                    print("OK:", lib_line)
        return lib_set


    def clean_files(self, stage):
        key_list = []
        if(stage == "qf"):
            key_list = self.file_dict["qf_clean_list"]


        for item in key_list:
            if(os.path.exists(self.file_dict[item])):
                os.remove(self.file_dict[item])

    def get_file_dict(self):
        return self.file_dict
    
    def create_split_filepaths(self, header, full_count, file_base_path, ext):
        full_count = int(full_count)
        for i in range (0, full_count):
            self.file_dict[header + "_" + str(i)] = file_base_path + "_" + str(i) + ext

    def __init__(self, config_dict, dir_dict):
        print("in file:", config_dict["pair_1"])
        #time.sleep(3)
        self.dir_dict = dir_dict
        self.config_dict = config_dict
        self.file_dict = dict()
        self.file_dict["raw_s"] = self.config_dict["single"]
        self.file_dict["raw_p1"] = self.config_dict["pair_1"]
        self.file_dict["raw_p2"] = self.config_dict["pair_2"]
        self.file_dict["qf_sort_p1"] = os.path.join(self.dir_dict["qf_sort"], "pair_1_sorted.fastq")
        self.file_dict["qf_sort_p2"] = os.path.join(self.dir_dict["qf_sort"], "pair_2_sorted.fastq") 
        self.file_dict["qf_adapt_s"] = os.path.join(self.dir_dict["qf_adapt"], "s_no_adapt.fastq")
        self.file_dict["qf_adapt_p1"] = os.path.join(self.dir_dict["qf_adapt"], "p1_no_adapt.fastq")
        self.file_dict["qf_adapt_p2"] = os.path.join(self.dir_dict["qf_adapt"], "p2_no_adapt.fastq")
        self.file_dict["qf_tags_p1"] = os.path.join(self.dir_dict["qf_tags"], "p1_no_tags.fastq")
        self.file_dict["qf_tags_p2"] = os.path.join(self.dir_dict["qf_tags"], "p2_no_tags.fastq")
        self.file_dict["qf_tags_s"] = os.path.join(self.dir_dict["qf_tags"], "s_no_tags.fastq")
        self.file_dict["qf_merge_p1"] = os.path.join(self.dir_dict["qf_merge"], "p1_no_merge.fastq")
        self.file_dict["qf_merge_p2"] = os.path.join(self.dir_dict["qf_merge"], "p2_no_merge.fastq")
        self.file_dict["qf_merge_s"] = os.path.join(self.dir_dict["qf_merge"], "s_merge.fastq")
        self.file_dict["qf_merge_s2"] = os.path.join(self.dir_dict["qf_merge"], "all_s.fastq")
        self.file_dict["qf_hq_s"] = os.path.join(self.dir_dict["qf_export"], "s_hq.fastq")
        self.file_dict["qf_hq_p1"] = os.path.join(self.dir_dict["qf_export"], "p1_hq.fastq")
        self.file_dict["qf_hq_p2"] = os.path.join(self.dir_dict["qf_export"], "p2_hq.fastq")
        self.file_dict["qf_o_p1"] = os.path.join(self.dir_dict["qf_export"], "p1_match.fastq")
        self.file_dict["qf_o_p2"] = os.path.join(self.dir_dict["qf_export"], "p2_match.fastq")
        self.file_dict["qf_o_s"] = os.path.join(self.dir_dict["qf_export"], "s_dupes.fastq")
        self.file_dict["qf_u_p1"] = os.path.join(self.dir_dict["qf_export"], "p1_unique.fastq")
        self.file_dict["qf_u_p2"] = os.path.join(self.dir_dict["qf_export"], "p2_unique.fastq")
        self.file_dict["qf_u_s"] = os.path.join(self.dir_dict["qf_export"], "s_unique.fastq")
        self.file_dict["qf_clstr_s"] = os.path.join(self.dir_dict["qf_export"], "s_unique.fastq.clstr")
        self.file_dict["qf_clstr_p1"] = os.path.join(self.dir_dict["qf_export"], "p1_unique.fastq.clstr")
        self.file_dict["qf_clstr_p2"] = os.path.join(self.dir_dict["qf_export"], "p2_unique.fastq.clstr")
        
        self.file_dict["qf_clean_list"] = [
            "qf_sort_p1", "qf_sort_p2", 
            "qf_adapt_s", "qf_adapt_p1", "qf_adapt_p2",
            "qf_tags_s", "qf_tags_p1", "qf_tags_p2",
            "qf_merge_s", "qf_merge_s2", "qf_merge_p1", "qf_merge_p2"                                           
        ]
        self.file_dict["qf_job"] = os.path.join(self.dir_dict["qf"], "qf_job.sh")
        for item in self.config_dict["Host_IDs"]:
            if(item == "none"):
                break
            else:
                self.file_dict[item + "_host_p1"] = os.path.join(self.dir_dict[item + "_host_export"], "p1_host.fastq")
                self.file_dict[item + "_host_p2"] = os.path.join(self.dir_dict[item + "_host_export"], "p2_host.fastq")
                self.file_dict[item + "_host_s"] = os.path.join(self.dir_dict[item + "_host_export"], "s_host.fastq")
                self.file_dict[item + "_host_o"] = os.path.join(self.dir_dict[item + "_host_export"], "o_host.fastq")
                self.file_dict[item + "_host_u"] = os.path.join(self.dir_dict[item + "_host_export"], "u_host.fastq")
                self.file_dict[item + "_no_host_p1"] = os.path.join(self.dir_dict[item + "_host_export"], "p1_no_host.fastq")
                self.file_dict[item + "_no_host_p2"] = os.path.join(self.dir_dict[item + "_host_export"], "p2_no_host.fastq")
                self.file_dict[item + "_no_host_s"] = os.path.join(self.dir_dict[item + "_host_export"], "s_no_host.fastq")
                self.file_dict[item + "_no_host_o"] = os.path.join(self.dir_dict[item + "_host_export"], "o_no_host.fastq")
                self.file_dict[item + "_no_host_u"] = os.path.join(self.dir_dict[item + "_host_export"], "u_no_host.fastq")
                self.file_dict[item + "_no_host_s_sam"] = os.path.join(self.dir_dict[item + "_host_scan"], "s_no_host.sam")
                self.file_dict[item + "_no_host_p_sam"] = os.path.join(self.dir_dict[item + "_host_scan"], "p_no_host.sam")
                self.file_dict[item + "_job"] = os.path.join(self.dir_dict["host_jobs"])
                self.file_dict["final_no_host_p1"] = os.path.join(self.dir_dict["main_host_export"], "p1_no_host.fastq")
                self.file_dict["final_no_host_p2"] = os.path.join(self.dir_dict["main_host_export"], "p2_no_host.fastq")
                self.file_dict["final_no_host_s"] = os.path.join(self.dir_dict["main_host_export"], "s_no_host.fastq")

                
                

            self.file_dict[item + "_host_clean_list"] = [
                item + "_host_p1", item + "_host_p2", item + "_host_s", item + "_host_o", item + "_host_u",
                item + "_no_host_p1", item + "_no_host_p2", item + "_no_host_s", item + "_no_host_o", item + "_no_host_u",
                item + "_no_host_s_sam", item + "_no_host_p_sam" 
            ]
        
        self.file_dict["vec_s_sam"] = os.path.join(self.dir_dict["vec_scan"], "s_no_vec.sam")
        self.file_dict["vec_s_bam"] = os.path.join(self.dir_dict["vec_scan"], "s_no_vec.bam")
        self.file_dict["no_vec_s"] = os.path.join(self.dir_dict["vec_export"], "s_no_vec.fastq")
        self.file_dict["vec_s"] = os.path.join(self.dir_dict["vec_export"], "s_vec.fastq") 
        self.file_dict["vec_p_sam"] = os.path.join(self.dir_dict["vec_scan"], "p_no_vec.sam")
        self.file_dict["vec_p1"] = os.path.join(self.dir_dict["vec_export"], "p1_vec.fastq")
        self.file_dict["vec_p2"] = os.path.join(self.dir_dict["vec_export"], "p2_vec.fastq")
        self.file_dict["no_vec_p1"] = os.path.join(self.dir_dict["vec_export"], "p1_no_vec.fastq")
        self.file_dict["no_vec_p2"] = os.path.join(self.dir_dict["vec_export"], "p2_no_vec.fastq")

        self.file_dict["vec_clean_list"] = [
            "vec_s_sam", "vec_s_bam", "vec_p_sam"
        ]

        self.file_dict["vec_job"] = os.path.join(self.dir_dict["vec"], "vector_job.sh")    
        
        self.file_dict["rRNA_s_fa"] = os.path.join(self.dir_dict["rRNA_convert"], "s.fasta")
        self.file_dict["rRNA_p1_fa"] = os.path.join(self.dir_dict["rRNA_convert"], "p1.fasta")
        self.file_dict["rRNA_p2_fa"] = os.path.join(self.dir_dict["rRNA_convert"], "p2.fasta")
        
        
        self.file_dict["rRNA_bnap_s"] = os.path.join(self.dir_dict["rRNA_bnap"], "s")
        self.file_dict["rRNA_bnap_p1"] = os.path.join(self.dir_dict["rRNA_bnap"], "p1")
        self.file_dict["rRNA_bnap_p2"] = os.path.join(self.dir_dict["rRNA_bnap"], "p2")

        self.file_dict["rRNA_bnap_job_s"] = os.path.join(self.dir_dict["rRNA_bnap"], "s_bnap.sh")
        self.file_dict["rRNA_bnap_job_p1"] = os.path.join(self.dir_dict["rRNA_bnap"], "p1_bnap.sh")
        self.file_dict["rRNA_bnap_job_p2"] = os.path.join(self.dir_dict["rRNA_bnap"], "p2_bnap.sh")

        self.file_dict["rRNA_bnap_other_s"] = os.path.join(self.dir_dict["rRNA_other"], "other_bnap_s.fastq")
        self.file_dict["rRNA_bnap_other_p1"] = os.path.join(self.dir_dict["rRNA_other"], "other_bnap_p1.fastq")
        self.file_dict["rRNA_bnap_other_p2"] = os.path.join(self.dir_dict["rRNA_other"], "other_bnap_p2.fastq")

        self.file_dict["rRNA_bnap_mRNA_s"] = os.path.join(self.dir_dict["rRNA_mRNA"], "mRNA_bnap_s.fastq")
        self.file_dict["rRNA_bnap_mRNA_p1"] = os.path.join(self.dir_dict["rRNA_mRNA"], "mRNA_bnap_p1.fastq")
        self.file_dict["rRNA_bnap_mRNA_p2"] = os.path.join(self.dir_dict["rRNA_mRNA"], "mRNA_bnap_p2.fastq")

        self.file_dict["rRNA_inf_job_s"] = os.path.join(self.dir_dict["rRNA_jobs"], "inf_s")
        self.file_dict["rRNA_inf_job_p1"] = os.path.join(self.dir_dict["rRNA_jobs"], "inf_p1")
        self.file_dict["rRNA_inf_job_p2"] = os.path.join(self.dir_dict["rRNA_jobs"], "inf_p2")
        
        self.file_dict["rRNA_pp_p_job"] = os.path.join(self.dir_dict["rRNA_jobs"], "rRNA_pp_p")
        self.file_dict["rRNA_pp_s_job"] = os.path.join(self.dir_dict["rRNA_jobs"], "rRNA_pp_s")

        self.file_dict["rRNA_inf_in_s"] = os.path.join(self.dir_dict["rRNA_split"], "inf_s")
        self.file_dict["rRNA_inf_in_p1"] = os.path.join(self.dir_dict["rRNA_split"], "inf_p1")
        self.file_dict["rRNA_inf_in_p2"] = os.path.join(self.dir_dict["rRNA_split"], "inf_p2")

        self.file_dict["rRNA_inf_out_s"] = os.path.join(self.dir_dict["rRNA_inf"], "s")
        self.file_dict["rRNA_inf_out_p1"] = os.path.join(self.dir_dict["rRNA_inf"], "p1")
        self.file_dict["rRNA_inf_out_p2"] = os.path.join(self.dir_dict["rRNA_inf"], "p2")

        self.file_dict["rRNA_bnap_all_s"] = os.path.join(self.dir_dict["rRNA_bnap"], "all_s.bnap_out")
        self.file_dict["rRNA_bnap_all_p1"] = os.path.join(self.dir_dict["rRNA_bnap"], "all_p1.bnap_out")
        self.file_dict["rRNA_bnap_all_p2"] = os.path.join(self.dir_dict["rRNA_bnap"], "all_p2.bnap_out")
        
        self.file_dict["rRNA_inf_all_s"] = os.path.join(self.dir_dict["rRNA_inf"], "all_s.inf_out")
        self.file_dict["rRNA_inf_all_p1"] = os.path.join(self.dir_dict["rRNA_inf"], "all_p1.inf_out")
        self.file_dict["rRNA_inf_all_p2"] = os.path.join(self.dir_dict["rRNA_inf"], "all_p2.inf_out")

        self.file_dict["rRNA_mRNA_s_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all_s_mRNA.fastq")
        self.file_dict["rRNA_mRNA_p1_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all_p1_mRNA.fastq")
        self.file_dict["rRNA_mRNA_p2_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all_p2_mRNA.fastq")
        
        self.file_dict["rRNA_other_s_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all_s_other.fastq")
        self.file_dict["rRNA_other_p1_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all_p1_other.fastq")
        self.file_dict["rRNA_other_p2_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all_p2_other.fastq")

        self.file_dict["repop_job"] = os.path.join(self.dir_dict["repop"], "repop_job.sh")
        self.file_dict["repop_s"] = os.path.join(self.dir_dict["repop_export"], "s.fastq")
        self.file_dict["repop_other_s"] = os.path.join(self.dir_dict["repop_export"], "s_other.fastq")
        self.file_dict["repop_s_extra"] = os.path.join(self.dir_dict["repop_export"], "s_extra.fastq")
        self.file_dict["repop_other_s_extra"] = os.path.join(self.dir_dict["repop_export"], "s_other_extra.fastq")

        self.file_dict["repop_p1"] = os.path.join(self.dir_dict["repop_export"], "p1.fastq")
        self.file_dict["repop_other_p1"] = os.path.join(self.dir_dict["repop_export"], "p1_other.fastq")

        self.file_dict["repop_p2"] = os.path.join(self.dir_dict["repop_export"], "p2.fastq")
        self.file_dict["repop_other_p2"] = os.path.join(self.dir_dict["repop_export"], "p2_other.fastq")

        self.file_dict["contigs_job"] = os.path.join(self.dir_dict["contigs"], "contigs.sh")
        self.file_dict["contigs_transcripts"] = os.path.join(self.dir_dict["contigs_spades"], "transcripts.fasta")
        self.file_dict["contigs_og_fa"] = os.path.join(self.dir_dict["contigs_spades"], "contigs.fasta")
        self.file_dict["contigs_split"] = os.path.join(self.dir_dict["contigs_mgm"], "disassembled_contigs.fasta")
        self.file_dict["contigs_gene_report"] = os.path.join(self.dir_dict["contigs_mgm"], "gene_report.txt")
        self.file_dict["contigs_spades_done"] = os.path.join(self.dir_dict["contigs_spades"], "spades_done")
        
        self.file_dict["contigs_out_fa"] = os.path.join(self.dir_dict["contigs_export"], "contigs.fasta")
        self.file_dict["contigs_map"] = os.path.join(self.dir_dict["contigs_export"], "contigs_map.tsv")
        self.file_dict["contigs_idx"] = os.path.join(self.dir_dict["contigs_export"], "contigs_bt2_idx")

        self.file_dict["contigs_p_sam"] = os.path.join(self.dir_dict["contigs_BT2"], "contigs_p.sam")
        self.file_dict["contigs_s_sam"] = os.path.join(self.dir_dict["contigs_BT2"], "contigs_s.sam")
        self.file_dict["contigs_p1"] = os.path.join(self.dir_dict["contigs_export"], "p1.fastq")
        self.file_dict["contigs_p2"] = os.path.join(self.dir_dict["contigs_export"], "p2.fastq")
        self.file_dict["contigs_s"] = os.path.join(self.dir_dict["contigs_export"], "s.fastq")
        
        self.file_dict["ga_ps_k2_c_job"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_k2_c.sh")
        self.file_dict["ga_ps_k2_s_job"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_k2_s.sh")
        self.file_dict["ga_ps_k2_p1_job"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_k2_p1.sh")
        self.file_dict["ga_ps_k2_p2_job"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_k2+p2.sh")
        
        self.file_dict["ga_ps_make_job"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_make.sh")
        self.file_dict["ga_lib_status"] = os.path.join(self.dir_dict["GA_ps_export"], "lib_status.txt")
        self.file_dict["ga_lib_reject"] = os.path.join(self.dir_dict["GA_ps_export"], "lib_reject.txt")
        self.file_dict["ga_lib_list"] = os.path.join(self.dir_dict["GA_ps_data"], "lib_list.txt")
        self.file_dict["ga_ps_k2_report_c"] = os.path.join(self.dir_dict["GA_ps_k2"], "kraken2_report_c.txt")
        self.file_dict["ga_ps_k2_report_s"] = os.path.join(self.dir_dict["GA_ps_k2"], "kraken2_report_s.txt")
        self.file_dict["ga_ps_k2_report_p"] = os.path.join(self.dir_dict["GA_ps_k2"], "kraken2_report_p.txt")
        self.file_dict["ga_ps_k2_report_all"] = os.path.join(self.dir_dict["GA_ps_k2"], "kraken2_report_all.txt")
        #self.file_dict["ga_ps_wevote"] = os.path.join(self.dir_dict["GA_ps_wevote"], "taxa_report.tsv")

        self.file_dict["ga_split_s"] = os.path.join(self.dir_dict["GA_split"], "s")
        self.file_dict["ga_split_p1"] = os.path.join(self.dir_dict["GA_split"], "p1")
        self.file_dict["ga_split_p2"] = os.path.join(self.dir_dict["GA_split"], "p2")

        self.file_dict["ga_bt2_p_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "p.sam")
        self.file_dict["ga_bt2_s_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "s.sam")
        self.file_dict["ga_bt2_c_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "c.sam")

        self.file_dict["ga_bt2_s_hits_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "s_hits.sam")
        self.file_dict["ga_bt2_s_miss_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "s_miss.sam")
        self.file_dict["ga_bt2_p_hits_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "p_hits.sam")
        self.file_dict["ga_bt2_p_miss_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "p_miss.sam")
        self.file_dict["ga_bt2_c_hits_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "c_hits.sam")
        self.file_dict["ga_bt2_c_miss_sam"] = os.path.join(self.dir_dict["GA_BT2_run"], "c_miss.sam")

        self.file_dict["ga_bt2_gene_map_c"] = os.path.join(self.dir_dict["GA_BT2_export_maps"], "c_gene_map.tsv")
        self.file_dict["ga_bt2_gene_map_p"] = os.path.join(self.dir_dict["GA_BT2_export_maps"], "p_gene_map.tsv")
        self.file_dict["ga_bt2_gene_map_s"] = os.path.join(self.dir_dict["GA_BT2_export_maps"], "s_gene_map.tsv")
        self.file_dict["ga_bt2_genes_c"] = os.path.join(self.dir_dict["GA_BT2_export_genes"], "genes_c.fasta")
        self.file_dict["ga_bt2_genes_p"] = os.path.join(self.dir_dict["GA_BT2_export_genes"], "genes_p.fasta")
        self.file_dict["ga_bt2_genes_s"] = os.path.join(self.dir_dict["GA_BT2_export_genes"], "genes_s.fasta")
        self.file_dict["ga_bt2_prot_s"] = os.path.join(self.dir_dict["GA_BT2_export_prot"], "prot_s.faa")
        self.file_dict["ga_bt2_prot_c"] = os.path.join(self.dir_dict["GA_BT2_export_prot"], "prot_c.faa")
        self.file_dict["ga_bt2_prot_p"] = os.path.join(self.dir_dict["GA_BT2_export_prot"], "prot_p.faa")
        self.file_dict["genes_full"] = os.path.join(self.dir_dict["GA_BT2_export_prot"], "genes_annot.fasta")

        self.file_dict["ga_rem_p1"] = os.path.join(self.dir_dict["GA_BT2_export"], "p1.fasta")
        self.file_dict["ga_rem_p2"] = os.path.join(self.dir_dict["GA_BT2_export"], "p2.fasta")
        self.file_dict["ga_rem_c"] = os.path.join(self.dir_dict["GA_BT2_export"], "c.fasta")
        self.file_dict["ga_rem_s"] = os.path.join(self.dir_dict["GA_BT2_export"], "s.fasta")
        
        self.file_dict["ga_dmd_p1_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_p1.sh")
        self.file_dict["ga_dmd_p2_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_p2.sh")
        self.file_dict["ga_dmd_s_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_s.sh")
        self.file_dict["ga_dmd_c_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_c.sh")
        

        self.file_dict["ga_dmd_p1_dmdout"] = os.path.join(self.dir_dict["GA_DMD_data"], "p1.dmdout")
        self.file_dict["ga_dmd_p2_dmdout"] = os.path.join(self.dir_dict["GA_DMD_data"], "p2.dmdout")
        self.file_dict["ga_dmd_s_dmdout"] = os.path.join(self.dir_dict["GA_DMD_data"], "s.dmdout")
        self.file_dict["ga_dmd_c_dmdout"] = os.path.join(self.dir_dict["GA_DMD_data"], "c.dmdout")

        self.file_dict["ga_dmd_prot_map_c"] = os.path.join(self.dir_dict["GA_DMD_export_maps"], "c_prot_map.tsv")
        self.file_dict["ga_dmd_prot_map_p1"] = os.path.join(self.dir_dict["GA_DMD_export_maps"], "p1_prot_map.tsv")
        self.file_dict["ga_dmd_prot_map_p2"] = os.path.join(self.dir_dict["GA_DMD_export_maps"], "p2_prot_map.tsv")
        self.file_dict["ga_dmd_prot_map_s"] = os.path.join(self.dir_dict["GA_DMD_export_maps"], "s_prot_map.tsv")

        #self.file_dict["ga_dmd_prot_map"] = os.path.join(self.dir_dict["GA_DMD_export"], "dmd_prot_map.tsv")
        #self.file_dict["ga_dmd_prot"] = os.path.join(self.dir_dict["GA_DMD_export"], "dmd_prot.faa")
        self.file_dict["ga_dmd_prot_c"] = os.path.join(self.dir_dict["GA_DMD_export_prot"], "prot_c.faa")
        self.file_dict["ga_dmd_prot_p1"] = os.path.join(self.dir_dict["GA_DMD_export_prot"], "prot_p1.faa")
        self.file_dict["ga_dmd_prot_p2"] = os.path.join(self.dir_dict["GA_DMD_export_prot"], "prot_p2.faa")
        self.file_dict["ga_dmd_prot_s"] = os.path.join(self.dir_dict["GA_DMD_export_prot"], "prot_s.faa")


        self.file_dict["ga_dmd_rem_p1"] = os.path.join(self.dir_dict["GA_DMD_export"], "p1_rem.fasta")
        self.file_dict["ga_dmd_rem_p2"] = os.path.join(self.dir_dict["GA_DMD_export"], "p2_rem.fasta")
        self.file_dict["ga_dmd_rem_s"] = os.path.join(self.dir_dict["GA_DMD_export"], "s_rem.fasta")
        self.file_dict["ga_dmd_rem_c"] = os.path.join(self.dir_dict["GA_DMD_export"], "c_rem.fasta")

        self.file_dict["ga_dmd_pp_p1_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_pp_p1.sh")
        self.file_dict["ga_dmd_pp_p2_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_pp_p2.sh")
        self.file_dict["ga_dmd_pp_s_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_pp_s.sh")
        self.file_dict["ga_dmd_pp_c_job"] = os.path.join(self.dir_dict["GA_DMD_jobs"], "dmd_pp_c.sh")

        self.file_dict["ga_fm_job"] = os.path.join(self.dir_dict["GA_FM"], "GA_FM.sh")
        self.file_dict["ga_fm_gene_map"] = os.path.join(self.dir_dict["GA_FM_export"], "final_gene_map.tsv")
        self.file_dict["ga_fm_bt2_genes"] = os.path.join(self.dir_dict["GA_FM_export"], "bt2_genes.fasta")
        self.file_dict["ga_fm_dmd_prot"] = os.path.join(self.dir_dict["GA_FM_export"], "dmd_proteins.faa")
        self.file_dict["ga_fm_all_prot"] = os.path.join(self.dir_dict["GA_FM_export"], "all_proteins.faa")
        

        self.file_dict["ta_report"] = os.path.join(self.dir_dict["TA_export"], "taxa_report.tsv")
        self.file_dict["ta_job"] = os.path.join(self.dir_dict["TA_jobs"], "TA.sh")


        self.file_dict["ec_job"] = os.path.join(self.dir_dict["EC_jobs"], "EC.sh")
        self.file_dict["ec_report"] = os.path.join(self.dir_dict["EC_data"], "tmp", "DL_prediction_results.txt")
        self.file_dict["ec_final_report"] = os.path.join(self.dir_dict["EC_data"], "DeepECv2_result.txt")


        self.file_dict["out_gene_map"] = os.path.join(self.dir_dict["out_export"], "gene_map.tsv")
        self.file_dict["out_rpkm"] = os.path.join(self.dir_dict["out_export"], "RPKM_table.tsv")
        self.file_dict["out_cytoscape"] = os.path.join(self.dir_dict["out_export"], "Cytoscape_network.tsv")
        self.file_dict["out_heatmap_rpkm"] = os.path.join(self.dir_dict["out_export"], "EC_heatmap_RPKM.tsv")
        self.file_dict["out_taxa_report"] = os.path.join(self.dir_dict["out_export"], "taxa_classification.tsv")

        for item in self.config_dict["Host_IDs"]:

            #self.file_dict[item + "_u_host_s"] = os.path.join(self.dir_dict["out_data"], item + "_u_s.fastq")
            #self.file_dict[item + "_u_host_p1"] = os.path.join(self.dir_dict["out_data"], item + "_u_p1.fastq")
            #self.file_dict[item + "_u_host_p2"] = os.path.join(self.dir_dict[item + "_full_host"], item + "_u_p2.fastq")
            self.file_dict[item + "_full_host_s"] = os.path.join(self.dir_dict[item + "_full_host"], item + "_s_full_host.fastq")
            self.file_dict[item + "_full_host_p1"] = os.path.join(self.dir_dict[item + "_full_host"], item +"_p1_full_host.fastq")
            self.file_dict[item + "_full_host_p2"] = os.path.join(self.dir_dict[item + "_full_host"], item + "_p2_full_host.fastq")

        self.file_dict["out_full_vec_s"] = os.path.join(self.dir_dict["out_data"], "full_s_vec.fastq")
        self.file_dict["out_full_vec_p1"] = os.path.join(self.dir_dict["out_data"], "full_p1_vec.fastq")
        self.file_dict["out_full_vec_p2"] = os.path.join(self.dir_dict["out_data"], "full_p2_vec.fastq")
        self.file_dict["out_contig_stats"] = os.path.join(self.dir_dict["out_export"], "contig_stats.txt")
        self.file_dict["out_read_count"] = os.path.join(self.dir_dict["out_export"], "read_counts.tsv")
        self.file_dict["out_taxa_groupby"] = os.path.join(self.dir_dict["out_export"], "taxa_summary.tsv")
        self.file_dict["out_rem_s"] = os.path.join(self.dir_dict["out_export"], "unannotated_s.fasta")
        self.file_dict["out_rem_c"] = os.path.join(self.dir_dict["out_export"], "unannotated_c.fasta")
        self.file_dict["out_rem_p1"] = os.path.join(self.dir_dict["out_export"], "unannotated_p1.fasta")
        self.file_dict["out_rem_p2"] = os.path.join(self.dir_dict["out_export"], "unannotated_p2.fasta")
        
        self.file_dict["out_copy_rem_job"] = os.path.join(self.dir_dict["out_jobs"], "copy_rem_reads.sh")
        self.file_dict["out_rpkm_job"] = os.path.join(self.dir_dict["out_jobs"], "rpkm.sh")
        self.file_dict["out_repop_job"] = os.path.join(self.dir_dict["out_jobs"], "repop.sh")
        self.file_dict["out_per_read_job"] = os.path.join(self.dir_dict["out_jobs"], "per_read.sh")
        self.file_dict["out_taxa_job"] = os.path.join(self.dir_dict["out_jobs"], "taxa.sh")
        self.file_dict["out_contig_stats_job"] = os.path.join(self.dir_dict["out_jobs"], "contig_stats.sh")
        self.file_dict["out_read_count_job"] = os.path.join(self.dir_dict["out_jobs"], "read_count.sh")
        self.file_dict["out_taxa_groupby_job"] = os.path.join(self.dir_dict["out_jobs"], "groupby_taxa.sh")
        self.file_dict["out_ec_heatmap_job"] = os.path.join(self.dir_dict["out_jobs"], "ec_heatmap.sh")