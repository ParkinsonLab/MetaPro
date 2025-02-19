
import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


        

class mpro_file_handler:
    #the file-interconnect. 
    


    def clean_files(self, stage):
        key_list = []
        if(stage == "qf"):
            key_list = self.file_dict["qf_clean_list"]


        for item in key_list:
            if(os.path.exists(self.file_dict[item])):
                os.remove(self.file_dict[item])

    def get_file_dict(self):
        return self.file_dict

    def __init__(self, config_dict, dir_dict):
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
        
        self.file_dict["qf_clean_list"] = [
            "qf_sort_p1", "qf_sort_p2", 
            "qf_adapt_s", "qf_adapt_p1", "qf_adapt_p2",
            "qf_tags_s", "qf_tags_p1", "qf_tags_p2",
            "qf_merge_s", "qf_merge_s2", "qf_merge_p1", "qf_merge_p2"                                           
        ]


        self.file_dict["no_host_s_sam"] = os.path.join(self.dir_dict["host_scan"], "s_no_host.sam")
        self.file_dict["no_host_s_bam"] = os.path.join(self.dir_dict["host_scan"], "s_no_host.bam")
        self.file_dict["no_host_s"] = os.path.join(self.dir_dict["host_export"], "s_no_host.fastq")
        self.file_dict["host_s"] = os.path.join(self.dir_dict["host_export"], "s_host_only.fastq")
        self.file_dict["no_host_p_sam"] = os.path.join(self.dir_dict["host_scan"], "p_no_host.sam")
        self.file_dict["no_host_p1"] = os.path.join(self.dir_dict["host_export"], "p1_no_host.fastq")
        self.file_dict["no_host_p2"] = os.path.join(self.dir_dict["host_export"], "p2_no_host.fastq")
        self.file_dict["host_p1"] = os.path.join(self.dir_dict["host_export"], "p1_host.fastq")
        self.file_dict["host_p2"] = os.path.join(self.dir_dict["host_export"], "p2_host.fastq")

        self.file_dict["host_clean_list"] = [
            "no_host_s_sam", "no_host_s_bam", "no_host_p_sam"
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

        self.file_dict["rRNA_split_s"] = os.path.join(self.dir_dict["rRNA_split"], "s_split") #note: the splits are guides. not the actual file
        self.file_dict["rRNA_split_p1"] = os.path.join(self.dir_dict["rRNA_split"], "p1_split")
        self.file_dict["rRNA_split_p2"] = os.path.join(self.dir_dict["rRNA_split"], "p2_split")

        

        self.file_dict["rRNA_barrnap_s"] = os.path.join(self.dir_dict["rRNA_barrnap"], "s")
        self.file_dict["rRNA_barrnap_p1"] = os.path.join(self.dir_dict["rRNA_barrnap"], "p1")
        self.file_dict["rRNA_barrnap_p2"] = os.path.join(self.dir_dict["rRNA_barrnap_p2"], "p2")

        self.file_dict["rRNA_barrnap_other_s"] = os.path.join(self.dir_dict["rRNA_other"], "s_barrnap")
        self.file_dict["rRNA_barrnap_other_p1"] = os.path.join(self.dir_dict["rRNA_other"], "p1_barrnap")
        self.file_dict["rRNA_barrnap_other_p2"] = os.path.join(self.dir_dict["rRNA_other"], "p2_barrnap")

        self.file_dict["rRNA_mRNA_barrnap_s"] = os.path.join(self.dir_dict["rRNA_mRNA"], "barrnap", "s")
        self.file_dict["rRNA_mRNA_barrnap_p1"] = os.path.join(self.dir_dict["rRNA_mRNA"], "barrnap", "p1")
        self.file_dict["rRNA_mRNA_barrnap_p2"] = os.path.join(self.dir_dict["rRNA_mRNA"], "barrnap", "p2")

        self.file_dict["rRNA_inf_out_s"] = os.path.join(self.dir_dict["rRNA_inf"], "s")
        self.file_dict["rRNA_inf_out_p1"] = os.path.join(self.dir_dict["rRNA_inf"], "p1")
        self.file_dict["rRNA_inf_out_p2"] = os.path.join(self.dir_dict["rRNA_inf"], "p2")

        self.file_dict["rRNA_barrnap_all_s"] = os.path.join(self.dir_dict["rRNA_barrnap"], "all_s.barrnap_out")
        self.file_dict["rRNA_barrnap_all_p1"] = os.path.join(self.dir_dict["rRNA_barrnap"], "all_p1.barrnap_out")
        self.file_dict["rRNA_barrnap_all_p2"] = os.path.join(self.dir_dict["rRNA_barrnap"], "all_p2.barrnap_out")
        
        self.file_dict["rRNA_inf_all_s"] = os.path.join(self.dir_dict["rRNA_inf"], "all_s.inf_out")
        self.file_dict["rRNA_inf_all_p1"] = os.path.join(self.dir_dict["rRNA_inf"], "all_p1.inf_out")
        self.file_dict["rRNA_inf_all_p2"] = os.path.join(self.dir_dict["rRNA_inf"], "all_p2.inf_out")

        self.file_dict["rRNA_mRNA_s_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all", "s_mRNA.fastq")
        self.file_dict["rRNA_mRNA_p1_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all", "p1_mRNA.fastq")
        self.file_dict["rRNA_mRNA_p2_fq"] = os.path.join(self.dir_dict["rRNA_mRNA"], "all", "p2_mRNA.fastq")
        
        self.file_dict["rRNA_other_s_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all", "s_other.fastq")
        self.file_dict["rRNA_other_p1_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all", "p1_other.fastq")
        self.file_dict["rRNA_other_p2_fq"] = os.path.join(self.dir_dict["rRNA_other"], "all", "p2_other.fastq")

        self.file_dict["repop_clstr_s"] = os.path.join(self.dir_dict["repop_data"], "s.clstr")
        self.file_dict["repop_s"] = os.path.join(self.dir_dict["repop_export"], "s.fastq")
        self.file_dict["repop_other_s"] = os.path.join(self.dir_dict["repop_export"], "s_other.fastq")

        