import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


        

class mpro_marker:
    
    def issue_marker(self, tag, marker_path):
    
        self.marker_dict[tag] = marker_path
    
    
    def issue_split_markers(self, header, count, location):
        #used for creating markers where there are split data.
        #auto-creates the marker to be used.  
        #also writes it to the main registry <marker_dict>
        for i in range(0, count):
            marker_name = header + "_" + str(i)
            #print("issuing:", marker_name)
            self.marker_dict[marker_name] = os.path.join(self.dir_dict[location], marker_name)
            

    def place_marker(self, tag):
        if(not os.path.exists(self.marker_dict[tag])):
            with open(self.marker_dict[tag], "w") as marker:
                pass
        
    def check_marker_list(self, marker_list):
        status = True
        for item in marker_list:
            print(dt.today(), "checking:", item)
            if(os.path.exists(item)):
                status = True
                print("OK:", item)
            else:
                print(dt.today(), "failed marker: ", item)
                return False
        print(dt.today(), "everything's fine")
        return status

    def check_marker(self, tag):
        print(dt.today(), "checking:", tag)
        if(os.path.exists(self.marker_dict[tag])):
            print(dt.today(), "skipping: ", self.marker_dict[tag])
            return False
        else:
            return True


    def __init__ (self, config_dict, dir_dict):
        self.dir_dict = dir_dict
        self.config_dict = config_dict
        self.bypass_log = config_dict["bypass_log"]
        self.m_name = dict()
        self.m_name["qf"] = "qf_marker"
        self.m_name["host"] = "host_marker"
        self.m_name["vec"] = "vec_marker"
        self.m_name["rRNA"] = "rRNA_marker"
        self.m_name["repop"] = "repop_marker"
        self.m_name["contigs"] = "contigs_marker"
        self.m_name["GA_pre_scan"] = "GA_pre_scan_marker"
        self.m_name["GA_BT2"] = "GA_BT2_marker"
        self.m_name["GA_DMD"] = "GA_DMD_marker"
        self.m_name["TA"] = "TA_marker"
        self.m_name["rRNA_barrnap_s"] = "rRNA_barrnap_s_marker"
        self.m_name["rRNA_barrnap_p1"] = "rRNA_barrnap_p1_marker"
        self.m_name["rRNA_barrnap_p2"] = "rrNA_barrnap_p2_marker"
        self.m_name["rRNA_inf_s"] = "rRNA_inf_s"
        self.m_name["rRNA_inf_p1"] = "rRNA_inf_p1"
        self.m_name["rRNA_inf_p2"] = "rRNA_inf_p2"

        self.marker_dict = dict()
        self.marker_dict["qf"] = os.path.join(self.dir_dict["qf"], self.m_name["qf"])
        for item in self.config_dict["host_IDs"]:
            self.marker_dict[item + "_host"] = os.path.join(self.dir_dict["main_host"], item + "_mkr")
        self.marker_dict["vec"] = os.path.join(self.dir_dict["vec"], self.m_name["vec"])
        self.marker_dict["rRNA"] = os.path.join(self.dir_dict["rRNA"], self.m_name["rRNA"])
        self.marker_dict["rRNA_bnap_s"] = os.path.join(self.dir_dict["rRNA"], "rRNA_bnap_s")
        self.marker_dict["rRNA_bnap_p1"] = os.path.join(self.dir_dict["rRNA"], "rRNA_bnap_p1")
        self.marker_dict["rRNA_bnap_p2"] = os.path.join(self.dir_dict["rRNA"], "rRNA_bnap_p2")
        
        self.marker_dict["rRNA_split_s"] = os.path.join(self.dir_dict["rRNA"], "split_s_marker")
        self.marker_dict["rRNA_split_p1"] = os.path.join(self.dir_dict["rRNA"], "split_p1_marker")
        self.marker_dict["rRNA_split_p2"] = os.path.join(self.dir_dict["rRNA"], "split_p2_marker")
        self.marker_dict["rRNA_inf_pp_paired"] = os.path.join(self.dir_dict["rRNA"], "inf_pp_paired_marker")
        self.marker_dict["rRNA_inf_pp_s"] = os.path.join(self.dir_dict["rRNA"], "inf_pp_s_marker")
        self.marker_dict["repop"] = os.path.join(self.dir_dict["repop"], self.m_name["repop"])
        self.marker_dict["contigs"] = os.path.join(self.dir_dict["contigs"], self.m_name["contigs"])
        self.marker_dict["GA_ps"] = os.path.join(self.dir_dict["GA_ps"], "GA_ps")
        self.marker_dict["ga_ps_s"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_s_marker")
        self.marker_dict["ga_ps_c"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_c_marker")
        self.marker_dict["ga_ps_p"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_p_marker")
        self.marker_dict["ga_ps_make"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_assemble")
        self.marker_dict["GA_BT2"] = os.path.join(self.dir_dict["GA_BT2"], "GA_BT2")
        self.marker_dict["GA_BT2_PP"] = os.path.join(self.dir_dict["GA_BT2"], "GA_BT2_PP")
        
        self.marker_dict["GA_DMD"] = os.path.join(self.dir_dict["GA_DMD"], "GA_DMD")
        self.marker_dict["GA_DMD_p1"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_p1")
        self.marker_dict["GA_DMD_p2"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_p2")
        self.marker_dict["GA_DMD_s"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_s")
        self.marker_dict["GA_DMD_c"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_c")
        self.marker_dict["GA_DMD_pp"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_pp")
        self.marker_dict["GA_DMD_pp_p1"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_pp_p1")
        self.marker_dict["GA_DMD_pp_p2"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_pp_p2")
        self.marker_dict["GA_DMD_pp_c"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_pp_c")
        self.marker_dict["GA_DMD_pp_s"] = os.path.join(self.dir_dict["GA_DMD_mkrs"], "GA_DMD_pp_s")
        
        

