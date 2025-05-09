import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


        

class mpro_marker:
    def issue_split_markers(self, header, count, location):
        #used for creating markers where there are split data.
        #auto-creates the marker to be used.  
        #also writes it to the main registry <marker_dict>
        for i in range(0, count):
            marker_name = header + "_" + str(i)
            self.marker_dict[marker_name] = os.path.join(self.dir_dict[location], marker_name)
            

    def place_marker(self, tag):
        if(not os.path.exist(self.marker_dict[tag])):
            marker = open(self.marker_dict[tag])
            marker.close()
        with open(self.bypass_log, "a") as log:
            log.write(self.m_name[tag] + "\n")

    def check_marker(self, tag):
        print(dt.today(), "checking:", tag)
        if(os.path.exists(tag)):
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
        self.m_name["GA_BWA"] = "GA_BWA_marker"
        self.m_name["GA_DMD"] = "GA_DMD_marker"
        self.m_name["TA"] = "TA_marker"
        self.m_name["rRNA_barrnap_s"] = "rRNA_barrnap_s_marker"
        self.m_name["rRNA_barrnap_p1"] = "rRNA_barrnap_p1_marker"
        self.m_name["rRNA_barrnap_p2"] = "rrNA_barrnap_p2_marker"
        self.m_name["rRNA_inf_s"] = "rRNA_inf_s_marker"
        self.m_name["rRNA_inf_p1"] = "rRNA_inf_p1_marker"
        self.m_name["rRNA_inf_p2"] = "rRNA_inf_p2_marker"

        self.marker_dict = dict()
        self.marker_dict["qf"] = os.path.join(self.dir_dict["qf"], self.m_name["qf"])
        for item in self.config_dict["host_IDs"]:
            self.marker_dict[item + "_host"] = os.path.join(self.dir_dict["main_host"], item + "_mkr")
        self.marker_dict["vec"] = os.path.join(self.dir_dict["vec"], self.m_name["vec"])
        self.marker_dict["rRNA"] = os.path.join(self.dir_dict["rRNA"], self.m_name["rRNA"])
        self.marker_dict["rRNA_split_s"] = os.path.join(self.dir_dict["rRNA"], "split_s_marker")
        self.marker_dict["rRNA_split_p1"] = os.path.join(self.dir_dict["rRNA"], "split_p1_marker")
        self.marker_dict["rrNA_split_p2"] = os.path.join(self.dir_dict["rRNA"], "split_p2_marker")
        self.marker_dict["rRNA_inf_pp_paired"] = os.path.join(self.dir_dict["rRNA"], "inf_pp_paired_marker")
        self.marker_dict["rRNA_inf_pp_s"] = os.path.join(self.dir_dict["rRNA"], "inf_pp_s_marker")
        self.marker_dict["repop"] = os.path.join(self.dir_dict["repop"], self.m_name["repop"])
        self.marker_dict["contigs"] = os.path.join(self.dir_dict["contigs"], self.m_name["contigs"])
        self.marker_dict["ga_ps_s"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_s_marker")
        self.marker_dict["ga_ps_c"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_c_marker")
        self.marker_dict["ga_ps_p"] = os.path.join(self.dir_dict["GA_ps"], "ga_ps_p_marker")
        
        


