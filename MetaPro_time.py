import os
import sys
from datetime import datetime as dt
import math
import time
from configparser import ConfigParser, ExtendedInterpolation


class mpro_timing:
    def compare_time(self, label):
        if(label == "qc"):
            return self.qc_end - self.qc_start
        elif(label == "host"):
            return self.host_end - self.host_start
        elif(label == "vec"):
            return self.vec_end - self.vec_start
        elif(label == "rRNA"):
            return self.rRNA_end - self.rRNA_start
        elif(label == "contigs"):
            return self.contigs_end - self.contigs_start
        elif(label == "GA_BWA"):
            return self.GA_BWA_end - self.GA_BWA_start
        elif(label == "GA_DMD"):
            return self.GA_DMD_end - self.GA_DMD_start
        elif(label == "TA"):
            return self.TA_end - self.TA_start
        elif(label == "EC"):
            return self.EC_end - self.EC_start
        elif(label == "out"):
            return self.out_end - self.out_start


    def measure_time(self, label, start_or_end):
        if(start_or_end == "start"):
            if(label == "qc"):
                self.qc_start = time.time()
            elif(label == "host"):
                self.host_start = time.time()
            elif(label == "vec"):
                self.vec_start = time.time()
            elif(label == "rRNA"):
                self.rRNA_start = time.time()
            elif(label == "contigs"):
                self.contigs_start = time.time()
            elif(label == "repop"):
                self.repop_start = time.time()
            elif(label == "GA_BWA"):
                self.GA_BWA_start = time.time()
            elif(label == "GA_BLAT"):
                self.GA_BLAT_start = time.time()
            elif(label == "GA_DMD"):
                self.GA_DMD_start = time.time()
            elif(label == "TA"):
                self.TA_start = time.time()
            elif(label == "EC"):
                self.EC_start = time.time()
            elif(label == "out"):
                self.out_start = time.time()

        elif(start_or_end == "end"):
            if(label == "qc"):
                self.qc_end = time.time()
            elif(label == "host"):
                self.host_end = time.time()
            elif(label == "vec"):
                self.vec_end = time.time()
            elif(label == "rRNA"):
                self.rRNA_end = time.time()
            elif(label == "contigs"):
                self.contigs_end = time.time()
            elif(label == "repop"):
                self.repop_end = time.time()
            elif(label == "GA_BWA"):
                self.GA_BWA_end = time.time()
            elif(label == "GA_BLAT"):
                self.GA_BLAT_end = time.time()
            elif(label == "GA_DMD"):
                self.GA_DMD_end = time.time()
            elif(label == "TA"):
                self.TA_end = time.time()
            elif(label == "EC"):
                self.EC_end = time.time()
            elif(label == "out"):
                self.out_end = time.time()
            
            

    def __init__(self):
        self.time_now = dt.today()
        #timing vars
        self.start_time     = time.time()
        self.end_time       = time.time()
        self.qc_start       = time.time()
        self.qc_end         = time.time()        
        self.host_start     = time.time()
        self.host_end       = time.time()        
        self.vec_start      = time.time()
        self.vec_end        = time.time()
        self.rRNA_start     = time.time()
        self.rRNA_end       = time.time()    
        self.repop_start    = time.time()
        self.repop_end      = time.time()
        self.contigs_start  = time.time()
        self.contigs_end    = time.time()
        self.GA_BWA_start   = time.time()
        self.GA_BWA_end     = time.time()
        self.GA_BLAT_start  = time.time()
        self.GA_BLAT_end    = time.time()        
        self.GA_DMD_start   = time.time()
        self.GA_DMD_end     = time.time()
        self.TA_start       = time.time()
        self.TA_end         = time.time()        
        self.EC_start       = time.time()
        self.EC_end         = time.time()
        self.out_start      = time.time()
        self.out_end        = time.time()
