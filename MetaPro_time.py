
import time

class metapro_time():
        
    def __init__ (self):
        #timing vars
        self.start_time                     = time.time()
        self.end_time                       = 0
        self.quality_start                  = 0
        self.quality_end                    = 0
        self.cleanup_quality_start          = 0
        self.cleanup_quality_end            = 0
        
        self.host_start                     = 0
        self.host_end                       = 0
        self.cleanup_host_start             = 0
        self.cleanup_host_end               = 0
        
        self.vector_start                   = 0
        self.vector_end                     = 0
        self.cleanup_vector_start           = 0
        self.cleanup_vector_end             = 0
        
        self.rRNA_filter_start              = 0  
        self.rRNA_filter_end                = 0
        self.cleanup_rRNA_filter_start      = 0
        self.cleanup_rRNA_filter_end        = 0   
        
        self.repop_start                    = 0
        self.repop_end                      = 0
        self.cleanup_repop_start            = 0
        self.cleanup_repop_end              = 0
        
        self.assemble_contigs_start         = 0
        self.assemble_contigs_end           = 0
        self.cleanup_assemble_contigs_start = 0
        self.cleanup_assemble_contigs_end   = 0
        
        self.destroy_contigs_start          = 0
        self.destroy_contigs_end            = 0
        self.cleanup_destroy_contigs_start  = 0 
        self.cleanup_destroy_contigs_end    = 0
        
        self.GA_BWA_start                   = 0
        self.GA_BWA_end                     = 0
        self.cleanup_GA_BWA_start           = 0
        self.cleanup_GA_BWA_end             = 0
        
        self.GA_BLAT_start                  = 0
        self.GA_BLAT_end                    = 0
        self.cleanup_GA_BLAT_start          = 0
        self.cleanup_GA_BLAT_end            = 0
        
        self.GA_DIAMOND_start               = 0
        self.GA_DIAMOND_end                 = 0
        self.cleanup_GA_DIAMOND_start       = 0
        self.cleanup_GA_DIAMOND_end         = 0
        
        self.TA_start                       = 0
        self.TA_end                         = 0
        self.cleanup_TA_start               = 0
        self.cleanup_TA_end                 = 0
        
        self.EC_start                       = 0
        self.EC_end                         = 0

        self.EC_DETECT_start                = 0  
        self.EC_DETECT_end                  = 0
        
        self.EC_PRIAM_start                 = 0
        self.EC_PRIAM_end                   = 0
        
        self.EC_DIAMOND_start               = 0
        self.EC_DIAMOND_end                 = 0
        
        self.cleanup_EC_start               = 0
        self.cleanup_EC_end                 = 0
        
        self.Cytoscape_start                = 0
        self.Cytoscape_end                  = 0
        self.cleanup_cytoscape_start        = 0
        self.cleanup_cytoscape_end          = 0
        
