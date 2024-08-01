#------------------------------------------------------------------
# MetaPro GA commands
# Aug 01, 2024: refactor/cleanup/reorg reasons. 

class mt_GA_commands:
    def __init__(self, no_host, Config_path, args_pack, Quality_score=33, tutorial_keyword = None, sequence_path_1=None, sequence_path_2=None, sequence_single=None, sequence_contigs = None):



    def create_GA_pre_scan_command(self, stage_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        final_folder    = os.path.join(subfolder, "final_results")
        data_folder     = os.path.join(subfolder, "data")
        dest_folder     = os.path.join(data_folder, "4_libs")
        jobs_folder     = os.path.join(subfolder, "jobs")
        
        self.make_folder(data_folder)
        self.make_folder(dest_folder)
        self.make_folder(jobs_folder)
        self.make_folder(final_folder)
        
        ga_get_lib = ">&2 echo GA pre-scan get libs | "
        ga_get_lib += self.tool_path_obj.Python + " "
        ga_get_lib += self.tool_path_obj.GA_pre_scan_get_lib + " "
        ga_get_lib += os.path.join(data_folder, "3_wevote", "taxonomic_classifications.tsv") + " "
        ga_get_lib += self.tool_path_obj.taxid_tree + " "
        ga_get_lib += self.tool_path_obj.nodes + " "
        ga_get_lib += os.path.join(dest_folder, "lib_list.txt") + " "
        ga_get_lib += os.path.join(dest_folder, "lib_reject.txt") + " "
        ga_get_lib += self.tool_path_obj.source_taxa_DB + " "
        ga_get_lib += str(self.tool_path_obj.taxa_exist_cutoff)
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        
        return [ga_get_lib + " && " + make_marker]
        
        
    def create_GA_pre_scan_assemble_lib_command(self, stage_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        final_folder    = os.path.join(subfolder, "final_results")
        data_folder     = os.path.join(subfolder, "data")
        dest_folder     = os.path.join(data_folder, "4_libs")
        jobs_folder     = os.path.join(subfolder, "jobs")
        
        self.make_folder(data_folder)
        self.make_folder(dest_folder)
        self.make_folder(jobs_folder)
        self.make_folder(final_folder)
    
        assemble_lib = ">&2 echo GA assemble libs | " 
        assemble_lib += self.tool_path_obj.Python + " " 
        assemble_lib += self.tool_path_obj.GA_pre_scan_assemble_lib + " "
        assemble_lib += os.path.join(dest_folder, "lib_list.txt") + " " 
        assemble_lib += self.tool_path_obj.source_taxa_DB +  " "
        assemble_lib += final_folder +  " " 
        assemble_lib += "all"
        
        #index_lib = "for i in $(ls " + final_folder + ");" + " "
        #index_lib += "do " + self.tool_path_obj.BWA + " index" + " "
        #index_lib += final_folder + "/$i; done" 
        
        make_marker = "touch" + " "
        make_marker += os.path.join(jobs_folder, marker_file)
        
        self.tool_path_obj.DNA_DB = final_folder
        
        #return [assemble_lib + " && " + index_lib + " && " + make_marker]
        return [assemble_lib + " && " + make_marker]
        
 
    def create_split_ga_fastq_data_command(self, stage_name, dependency_stage_name, category, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        final_folder    = os.path.join(subfolder, "final_results")
        data_folder     = os.path.join(subfolder, "data")
        split_folder    = os.path.join(final_folder, category)#, os.path.join(data_folder, "0_read_split", category)
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        
        jobs_folder     = os.path.join(data_folder, "jobs")
        self.make_folder(subfolder)
        self.make_folder(final_folder)
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
        final_folder    = os.path.join(subfolder, "final_results")
        split_folder    = os.path.join(final_folder, category)#os.path.join(data_folder, "0_read_split", category)
        dep_folder      = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
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
        bwa_job += self.tool_path_obj.BWA + " mem -t " + self.threads_str + " "
        bwa_job += ref_path + " "
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
        pp_folder       = os.path.join(data_folder, "2_bwa_pp")
        final_folder    = os.path.join(subfolder, "final_results")
        dep_loc         = os.path.join(self.Output_Path, dependency_stage_name, "final_results")
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        self.make_folder(pp_folder)
        
        reads_in    = query_file
        bwa_in      = os.path.join(bwa_folder, sample_root_name + "_" + ref_tag + ".sam")
        reads_out = ""
        if(self.tool_path_obj.GA_DB_mode == "multi"):
            print(dt.today(), "BWA_pp running in split-mode")
            reads_out   = os.path.join(pp_folder, sample_root_name + "_" + ref_tag + ".fasta")
        else:
            print(dt.today(), "BWA_pp running in single-mode")
            reads_out = os.path.join(final_folder, sample_root_name + "_" + ref_tag + ".fasta")
        

        map_read_bwa = ">&2 echo " + str(dt.today()) + " GA BWA PP generic: " + sample_root_name + " | "
        map_read_bwa += self.tool_path_obj.Python + " "
        map_read_bwa += self.tool_path_obj.Map_reads_gene_BWA + " "
        map_read_bwa += str(self.tool_path_obj.BWA_cigar_cutoff) + " "
        map_read_bwa += ref_path + " "
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
        ]

        return COMMANDS_Annotate_BWA


 

    def create_BWA_copy_contig_map_command(self, stage_name, dependency_stage_name, marker_file):
        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        bwa_folder      = os.path.join(data_folder, "1_bwa")
        pp_folder       = os.path.join(data_folder, "2_bwa_pp")
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

    def create_merge_BWA_fasta_command(self, stage_name, query_file, marker_file):
        sample_root_name = os.path.basename(query_file)
        sample_root_name = os.path.splitext(sample_root_name)[0]

        subfolder       = os.path.join(self.Output_Path, stage_name)
        data_folder     = os.path.join(subfolder, "data")
        bwa_folder      = os.path.join(data_folder, "1_bwa")
        split_folder    = os.path.join(data_folder, "0_read_split")
        pp_folder       = os.path.join(data_folder, "2_bwa_pp")
        final_folder    = os.path.join(subfolder, "final_results")
        
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        self.make_folder(pp_folder)

        merge_bwa_fastas = ">&2 echo " + str(dt.today()) + " GA BWA merge leftover reads " + sample_root_name + " | "
        merge_bwa_fastas += self.tool_path_obj.Python + " "
        merge_bwa_fastas += self.tool_path_obj.GA_merge_fasta + " "
        merge_bwa_fastas += pp_folder + " " 
        merge_bwa_fastas += sample_root_name + " " 
        merge_bwa_fastas += final_folder

        make_marker = ">&2 echo merge BWA leftover fastas: " + marker_file + " | " 
        make_marker += "touch" + " " 
        make_marker += os.path.join(jobs_folder, marker_file)

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
        blat_command += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
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
        blat_pp += self.tool_path_obj.Python + " "
        blat_pp += self.tool_path_obj.Map_reads_gene_BLAT + " "
        blat_pp += str(self.tool_path_obj.BLAT_identity_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_length_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_score_cutoff) + " "
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
        blat_pp += self.tool_path_obj.Python + " "
        blat_pp += self.tool_path_obj.Map_reads_gene_BLAT + " "
        blat_pp += str(self.tool_path_obj.BLAT_identity_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_length_cutoff) + " "
        blat_pp += str(self.tool_path_obj.BLAT_score_cutoff) + " "
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
        merge_blat_fastas += self.tool_path_obj.Python + " "
        merge_blat_fastas += self.tool_path_obj.GA_merge_fasta + " "
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
        diamond_annotate += self.tool_path_obj.DIAMOND
        diamond_annotate += " blastx -p " + self.threads_str
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
        dep_0_path      = os.path.join(self.Output_Path, dep_0_name, "final_results")   #assemble-contigs
        dep_1_path      = os.path.join(self.Output_Path, dep_1_name, "final_results")   #bwa
        dep_2_path      = os.path.join(self.Output_Path, dep_2_name, "final_results")   #blat
        dep_3_path      = os.path.join(self.Output_Path, dep_3_name, "final_results")   #dmd
        jobs_folder     = os.path.join(data_folder, "jobs")
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        self.make_folder(jobs_folder)
        
        final_merge_fastq = self.tool_path_obj.Python + " "
        final_merge_fastq += self.tool_path_obj.GA_final_merge_fasta + " "
        final_merge_fastq += dep_0_path + " "
        final_merge_fastq += dep_3_path + " "
        final_merge_fastq += self.read_mode + " "
        final_merge_fastq += final_folder
        
        final_merge_proteins = self.tool_path_obj.Python + " "
        final_merge_proteins += self.tool_path_obj.GA_final_merge_proteins + " "
        final_merge_proteins += dep_1_path + " "
        final_merge_proteins += dep_2_path + " "
        final_merge_proteins += dep_3_path + " "
        final_merge_proteins += final_folder
        
        final_merge_maps = self.tool_path_obj.Python + " "
        final_merge_maps += self.tool_path_obj.GA_final_merge_maps + " "
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