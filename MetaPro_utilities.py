#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

#April 12, 2021
#---------------------------------------------
#MetaPro_utilities.py
#This code houses the various helper functions MetaPro uses to coordinate the multi-threaded traffic.

import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import subprocess as sp
import numpy as np
import re

class mp_seq_handler:
    #internalized because we need the number of files generated + saves another file from being generated
    def __init__(self):
        self.file_count = 0


    


    def fasta_to_dataframe_simple(fasta_file):
        # Read the file content
        with open(fasta_file, 'r') as file:
            content = file.read()
        
        # Use regex to find all entries, capturing just the identifier and sequence
        pattern = r'>([^\s]+)[^\n]*\n([^>]+)'
        
        # Find all matches
        matches = re.findall(pattern, content + '>')
        
        # Convert matches to DataFrame
        df = pd.DataFrame(matches, columns=['identifier', 'sequence'])
        
        # Clean up the sequences
        df['sequence'] = df['sequence'].str.replace('\n', '')
        
        return df

    def fastq_to_dataframe_efficient(self, fastq_file):
        # Read all lines at once
        with open(fastq_file, 'r') as file:
            lines = file.readlines()
        
        # Convert to numpy array for efficient slicing
        lines_array = np.array([line.strip() for line in lines])
        
        # Calculate total number of entries
        num_entries = len(lines_array) // 4
        
        # Extract components using vectorized operations
        identifiers = lines_array[0::4]  # Every 4th line starting from 0
        sequences = lines_array[1::4]    # Every 4th line starting from 1
        quality_scores = lines_array[3::4]  # Every 4th line starting from 3
        
        # Remove the @ prefix from identifiers
        identifiers = np.char.replace(identifiers.astype(str), '@', '', count=1)
        
        # Create DataFrame
        df = pd.DataFrame({
            'ID': identifiers,
            'seq': sequences,
            'qual': quality_scores
        })
        
        return df    
        
    def fastq_to_fasta(self, fastq_file, fasta_file):
        with open(fastq_file, "r") as fin:
            with open(fasta_file, 'w') as fout:
                while True:
                    header = fin.readline().strip()
                    if not header:
                        break
                    seq = fin.readline().strip()
                    fin.readline()  # Skip '+' line
                    fin.readline()  # Skip quality line
                    fout.write(f">{header[1:]}\n{seq}\n")
                    
    def get_split_count(self, split_dir, header):
        list_of_files = os.listdir(os.path.abspath(split_dir))
        split_count = 0
        for item in list_of_files:
            if(item.startswith(header)):
                split_count += 1
        return split_count

    def split_fastq(self, file_name_in, file_name_out, chunks, export_mode):
        print(dt.today(), "FASTQ file name in:", file_name_in)
        #FASTQ has 4 lines per entry.
        
        fastq_df = self.fastq_to_dataframe_efficient(file_name_in)

        print(fastq_df)
        #At this point, we've already got the number of reads.
        #chunks = m.ceil(len(fastq_df) / split_count) #how many sequences each split file will have
        #print("total df length:", len(fastq_df))
        print("chunk size:", chunks)
        if(chunks < 1):
            print(dt.today(), "chunks is set to a default of 90000.  originally:", chunks)
            chunks = 90000
        #if(chunks < 1):
        #    print("split count too large. not enough info to split")
        #    chunks = 1
        #    print("new split count:", split_count)
            
        #for i in range(0, split_count):
        if(export_mode == "fastq"):
            index_count = 0
            while(True):
                print("working on segment :", index_count +1, "of fastq splitter")
                #fancy naming
                new_file_name = file_name_out + "_" + str(index_count) + ".fastq"
                
                
                #split file by selective selection, and writing
                start_index = int(index_count * chunks)
                end_index = int(((index_count+1) * chunks))
                #if(chunks == 1):
                #    end_index += 1 #override on splits that only have 1
                index_count += 1
                fastq_df["ID"] = fastq_df["ID"].apply(lambda x: "@" + x)
                if not(fastq_df.iloc[start_index:end_index, :].empty):
                    fastq_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize = chunks, mode = "w+", index=False, sep='\n', header=False, quoting = 3)
                else:
                    #print("empty frame detected.  no sense in running the rest of the fastq splitter")
                    break
        else:
            #export as fasta.  Save a step.
            index_count = 0
            file_count = 0
            while(True):
                new_file_name = file_name_out + "_" + str(index_count) + ".fasta"
                start_index = int(index_count * chunks)
                end_index = int(((index_count+1) * chunks))
                index_count += 1
                subselect_df = fastq_df[["ID", "seq"]]
                subselect_df = subselect_df.iloc[start_index:end_index, :]
                subselect_df["ID"] = subselect_df["ID"].apply(lambda x: ">" + x)
                if(not subselect_df.empty):
                    print(dt.today(), "exporting FASTQ split to FASTA:", new_file_name)
                    #subselect_df.to_csv(new_file_name, mode = "w+", index = False, sep = "\n", header = False, quoting=3)
                    np.savetxt(new_file_name, subselect_df.values, delimiter = "\n", fmt='%s')
                    file_count += 1
                else:
                    break


        return (file_count)
    

    def split_fasta(self, file_name_in, file_name_out, chunks):#split_count = 4):
        #modded to take in fixed chunks
        fasta_df = self.fasta_to_dataframe_pandas(file_name_in)
        fasta_df.columns = ["row"]
        #There's apparently a possibility for NaNs to be introduced in the raw fasta.  We have to strip it before we process (from DIAMOND proteins.faa)
        fasta_df.dropna(inplace=True)
        new_df = pd.DataFrame(fasta_df.loc[fasta_df.row.str.contains('>')])  # grab all the IDs
        new_df.columns = ["names"]
        new_data_df = fasta_df.loc[~fasta_df.row.str.contains('>')]  # grab the data
        new_data_df.columns = ["data"]
        fasta_df = new_df.join(new_data_df, how='outer')  # join them into a 2-col DF
        fasta_df["names"] = fasta_df.fillna(method='ffill')  # fill in the blank spaces in the name section
        fasta_df.dropna(inplace=True)  # remove all rows with no sequences
        fasta_df.index = fasta_df.groupby('names').cumcount()  # index it for transform
        temp_columns = fasta_df.index  # save the index names for later
        fasta_df = fasta_df.pivot(values='data', columns='names')  # pivot
        fasta_df = fasta_df.T  # transpose
        fasta_df["seq"] = fasta_df[fasta_df.columns[:]].apply(lambda x: "".join(x.dropna()), axis=1)  # consolidate all cols into a single sequence
        fasta_df.drop(temp_columns, axis=1, inplace=True)
        # At this point, we've already got the number of reads.
        #chunks = m.ceil(len(fasta_df) / split_count)  # how many sequences each split file will have
        print("total df length:", len(fasta_df))
        print("chunk size:", chunks)
        if (chunks < 1):
            print("split count too large. not enough info to split")
            chunks = 1
            


        #for i in range(0, split_count):
        index_count = 0
        while(True):
            print("working on segment :", index_count + 1, "of FASTA splitter" )
            # fancy naming
            new_file_name = file_name_out + "_" + str(index_count) + ".fasta"
            
            # split file by selective selection, and writing
            start_index = int(index_count * chunks)
            end_index = int(((index_count + 1) * chunks))
            # if(chunks == 1):
            #    end_index += 1 #override on splits that only have 1
            index_count += 1
            if not (fasta_df.iloc[start_index:end_index, :].empty):
                fasta_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize=chunks, mode="w+", sep='\n', header=False)
            else:
                print("empty frame detected.  no sense in running the rest of FASTA splitter")
                break
        return (index_count + 1)

class mp_util:
    def __init__(self, config_dict, dir_dict): #, config_obj):
        self.mp_store = []
        self.output_folder_path = dir_dict["main"]
        #self.paths = config_obj
        self.bypass_log_name = os.path.basename(config_dict["bypass_log"])
        self.mpsh = mp_seq_handler()

    def get_query_size_gb(self, query_file):
        """Get query file size in GB"""
        try:
            file_size_bytes = os.path.getsize(query_file)
            file_size_gb = file_size_bytes / (1024**3)  # Convert bytes to GB
            return file_size_gb
        except FileNotFoundError:
            print(f"Warning: Query file {query_file} not found. Using default size.")
            return 7.0  # Default fallback
        except Exception as e:
            print(f"Error getting file size for {query_file}: {e}. Using default size.")
            return 7.0  # Default fallback

    def determine_dmd_mem_limit(self, query_file):
        mem = psu.virtual_memory()
        total_mem = mem.total/(1024*1024*1024)
        
        # Get actual query file size
        query_file_size_gb = self.get_query_size_gb(query_file)
        
        # Sequential processing - use 80% of total RAM for single job
        safety_buffer = 0.8
        available_memory = total_mem * safety_buffer
        
        # Account for actual DIAMOND memory usage (not full 373GB database)
        nr_index_memory = 15   # GB (database metadata/index structures only)
        working_memory = 10    # GB (alignment buffers, etc.)
        query_memory = query_file_size_gb + 2  # Add 2GB overhead
        base_memory = nr_index_memory + query_memory + working_memory
        
        # Available memory for block-size after accounting for base usage
        available_for_blocks = available_memory - base_memory
        
        print(f"Total RAM: {total_mem:.1f}GB, Available for blocks: {available_for_blocks:.1f}GB")
        
        # Set block-size based on available memory
        if available_for_blocks > 100:     # 400GB+ systems
            block_size = min(40.0, available_for_blocks / 3)
        elif available_for_blocks > 50:    # 200GB+ systems
            block_size = min(25.0, available_for_blocks / 2)
        elif available_for_blocks > 20:    # 180GB systems
            block_size = min(15.0, available_for_blocks / 2)
        elif available_for_blocks > 5:     # Lower memory systems
            block_size = min(8.0, available_for_blocks / 2)
        elif available_for_blocks > 0:     # Minimal systems
            block_size = min(4.0, available_for_blocks)
        else:
            print("WARNING: Very low available memory for DIAMOND processing.")
            block_size = 2.0  # Conservative fallback
        return block_size

    def concatenate_files_efficient(self, file_list, output_file):
        with open(output_file, 'w') as outfile:
            for filename in file_list:
                with open(filename, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

    def mem_checker(self, threshold):
        #threshold is a percentage for available memory.  
        
        mem = psu.virtual_memory()
        available_mem = mem.available
        total_mem = mem.total
        
        available_pct = 100 * available_mem / total_mem
        
        if(float(available_pct) <= float(threshold)):
            return False
        else:
            return True
            
    def mem_footprint_checker(self, count, mem_footprint):
        #assumes some pre-determinted footprint for the program in question.
        #then it doesn't let the program run
        mem = psu.virtual_memory()
        total_mem = mem.total/(1024*1024*1024)

        occupied_mem = count * mem_footprint

        max_limit = total_mem * 0.9
        if(occupied_mem > max_limit):
            return False
        else:
            return True

            
    def make_folder(self, folder_path):
        if not (os.path.exists(folder_path)):
            os.makedirs(folder_path)
            
    def delete_folder_simple(self, folder_path):
        if(os.path.exists(folder_path)):
            print(dt.today(), "Deleting:", folder_path)
            shutil.rmtree(folder_path)
            print(dt.today(), "finished deleting:", folder_path)

    def delete_folder(self, folder_path):
        if (os.path.exists(os.path.join(folder_path, "data"))):
            print("deleting", os.path.join(folder_path, "data"))
            shutil.rmtree(os.path.join(folder_path, "data"))
        else:
            print("can't delete folder: doesn't exist:", folder_path)
            
    def compress_folder(self, folder_path):
        zip_loc = os.path.join(folder_path, "data")
        z = zipfile.ZipFile(folder_path + "_data.zip", "a", zipfile.ZIP_DEFLATED)
        print("compressing interim files:", folder_path)
        for root, dirs, files in os.walk(zip_loc):
            #print("root:", root)
            #print("dirs:", dirs)
            #print("files:", files)
            #print("===============================")
            for file in files:
                z.write(os.path.join(root, file))
        z.close()
            
    def write_to_bypass_log(self, folder_path, message):
        bypass_log_path = os.path.join(folder_path, self.bypass_log_name)
        with open(bypass_log_path, "a") as bypass_log:
            bypass_log.write("\n")
            new_message = message + "\n"
            bypass_log.write(new_message)
            


    def check_bypass_log(self, folder_path, message):
        stop_message = "stop_" + str(message)
        bypass_keys_list = list()
        bypass_log_path = os.path.join(folder_path, self.bypass_log_name)
        if(os.path.exists(bypass_log_path)):
            with open(bypass_log_path, "r") as bypass_log:
                for line in bypass_log:
                    bypass_key = line.strip("\n")
                    bypass_keys_list.append(bypass_key)
            
            if(stop_message in bypass_keys_list):
                print(dt.today(), "stopping at:", message)
                print("to continue, remove:", stop_message, "from the bypass_log")
                sys.exit("brakes engaged")
            
            elif(message in bypass_keys_list):
                print(dt.today(), "bypassing:", message)
                return False
            
            else:
                print(dt.today(), "running:", message) 
                return True
        else:
            open(bypass_log_path, "a").close()
            print(dt.today(), "no bypass log.  running:", message)
            return True
            
    def conditional_write_to_bypass_log(self, label, stage_folder, file_name): 
        #convenience for checking if a file exists, and writing to the bypass log
        if self.check_bypass_log (self.output_folder_path, label):
            file_path = os.path.join(self.output_folder_path, stage_folder, file_name)
            if(os.path.exists(file_path)):
                self.write_to_bypass_log(self.output_folder_path, label)
        

    # Used to determine quality encoding of fastq sequences.
    # Assumes Phred+64 unless there is a character within the first 10000 reads with encoding in the Phred+33 range.
    def check_code(self, segment):
        encoding = 64
        for item in segment:
            if(ord(item) < 64):
                encoding = 33
                break
        return encoding

    def determine_encoding(self, fastq_path):
        #import the first 10k lines, then check the quality scores.
        #if the quality score symbols are below 76, it's phred33.  
        print("using:", [fastq_path])
        fastq_df = self.mpsh.fastq_to_dataframe_efficient(fastq_path)
        quality_encoding = fastq_df["qual"].apply(lambda x: self.check_code(x)).mean() #condense into a single number.
        if(quality_encoding == 64): #all must be 64 or else it's 33
            quality_encoding = 64
        else:
            quality_encoding =  33
        return str(quality_encoding)


    # handles where to kill the pipeline, due to the prev step behaving badly
    # logic is:  if the files inside the dep_path (or dep job label shortcut to the final_results)
    #            are empty, then there's an error.  kill the pipeline 
    def check_where_kill(self, dep_job_label=None, dep_path=None):
        if dep_job_label is None:
            if dep_path is None:
                return True
            else:
                dep_job_path = dep_path
        else:
            dep_job_path = os.path.join(dep_job_label, "final_results")

        file_list = os.listdir(dep_job_path)
        if len(file_list) > 0:
            for item in file_list:
                file_check_path = os.path.join(dep_job_path, item)
                if (os.path.getsize(file_check_path)) == 0:
                    print("empty file detected: rerunning stage")
                    sys.exit("bad dep")
            # run the job, silently
            return True
        else:
            print("stopping the pipeline.  dependencies don't exist")
            sys.exit("no dep")


    # handles where to auto-resume the pipeline on a subsequent run
    # label: used as a shorthand for paths we expect
    # full path: a bypass for when we want to use it for detecting a location that doesn't fall into the normal format (final_results)
    # dep: for checking if the job's dependencies are satisfied-> meant to point to the last stage's "final_results"
    # logic is: if the full_path has no files (or the job label shortcut to final_results)
    #           and the dependencies are ok, start the stage
    #Aug 19, 2019: There's a tweak to this:  DIAMOND will generate zero-size files, due to no-matches
    #it's allowable.
    def check_where_resume(self, job_label=None, full_path=None, dep_job_path=None, file_check_bypass = False):
        if(not file_check_bypass):
            self.check_where_kill(dep_job_path)
        if job_label:
            job_path = os.path.join(job_label, "final_results")
        else:
            job_path = full_path

        print("looking at:", job_path)

        if os.path.exists(job_path):
            file_list = os.listdir(job_path)
            if(not file_check_bypass):
                if len(file_list) > 0:
                    for item in file_list:
                        file_check_path = os.path.join(job_path, item)
                        if (os.path.getsize(file_check_path)) == 0:
                            print("empty file detected: rerunning stage")
                            return False
                    print("bypassing!")
                    return True
                else:
                    print("no files detected: running")
                    return False
            else:
                print("bypassing for special reasons")
                return True
        else:
            print("doesn't exist: running")
            return False

    def make_script(self, job_path, command_list):
        
        # create the pbs job, and launch items
        # job path: name and location of where to write the script.
        # command list:  list of command statements for writing
        # returns back the job ID given from sbatch

        # docker mode: single cpu
        # no ID, no sbatch.  just run the command
        
        shell_script_full_path = job_path
        print("using:", shell_script_full_path)
        #time.sleep(3)
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


    def run_subjob_simple(self, job_location, commands):
        #just launches a job.  no multi-process.
        process = mp.Process(
            target=self.make_script,
            args=(job_location, commands)
        )
        process.start()
        process.join()

<<<<<<< HEAD
    def run_subjob_with_mp_store(self, job_location, commands):
=======
    def run_subjob_with_mp_store(self, job_location, job_label, commands):
>>>>>>> origin/feature/2025_tool_update
        #just launches a job.  no multi-process.
        process = mp.Process(
            target=self.make_script,
            args=(job_location, commands)
        )
        process.start()
        self.mp_store.append(process)

    def launch_only_simple(self, commands):
        process = mp.Process(
            target=self.launch_only,
            args=(commands, len(commands))
        )
        process.start()
        process.join()
        
<<<<<<< HEAD
    def subdivide_and_launch(self, job_delay, mem_threshold, job_limit, job_file, commands):
=======
    def subdivide_and_launch(self, job_delay, mem_threshold, job_limit, job_location, job_label, commands):
>>>>>>> origin/feature/2025_tool_update
        #just launches a job.  no multi-process.
        #Jan 25, 2022: now adding job controls.
        job_counter = 0
        for item in commands:
            job_location = os.path.dirname(job_file)
            job_name = os.path.basename(job_file)
            full_job_name = "job_" + str(job_counter) + "_" + job_name
            full_job_path = os.path.join(job_location, full_job_name)
            job_counter += 1
            job_submitted = False
            while(not job_submitted):
                if(len(self.mp_store) < job_limit):
                    if(self.mem_checker(mem_threshold)):

                        process = mp.Process(
                            target=self.make_script,
                            args=(full_job_path, [item])
                        )
                        process.start()
                        self.mp_store.append(process)
                        print(dt.today(), full_job_path, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                        job_submitted = True
                    else:
                        time.sleep(job_delay)
                else:
                    self.wait_for_mp_store()
            time.sleep(job_delay)
        #final wait for everything to be done
        self.wait_for_mp_store()
                
        
    def launch_only_with_hold(self, mem_threshold, job_limit, job_delay, job_name, command):
        #launch a job in launch-only mode
        job_submitted = False
        while(not job_submitted):
            if(len(self.mp_store) < job_limit):
                if(self.mem_checker(mem_threshold)):
                    process = mp.Process(
                        target = self.launch_only,
                        args = (command, len(command))
                    )
                    process.start()
                    self.mp_store.append(process)
                    print(dt.today(), job_name, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    job_submitted = True
                else:
                    #print(dt.today(), job_name, "Pausing. mem limit reached:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    time.sleep(job_delay)
            else:
                print(dt.today(), "job limit reached.  waiting for queue to flush")
                self.wait_for_mp_store()
        time.sleep(job_delay)
        

    def run_subjob_with_hold(self, mem_threshold, job_limit, job_delay, job_location, command):
        #launch a job in launch-with-create mode
        job_submitted = False
        while(not job_submitted):
                
            if(len(self.mp_store) < job_limit):
                if(self.mem_checker(mem_threshold)):
                    process = mp.Process(
                        target = self.make_script,
                        args = (job_location, command)
                    )
                    process.start()
                    self.mp_store.append(process)
                    print(dt.today(), job_location, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    job_submitted = True
                else:
                    #print(dt.today(), job_name, "Pausing. mem limit reached:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    time.sleep(job_delay)
            else:
                print(dt.today(), "job limit reached.  waiting for queue to flush")
                self.wait_for_mp_store()
        #final wait
        self.wait_for_mp_store()

    def run_subjob_no_final(self, mem_threshold, job_limit, job_delay, job_location, command):
        #launch a job in launch-with-create mode
        job_submitted = False
        while(not job_submitted):
                
            if(len(self.mp_store) < job_limit):
                if(self.mem_checker(mem_threshold)):
                    process = mp.Process(
                        target = self.make_script,
                        args = (job_location, command)
                    )
                    process.start()
                    self.mp_store.append(process)
                    print(dt.today(), job_location, "job submitted.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    job_submitted = True
                else:
                    #print(dt.today(), job_name, "Pausing. mem limit reached:", psu.virtual_memory().available/(1024*1024*1000), "GB", end='\r')
                    time.sleep(job_delay)
            else:
                print(dt.today(), "job limit reached.  waiting for queue to flush")
                self.wait_for_mp_store()
        
    def run_subjob_with_mem_footprint(self, mem_footprint, job_limit, job_location, command):
        #launch a job in launch-with-create mode
        #this controller won't be optimized for the system. It's made to keep the node from exploding.
        job_submitted = False
        
        while(not job_submitted):

            if(len(self.mp_store) < job_limit):    
                if(self.mem_footprint_checker(len(self.mp_store), mem_footprint)):
                    process = mp.Process(
                        target = self.make_script,
                        args = (job_location, command)
                    )
                    process.start()
                    self.mp_store.append(process)
                    job_submitted = True
                    
                    print(dt.today(), job_location, "job submitted.  mem:", len(self.mp_store) * mem_footprint, "GB", end='\r')
                else:
                    print(dt.today(), "job limit reached.  waiting for queue to flush")
                    self.wait_for_mp_store()
            else:
                print(dt.today(), "job limit reached.  waiting for queue to flush")
                self.wait_for_mp_store()    
                

    #check if all jobs ran
    def check_all_job_markers(self, job_marker_list, final_folder_checklist):
        time.sleep(2)
        #if it's already been created, that means the job was killed.
        if(os.path.exists(final_folder_checklist)):
            print(dt.today(), final_folder_checklist, "exists: adding to it")
            #open it, import it.
            with open(final_folder_checklist, "r") as old_list:
                for line in old_list:
                    cleaned_line = line.strip("\n")
                    job_marker_list.append(cleaned_line)
            #then overwrite it
            with open(final_folder_checklist, "w") as checklist:
                for item in job_marker_list:
                    checklist.write(item +"\n")
                    
                for item in job_marker_list:
                    if(not os.path.exists(item)):
                        print(dt.today(), item, "not found.  kill the pipe.  restart this stage")
                        sys.exit("not all jobs completed")
            
        else:
            
            with open(final_folder_checklist, "w") as checklist:
                for item in job_marker_list:
                    checklist.write(item +"\n")
                    
                for item in job_marker_list:
                    if(not os.path.exists(item)):
                        print(dt.today(), item, "not found.  kill the pipe.  restart this stage")
                        sys.exit("not all jobs completed")

    def limited_wait_for_mp_store(self, limit):
        print(dt.today(), "limit for mp_store:", limit)
        if(len(self.mp_store) >= limit):
            print(dt.today(), "mp_store limit reached. pausing to flush")
            for item in self.mp_store:
                item.join()
            self.mp_store[:] = []
            print(dt.today(), "mp_store flushed. continuing")
    
    
    def wait_for_mp_store(self):
        print(dt.today(), "closing down processes: ", len(self.mp_store))
        count = 0
        for item in self.mp_store:
            
            print(dt.today(), "closed down: " + str(count) +  "/" + str(len(self.mp_store)) + "            ", end = "\r") 
            count += 1
            item.join()
        self.mp_store[:] = []

    def clean_or_compress(self, analysis_path, keep_all, keep_stage):
        if(keep_all == "no" and keep_stage == "no"):
            self.delete_folder(analysis_path)
        elif(keep_all == "compress" or keep_stage == "compress"):
            self.compress_folder(analysis_path)
            self.delete_folder(analysis_path)


    def launch_stage_simple(self, job_path, command_list, keep_all, keep_job):
        #wrapper for simple job launches (quality, host)
        #cleanup_job_start = 0
        #cleanup_job_end = 0
        print("job path:", job_path)
        job_label = os.path.basename(job_path)

        if self.check_bypass_log(self.output_folder_path, job_label):
            print(dt.today(), "NEW CHECK running:", job_label)
            self.run_subjob_simple(job_path, command_list)
            
            #self.write_to_bypass_log(self.output_folder_path, job_label)
            #cleanup_job_start = time.time()
            self.clean_or_compress(job_path, keep_all, keep_job)
            #cleanup_job_end = time.time()    
        else:
            print(dt.today(), "skipping job:", job_label)

        #return cleanup_job_start, cleanup_job_end