#this isn't the final design.  this is just to assemble the pieces to test the efficiency of the new GA lib-maker
#new thing: copies all files, and their index files.  
#there are a handful of exceptions
#actually, this file should never be used. we shouldn't copy.  this should be code inside the pipe to just generate commands using the names_list

#note: only the class-level split has a size problem. everyone else is under 4.5GB



import os
import sys
import time 
from datetime import datetime as dt
import shutil as sh

def import_lib_names(names_file):
    files_list = []
    with open(names_file, "r") as names_in:
        for line in names_in:
            if("can't" in line):
                continue
            else:
                
                cleaned_line = line.strip("\n")
                
                src_path = cleaned_line.split("|")[3]
                if(os.path.exists(src_path)):
                    files_list.append(src_path)
    return files_list
    
def copy_all_files(root_name, lib_dir, dest_path):
    amb_path = os.path.join(lib_dir, root_name + ".amb")
    ann_path = os.path.join(lib_dir, root_name + ".ann")
    bwt_path = os.path.join(lib_dir, root_name + ".bwt")
    pac_path = os.path.join(lib_dir, root_name + ".pac")
    sa_path = os.path.join(lib_dir, root_name + ".sa")
    root_path = os.path.join(lib_dir, root_name)

    sh.copy(root_path, dest_path)
    sh.copy(amb_path, dest_path)
    sh.copy(ann_path, dest_path)
    sh.copy(bwt_path, dest_path)
    sh.copy(pac_path, dest_path)
    sh.copy(sa_path, dest_path)
    
def copy_fasta(root_name, lib_dir, dest_path):
    root_path = os.path.join(lib_dir, root_name)
    sh.copy(root_path, dest_path)
        
    

if __name__ == "__main__":
    names_file = sys.argv[1]
    lib_dir = sys.argv[2]
    names_list = import_lib_names(names_file)
    dest_dir = sys.argv[3]
    mode = sys.argv[4]
    
    for item in names_list:
        
        if(item == "1236.fasta" or item  == "1760.fasta" or item == "28211.fasta"):
            root_name = item.split(".")[0]
            for i in range(0, 3):
                real_name = root_name + "_" + str(i) + ".fasta"
                if(mode == "all"):
                    copy_all_files(real_name, lib_dir, dest_dir)
                else:
                    copy_fasta(real_name, lib_dir, dest_dir)
        
        elif(item == "91061.fasta"):
            for i in range(0, 2):
                root_name = item.split(".")[0]
                real_name = root_name + "_" + str(i) + ".fasta"
                if(mode == "all"):
                    copy_all_files(real_name, lib_dir, dest_dir)
                else:
                    copy_fasta(real_name, lib_dir, dest_dir)
        else:
            if(mode == "all"):
                copy_all_files(item, lib_dir, dest_dir)
            else:
                copy_fasta(item, lib_dir, dest_dir)
    