#this code is entirely built to cater to centrifuge's shortcomings.
#it'll change names of ">NODE" to "_NODE", and back

import os
import sys


if __name__ == "__main__":
    contig_in_file = sys.argv[1]
    contig_out_file = sys.argv[2]
    
    new_write_list = []
    
    with open(contig_in_file, "r") as contig_in:
        for line in contig_in:
            if(">NODE" in line):
                new_header = line.replace(">NODE", "_NODE")
                print("old:", line)
                print("new:", new_header)
                new_write_list.append(new_header)
            elif("_NODE" in line):    
                new_header = line.replace("_NODE", ">NODE")
                print("swap old:", line)
                print("swap new:", new_header)
                new_write_list.append(new_header)
            else:
                new_write_list.append(line)
                
                
    with open(contig_out_file, "w") as contig_out:
        for line in new_write_list:
            contig_out.write(line)