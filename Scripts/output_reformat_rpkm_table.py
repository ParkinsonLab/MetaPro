import os
import sys

#This code just flattens the ECs so that they'll duplicate the row info if there's more than 1 EC per gene

if __name__ == "__main__":
    
    input_rpkm = sys.argv[1]
    output_rpkm = sys.argv[2]
    
    with open(output_rpkm, "w") as new_file:
        with open(input_rpkm, "r") as rpkm_file:
            for line in rpkm_file:
                list_line = line.split("\t")
                if(list_line[0] == "GeneID"):
                    new_file.write(line)
                else:
                    the_first_part = list_line[0:3]
                    ec_part = list_line[3]
                    the_rest = list_line[4:]
                    
                    ec_list = ec_part.split("|")
                    
                    new_string = ""
                    for entry in the_first_part:
                        new_string += entry + "\t"
                    print(new_string)
                    
                    for item in ec_list:
                        new_suffix = item + "\t"
                        
                        for entry in the_rest:
                            new_suffix += entry + "\t"
                        new_suffix = new_suffix[:-1]

                        new_file.write(new_string + new_suffix)
                        
