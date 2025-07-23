import os
import sys
import time
from datetime import datetime as dt

#Jul 2025: due to some weird name mangling with BT2, the contig headers needed reformatting.
#This just removes the extra ">" found in the name from mgm1


if __name__ == "__main__":
    contigs_in_file = sys.argv[1]
    contigs_out_file = sys.argv[2]
    with open(contigs_out_file, "w") as contigs_out:
        with open(contigs_in_file, "r") as contigs_in:
            for line in contigs_in:
                if(line.startswith(">")):
                    header = str(line.strip("\n"))
                    header = str(header[1:])
                    header = header.replace(">", "|")
                    header = ">" + header + "\n"
                    
                    contigs_out.write(header)
                else:
                    contigs_out.write(line)