#splits a fastq and converts it into fasta segments.
#used specifically for rRNA, to reduce the number of files

import os
import sys


def import_fastq(fastq_file):
    fastq_dict = dict()
    with open(fastq_file, "r") as fastq_in:
        for line in fastq_in:
            if(line.startswith("@")):
                fastq_ID = line.strip("\n").split(" ")[0]


if __name__ == "__main__":
    fastq_in_file = sys.argv[1]
    fasta_out_dir = sys.argv[2]
    read_split_count = sys.argv[3]

