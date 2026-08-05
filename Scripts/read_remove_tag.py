#removes tags from the read IDs

import sys
import time
import pandas as pd
import os
import numpy as np
from datetime import datetime as dt
def eat_tag(x):
    new_name = x
    if (x.endswith("/1")): 
        new_name = x.split("/")[0]
    elif(x.endswith("/2")):
        new_name = x.split("/")[0]
        
    return new_name

def remove_tag(input_file, output_file):
    df = pd.read_csv(input_file, header=None, names=[None], sep=r'\n', engine = "python", quoting=3, skip_blank_lines = False)
    end_read_time = time.perf_counter()
    
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))
    df.columns = ["ID", "sequences", "junk", "quality"]
    df["ID"] = df["ID"].apply(lambda x: eat_tag(x))
    np.savetxt(output_file, df.values.flatten(), fmt='%s')
    #df.to_csv(output_file, sep='\n', mode = 'w+', header=False, index=False, quoting = 3)
    
if __name__ == "__main__":
    scinet_flag = False
    if(len(sys.argv) != 3):
        print("something wrong with the args")
        sys.exit()
    else:
        input_fastq_path = os.path.abspath(sys.argv[1])      # expects a path to the fastq file, including the name
        export_path = os.path.abspath(sys.argv[2])           # expects a place to 

        print(dt.today(), "remove tag in:", input_fastq_path)
        print(dt.today(), "remove tag out:", export_path)

        if(not os.path.exists(input_fastq_path)):
            print(dt.today(), input_fastq_path, "does not exist. skipping")
            sys.exit()
        elif(os.path.getsize(input_fastq_path) == 0):
            print(dt.today(), input_fastq_path, "empty. skipping")
            sys.exit()
    
    remove_tag(input_fastq_path, export_path)
