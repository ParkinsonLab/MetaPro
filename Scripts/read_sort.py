import sys
import time
import pandas as pd
import os
from datetime import datetime as dt
import numpy as np

# This code started off as a testbed to develop a method to import FASTQs effectively into pandas.
# It simply imports the FASTQ file, sorts by the ID, and exports it back. 

def sort_and_export(input_file, output_file):

    start_total_call = time.perf_counter()
    print("input file:", input_file)
    df = pd.read_csv(input_file, header=None, names=[None], sep=r'\n', engine = "python", quoting=3, skip_blank_lines = False)
    end_read_time = time.perf_counter()
    
    # reshaping, aka: unflattening the 1-D array into something meaningful to us
    # we've chosen to leave the @ in the ID.  we didn't feel a strong need to strip it away.
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))
    df.columns = ["ID", "sequences", "junk", "quality"]
    df["ID"] = df["ID"].apply(lambda x: x.split(" ")[0])
    df["ID"] = df["ID"].apply(lambda x: x.replace(".", "_"))
    df = df.sort_values(by=['ID'])
    
        
    end_df_time = time.perf_counter()
    np.savetxt(output_file, df.values.flatten(), fmt='%s')
    #df.to_csv(output_file, sep='\n', mode = 'w+', header=False, index=False, quoting = 3)
    end_total_call = time.perf_counter()
    #-------------------------------------------------------------------------------------------
    
    print("import csv time:", end_read_time - start_total_call, "s")
    print("dataframe interpret time:", end_df_time - end_read_time, "s")
    print("write time:", end_total_call - end_df_time, "s")
    print("total runtime:", end_total_call - start_total_call, "s")

    
if __name__ == "__main__":
    print(dt.today(), "START TIME")
    input_fastq_path = os.path.abspath(sys.argv[1])      # expects a path to the fastq file, including the name
    export_path = os.path.abspath(sys.argv[2])           # expects a place to dump the new fastq, including the name
    
    
    sort_and_export(input_fastq_path, export_path)
    print(dt.today(), "DONE")