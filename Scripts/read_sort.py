import sys
import time
import pandas as pd

# This code started off as a testbed to develop a method to import FASTQs effectively into pandas.
# It simply imports the FASTQ file, sorts by the ID, and exports it back. 

def sort_and_export(input_file, output_file, direction):

    start_total_call = time.clock()
    print("input file:", input_file)
    df = pd.read_csv(input_file, header=None, names=[None], sep='\n', quoting=3, skip_blank_lines = False)
    end_read_time = time.clock()
    
    # reshaping, aka: unflattening the 1-D array into something meaningful to us
    # we've chosen to leave the @ in the ID.  we didn't feel a strong need to strip it away.
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))
    df.columns = ["ID", "sequences", "junk", "quality"]
    df["ID"] = df["ID"].apply(lambda x: x.split(" ")[0])
    df["ID"] = df["ID"].apply(lambda x: x.replace(".", "_"))
    df = df.sort_values(by=['ID'])
    
        
    end_df_time = time.clock()
    df.to_csv(output_file, sep='\n', mode = 'w+', header=False, index=False, quoting = 3)
    end_total_call = time.clock()
    #-------------------------------------------------------------------------------------------
    
    print("import csv time:", end_read_time - start_total_call, "s")
    print("dataframe interpret time:", end_df_time - end_read_time, "s")
    print("write time:", end_total_call - end_df_time, "s")
    print("total runtime:", end_total_call - start_total_call, "s")

    
if __name__ == "__main__":
    scinet_flag = False
    if(len(sys.argv) != 4):
        print("something wrong with the args")
        sys.exit()
    else:
        input_fastq_path = sys.argv[1]      # expects a path to the fastq file, including the name
        export_path = sys.argv[2]           # expects a place to dump the new fastq, including the name
        direction = sys.argv[3]
    
    sort_and_export(input_fastq_path, export_path, direction)
