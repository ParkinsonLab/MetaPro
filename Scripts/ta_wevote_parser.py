import os
import sys
import pandas as pd
#Sifts through the WEVOTE results in a neat way. 
#We only want the first and last rows. 
#Avoids using awk, so it's portable to other OSes, should we want that functionality.
#Mar 15, 2021: this needs to handle cases without a contig map (meaning no tool_0, 1, 2)

if __name__ == "__main__":
    wevote_file = sys.argv[1]
    output_file = sys.argv[2]
    
    
    
    wevote_df = pd.read_csv(wevote_file, sep = "\t", header = None)
    try:
        wevote_df.columns = ["read_id", "tool_count", "tools_can_classify", "agreed", "score", "junk", "tool_3", "tool_4", "assignment"]
    except ValueError:
        wevote_df.columns = ["read_id", "tool_count", "tools_can_classify", "agreed", "score", "junk", "tool_0", "tool_1", "tool_2", "tool_3", "tool_4", "assignment"]
    wevote_df = wevote_df[["read_id", "assignment"]]
    wevote_df.to_csv(output_file, sep = "\t", index = False, header = None)