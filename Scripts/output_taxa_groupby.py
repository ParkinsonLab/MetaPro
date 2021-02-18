import pandas as pd
import sys
import os
#does a simple groupby and sums up the taxa found
if __name__ == "__main__":

    taxa_class_file = sys.argv[1]
    taxa_df = pd.read_csv(taxa_class_file, sep = "\t", error_bad_lines = False, quoting = 3)
    taxa_df.columns = ["class", "read", "taxa"]
    taxa_df["count"] = 1
    taxa_df = taxa_df.groupby("taxa", as_index = False).sum()
    #print(taxa_df)
    taxa_df.to_csv(sys.argv[2], sep = "\t", index = False)#, header = False)