#adds a C or a U depending on the taxa.  it apparently needs this label for Krona to make a table.
#

import os
import pandas as pd
import sys

def apply_record(x):
    if(x == 0):
        return "U"
    else:
        return "C"

if __name__ == "__main__":
    taxa_in = sys.argv[1]
    taxa_out = sys.argv[2]
    
    taxa_df = pd.read_csv(taxa_in, sep = "\t", header = None)
    taxa_df.columns = ["read", "taxa"]
    taxa_df["record"] = taxa_df["taxa"].apply(lambda x: apply_record(x))
    
    taxa_df = taxa_df[["record", "read", "taxa"]]
    print(taxa_df)
    
    taxa_df.to_csv(taxa_out, sep = "\t", header = False, index = False)