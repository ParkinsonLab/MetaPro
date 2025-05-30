import os
import sys
import pandas as pd
import numpy as np

# what this code does:
# takes in 2 pair fastq files, splits them up into 3: matching IDs (in both), and what doesn't match
# to scale this, break up the inputs, and build a reduction. 
# we expect this code to be called multiple times

def filter_for_orphans(pair_0_path_i, pair_1_path_i, orphans_path_i, pair_0_path_o, pair_1_path_o, unique_path_o):
    pre_df_0 = pd.read_csv(pair_0_path_i, header=None, names=[None], sep = r'\n', engine="python", skip_blank_lines = False, quoting=3)
    pre_df_1 = pd.read_csv(pair_1_path_i, header=None, names=[None], sep = r'\n', engine="python",  skip_blank_lines = False, quoting=3)
    df_0 = pd.DataFrame(pre_df_0.values.reshape(int(len(pre_df_0)/4), 4))
    df_1 = pd.DataFrame(pre_df_1.values.reshape(int(len(pre_df_1)/4), 4))
    
    df_0.columns = ["ID", "seq", "junk", "quality"]
    df_1.columns = ["ID", "seq", "junk", "quality"]
    #df_0["ID"] = df_0["ID"].apply(lambda x: x.split(" ")[0]) #There may be a space in the ID, left in some types of FASTQs.  
    #df_1["ID"] = df_1["ID"].apply(lambda x: x.split(" ")[0]) #we remove them to that the ID-matching will work
    common = df_0.merge(df_1, on=["ID"])
    
    #stuff that belongs go here
    #df_0[df_0.ID.isin(common.ID)].to_csv(pair_0_path_o, sep = '\n', mode = 'w+', header = False, index = False, quoting = 3)
    #df_1[df_1.ID.isin(common.ID)].to_csv(pair_1_path_o, sep = '\n', mode = 'w+', header = False, index = False, quoting = 3)
    df_0_yes = df_0[df_0.ID.isin(common.ID)]
    df_1_yes = df_1[df_1.ID.isin(common.ID)]

    np.savetxt(pair_0_path_o, df_0_yes.values, delimiter = "\n", fmt='%s')
    np.savetxt(pair_1_path_o, df_1_yes.values, delimiter = "\n", fmt='%s')
    
    #stuff that doesn't belong go to another pile
    # First dataframe (write mode)
    data_0 = df_0[~df_0.ID.isin(common.ID)]
    np.savetxt(unique_path_o, data_0.values, delimiter='\n', fmt='%s')

    # Second dataframe (append mode)
    data_1 = df_1[~df_1.ID.isin(common.ID)]
    with open(unique_path_o, 'a') as f:
        np.savetxt(f, data_1.values, delimiter='\n', fmt='%s')
    
    #There's some situations where there's no orphans generated from the previous steps.
    # feb 19, 2025: this code should work if it's opening the same file it's appending to.  
    if(os.path.exists(orphans_path_i)):
        orphans_i_file = pd.read_csv(orphans_path_i, header=None, names=[None], sep = r'\n', engine="python",  skip_blank_lines = False, quoting = 3)
        orphans_df = pd.DataFrame(orphans_i_file.values.reshape(int(len(orphans_i_file)/4), 4))
        orphans_df.columns = ["ID", "seq", "junk", "quality"]
        np.savetxt(unique_path_o, orphans_df.values.flatten(), fmt='%s')
        #orphans_df.to_csv(unique_path_o, sep = '\n', mode = 'a', header=False, index=False, quoting = 3)

if __name__ == "__main__":
    if(len(sys.argv) < 7):
        print("Too few input arguements.  Not filtering for orphans")
    elif(len(sys.argv) > 7):
        print("Too many input arguments.  Not filtering for orphans")
    else:
        pair_0_path_i = sys.argv[1]
        pair_1_path_i = sys.argv[2]
        orphans_path_i = sys.argv[3]
        pair_0_path_o = sys.argv[4]
        pair_1_path_o = sys.argv[5]
        unique_path_o = sys.argv[6]
        print("pair 0 in:", pair_0_path_i)
        print("pair 1 in:", pair_1_path_i)
        print("pair 0 out:", pair_0_path_o)
        print("pair 1 out:", pair_1_path_o)
        filter_for_orphans(pair_0_path_i, pair_1_path_i, orphans_path_i, pair_0_path_o, pair_1_path_o, unique_path_o)
        
