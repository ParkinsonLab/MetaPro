import pandas as pd
import sys
import os
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

#Sept 18, 2018:  This is not the greatest of scripts.  It's a series of DF manipulations, that produces what we need.  
#Doesn't have to be pretty, just readable.  We don't expect to change anything here.  It manipulates 3 files, and spits out
#an annotation of the enzymes we collected in the run

if __name__ == "__main__":
    pd.options.mode.chained_assignment = None  # default='warn' -> program makes too much noise.  this quiets it
    plt.ioff()
    
    ec_pathway_file = sys.argv[1]   #path->EC file  -> EC_pathway.txt
    rpkm_table_file = sys.argv[2]   #EC-> gene file -> RPKM_table.tsv
    pathway_superpathway_file = sys.argv[3] #path->superpath file -> pathway_to_superpathway.csv
    output_dir = sys.argv[4]
    if(not output_dir.endswith("/")):
        output_dir += "/"
    
    
    #imports the path->ec file
    path_df = pd.read_csv(ec_pathway_file, header=None, names=None, sep = '\t', skip_blank_lines = False)
    path_df.columns = (["path", "EC"])
    path_df["path"] = path_df["path"].apply(lambda x: x.split(":")[1])  
    path_df.index = path_df.groupby(["path"]).cumcount()
    temp_columns = path_df.index
    path_df = path_df.pivot(values='EC', columns = 'path')
    path_df = path_df.T
    path_df["EC_list"] = path_df[path_df.columns[:]].apply(lambda x: ",".join(x.dropna()), axis=1)
    path_df.drop(temp_columns, axis=1, inplace=True)
    
    
    #imports the path->superpath file
    super_df = pd.read_csv(pathway_superpathway_file, sep = ',', skip_blank_lines = False)
    super_df.drop(columns = ["Pathway"], inplace=True)
    
    #merge the path->EC and path->superpath dfs together, and put all the enzymes into categories
    super_df = super_df.join(path_df, how = 'left', on = 'Pathway ID')
    super_df.drop(columns = ["Pathway ID"], inplace = True)
    super_df.index = super_df.groupby("Superpathway").cumcount()
    super_df = super_df.pivot(values = "EC_list", columns = 'Superpathway')
    super_df = super_df.T
    super_temp_columns = super_df.columns
    super_df["combined"] = super_df[super_df.columns[:]].apply(lambda x: ",".join(x.dropna()), axis=1)
    super_df.drop(super_temp_columns, axis=1, inplace=True)
    super_enzyme_df = pd.DataFrame(super_df.combined.str.split(',').tolist(), index = super_df.index)
    superpath_list = list(super_enzyme_df.index)
    super_enzyme_df = super_enzyme_df.T
    
    #sort the ECs we've got and only retain the unique ones (for each superpath) 
    #pull each column from the DF (1 per superpath, so the loop isn't that bad), get their unique, and glue back into a new DF
    count = 0
    cleaned_enzyme_df = None
    for item in superpath_list:        
        new_df = pd.DataFrame(super_enzyme_df[item])
        new_df = pd.DataFrame(new_df[item].unique())
        new_name = item + "_decoupled.csv"
        if(count == 0):
            cleaned_enzyme_df = new_df
            cleaned_enzyme_df.columns = [item]
        else:
            cleaned_enzyme_df[item] = new_df
        count += 1
    super_enzyme_df = cleaned_enzyme_df
    super_enzyme_df = super_enzyme_df.T
    super_enzyme_df["count"] = super_enzyme_df.count(axis = 1) #get the total for each category (and thus total for all superpaths)
    #we have superpathway -> enzyme now, and they're unique per superpath.  gets us how many enzymes make up the complete set in a superpath
    
    #now we get the number of enzymes, annotated against the list to figure out what we have

    #This df now contains all genes -> enzymes, listed with their superpathway
    rpkm_df = pd.read_csv(rpkm_table_file, sep = '\t', skip_blank_lines = False)
    actual_read_df = rpkm_df[["GeneID", "EC#"]]
    actual_read_df["EC#"] = "ec:" + actual_read_df["EC#"]
    count = 0
    just_matched_enzymes_df = None
    #same logic as above, just on a different DF (pull each col, unique, glue to new df)
    for item in cleaned_enzyme_df.columns.values:      
        actual_read_df[item] = actual_read_df["EC#"].isin(cleaned_enzyme_df[item])
        new_df = pd.DataFrame(actual_read_df["EC#"].loc[actual_read_df["EC#"].isin(cleaned_enzyme_df[item])].unique())
        if(count == 0):
            just_matched_enzymes_df = new_df       
            just_matched_enzymes_df.columns = [item]
        else:
            
            just_matched_enzymes_df[item] = new_df
        count += 1
        
    just_matched_enzymes_df = just_matched_enzymes_df.T
    just_matched_enzymes_df["count"] = just_matched_enzymes_df.count(axis=1)
    just_matched_enzymes_df["total for category"] = super_enzyme_df["count"]
    just_matched_enzymes_df["% coverage"] = round((just_matched_enzymes_df["count"] / just_matched_enzymes_df["total for category"]) * 100, 2)
    just_matched_enzymes_df = just_matched_enzymes_df.T
    just_matched_enzymes_df.to_csv(output_dir + "EC_coverage.csv", mode = "w")
    #We've got our % coverage now.  
    
    #this takes the rpkm df and transforms it to get our heatmap, with help from cleaned_enzyme_df
    fixed_rpkm_df = rpkm_df
    fixed_rpkm_df["EC#"] = "ec:" + fixed_rpkm_df["EC#"]
    fixed_rpkm_df.drop(columns = ["Length", "RPKM", "#Reads", "Other"], inplace = True)  #removes the cols we don't need.  Also removes "Other"
    count = 0
    selected_heatmap_df = None
    for item in cleaned_enzyme_df.columns.values:
        new_df = pd.DataFrame(fixed_rpkm_df.loc[fixed_rpkm_df["EC#"].isin(cleaned_enzyme_df[item])])
        new_df["Superpathway"] = item
        if(count == 0):
            selected_heatmap_df = new_df
        else:
            selected_heatmap_df = pd.concat([selected_heatmap_df, new_df])
        count+= 1
    
    selected_heatmap_df = selected_heatmap_df.groupby(["Superpathway"]).sum()           #collapse the rows, grouped by the superpath`
    plt.subplots(figsize = (25,10))                                                     #set the image size    
    heatmap = sns.heatmap(selected_heatmap_df)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation = 40, ha = "right")     #make the labels pretty
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation = 40, ha = "right")
    heatmap.figure.savefig(output_dir + "enzyme_superpathway_heatmap.jpg")              #export it
    selected_heatmap_df.to_csv(output_dir + "enzyme_superpathway_heatmap.csv", mode="w")
    
    