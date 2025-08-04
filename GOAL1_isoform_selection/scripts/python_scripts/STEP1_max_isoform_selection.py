"""
The purpose of this script is to identify the longest protein isoforms that
will be used for orthogrouping.

Proteins will be grouped by GeneID and the longest protein in that group will
be designated as the longest isoform. For keeping records, these longest isoforms
are written to disk by the value in the Protein Product column of the selected row.
"""

import os as os
import pandas as pd

#Specify working directory if necessary
os.chdir("/media/damian/Memorex USB/comp_mutualism")

file_list = os.listdir("./input_data/genome_metadata/")
txt_file_list = []
for item in file_list:
    if "txt" in item:
        txt_file_list.append(item)
    else:
        print "Not appending " + item
print txt_file_list
print len(txt_file_list)

#os.makedirs("./processed_data/longest_protein/")
for item in txt_file_list:
    print item
    print type(item)
    df = pd.read_csv("./input_data/genome_metadata/" + item, sep = "\t")
    max_index = df.groupby("GeneID")["Length"].idxmax()
    print item + " max_index: " + str(len(max_index))
    print item + " set(GeneID): "+ str(len(set(df["GeneID"])))
    max_index_df = pd.DataFrame(max_index)
    max_index_list = []
    for idx in max_index_df.Length:
        max_index_list.append(idx)
    max_df = df[df.index.isin(max_index_list)]
    longest_protein = []
    for product in max_df["Protein product"]:
        longest_protein.append(product)
    print item + " len(longest_protein): " + str(len(longest_protein))
    with open("./processed_data/longest_protein/" + item + "p", "w") as f:
        f.write("\n".join(longest_protein))