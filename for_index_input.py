

import numpy as np
import pandas as pd
import os
import sys

'''
this program reads sciMET oligos csv file containing scimet tags
this extracts sciMET name and tags to a new file
this then copies the file into yet another file that is formatted
for input to validating_edit_metric_tags.py

run this file as: python3 for_index_input.py path/to/csv/file \
    path/where/to/save
'''
# print(sys.argv[1].split("/"))
file_path_split = sys.argv[1].split("/")
file_path = "/" + "/".join(file_path_split[1:len(file_path_split) - 1]) + "/"
# print(file_path)

csv_file = sys.argv[1]
df = pd.read_csv(csv_file, sep=",")
# print(df)
pd.set_option("display.max_rows", 400) # this sets up the number of rows displayed

df = df[["Name", "Index"]]
# print(df)
df = df.dropna()
# print(df)
# print(np.arange(135, 151))
df = df.drop(np.arange(135, 151))
# print(df)

extract_file = file_path + "existing_indexes.txt"

## be careful with the path given -> change it according to your system
with open(extract_file, 'w') as f:
    f.write(df.to_string(header = True, index = False))

output_file = sys.argv[2]
with open(extract_file, "r") as rf:
    with open(output_file, "w") as wf:
        for line in rf:
            line = line.strip()
            if line.startswith("#"):
                wf.write(line + "\n")
            elif line.startswith("Name"):
                line_split = line.split()
                new_line = ":".join(line_split)
                wf.write("[" + new_line + "] \n")
            else:
                line = line.strip()
                line_split = line.split()
                ## try and accept accounts for blank lines here
                try:
                    line_split[1] = line_split[1].upper()
                except:
                    pass
                new_line = ":".join(line_split)
                wf.write(new_line + "\n")
