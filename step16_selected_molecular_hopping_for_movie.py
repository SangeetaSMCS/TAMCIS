import numpy as np
import json,os 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter
import matplotlib.pylab as pl


def CountFrequency(my_list):
 
    # Creating an empty dictionary
    freq = {}
    for item in my_list:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1
 
    sorted_freq = sorted(freq.items(), key=lambda x:x[1], reverse=True)
    converted_dict = dict(sorted_freq)
    return converted_dict
def convert_list_into_string(python_list):
    str_1 = ""
    num_elem = len(python_list)
    for i in range(num_elem):
        if(i != (num_elem-1)):
            str_1 = str_1 + python_list[i] + " "
        else:
            str_1 = str_1 + python_list[i]
    return str_1

def write_tcl_file(out, string_to_write,idx):
    print_string = 'set ' + 'cluster_list_'+str(idx) + " {" + string_to_write\
                                                + "}"
    print(print_string, file=out, end='\n')
    print_string = 'set sel_string ' + "\"segname " + string_to_write\
                                              + "\""
    print(print_string, file=out, end='\n')

input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
categories = dict_TAMCIS["prominent_class"]
len_promt_class = len(categories)

with open("VMD_color_prominent_class.json","r") as f:
    color_pclass_in_VMD = json.load(f)

###################################################################################################################################
#-----------------------------------------USER INPUT JSON -------------------------------------------------------------------------------
path = "/home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/"+peptide+"/pdb_file/500ns"
all_segname_json = path+"/segment_names_"+peptide+"_500ns_f1.json"
tag_fiber = "fiber1"


#-----------------------------------------------------------------------------------------------------------------------------


#################################################################################################################################
with open(all_segname_json) as f:
    all_segname_data = json.load(f)
all_segname_list = all_segname_data["peptide_segnames"]


# Read the CSV file into a DataFrame
csv_file = dict_TAMCIS["molecular_hopping_csv"]
df = pd.read_csv(csv_file)

drop_columns = []
for cele in df.columns:
    if('.' in cele):
        drop_columns.append(cele)
    elif(cele not in all_segname_list):
        drop_columns.append(cele)


#print(drop_columns)
df.drop(drop_columns, inplace=True, axis=1)

len_segname = len(df.columns)

# Select the last two rows
last_two_rows = df.iloc[-10:]

# Calculate the sum of the last two rows for each column
sum_last_two_rows = last_two_rows.sum()

# Sort the columns based on the sum of the last two rows
df = df[sum_last_two_rows.sort_values().index]



select_indices = list(range(len_segname))

df = df.iloc[:,select_indices]
print(df)

df_trans = df.transpose()
#fig, ax = plt.subplots()
'''
df = df.fillna(-1)
for index, row in df.iterrows():
    dict_prompt_class = {}
    for col_name, value in zip(df.columns, row):
        value_int = int(value)
        if(value_int not in dict_prompt_class.keys()):
            dict_prompt_class[value_int] = []
        dict_prompt_class[value_int].append(col_name)
    print(index,dict_prompt_class)
    frame_dir = "frame_"+str(index)
    if(os.path.isdir(frame_dir)==False):
        os.mkdir(frame_dir)
    tcl_file =frame_dir+"/"+peptide+"_category_segname_"+tag+".tcl"
    with open(tcl_file,"w+") as f:
        for idx in range(len_promt_class):
            c_idx = idx+1
            if(c_idx in dict_prompt_class.keys()):
                str_test = convert_list_into_string(dict_prompt_class[c_idx])
                write_tcl_file(f, str_test,idx)
            else:
                empty_list = []
                str_test = convert_list_into_string(empty_list)
                write_tcl_file(f, str_test,idx)
'''
# Clusterwise tcl file make
df = df.fillna(-1)

# row : segnames
# column : time
# Iter over each time and get segnames
for index, row in df.iterrows():
    print(":::",index,row)
    dict_prompt_class = {} # store segnames wrt the class at that time
    #print(index,row)
    for col_name, value in zip(df.columns, row):
        value_int = int(value)
        if(value_int not in dict_prompt_class.keys()):
            dict_prompt_class[value_int] = []
        dict_prompt_class[value_int].append(col_name)
        #print(col_name, value)
    print(index,dict_prompt_class)
    frame_dir = "frame_"+str(index)
    if(os.path.isdir(frame_dir)==False):
        os.mkdir(frame_dir)
    tcl_file =frame_dir+"/"+peptide+"_category_segname_"+tag+".tcl"
    with open(tcl_file,"w+") as f:
        for idx in range(len_promt_class):
            c_idx = idx+1
            if(c_idx in dict_prompt_class.keys()):
                str_test = convert_list_into_string(dict_prompt_class[c_idx])
                write_tcl_file(f, str_test,idx)
            else:
                empty_list = []
                str_test = convert_list_into_string(empty_list)
                write_tcl_file(f, str_test,idx)
    # ---------------------------------------------------------------
    for idx in range(len_promt_class):
        c_idx = idx+1
        tcl_file =frame_dir+"/"+peptide+"_category_segname_"+tag+"_"+str(c_idx)+".tcl"
        with open(tcl_file,"w+") as f:
            idx_vmd = "seg"
            if(c_idx in dict_prompt_class.keys()):
                str_test = convert_list_into_string(dict_prompt_class[c_idx])
                write_tcl_file(f, str_test,idx_vmd)
            else:
                empty_list = []
                str_test = convert_list_into_string(empty_list)
                write_tcl_file(f, str_test,idx_vmd)
            print(c_idx)
            print(color_pclass_in_VMD)
            VMD_cid = color_pclass_in_VMD[str(c_idx)]
            print("set color_id "+str(VMD_cid[0]),file=f)



