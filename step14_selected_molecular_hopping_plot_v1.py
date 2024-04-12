import numpy as np
import json 
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


input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
categories = dict_TAMCIS["prominent_class"]


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


print(len_segname)



#-----------------------------USER CHOICE------------------------------------------------------
#select_indices = list(range(len_segname))
#select_indices = list(range(235))
#select_indices = list(range(5,10))+list(range(15,25))+list(range(115,150))+list(range(212,216)) + list(range(218,220)) + list(range(225,232)) #iii
#select_indices = list(range(225,232))  
#select_indices = list(range(0,8))+list(range(42,50))+list(range(200,235))
select_indices = list(range(0,8))+list(range(42,50))+list(range(135,150))+list(range(250,280))+list(range(430,len_segname)) #lll
#select_indices = list(range(0,235))
#select_indices = list(range(20))+list(range(290,300))+list(range(390,400))+list(range(480,499))
#------------------------------------------------------------------------------------------------


df = df.iloc[:,select_indices]
print(df)
len_lab = 9
df_trans = df.transpose()
#fig, ax = plt.subplots()
fig = plt.figure(figsize=(7,3.2))
ax = fig.add_subplot(111)


#color_list = ['pink', 'violet','purple', 'blue','green', 'yellow','orange','red','brown']

col_dict = {}
for k in categories.keys():
    col_dict.update({categories[k][2]:categories[k][1]})

cmap = ListedColormap([col_dict[x] for x in col_dict.keys()])
# Let's also define the description of each category : 1 (blue) is Sea; 2 (red) is burnt, etc... Order should be respected here ! Or using another dict maybe could help.
labels = np.array([str(k) for k in categories.keys()])
len_lab = len(labels)

#res1 = sns.heatmap(df_trans,cbar_kws={'label': 'Subordinate $\Gamma$ Clusters','orientation': 'vertical'},xticklabels=50,yticklabels=2,fmt='' ,cmap=cmap, ax=ax)
res1 = sns.heatmap(df_trans,xticklabels=100,yticklabels=1,fmt='' ,cmap=cmap, ax=ax)


font_size  = 20
font_size_lab = 12
# Set fontsize for heatmap ticks
ax.set_xticklabels(ax.get_xticklabels(), fontsize=font_size)
#ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)


# Set y-axis tick labels to column indices at a gap of 5
tick_indices = range(0, len(df.columns), 10)
ax.set_yticks(tick_indices)
ax.set_yticklabels(tick_indices, fontsize=font_size)

ax.set_xlabel('t(ns)', fontsize=font_size_lab)
ax.set_ylabel('Molecular index', fontsize=font_size_lab)

cbar = res1.collections[0].colorbar
cgap = 1.*(len_lab-1)/len_lab
cstart = 1. +(cgap/2.)
cbar_label_pos = [(i*cgap + cstart) for i in range(len_lab)]
print(cbar_label_pos)
cbar.set_ticks(cbar_label_pos)
cbar.set_ticklabels(labels, fontsize=font_size_lab)

fig_name = peptide+"_temporal_hopping_for_"+tag_fiber+"_"+tag+"_v1.png"
plt.savefig(fig_name, bbox_inches='tight', dpi=300)
#plt.close()
plt.show()



