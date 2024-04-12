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
all_segname_json = dict_TAMCIS["fiberpath_seganme_json"]


# Read the CSV file into a DataFrame
csv_file = dict_TAMCIS["molecular_hopping_csv"]
df = pd.read_csv(csv_file)

drop_columns = []
for cele in df.columns:
    if('.' in cele):
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

'''
# Get the last row values
last_row_values = df.iloc[-1]

# Sort the columns based on the last row values
df = df[sorted(df.columns, key=lambda col: last_row_values[col])]
'''

select_indices = list(range(len_segname))
#select_indices = list(range(235))
#select_indices = list(range(5,10))+list(range(15,25))+list(range(115,150))+list(range(212,216)) + list(range(218,220)) + list(range(225,232)) #iii
#select_indices = list(range(225,232))  
#select_indices = list(range(0,8))+list(range(42,50))+list(range(200,235))
#select_indices = list(range(0,8))+list(range(42,50))+list(range(135,143))+list(range(350,375))+list(range(485,499)) #lll
#select_indices = list(range(0,235))
#select_indices = list(range(20))+list(range(290,300))+list(range(390,400))+list(range(480,499))
df = df.iloc[:,select_indices]
print(df)
len_lab = len(categories)
df_trans = df.transpose()
#fig, ax = plt.subplots()
fig = plt.figure(figsize=(7,3.2))
ax = fig.add_subplot(111)


#color_list = ['pink', 'violet','purple', 'blue','green', 'yellow','orange','red','brown']
#color_list = [categories[k][1] for k in categories.keys()]
'''
color_list = []
for k in categories.keys():
    print(k,categories[k][1])
    color_list.append(categories[k][1])

col_dict = {}
for i in range(len(color_list)):
    col_dict.update({i:color_list[i]})
'''
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

fig_name = peptide+"_temporal_hopping_"+tag+".png"
plt.savefig(fig_name, bbox_inches='tight', dpi=300)
#plt.close()
plt.show()

'''
column_names = df.columns.tolist()
dict_chosen_mol = {}
dict_chosen_mol["peptide_segnames"] = column_names

json_file= peptide+"_chosen_mol_fiber.json"
with open(json_file,"w+") as f:
    json.dump(dict_chosen_mol,f)
print(len(column_names))




# Select two rows
df = df.fillna(0)
selected_rows = df.iloc[-1,:]  # Selects rows 1 and 2 (0-based index)
row_0 = selected_rows.tolist()


peptide = "lll"
#categories = {"$N_{e2e}^{ppH}$":[[1,0,0,0],"orange",1],"$N_{sp}^C$":[[0,0,0,1],"red",2],"$N_{e2e}^{ppH}-N_{sc}^C$":[[1,0,0,1],"green",3],"$N_{m2m}^{ppH}-N_{sc}^C$":[[0,0,1,1],"gray",4],"$N_{e2e}^{ppH}-N_{m2m}^{ppH}-N_{sc}^C$":[[1,0,1,1],"blue",5],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}-N_{sc}^C$":[[1,1,0,1],"cyan",6],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}-N_{m2m}^{ppH}-N_{sc}^C$":[[1,1,1,1],"purple",7]}
categories = {"E2":[[1,0,0,0],"orange",1],"SP":[[0,0,0,1],"red",2],"E2SP":[[1,0,0,1],"green",3],"M2SP":[[0,0,1,1],"gray",4],"E2M2SP":[[1,0,1,1],"blue",5],"E2EMSP":[[1,1,0,1],"cyan",6],"AHBSP":[[1,1,1,1],"purple",7]}


# Read the CSV file into a DataFrame
csv_file = peptide+"_temporal_hopping.csv"
df = pd.read_csv(csv_file)
print(df)

# Select the last two rows
last_two_rows = df.iloc[-10:]

# Calculate the sum of the last two rows for each column
sum_last_two_rows = last_two_rows.sum()

# Sort the columns based on the sum of the last two rows
df = df[sum_last_two_rows.sort_values().index]

selected_rows = df.iloc[-1,:]  # Selects rows 1 and 2 (0-based index)
row_1 = selected_rows.tolist()

dict_cnt_r0 = CountFrequency(row_0)
dict_cnt_r1 = CountFrequency(row_1)

print(dict_cnt_r0)

x = []
y = []
for i in range(1,8):
    i_flt = float(i)
    if(i_flt in dict_cnt_r0.keys()):
        value = dict_cnt_r0[i_flt]
        x.append(i)
        y.append(100.*value/235)
    else:
        x.append(i)
        y.append(0/235)

plt.plot(x,y,marker='o')

x1 = []
y1 = []
for i in range(1,8):
    i_flt = float(i)
    if(i_flt in dict_cnt_r1.keys()):
        value = dict_cnt_r1[i_flt]
        x1.append(i)
        y1.append(100.*value/235)
    else:
        x1.append(i)
        y1.append(0/235)

plt.plot(x1,y1,marker='s')

plt.show()

cater = []
bar_colors = []
for ctg,value in categories.items():
    cater.append(ctg)
    bar_colors.append(value[1])

# Generate different colors for each bar

# Create a bar plot with different colors for each bar
#plt.bar(categories, values, color=bar_colors)

X_axis = np.arange(len(cater))
# Create a bar plot with two sets of data
plt.bar(X_axis - 0.2, y, label='III', color=bar_colors, width=0.4,edgecolor="white",hatch='///')
plt.bar(X_axis + 0.2, y1, label='LLL', color=bar_colors, width=0.4,edgecolor="white",hatch='.')




plt.xticks(X_axis,cater,rotation=90,fontsize = 20)
plt.yticks(fontsize = 15)

# Add labels and title
plt.ylabel('% frequency',fontsize = 12)
#plt.ylabel('Y-axis')
#plt.title('Bar Plot with Two Sets of Data')

# Display the legend
plt.legend()

# Display the plot
#plt.savefig("comparism_lll_iii_fiber.png",bbox_inches='tight', dpi=300)
plt.show()

'''
'''
# Count occurrences of similar values in each column
value_counts = selected_rows.apply(lambda x: x.value_counts())

# Plot the results as a line plot
value_counts.T.plot(kind='line', marker='o')

# Display the plot
plt.xlabel('Columns')
plt.ylabel('Count')
plt.title('Count of Similar Values in Selected Rows')
plt.show()
'''

