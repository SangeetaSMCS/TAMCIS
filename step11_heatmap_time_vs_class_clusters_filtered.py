import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
import sort_cluster_name_based_binned_v1 as sortc

# ------------------USER INPUT-----------------------------------------------------------------------------------------------
percent_population = 15
cutoff_num_occurences = 10
cutoff_num_mol_per_frame = 10


#---------------------------------------------------------------------------------------------------
input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"] 
cluster_json = dict_TAMCIS["cluster_alphanum_file"]
system_size = dict_TAMCIS["system_size"]
param_list = dict_TAMCIS["parameters_list"]


# Provide prominent class which is simple combination of order parameter with color code 
#prominent_categories = {"E2":[[1,0,0,0],"orange",1],"SP":[[0,0,0,1],"red",2],"E2SP":[[1,0,0,1],"green",3],"M2SP":[[0,0,1,1],"gray",4],"E2M2SP":[[1,0,1,1],"blue",5],"E2EMSP":[[1,1,0,1],"cyan",6],"AHBSP":[[1,1,1,1],"purple",7]}
json_file_sorted_class = peptide+"_sorted_dict_classes_gamma_cluster.json"
with open(json_file_sorted_class,"r") as f:
    sorted_class_info = json.load(f)

prominent_categories = {}
for prmt_class in sorted_class_info:
    list_prmt = [sorted_class_info[prmt_class]["weight"],sorted_class_info[prmt_class]["color"],sorted_class_info[prmt_class]["weight_factor"]]
    prominent_categories[prmt_class] = list_prmt


#categories = {"$N_{e2e}^{ppH}$":[[1,0,0,0],"orange"],"$N_{sp}^C$":[[0,0,0,1],"red"],"$N_{e2e}^{ppH}-N_{sc}^C$":[[1,0,0,1],"green"],"$N_{m2m}^{ppH}-N_{sc}^C$":[[0,0,1,1],"gray"],"$N_{e2e}^{ppH}-N_{m2m}^{ppH}$\n$-N_{sc}^C$":[[1,0,1,1],"blue"],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}$\n$-N_{sc}^C$":[[1,1,0,1],"cyan"],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}$\n$-N_{m2m}^{ppH}-N_{sc}^C$":[[1,1,1,1],"purple"]}

categories = prominent_categories


# Input for sorting and param list order in sorting
#sorted_seq=['0','L','M1','M2','H']
#param_list = [0,1]


sorted_seq = list(dict_TAMCIS["binning"][param_list[0]].keys())
print(sorted_seq)


#sorted_seq = ["0","L","M","H"]
#param_list_sort_seq = [0,1,2,3]
param_list_sort_seq = [i for i in range(len(sorted_seq))]

# ------------------------------------------------------------
with open(cluster_json) as f:
    cluster_data = json.load(f)

cluster_key_list = list(cluster_data.keys())

list_cluster_cid_and_name = []

for ckey in cluster_key_list:
    time_list = list(cluster_data[ckey].keys())
    #print(ckey,cluster_data[ckey][time_list[0]]["name"],len(time_list),time_list)
    list_cluster_cid_and_name.append([ckey,cluster_data[ckey][time_list[0]]["name"]])

sorted_list_cluster_cid_and_name = sortc.sort_cluster_name_based_bin(list_cluster_cid_and_name,sorted_seq,param_list_sort_seq)

x = []
y = []
z = []

sorted_list_timp_cluster_cid_and_name = []

dict_num_cluster_t = {}
for t in range(num_frame):
    t_str = str(float(t))
    num_cluster_per_time = 0
    for ckey_list in  sorted_list_cluster_cid_and_name:
        cid = ckey_list[0]
        cname = ckey_list[1]
        cluster = cluster_data[cid]
        time_list = list(cluster.keys())
        num_occurences = len(time_list)
        list_t_count_non_stray = []
        check_atleast_cutoff_num_mol = False
        if(num_occurences >= cutoff_num_occurences):
            for ct in time_list:
                cnt_non_stray = len(cluster[ct]["non_stray"])
                if(cnt_non_stray >= cutoff_num_mol_per_frame):
                    check_atleast_cutoff_num_mol = True
                    break
        

        #if(num_occurences >= cutoff_num_occurences):
        if(t_str in time_list and check_atleast_cutoff_num_mol == True):
            num_non_stray = len(cluster[t_str]["non_stray"])
            if(num_non_stray != 0):
                x.append(t) 
                y.append(cname)
                z.append(num_non_stray)
                num_cluster_per_time = num_cluster_per_time + 1
            else:
                x.append(t)
                y.append(cname)
                z.append(float('nan'))
        else:
            x.append(t)
            y.append(cname)
            z.append(float('nan'))            
        if(t == 0 and check_atleast_cutoff_num_mol == True):
            sorted_list_timp_cluster_cid_and_name.append(ckey_list)
    dict_num_cluster_t.update({t:num_cluster_per_time})


# Create a DataFrame from the lists
df = pd.DataFrame({'x': x, 'y': y, 'z': z})

y_tickslabel = [c[1] for c in sorted_list_timp_cluster_cid_and_name]

# Pivot the DataFrame to get a suitable format for sns.heatmap
#heatmap_data = df.pivot('y', 'x', 'z')
heatmap_data = df.pivot(index='y', columns='x', values='z')
'''
for index, row in heatmap_data.iterrows():
    print(row.name)
    for i in row:
        if(i > 0):
            print(i,end=',')
    print()
'''
heatmap_data = heatmap_data.transpose()
heatmap_data = heatmap_data[y_tickslabel]

heatmap_data = heatmap_data.fillna(0)
#categories = {"e2e-m2m":[1,0,1,1],"e2e-e2m":[1,1,0,1],"e2e-e2m-m2m":[1,1,1,1],"e2e-sc":[1,0,0,1],"$N_{sp}^C":[0,0,0,1],"e2e":[1,0,0,0],"m2m-sc":[0,0,1,1]}

#categories = {"$\mathbf{N_{e2e}^{ppH}}$":[1,0,0,0],"$\mathbf{N_{sp}^C}$":[0,0,0,1],"$\mathbf{N_{e2e}^{ppH}-N_{sc}^C}$":[1,0,0,1],"$\mathbf{N_{m2m}^{ppH}-N_{sc}^C}$":[0,0,1,1],"$\mathbf{N_{e2e}^{ppH}-N_{m2m}^{ppH}}$\n$\mathbf{-N_{sc}^C}$":[1,0,1,1],"$\mathbf{N_{e2e}^{ppH}-N_{e2m}^{ppH}}$\n$\mathbf{-N_{sc}^C}$":[1,1,0,1],"$\mathbf{N_{e2e}^{ppH}-N_{e2m}^{ppH}}$\n$\mathbf{-N_{m2m}^{ppH}-N_{sc}^C}$":[1,1,1,1]}

order_categories_list = list(categories.keys())
categories_clist = {}
categories_color_list = {}
#print(heatmap_data.columns)
for col_name in heatmap_data.columns:
    col_id = []
    for i in range(len(col_name)):
        str_i = col_name[i]
        if(str_i == "0"):
            col_id.append(0)
        else:
            col_id.append(1)
    print(col_name,col_id)
    for key,value in categories.items():
        if(col_id == value[0]):
            if(key not in categories_clist.keys()):
                categories_clist[key] = []
            categories_clist[key].append(col_name)
            categories_color_list[col_name]=value[1]
            

print(categories_clist)



nw_heatmap_data = pd.DataFrame()
selected_categories = []

chk_colored_xticks = True


c = 0
for ckey in categories_clist.keys():
    list_of_clusters = categories_clist[ckey]
    for clr in list_of_clusters:
        nw_heatmap_data[clr] = heatmap_data[clr]
    c += len(list_of_clusters)
    selected_categories.append(c)

selected_categories.pop()


cmap = 'magma_r'
#cmap = "viridis_r"
column_names = nw_heatmap_data.columns.tolist()
nw_heatmap_data = 100*nw_heatmap_data/system_size
fig, ax = plt.subplots(1, 1, figsize=(6, 2))  # Increase figsize if needed
im = sns.heatmap(nw_heatmap_data,xticklabels=column_names,  cmap=cmap,vmin=0,vmax=percent_population, cbar=True, ax=ax)
#im = sns.heatmap(nw_heatmap_data,xticklabels=column_names,  cmap=cmap, cbar=True, ax=ax)
fontsize = 12
fontsize_label = 10

ax.invert_yaxis()
yticks_gap = 200
yticks_values = np.arange(0, num_frame+1, yticks_gap)
ax.set_yticks(yticks_values)
ax.set_yticklabels(yticks_values,fontsize=fontsize)
ax.set_xlabel('$\Gamma$ Cluster',fontsize=fontsize_label)
ax.set_ylabel('t(ns)',fontsize=fontsize_label)


ax.set_xticklabels(column_names, rotation=90, ha='center',fontsize=fontsize)

if(chk_colored_xticks == True):
    # Loop through each tick label and set its color
    for i, tick_label in enumerate(ax.get_xticklabels()):
        xtick_color = categories_color_list[tick_label.get_text()]
        tick_label.set_color(xtick_color)


colorbar = ax.collections[0].colorbar
colorbar.set_label('% Population (per t)', labelpad=2, fontsize=fontsize_label)  # Set the font size for colorbar label

# Set the font size for colorbar tick labels
colorbar.ax.yaxis.set_tick_params(labelsize=fontsize)



# Add vertical partition lines for selected categories
for category in selected_categories:
    ax.axvline(category, color='black', linewidth=0.5, linestyle='dashed')  # Adjust linewidth and linestyle as needed


fig_name = peptide+"_CI_cluster_heatmap_time_evolution_"+tag+"_v2.png"
plt.savefig(fig_name, bbox_inches='tight', dpi=300)
plt.show()



# Reinitialization of weight for classes
imp_prominent_categories = {}
set_imp_gamma_cluster_filtered = set(y_tickslabel)
cnt = 1
for prmt_class in sorted_class_info:
    list_gamma_clusters = [gclr[1] for gclr in sorted_class_info[prmt_class]["sorted_gamma_cluster"]]
    intersection_list = list(set_imp_gamma_cluster_filtered.intersection(set(list_gamma_clusters)))
    print(intersection_list)
    if(len(intersection_list)>0):
        list_prmt = [sorted_class_info[prmt_class]["weight"],sorted_class_info[prmt_class]["color"],cnt]
        imp_prominent_categories[prmt_class] = list_prmt
        cnt += 1


dict_TAMCIS["prominent_class"] = imp_prominent_categories
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)


