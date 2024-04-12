import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
import sort_cluster_name_based_binned_v1 as sortc
import itertools

#-----------------------Relative percentage population for heatmap ----------------------------------------------------------------------------
dict_population_norm_bar = {"E2.EM.M2.SP":10,"E2.EM.M2.SP'":0.15,"E2.EM.M2'.SP":15,"E2.EM.M2'.SP'":0.6,"E2.EM'.M2.SP":15,"E2.EM'.M2'.SP":20,"E2.EM'.M2'.SP'":3,"E2.EM'.M2;.SP":3,"E2'.EM'.M2.SP":3,"E2'.EM'.M2'.SP":20}

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
# JSON file containing the infomation of all gamma  clusters under all possible classes
json_file_sorted_class = peptide+"_sorted_dict_classes_gamma_cluster.json"
with open(json_file_sorted_class,"r") as f:
    sorted_class_info = json.load(f)

prominent_categories = {}
for prmt_class in sorted_class_info:
    #e.g.  list_pmt :[[1,0,0,0],"orange",1] 
    list_prmt = [sorted_class_info[prmt_class]["weight"],sorted_class_info[prmt_class]["color"],sorted_class_info[prmt_class]["weight_factor"]]
    prominent_categories[prmt_class] = list_prmt


categories = prominent_categories
#------------------------------------------------
cutoff_num_occurences = 1
cutoff_num_mol_per_frame = 1

# Input for sorting and param list order in sorting
#sorted_seq=['0','L','M1','M2','H']
#param_list = [0,1]
#param_list = [0,1,2,3]
#sorted_seq=['0','L','M','H']

# sorted sequence for bin name ['0','L','M','H']
sorted_seq = list(dict_TAMCIS["binning"][param_list[0]].keys())
print(sorted_seq)


#sorted_seq = ["0","L","M","H"]
#param_list_sort_seq = [0,1,2,3]
param_list_sort_seq = [i for i in range(len(sorted_seq))]

# Gamma clusters information ------------------------------------------------------------
with open(cluster_json) as f:
    cluster_data = json.load(f)

cluster_key_list = list(cluster_data.keys())

list_cluster_cid_and_name = []

for ckey in cluster_key_list:
    time_list = list(cluster_data[ckey].keys())
    #print(ckey,cluster_data[ckey][time_list[0]]["name"],len(time_list),time_list)
    list_cluster_cid_and_name.append([ckey,cluster_data[ckey][time_list[0]]["name"]])

sorted_list_cluster_cid_and_name = sortc.sort_cluster_name_based_bin(list_cluster_cid_and_name,sorted_seq,param_list_sort_seq)

x = [] # store time 
y = [] # gamma cluster name
z = [] # number of molecules present at particular time in that cluster

sorted_list_timp_cluster_cid_and_name = []

dict_num_cluster_t = {}
for t in range(num_frame):
    t_str = str(float(t))
    num_cluster_per_time = 0
    for ckey_list in sorted_list_cluster_cid_and_name:
        cid = ckey_list[0]
        cname = ckey_list[1]
        cluster = cluster_data[cid]
        time_list = list(cluster.keys())
        num_occurences = len(time_list)
        list_t_count_non_stray = []
        check_atleast_cutoff_num_mol = False # 
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
                x.append(t) # store time 
                y.append(cname) # gamma cluster name
                z.append(num_non_stray) # number of molecules present at particular time in that cluster
                num_cluster_per_time = num_cluster_per_time + 1
            else:
                x.append(t)
                y.append(cname)
                z.append(float('nan')) # setting nan if no molecules present for getting the heatmap color white (means no data present)
        else: # check condition not satisfied
            x.append(t)
            y.append(cname)
            z.append(float('nan'))            
        if(t == 0 and check_atleast_cutoff_num_mol == True):
            sorted_list_timp_cluster_cid_and_name.append(ckey_list)
    dict_num_cluster_t.update({t:num_cluster_per_time})


# Create a DataFrame from the lists
# x : time 
# y : gamma cluster name
# z : number of molecules present at particular time in that cluster
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
# 
heatmap_data = heatmap_data.transpose()
heatmap_data = heatmap_data[y_tickslabel]

heatmap_data = heatmap_data.fillna(0) # filling nan with zero value

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
selected_categories_class_name = []

chk_colored_xticks = True

c = 0
for ckey in categories_clist.keys():
    list_of_clusters = categories_clist[ckey]
    print(ckey)
    for clr in list_of_clusters:
        nw_heatmap_data[clr] = heatmap_data[clr]
    c += len(list_of_clusters)
    selected_categories.append(c)
    selected_categories_class_name.append(ckey)

selected_categories.pop()
#selected_categories_class_name.append("AHBSP")


# selective_vmax can be False or the % population for normalizing % for comparison purpose   
def heatmap_class(classwise_heapmap_data,tag_png,selective_vmax):
    cmap = 'magma_r'
    fig, ax = plt.subplots(1, 1, figsize=(6, 2))  # Increase figsize if needed
    column_names = classwise_heapmap_data.columns.tolist()
    #im = sns.heatmap(classwise_heapmap_data,xticklabels=column_names,cmap=cmap,vmin=0,vmax=15, cbar=True, ax=ax)
    if(selective_vmax != False):
        im = sns.heatmap(classwise_heapmap_data,xticklabels=column_names,cmap=cmap,vmin=0,vmax=selective_vmax, cbar=True, ax=ax)
    else:
        im = sns.heatmap(classwise_heapmap_data,xticklabels=column_names,cmap=cmap, cbar=True, ax=ax)
    fontsize = 12
    fontsize_label = 10

    ax.invert_yaxis()
    yticks_gap = 100
    yticks_values = np.arange(0, num_frame+1, yticks_gap)
    ax.set_yticks(yticks_values)
    ax.set_yticklabels(yticks_values,fontsize=fontsize)
    ax.set_xlabel('$\Gamma$ Cluster',fontsize=fontsize_label)
    ax.set_ylabel('t(ns)',fontsize=fontsize_label)
    ax.set_xticklabels(column_names, rotation=90, ha='center',fontsize=fontsize)
    colorbar = ax.collections[0].colorbar
    colorbar.set_label('% Population', labelpad=2, fontsize=fontsize_label)  # Set the font size for colorbar label

    # Set the font size for colorbar tick labels
    colorbar.ax.yaxis.set_tick_params(labelsize=fontsize)
    plt.title(peptide.upper()+":"+tag_png)
    fig_name = peptide+"_CI_cluster_heatmap_time_evolution_"+tag+"_"+tag_png+"_v4.png"
    plt.savefig(fig_name, bbox_inches='tight', dpi=300)
    plt.show()



i_start = 0
gen_cnt = 0

print("start==========================")
for category in selected_categories:
    #ax.axvline(category, color='black', linewidth=0.5, linestyle='dashed')  # Adjust linewidth and linestyle as needed
    print(category,"+=============+")
    i_end = category
    classwise_heapmap_data = nw_heatmap_data.iloc[:,i_start:i_end]
    classwise_heapmap_data = 100*classwise_heapmap_data/system_size
    print(classwise_heapmap_data)
    png_tag = selected_categories_class_name[gen_cnt]
    if(png_tag in dict_population_norm_bar.keys()):
        selective_population_vmax = dict_population_norm_bar[png_tag]
    else:
        selective_population_vmax = False
    heatmap_class(classwise_heapmap_data,png_tag,selective_population_vmax)
    gen_cnt += 1
    i_start = i_end

i_end = -1
classwise_heapmap_data = nw_heatmap_data.iloc[:,i_start:i_end]
classwise_heapmap_data = 100*classwise_heapmap_data/system_size
print(classwise_heapmap_data)
png_tag = selected_categories_class_name[gen_cnt]
if(png_tag in dict_population_norm_bar.keys()):
    selective_population_vmax = dict_population_norm_bar[png_tag]
else:
    selective_population_vmax = False
print(png_tag,"=====================")
heatmap_class(classwise_heapmap_data,png_tag,selective_population_vmax)

