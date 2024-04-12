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

def convert_list_into_string(python_list):
    str_1 = ""
    num_elem = len(python_list)
    for i in range(num_elem):
        if(i != (num_elem-1)):
            str_1 = str_1 + python_list[i] + " "
        else:
            str_1 = str_1 + python_list[i]
    return str_1

def write_tcl_file(out, string_to_write):
    print_string = 'set ' + 'segment_list ' + "{" + string_to_write\
                                                + "}"
    print(print_string, file=out, end='\n')
    print_string = 'set sel_string ' + "\"segname " + string_to_write\
                                              + "\""
    print(print_string, file=out, end='\n')

###########################################################################################
# Input files : CI cluster json and all segnames present in the system


input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
categories = dict_TAMCIS["prominent_class"]
all_segname_json = dict_TAMCIS["fiberpath_seganme_json"]
param_list = dict_TAMCIS["parameters_list"]





with open(input_file) as f:
    cluster_data = json.load(f)

with open(all_segname_json) as f:
    all_segname_data = json.load(f)
all_segname_list = all_segname_data["peptide_segnames"]

#-------------------------------------------------------------------------------------------
# determination of significant cluster based on cutoff which is each CI cluster must present atleast x time points of the trajectory , id cutoff_time=x  

cutoff_time = 1
significant_cluster_list = []
significant_cluster_idx_list = []
for cidx in cluster_data:
    if(cidx != '-1'):
        cluster = cluster_data[cidx]
        for t in range(num_frame):
            t_key_list = list(cluster.keys())
            t_str = str(float(t))
            if(t_str in t_key_list):
                t_cluster = cluster[t_str]
                len_non_stray = len(t_cluster['non_stray'])
                if(len_non_stray >= cutoff_time):
                    #print(t,t_cluster['name'],len_non_stray)
                    significant_cluster_list.append(t_cluster['name'])
                    significant_cluster_idx_list.append(cidx)
uniq_significant_cluster_list = list(np.unique(significant_cluster_list))
uniq_significant_cluster_list.append("others")
uniq_significant_cluster_idx_list = list(np.unique(significant_cluster_idx_list))
#uniq_significant_cluster_list.append("others")
print(uniq_significant_cluster_list,uniq_significant_cluster_idx_list)

#--------------------------------------------------------------------------------------------
# Sort according to the weight of CI cluster in accending order according to the param_list order 
first_param = dict_TAMCIS["parameters_list"][0]
sorted_seq=list(dict_TAMCIS["binning"][first_param].keys())
sign_cluster_list = []
weight_list = []
# getting the weight of gamma clusters in sign_cluster_list e.g. [['L000','cid',[1,0,0,1]]]
for cidx in uniq_significant_cluster_idx_list:
    imp_cluster = cluster_data[cidx]
    t_list = list(imp_cluster.keys())
    cluster_name = imp_cluster[t_list[0]]['name']
    len_cname = len(cluster_name)
    weight_cname = []
    for i in range(len_cname):
        char_pos = cluster_name[i]
        weight_cname.append(sorted_seq.index(char_pos))
    weight_list.append(weight_cname)
    sign_cluster_list.append([cluster_name,cidx,weight_cname])

#print(sign_cluster_list)
#order of parameter we want
#param_list = [0,1,2,3]

sorted_seq = list(dict_TAMCIS["binning"][param_list[0]].keys())
print(sorted_seq)


#sorted_seq = ["0","L","M","H"]
#param_list_sort_seq = [0,1,2,3]
param_list_sort_seq = [i for i in range(len(sorted_seq))]
sorted_weight_list = (sorted(weight_list,key=lambda x: [x[i] for i in param_list_sort_seq]))

sorted_uniq_significant_cluster_idx_list = []
for wgt in sorted_weight_list:
    for u_cidx in sign_cluster_list:
        if(wgt == u_cidx[2]):
            #sorted_uniq_significant_cluster_idx_list.append([u_cidx[1],wgt])
            sorted_uniq_significant_cluster_idx_list.append([u_cidx[1],u_cidx[0]])
#print(sorted_uniq_significant_cluster_idx_list)

sorted_uniq_significant_cluster_name_list = ['others']
sorted_uniq_significant_cluster_idx_rank_list = []

for s_imp_id in sorted_uniq_significant_cluster_idx_list:
    sorted_uniq_significant_cluster_name_list.append(s_imp_id[1])

#print("This is the list:",sorted_uniq_significant_cluster_idx_list)

list_sorted_uniq_sign_cluster_name = [i[1] for i in sorted_uniq_significant_cluster_idx_list]
#sorted_uniq_significant_cluster_idx_list = ['10']

#--------------------------------------------------------------------------------------------------------------------------
def catergory_weight(list_CI_clusters):
    dict_caterg_wgt_CI_cluster = {}
    for col_name in list_CI_clusters:
        col_id = [] # list of weight 4D cluster: [1,0,0,1] etc..
        for i in range(len(col_name)):
            str_i = col_name[i]
            if(str_i == "0"):
                col_id.append(0)
            else:
                col_id.append(1)
        dict_caterg_wgt_CI_cluster[col_name] = col_id
    return dict_caterg_wgt_CI_cluster


#categories = {"$N_{e2e}^{ppH}$":[[1,0,0,0],"orange",1],"$N_{sp}^C$":[[0,0,0,1],"red",2],"$N_{e2e}^{ppH}-N_{sc}^C$":[[1,0,0,1],"green",3],"$N_{m2m}^{ppH}-N_{sc}^C$":[[0,0,1,1],"gray",4],"$N_{e2e}^{ppH}-N_{m2m}^{ppH}$\n$-N_{sc}^C$":[[1,0,1,1],"blue",5],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}$\n$-N_{sc}^C$":[[1,1,0,1],"cyan",6],"$N_{e2e}^{ppH}-N_{e2m}^{ppH}$\n$-N_{m2m}^{ppH}-N_{sc}^C$":[[1,1,1,1],"purple",7]}

# user defined superlative catergory provide 
#categories = {"E2":[[1,0,0,0],"orange",1],"SP":[[0,0,0,1],"red",2],"E2SP":[[1,0,0,1],"green",3],"M2SP":[[0,0,1,1],"gray",4],"E2M2SP":[[1,0,1,1],"blue",5],"E2EMSP":[[1,1,0,1],"cyan",6],"AHBSP":[[1,1,1,1],"purple",7]}


# Extracting the weight of the CI clusters according to that of superlative catergory class [Must have correlation with caterories we want]
# Output: dictionary with CI cluster as key and value as weight list
dict_caterg_wgt_CI_cluster = catergory_weight(list_sorted_uniq_sign_cluster_name)

#Time list at which we want to examine the hopping characteristic 
time_list = list(range(0,num_frame,1))


#all_segname_list = all_segname_list[:50]
i_seg_gap = len(all_segname_list)
c = 0
avg_hop_count = 0
for i_seg in range(0,len(all_segname_list),i_seg_gap):
    segname_list = all_segname_list[i_seg:i_seg+i_seg_gap]
    segname_sign_cluster_info = {}
    segname_sign_cluster_info_fluctuation = {}
    for segname in segname_list:
        segname_sign_cluster_info[segname] = []
        segname_sign_cluster_info_fluctuation[segname] = {}
        count_hop = 0
        for t in time_list:
            t_str = str(float(t))
            # -1 is the initial flag for the segname not present in the  cluster at that time
            belong_in = np.nan
            for sub_sorted_uniq_significant_cluster_idx_list in sorted_uniq_significant_cluster_idx_list:
                #print(sub_sorted_uniq_significant_cluster_idx_list)
                imp_cidx,imp_cname = sub_sorted_uniq_significant_cluster_idx_list
                imp_cluster = cluster_data[imp_cidx]            
                t_key_list = list(imp_cluster.keys())
                cater_wgt = dict_caterg_wgt_CI_cluster[imp_cname]
                if(t_str in t_key_list):
                    t_imp_cluster = imp_cluster[t_str]['segid']
                    imp_cluster_name = imp_cluster[t_str]['name']
                    if(segname in t_imp_cluster):
                        #segname_sign_cluster_info[segname].append(imp_cluster_name)
                        #segname_sign_cluster_info[segname].append(int(imp_cidx))
                        # key = class name, value : [[1,0,0,0],"color",weight]
                        for key,value in categories.items():
                            # class weight == cluster weight
                            if(cater_wgt == value[0]):
                                belong_in = value[2]

                                # Calculation for fluctuations
                                if(key not in segname_sign_cluster_info_fluctuation[segname].keys()):
                                    segname_sign_cluster_info_fluctuation[segname][key]={}
                                if(imp_cname not in segname_sign_cluster_info_fluctuation[segname][key].keys()):
                                    segname_sign_cluster_info_fluctuation[segname][key][imp_cname] = 0
                                segname_sign_cluster_info_fluctuation[segname][key][imp_cname] = segname_sign_cluster_info_fluctuation[segname][key][imp_cname] + 1
                        if(belong_in == np.nan):
                            belong_in = 0
                        
                        # Calculation of per molecule hopping freq
                        if(t == 0):
                            prev_belong_in = belong_in
                        else:
                            if(prev_belong_in != belong_in):
                                count_hop = count_hop + 1
                                prev_belong_in = belong_in
                        #break
            segname_sign_cluster_info[segname].append(belong_in)
        #segname_sign_cluster_info_fluctuation[segname]["hopping_freq"] = count_hop
        avg_hop_count += count_hop

    df = pd.DataFrame(segname_sign_cluster_info)
    #df = df.iloc[:, df.iloc[-1].argsort()]

    # Save the selected DataFrame to a CSV file
    csv_file = peptide+"_temporal_hopping_"+tag+".csv"
    df.to_csv(csv_file, index=False)

    json_file_fluctuation = peptide+"_CI_cluster_hop_count_"+tag+".json"
    with open(json_file_fluctuation,"w+") as f:
        json.dump(segname_sign_cluster_info_fluctuation,f,indent=True)


'''
# average time spend for each molecules in 7 clusters 
count_dict_hop_average = {}
for seg_key in segname_sign_cluster_info_fluctuation.keys():
    for pclass_key in segname_sign_cluster_info_fluctuation[seg_key].keys():
        pclass_value = segname_sign_cluster_info_fluctuation[seg_key][pclass_key]
        #print(pclass_value)
        if(pclass_key not in count_dict_hop_average.keys()):
            count_dict_hop_average[pclass_key] = [0,0]
        for CI_key in pclass_value:
            count_dict_hop_average[pclass_key][0] += pclass_value[CI_key]
            count_dict_hop_average[pclass_key][1] += 1


for key in count_dict_hop_average.keys():
    avg_hop = count_dict_hop_average[key][0]/count_dict_hop_average[key][1]
    print(key,count_dict_hop_average[key],avg_hop)


print("average hopping : ",avg_hop_count/len(all_segname_list))  
'''
dict_TAMCIS["molecular_hopping_csv"] = csv_file
# json file storing molecular hopping count 
dict_TAMCIS["prominent"] = json_file_fluctuation
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)
