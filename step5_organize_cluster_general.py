import numpy as np
import json

def cluster_sort(cluster_dict,sorted_cluster_filename_json):
    list_key = cluster_dict.keys()
    cluster_len_total = []
    for key in list_key:
        x = cluster_dict[key]["time"]
        cluster_len_total.append(len(x))
    len_id = cluster_len_total.copy()

    sorted_list_key = [x for _,x in sorted(zip(len_id,list_key),reverse=True)]

    if('-1' in sorted_list_key):
        sorted_list_key.remove('-1')
    sorted_cluster_dict = {}
    for idx in range(len(sorted_list_key)):
        key = sorted_list_key[idx]
        sorted_key = str(idx)
        sorted_cluster_dict[sorted_key] = cluster_dict[key]
        sorted_cluster_dict[sorted_key]["length"] = len(sorted_cluster_dict[sorted_key]["segid"])
        #print(key,len(cluster_dict[key]["x"]),len(sorted_cluster_dict[sorted_key]["segid"]))
    if('-1' in cluster_dict):
        sorted_cluster_dict['-1'] = cluster_dict['-1']
        sorted_cluster_dict['-1']["length"] = len(sorted_cluster_dict['-1']["segid"])

    sorted_cluster_dict["number_cluster"] = len(sorted_list_key)
    with open(sorted_cluster_filename_json,'w+') as f:
        json.dump(sorted_cluster_dict,f,indent=8)

    return sorted_cluster_dict


def organize_cluster_json(cluster_filename,base_path):
    cluster = np.genfromtxt(base_path+"/"+cluster_filename,skip_header = 2,dtype=int)
    with open(base_path+"/"+cluster_filename) as f:
        cdata_linewise = f.readlines()

    # original file
    xyz_filename = cdata_linewise[0].split(' ')[1]
    #data = np.genfromtxt(base_path+"/"+xyz_filename,skip_header = 2)
    #data_str = np.genfromtxt(base_path+"/"+xyz_filename,skip_header = 2,dtype=str)
    #with open(base_path+"/"+xyz_filename) as f:
    #    data_linewise = f.readlines()
    data = np.genfromtxt(xyz_filename,skip_header = 2)
    data_str = np.genfromtxt(xyz_filename,skip_header = 2,dtype=str)
    with open(xyz_filename) as f:
        data_linewise = f.readlines()

    # Get the parameter list
    parameter_list = data_linewise[0].split(' ')
    print(parameter_list)

    labels = parameter_list[1:-1]
    print(labels)

    # Cluster key list
    param_list = ["segid"]
    param_list.extend(labels)
    param_list.append("time")
    print("parameter list :",labels,param_list)

    # name of the time independent sorted json file 
    sorted_cluster_filename_json = cluster_filename.replace("cluster_data","cluster_sorted_data")
    sorted_cluster_filename_json = base_path+"/"+sorted_cluster_filename_json.replace(".dat",".json")
    
    segid = data_str[:,0]
    x = data[:,1]
    
    c = cluster[:,2]
    data_c_idx = cluster[:,1]
    c_len = len(c)

    cluster_dict = {}
    for idx in range(len(data_c_idx)):
        key = str(c[idx])
        if(key  not in cluster_dict.keys()):
            cluster_dict[key] = {}
            for klist in param_list:
                cluster_dict[key][klist] = []
        for idx_klist in range(len(param_list)):
            klist = param_list[idx_klist]
            data_idx = data_c_idx[idx]
            if(idx_klist == 0):
                cluster_dict[key][klist].append(data_str[data_idx,idx_klist])
            else:
                cluster_dict[key][klist].append(data[data_idx,idx_klist])

    # Apply sorting on cluster dictionary with respect to cluster population 
    sorted_cluster_dict = cluster_sort(cluster_dict,sorted_cluster_filename_json)
    print("FILE CREATED :",sorted_cluster_filename_json)
    return sorted_cluster_filename_json

            
########################################################################################################################################
#  USER INPUT
#---------------------------------------------------------------------------------------------------------------------------------

input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

#--------------------------------------------------------------------------------------------------------------------
peptide = dict_TAMCIS["peptide"]
xyz_filename = dict_TAMCIS["xyz_raw_file"]

base_path = "."
cluster_filename = dict_TAMCIS["cluster_file"]
sorted_cluster_filename_json = organize_cluster_json(cluster_filename,base_path)

dict_TAMCIS["cluster_jsonfile"] = sorted_cluster_filename_json
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)

# ---------------------------------------------------------------------------------------------------------


    
