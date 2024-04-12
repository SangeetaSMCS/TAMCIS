import numpy as np 
import json 


input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
output_file = dict_TAMCIS["cluster_alphanum_file"]
with open(input_file) as f:
    cluster_data = json.load(f)

pair_cid_name_list = []
cname_list = []
for c_id in cluster_data.keys():
    time_list = list(cluster_data[c_id])
    cname_list.append(cluster_data[c_id][time_list[0]]['name'])
    pair_cid_name_list.append([c_id,cluster_data[c_id][time_list[0]]['name']])

uniq_cname_list = list(np.unique(cname_list))
merge_cluster_dict = {}
for cname in uniq_cname_list:
    if(cname not in merge_cluster_dict.keys()):
        merge_cluster_dict[cname] = []
    for cid_name in pair_cid_name_list:
        if(cid_name[1]==cname):
            merge_cluster_dict[cname].append(cid_name[0])
#print(pair_cid_name_list)
#print(merge_cluster_dict)



for merge_key in merge_cluster_dict.keys():
    list_merge = merge_cluster_dict[merge_key]
    #print("merge list =",list_merge)
    # taken the 1st cluster id in the merge list
    init_comb_cidx = list_merge[0]
    # Extarct information from original cluster data
    init_comb_cluster = cluster_data[init_comb_cidx]
    init_comb_time_list = list(init_comb_cluster.keys())
    param_list = list(init_comb_cluster[init_comb_time_list[0]].keys())
    #print("list:",list(init_comb_cluster['325.0'].keys()))
    #print(len(init_comb_time_list))
    for comb_cid in list_merge[1:]:
        print("COMMON CLUSTER name found, trying to fix it ",comb_cid)
        comb_cluster = cluster_data[comb_cid]
        comb_time_list = list(cluster_data[comb_cid].keys())
        for c_time in comb_time_list:
            if(c_time in init_comb_time_list):
                for param_name in param_list:
                    init_param = init_comb_cluster[c_time][param_name]
                    c_param = comb_cluster[c_time][param_name]
                    if(param_name == 'param'):
                        for op in init_param.keys():
                            init_value_op = init_param[op].split("-")
                            c_value_op = c_param[op].split("-")
                            combine_value_op = []
                            combine_value_op.extend(init_value_op)
                            combine_value_op.extend(c_value_op)
                            update_value_op = str(min(combine_value_op))+"-"+str(max(combine_value_op))
                            cluster_data[init_comb_cidx][c_time][param_name][op] = update_value_op
                    elif(param_name != 'name'):
                        #print(param_name)
                        cluster_data[init_comb_cidx][c_time][param_name].extend(c_param)
                        #print(c_time,"EXtend : ",cluster_data[init_comb_cidx][c_time][param_name],c_param)
                    #else:
                    #    print(c_param,init_param)
            else:
                cluster_data[init_comb_cidx][c_time] = cluster_data[comb_cid][c_time]
                #print(c_time,"NOT there",cluster_data[init_comb_cidx][c_time])
                init_comb_time_list = list(init_comb_cluster.keys())
                #print(len(init_comb_time_list),"\n\n")
        del cluster_data[comb_cid]


with open(output_file,"w+") as f:
    json.dump(cluster_data,f,indent=True)

print("OUTPUT file :",output_file)