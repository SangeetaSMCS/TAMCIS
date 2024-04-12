import numpy as np
import json
import os

def extract_from_tcl(tcl_filename):
    data=np.genfromtxt(tcl_filename,dtype=str)
    segname=data[2:]
    segname[0]=segname[0].replace('{','')
    segname[-1]=segname[-1].replace('}','')
    return segname

def membership_size_bining(list_member_size,dict_membership_binned,cluster_id,frame_id,cluster_seglist):
    #cluster_size_bin = [1,3,6,11,16,21,31,51,101]
    cluster_size_bin = [1,2,6,11,21,51,101,151]

    for population in list_member_size :
        freq_chk_flag = 0
        for idx1 in range(len(cluster_size_bin)-1):
            lower_bound = cluster_size_bin[idx1]
            upper_bound = cluster_size_bin[idx1+1]
            dict_key = str(lower_bound)+"_"+str(upper_bound-1)
            # Initiallization
            if(dict_key not in dict_membership_binned.keys()):
                dict_membership_binned[dict_key] = {}
            
            if(frame_id not in dict_membership_binned[dict_key].keys()):
                dict_membership_binned[dict_key][frame_id]={}


            if(population>= lower_bound and population < upper_bound and freq_chk_flag ==0):
                dict_key = str(lower_bound)+"_"+str(upper_bound-1)
                dict_membership_binned[dict_key][frame_id].update({cluster_id:cluster_seglist.tolist()})
                freq_chk_flag = 1

        dict_key = str(cluster_size_bin[-1])+"_above"
        if(dict_key not in dict_membership_binned.keys()):
                dict_membership_binned[dict_key] = {}
            
        if(frame_id not in dict_membership_binned[dict_key].keys()):
            dict_membership_binned[dict_key][frame_id]={}
        
        if(population >= upper_bound and freq_chk_flag ==0):
            dict_key = str(cluster_size_bin[-1])+"_above"
            dict_membership_binned[dict_key][frame_id].update({cluster_id:cluster_seglist.tolist()})
            freq_chk_flag = 1

cluster_size_bin = [1,2,5,10,20,50,100,150]
dict_cluster_size_binned = {}
num_cluster = 600

################################################ PROVIDE the num_frame and peptide name #############################################
num_frame = 500
peptide = "lll"
##########################################################################################################

dict_membership_binned = {}
for frame_id in range(num_frame):
    path_frame = "frame_"+str(frame_id)
    for cluster_id in range(num_cluster):
        cluster_file = "cluster_"+str(cluster_id)+".tcl"
        tcl_filename = path_frame+"/"+cluster_file
        if(os.path.isfile(tcl_filename) == True):
            cluster_seglist = extract_from_tcl(tcl_filename)
            cluster_len = len(cluster_seglist)
            list_member_size = [cluster_len]
            membership_size_bining(list_member_size,dict_membership_binned,cluster_id,frame_id,cluster_seglist)


for key in dict_membership_binned.keys():
    print(key,dict_membership_binned[key])
    dict_get = dict_membership_binned[key]
    with open(peptide+"_all_cluster_of_size_"+key+".json","w+") as f:
        json.dump(dict_get,f,indent=4)


