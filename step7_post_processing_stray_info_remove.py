
import json,os
from subprocess import list2cmdline
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import PercentFormatter
import math
import pandas as pd
import plotly.express as px
from pandas.plotting import parallel_coordinates
import matplotlib as mpl
import seaborn as sns
#%matplotlib inline


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



def alphanumeric_name_cluster(labels,parameter_bins,alpha_name):
    compt_MV_bin = {}
    for idx,oparam in enumerate(labels):
        idx_param = idx+1
        compt_MV_bin[oparam] = {}
        oparam_bin_value = parameter_bins[str(idx_param)]
        for idx_v,vbin in enumerate(oparam_bin_value):
            if(idx_v < len(oparam_bin_value)-1):
                compt_MV_bin[oparam][alpha_name[idx_v]] = [vbin,oparam_bin_value[idx_v+1]-1]
            else:
                compt_MV_bin[oparam][alpha_name[idx_v]] = [vbin,vbin+100]
    return compt_MV_bin


# BE VERY SPECIFIC ABOUT BINNING
def OP_human_readable_tag_v1(param,cp_min,cp_max,compt_MV_bin):
    '''compt_MV_bin = {'total_H_bond_contact':{"0":[0,0],"L":[1,2],"M":[3,4],"H":[5,30]},
                    'sc_wat_res1_contact':{"0":[0,0],"L":[1,9],"M":[10,19],"H":[20,60]},
                    'sc_wat_res2_contact':{"0":[0,0],"L":[1,9],"M":[10,19],"H":[20,60]},
                    'sc_wat_res3_contact':{"0":[0,0],"L":[1,9],"M":[10,19],"H":[20,60]},
                    }
    compt_MV_bin = {'pep_H_end_end_contact':{"0":[0,0],"L":[1,2],"M":[3,4],"H":[5,20]},
                    'pep_H_end_mid_contact':{"0":[0,0],"L":[1,2],"M":[3,4],"H":[5,20]},
                    'pep_H_mid_mid_contact':{"0":[0,0],"L":[1,2],"M":[3,4],"H":[5,20]},
                    'sidechain_pep_contact':{"0":[0,0],"L":[1,9],"M":[10,29],"H":[30,120]}
                    }
    
    compt_MV_bin = {'total_H_bond_contact':{"0":[0,0],"L":[1,3],"M1":[4,5],"M2":[6,7],"H":[8,12]},
                    'sidechain_pep_contact':{"0":[0,0],"L":[1,9],"M1":[10,24],"M2":[25,44],"H":[45,120]}
                    }'''
    current_param_bin = compt_MV_bin[param]
    OP_tag = ''
    for bin_key in current_param_bin:
        min_bin = current_param_bin[bin_key][0]
        max_bin = current_param_bin[bin_key][1]
        if(cp_min >=min_bin and cp_max <= max_bin):
            OP_tag = bin_key
    
    return OP_tag


def postprocess_remove_stray_molecules(filename,num_frame,parameter_bins,alphanumeric_name):
    with open(filename) as f:
        cluster_data = json.load(f)

    # cluster keys, eg. ['0','1',....]
    cluster_key_list = list(cluster_data.keys())
    cluster_key_list.remove("number_cluster")
    num_cluster = len(cluster_key_list)

    # parameter list : example: "total_H_bonds","sidechain_pep_contact", etc.
    parameter_list = list(cluster_data['0'].keys())
    labels = parameter_list[1:-2]
    
    print(parameter_list,labels,"-------------------------")
    
    compt_MV_bin = alphanumeric_name_cluster(labels,parameter_bins,alphanumeric_name)
    timewise_cluster_data = {}

    ## Collecting all the stray molecules and store framewise in a dictionary
    stray_json_file = peptide+"_all_cluster_of_size_1_1.json"
    stray_mol_dict = {}
    with open(stray_json_file) as f:
        stray_data = json.load(f)
    for t_stray in range(num_frame):
        stray_mol = stray_data[str(t_stray)]
        stray_mol_list = []
        for clr in stray_mol:
            stray_mol_list.append(*stray_mol[clr])
        stray_mol_dict.update({t_stray:stray_mol_list})
    
    ## Collecting timewise cluster segregation 
    for cidx in cluster_key_list:
        ctime = cluster_data[cidx]['time']
        csegid = cluster_data[cidx]['segid']
        clength = cluster_data[cidx]['length']
        cparam_dict = {}
        # getting the range of each OP value range for each classified cluster
        cluster_name = ""
        for param in labels:
            c_param = cluster_data[cidx][param]
            #print(cidx,param)
            cp_max = max(c_param)
            cp_min = min(c_param)
            OP_tag = OP_human_readable_tag_v1(param,cp_min,cp_max,compt_MV_bin)
            cluster_name = cluster_name + OP_tag
            cp_value = str(cp_min)+"-"+str(cp_max)
            cparam_dict.update({param:cp_value})
        timewise_cluster_data[cidx] = {}
        
        for idx,t in enumerate(ctime):
            if(t not in timewise_cluster_data[cidx].keys()):
                timewise_cluster_data[cidx][t]={}
                timewise_cluster_data[cidx][t]['segid']=[]
                timewise_cluster_data[cidx][t]['param'] = cparam_dict
                timewise_cluster_data[cidx][t]['name'] = cluster_name
                if(len(cluster_name) == 3):
                    print(cluster_name,cidx,cparam_dict)
                    print("ERROR in Cluster name assignment")
            timewise_cluster_data[cidx][t]['segid'].append(csegid[idx])
            

    
    for cidx in timewise_cluster_data.keys():
        if(cidx != '-1'):
            #print("CLUSTER ",cidx,"----------------------------------------")
            count_t = 0
            for t in timewise_cluster_data[cidx].keys():
                t_int = int(t)
                list_stray_mol = stray_mol_dict[t_int]
                #print(t_int)
                ct_segid = timewise_cluster_data[cidx][t]['segid']
                intersect_ct_stray  = list(set(ct_segid).intersection(set(list_stray_mol)))
                ct_segid_non_stray = list(set(ct_segid)-set(intersect_ct_stray))
                timewise_cluster_data[cidx][t]['non_stray'] = ct_segid_non_stray
                timewise_cluster_data[cidx][t]['stray'] = intersect_ct_stray
                #if(intersect_ct_stray != []):
                    #print("list :",ct_segid)
                    #print("non stray :",ct_segid_non_stray)
                    #print(cidx, t, len(ct_segid), len(intersect_ct_stray))
                count_t = count_t + 1
            #print("time count",count_t," length",cluster_data[cidx]['length'])
    output_json_file = peptide+"_cluster_alphanumeric_name_"+tag+".json"
    with open(output_json_file,"w+") as f:
        json.dump(timewise_cluster_data,f,indent=True)
    return output_json_file,stray_json_file,compt_MV_bin
    
    


################################## USER INPUT ###############################
#set a short tag for gamma cluster : to understand the orderparameter space
# Like for OP : 3 classified H-Bonds and sidechain-peptide : tag = chbsp 
# DON'T FORGET TO CHANGE this tag
tag = "chbsp"

# assumption that every order parameters has 4 bins
# alphanumeric name is the name given to each bin so that we acn understand the properties of gamma cluster from it name 
# If you want to increse bins then we have go back to step 3 repeat the process and also change this alphanumeric_name like for 5 bins say :
# The alphanumic name will be ['0','L','M1','M2','H']
alphanumeric_name = ['0','L','M','H']

################################################################################################


input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)






#--------------------------------------------------------------------------------------------------------------------
peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
cluster_filename = dict_TAMCIS["cluster_jsonfile"]
parameter_bins = dict_TAMCIS["manual_bin"]


output_json_file,stray_json_file,compt_MV_bin = postprocess_remove_stray_molecules(cluster_filename,num_frame,parameter_bins,alphanumeric_name)

#alphanumeric_name_cluster(labels,parameter_bins,alpha_name)
dict_TAMCIS["tag"] = tag
dict_TAMCIS["stray_json"] = stray_json_file
dict_TAMCIS["binning"] = compt_MV_bin
dict_TAMCIS["cluster_alphanum_file"] = output_json_file
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)




