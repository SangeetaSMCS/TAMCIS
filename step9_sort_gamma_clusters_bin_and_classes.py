import numpy as np
import json
import itertools
import matplotlib.pyplot as plt

def int_to_hex(value, max_value):
    """
    Converts an integer value to a hexadecimal color code.

    Args:
        value: The integer value (0-max_value).
        max_value: The maximum allowed integer value (default: 32).

    Returns:
        A string representing the hexadecimal color code (e.g., '#FFFFFF').

    Raises:
        ValueError: If the provided value is outside the valid range (0-max_value).
    """

    if not (0 <= value <= max_value):
        raise ValueError(f"Value must be between 0 and {max_value}")

    # Convert the integer value to a 3-digit hex string with leading zeros
    hex_code = '#' + format(value, '06x')
    print(hex_code)

    # Create the final hexadecimal color code by prepending '#'
    return hex_code

def int_to_hex_color(number):
    # Ensure the input number is within the valid range for RGB color values
    number = max(0, min(number*100000, 16777215))  # 16777215 is the maximum value for RGB

    # Convert the number to hexadecimal format
    hex_color = '#' + format(number, '06x')  # '06x' ensures the output is 6 characters long

    return hex_color


#--------------------------------------------------------------------------------------------------------------------------
# Assigning gamma cluster weigth for class 
def catergory_weight_v1(list_CI_clusters_with_other_info,idx_CI_name):
    dict_caterg_wgt_CI_cluster = {}
    for col_name_info in list_CI_clusters_with_other_info:
        col_id = [] # list of weight 4D cluster: [1,0,0,1] etc.
        col_name = col_name_info[idx_CI_name]
        for i in range(len(col_name)):
            str_i = col_name[i]
            if(str_i == "0"):
                col_id.append(0)
            else:
                col_id.append(1)
        col_name_info_update = []
        col_name_info_update.extend(col_name_info)
        col_name_info_update.append(col_id)
        dict_caterg_wgt_CI_cluster[col_name] = col_name_info_update
    return dict_caterg_wgt_CI_cluster


def sort_cluster_name_based_bin_v1(list_cluster_cid_and_name,sorted_seq,param_list):
    '''Logic : any further seperation for M1 M2,.. should be mentioned in the sorted_seg list and subindex M1,M2 ...much start from 1, 2 , not from 0 as 
    0 is seperate entity here which means absolute zero value'''
    #sorted_seq=['0','L','M1','M2','H']
    list_cluster_with_wgt = []
    weight_list = []
    for cname_cid in list_cluster_cid_and_name:
        cluster_name = cname_cid[1] 
        len_cname = len(cluster_name)
        cidx = cname_cid[0]
        weight_cname = []
        for i in range(len_cname):
            char_pos = cluster_name[i]
            if(char_pos.isdigit() == False or char_pos == '0'):
                if(i != len_cname-1 ):
                    char_pos1 = cluster_name[i+1]
                    if(char_pos1.isdigit() and char_pos1 != '0'):
                        char_pos = char_pos+char_pos1
                        i = i + 1
                weight_cname.append(sorted_seq.index(char_pos))
        weight_list.append(weight_cname)
        list_cluster_with_wgt.append([cidx,cluster_name,weight_cname])

    #print(list_cluster_with_wgt)
    #order of parameter we want
    #param_list = [0,1]
    sorted_weight_list = (sorted(weight_list,key=lambda x: [x[i] for i in param_list]))

    sorted_list_cluster_cid_and_name = []
    for wgt in sorted_weight_list:
        for u_cidx in list_cluster_with_wgt:
            if(wgt == u_cidx[2]):
                sorted_list_cluster_cid_and_name.append([u_cidx[0],u_cidx[1]])
    return sorted_list_cluster_cid_and_name

#############################################################################################
def get_position_value_in_combination(gamma_dim,gamma_name_compliment,wgt):
    prm_class_name = ""
    for i in range(gamma_dim):
        if(wgt[i] == 1):
            prm_class_name += gamma_name[i]
        else:
            prm_class_name += gamma_name_compliment[i]
    return prm_class_name





################### USER INPUT ##########################################################################################################################################################
#order of parameter we want
#gamma_name = ["A","B","C","D"]
#gamma_name_compliment = ["A'","B'","C'","D'"]

# name of the each interaction in set theory notation
gamma_name = ["E2.","EM.","M2.","SP"]
# name of the compliment of each interaction  in set theory notation
gamma_name_compliment = ["E2'.","EM'.","M2'.","SP'"]

# Color assign for visualization of specific important classes
dict_color_for_class = {"E2.EM.M2.SP":"purple","E2.EM.M2.SP'":"gray","E2.EM.M2'.SP":"cyan","E2.EM.M2'.SP'":"gray","E2.EM'.M2.SP":"blue","E2.EM'.M2'.SP":"green","E2.EM'.M2'.SP'":"yellow","E2'.EM'.M2.SP":"gray","E2'.EM'.M2'.SP":"red"}


# unimportant classes color
unimp_cls_color = "gray"
#######################################################################################################################################################################








input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
param_list = dict_TAMCIS["parameters_list"]


prominent_categories = {}
gamma_dim = len(gamma_name)
count = 1

# assigning weight representation of each class e.g. : class name E2.EM.M2'.SP will have weight [1,1,0,1]
for CI_comb_cnt in range(gamma_dim):
    # order of weight mentained as [0,0,0,0], [0,0,0,1], [0,0,1,0]....
    for x in itertools.combinations(range(gamma_dim), CI_comb_cnt ) :
        wgt = [ 1 if i in x else 0 for i in range(gamma_dim) ]
        prm_class_name = get_position_value_in_combination(gamma_dim,gamma_name_compliment,wgt)
        # Assigning colors 
        if(prm_class_name in dict_color_for_class.keys()):
            prm_class = dict_color_for_class[prm_class_name]
        else:
            prm_class = unimp_cls_color
        prominent_categories[prm_class_name] = [wgt,prm_class,count] 
        count += 1 # count is just id assigned in ascending order of weight 

# For [1,1,1,1] class 
wgt = [1]*gamma_dim
prm_class_name = get_position_value_in_combination(gamma_dim,gamma_name_compliment,wgt)
if(prm_class_name in dict_color_for_class.keys()):
    prm_class = dict_color_for_class[prm_class_name]
else:
    prm_class = unimp_cls_color
prominent_categories[prm_class_name] = [wgt,prm_class,count]

print(prominent_categories,"\nNumber of sets = ",len(prominent_categories),"\nNumber of proper set + 1 (analitically) = ",2**len(gamma_name))




# For gamma clusters and assigning gamma clusters to respective classes
with open(input_file) as f:
    cluster_data = json.load(f)

# list of all gamma clusters present
cluster_key_list = list(cluster_data.keys())
list_cluster_cid_and_name = []

for ckey in cluster_key_list:
    time_list = list(cluster_data[ckey].keys())
    #print(ckey,cluster_data[ckey][time_list[0]]["name"],len(time_list),time_list)
    list_cluster_cid_and_name.append([ckey,cluster_data[ckey][time_list[0]]["name"]])


idx_CI_name = 1
dict_cluster_cid_name_wgt = catergory_weight_v1(list_cluster_cid_and_name,idx_CI_name)

#print(dict_cluster_cid_name_wgt)



dict_gamma_clr_under_prmt_classes = {}
for pmt_class_key in prominent_categories.keys():
    pmt_class = prominent_categories[pmt_class_key]
    pmt_class_wgt = pmt_class[0]
    for gamma_clr_key in dict_cluster_cid_name_wgt.keys():
        gamma_clr = dict_cluster_cid_name_wgt[gamma_clr_key]
        gamma_clr_wgt = gamma_clr[2]
        if(pmt_class_wgt == gamma_clr_wgt):
            if(pmt_class_key not in dict_gamma_clr_under_prmt_classes.keys()):
                dict_gamma_clr_under_prmt_classes[pmt_class_key] = []
            dict_gamma_clr_under_prmt_classes[pmt_class_key].append(gamma_clr[:2])

#print(dict_gamma_clr_under_prmt_classes)

sorted_seq = list(dict_TAMCIS["binning"][param_list[0]].keys())
print(sorted_seq)


#sorted_seq = ["0","L","M","H"]
#param_list_sort_seq = [0,1,2,3]
param_list_sort_seq = [i for i in range(len(sorted_seq))]


for gamma_class_key in dict_gamma_clr_under_prmt_classes.keys():
    gamma_class = dict_gamma_clr_under_prmt_classes[gamma_class_key]
    sorted_list_cluster_cid_and_name = sort_cluster_name_based_bin_v1(gamma_class,sorted_seq,param_list_sort_seq)
    #print(gamma_class)
    #print(sorted_list_cluster_cid_and_name)
    dict_gamma_clr_under_prmt_classes[gamma_class_key] = {}
    dict_gamma_clr_under_prmt_classes[gamma_class_key]["sorted_gamma_cluster"] = sorted_list_cluster_cid_and_name
    dict_gamma_clr_under_prmt_classes[gamma_class_key]["weight"] = prominent_categories[gamma_class_key][0]
    dict_gamma_clr_under_prmt_classes[gamma_class_key]["color"] = prominent_categories[gamma_class_key][1]
    dict_gamma_clr_under_prmt_classes[gamma_class_key]["weight_factor"] = prominent_categories[gamma_class_key][2]


# eliminating stray class e.g. [0,0,0,0]
stray_class_name = ""
for elem in gamma_name_compliment:
    stray_class_name += elem
del dict_gamma_clr_under_prmt_classes[stray_class_name]
print(dict_gamma_clr_under_prmt_classes)



# Output file for sorted json file for classes and gamma clusters 
outjson_filename = peptide+"_sorted_dict_classes_gamma_cluster.json"
with open(outjson_filename,"w+") as f:
    json.dump(dict_gamma_clr_under_prmt_classes,f,indent=4)


dict_TAMCIS["class_name"] = gamma_name
dict_TAMCIS["class_compliment_name"] = gamma_name_compliment
dict_TAMCIS["sorted_classes_json"] =  outjson_filename
# Update the input_TAMCIS.json
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)




