import numpy as np
import json
import general_code as gc

''' PROVIDE THE MANUAL BINING
 IMPORTANT NOTE: EVENIF NOT BINNING REQUIRE, THEN ALSO YOU HAVE DO THIS STEP WITH EMPLY "list_idx rescale" AND "list_bins_idxwise" and ENLIST all the idx's in "idx_non_rescale_list" '''

'''CONDITION'''
## item >= bin_start and item < bin_end ###########
# EDIT ON THIS DICTIONARY
list_bins_idxwise = {1:[0,1,3,5],2:[0,1,3,5],3:[0,1,3,5],4:[0,1,10,30]}

input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)





#--------------------------------------------------------------------------------------------------------------------
peptide = dict_TAMCIS["peptide"]
xyz_filename = dict_TAMCIS["xyz_raw_file"]
# index of bining parameter sc : 2 and pep_water : 3 according to xyz file 
list_idx_rescale = list(list_bins_idxwise.keys())
#list_idx_rescale = {1,2}
# bining
#list_bins_idxwise = {4:[0,1,10,20,30,40], 5:[0,20,40,60,80,100]}
#list_bins_idxwise = {4:[0,1,10,20,30,40]}
#list_bins_idxwise = {4:[0,1,5,10,15,20,30,40]}
#list_bins_idxwise = {2:[0,1,10,20],3:[0,1,10,20],4:[0,1,10,20]}

#list_bins_idxwise = {1:[0,1,4,6,8], 2:[0,1,10,25,45]}
#list_bins_idxwise = {4:[0,1,10,20,30,40]}
# indexes of other parameters
idx_non_rescale_list = {}

xyz_renorm = "_BINNED.xyz"
xyz_file_renormalized = xyz_filename.replace(".xyz",xyz_renorm)

#------------------------------------------------------------------------------------------------------------------------

data = np.genfromtxt(xyz_filename,skip_header = 2)
data_str = np.genfromtxt(xyz_filename,skip_header = 2,dtype=str)
with open(xyz_filename) as f:
    data_linewise = f.readlines()

# Get the parameter list
#parameter_list = data_linewise[0].split(' ')
#print(parameter_list)

parameter_list = dict_TAMCIS["parameters_list"]
print(parameter_list)
# number parameters present in the input xyz file excluding time parameter
num_parameters = len(parameter_list)

ext_parameter_list = [peptide]
for elem in parameter_list:
    ext_parameter_list.append(elem)
ext_parameter_list.append("time")


dict_transformed_param = {}
# Looping over the rescaling indices and different rescaling index has different binning
for idx_rescale in list_bins_idxwise.keys():
    # repective index bin
    print(ext_parameter_list[idx_rescale],"------------")
    non_unif_bin = list_bins_idxwise[idx_rescale]
    len_bin = len(non_unif_bin)

    # rescale parameter
    rescale_param = data[:,idx_rescale]
    rescale_param_max = max(rescale_param)

    #transformed parameter according to the binning
    transformed_param = []
    for item in rescale_param:
        flag = 0
        for i in range(len_bin-1):
            bin_start = non_unif_bin[i]
            bin_end = non_unif_bin[i+1]
            if(item >= bin_start and item < bin_end):
                transformed_param.append(i)
                flag = 1
                print(item, " in the bin [",bin_start,",",bin_end,"] and idx = ",i )
                break
        if(flag == 0):
            transformed_param.append(len_bin-1)
            bin_start = non_unif_bin[-1]
            bin_end = "..."
            print(item, " in the bin [",bin_start,",",bin_end,"] and idx = ",len_bin-1)
        #Storing the transformed param in the dict
        dict_transformed_param[idx_rescale] = transformed_param


#for idx in range(len(transformed_param)):
#    print(rescale_param[idx],"---",transformed_param[idx])


parameters = []
for i in range(num_parameters):
    parameters.append([])
segname_list = data_str[:,0]
time = data[:,-1]
# ==========================================================================================
for idx in range(num_parameters):
    if(idx+1 in idx_non_rescale_list):
        parameters[idx]=data[:,idx+1]
    elif(idx+1 in list_idx_rescale):
        parameters[idx]=dict_transformed_param[idx+1]
#print(parameters)

first_info = "Original_filename "+xyz_filename+" \n"+data_linewise[0].replace("\n","")
gc.make_xyz_file_general(num_parameters,parameters,time,segname_list,xyz_file_renormalized,first_info) 
print("Modified file : ",xyz_file_renormalized)

dict_TAMCIS["manual_bin"] = list_bins_idxwise
dict_TAMCIS["xyz_binned_file"] = xyz_file_renormalized
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)





