import numpy as np
import json
import matplotlib.pyplot as plt
#%matplotlib inline
from mpl_toolkits.mplot3d import Axes3D
import general_code as gc





# Extract interaction from each molecule for a particular frame  
def calculate_interaction_for_list_of_mol(mol_list_cluster,base_interaction_path,frame_id,json_interaction_base,interaction_key):
    list_num_contact = []
    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+".json"
        data = gc.read_json(interaction_jsonfile)
        num_contact = data[frame_id][interaction_key]
        list_num_contact.append(num_contact)
    
    return list_num_contact

# Extract interaction from each molecule for a RANGE OF frames  
def calculate_interaction_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id,json_interaction_base,interaction_key):
    list_num_contact = []
    list_mol_time = []
    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+".json"
        #print(interaction_jsonfile)
        data = gc.read_json(interaction_jsonfile)
        for frame_id in list_frame_id:
            #print(molid,frame_id,interaction_key)
            num_contact = data[frame_id][interaction_key]
            list_num_contact.append(num_contact)
            list_mol_time.append([molid,frame_id])    
    return list_num_contact,list_mol_time

#  Extract interaction from each molecule for a RANGE OF frames for SIZE OF MEMBERSHIP OF CLUSTER
def calculate_member_cluster_size_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id,json_interaction_base):
    list_num_contact = []
    list_mol_time = []
    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+".json"
        #print(interaction_jsonfile)
        data = gc.read_json(interaction_jsonfile)
        for frame_id in list_frame_id:
            #print(molid,frame_id,interaction_key)
            num_contact = data["frame"][frame_id]
            list_num_contact.append(num_contact)
            list_mol_time.append([molid,frame_id])
    
    return list_num_contact,list_mol_time

# Extract interaction from each molecule for a RANGE OF frames specially for SALT BRIDGE CALCULATION 
def calculate_electrostatic_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id,json_interaction_base):
    list_num_contact = []
    list_mol_time = []
    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+".json"
        data = gc.read_json(interaction_jsonfile)
        for frame_id in list_frame_id:
            #print(molid,frame_id,interaction_key)
            num_contact = len(data[frame_id]["N_O_pair_atleast_one_within_4A"])+len(data[frame_id]["O_N_pair_atleast_one_within_4A"])
            list_num_contact.append(num_contact)
            list_mol_time.append([molid,frame_id])
    
    return list_num_contact,list_mol_time

# Extract interaction from each molecule for a RANDE OF frames  for CLASSIFIED H-bonding
def calculate_classification_pep_H_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4):

    list_mol_time = []
    pep_H_end_end_contact = []
    pep_H_end_middle_contact = []
    pep_H_middle_middle_contact = []
    total_H_bond = []
    
    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+"_classified.json"
        print(interaction_jsonfile)
        data = gc.read_json(interaction_jsonfile)
        for frame_id in list_frame_id:
            num_contact = data[frame_id][interaction_key1]
            pep_H_end_end_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key2]
            pep_H_end_middle_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key3]
            pep_H_middle_middle_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key4]
            total_H_bond.append(num_contact)
            list_mol_time.append([molid,frame_id])

    return pep_H_end_end_contact,pep_H_end_middle_contact,pep_H_middle_middle_contact,total_H_bond,list_mol_time

def calculate_classified_sc_wat_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4):

    list_mol_time = []
    sc_wat_res1_contact = []
    sc_wat_res2_contact = []
    sc_wat_res3_contact = []
    total_sc_wat_contact = []

    for molid in mol_list_cluster:
        #base_interaction_path = base_interaction_json+"/pep-pep-contact"
        interaction_jsonfile = base_interaction_path+"/"+json_interaction_base+molid+".json"
        print(interaction_jsonfile)
        data = gc.read_json(interaction_jsonfile)
        for frame_id in list_frame_id:
            num_contact = data[frame_id][interaction_key1]["num_contact_atoms"]
            sc_wat_res1_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key2]["num_contact_atoms"]
            sc_wat_res2_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key3]["num_contact_atoms"]
            sc_wat_res3_contact.append(num_contact)
            num_contact = data[frame_id][interaction_key4]
            total_sc_wat_contact.append(num_contact)
            list_mol_time.append([molid,frame_id])

    return sc_wat_res1_contact,sc_wat_res2_contact,sc_wat_res3_contact,total_sc_wat_contact,list_mol_time


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# data extraction for list of order parameter
def composite_OP_extraction(peptide,base_interaction_json,fiber_seganme_json,frame_end,parameters_list):

    data = gc.read_json(fiber_seganme_json)
    mol_list_cluster = data["peptide_segnames"]  

    #mol_list_cluster = mol_list_cluster[0:50]
    len_seglist = len(mol_list_cluster)
 
    # INCLUSION OF PARAMETER SPACE NEEDS INCLUSTION OF THE PARAMETER NAME HERE AND INCREASE THE NUMBER OF PARAMETERS
    # ==============================================================================================================================
    num_parameters = len(parameters_list)
    
    parameters = []

    for i in range(num_parameters):
        parameters.append([])

    # collecting frame id list
    list_frame_id_z3 = []
    list_frame_id_z4 = []
    time1 = []
    for f_idx in range(0,frame_end,1):
        frame_id = str(f_idx).zfill(4)
        list_frame_id_z4.append(frame_id)
        frame_id = str(f_idx).zfill(3)
        list_frame_id_z3.append(frame_id)
        time1.append(f_idx)
    
    time = []
    for count in range(len_seglist):
        time.extend(time1)
    
    clusterwise_interaction_dict = {}

    if(peptide == "fff"):
        list_frame_id_z = list_frame_id_z4
    else:
        list_frame_id_z = list_frame_id_z3

    # INCLUDE THE PARAMETER SPACE HERE
    # =======================================================================================================================================================================================================================================================================================================================================================================================
    # pep H : classifications : e-e, m-m , e-m
    param_name = "total_H_bond_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/pep-pep-H-bond/raw_data_all_frames/angle_filtered_data"
        json_interaction_base = "pep_pep_angle_filterd_data_"
        interaction_key1 = "num_end_to_end"
        interaction_key2 = "num_end_to_middle"
        interaction_key3 = "num_middle_to_middle"
        interaction_key4 = "num_H_bond"
        clusterwise_interaction_dict["pep_H_end_end_contact"],clusterwise_interaction_dict["pep_H_end_middle_contact"],clusterwise_interaction_dict["pep_H_middle_middle_contact"],clusterwise_interaction_dict["total_H_bond_contact"],list_mol_time = calculate_classification_pep_H_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id_z,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4)
    
    # pep H : classifications : e-e, m-m , e-m
    param_name = "pep_H_end_end_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/pep-pep-H-bond/raw_data_all_frames/angle_filtered_data"
        json_interaction_base = "pep_pep_angle_filterd_data_"
        interaction_key1 = "num_end_to_end"
        interaction_key2 = "num_end_to_middle"
        interaction_key3 = "num_middle_to_middle"
        interaction_key4 = "num_H_bond"
        clusterwise_interaction_dict["pep_H_end_end_contact"],clusterwise_interaction_dict["pep_H_end_middle_contact"],clusterwise_interaction_dict["pep_H_middle_middle_contact"],clusterwise_interaction_dict["total_H_bond_contact"],list_mol_time = calculate_classification_pep_H_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id_z,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4)
        
    #sidechain pep contact : list of numbers for each cluster
    param_name = "sidechain_pep_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/sidechain-pep-contact/raw_data_all_frames"
        json_interaction_base = "side_chain_protein_contacts_"
        interaction_key = "num_side_chain_protein_contacts"
        clusterwise_interaction_dict["sidechain_pep_contact"],list_mol_time = calculate_interaction_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id_z3,json_interaction_base,interaction_key)
    
    #sidechain water contact : residuewise contact 
    param_name = "total_sc_wat_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/sidechain-water-contact/raw_data_all_frames"
        json_interaction_base = "side_chain_water_contacts_"
        interaction_key1 = "1"
        interaction_key2 = "2"
        interaction_key3 = "3"
        interaction_key4 = "num_side_chain_water_contacts"
        clusterwise_interaction_dict["sc_wat_res1_contact"],clusterwise_interaction_dict["sc_wat_res2_contact"],clusterwise_interaction_dict["sc_wat_res3_contact"],clusterwise_interaction_dict["total_sc_wat_contact"],list_mol_time = calculate_classified_sc_wat_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id_z,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4)
    #classified sc wat
    param_name = "sc_wat_res1_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/sidechain-water-contact/raw_data_all_frames"
        json_interaction_base = "side_chain_water_contacts_"
        interaction_key1 = "1"
        interaction_key2 = "2"
        interaction_key3 = "3"
        interaction_key4 = "num_side_chain_water_contacts"
        clusterwise_interaction_dict["sc_wat_res1_contact"],clusterwise_interaction_dict["sc_wat_res2_contact"],clusterwise_interaction_dict["sc_wat_res3_contact"],clusterwise_interaction_dict["total_sc_wat_contact"],list_mol_time = calculate_classified_sc_wat_for_list_of_mol(mol_list_cluster,base_interaction_path,list_frame_id_z,json_interaction_base,interaction_key1,interaction_key2,interaction_key3,interaction_key4)

    #electrostatic/salt-bridge contact : list of numbers for each cluster
    param_name = "electrostatic_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/salt-bridge-interaction/raw_data_all_frames"
        json_interaction_base = "pep_pep_salt_bridge_"
        #interaction_key = "num_side_chain_protein_contacts"
        clusterwise_interaction_dict["electrostatic_contact"],list_mol_time = calculate_electrostatic_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id_z3,json_interaction_base)
        
    # pep water contact : 
    param_name = "pep_pep_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/pep-pep-contact/raw_data_all_frames"
        json_interaction_base = "pep_pep_contacts_"
        interaction_key = "num_atom_contacts"
        clusterwise_interaction_dict["pep_pep_contact"],list_mol_time = calculate_interaction_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id_z3,json_interaction_base,interaction_key)


    # pep water contact : 
    param_name = "pep_wat_contact"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/pep-water-contact/raw_data_all_frames"
        json_interaction_base = "pep_wat_contacts_"
        interaction_key = "num_atom_contacts"
        clusterwise_interaction_dict["pep_wat_contact"],list_mol_time = calculate_interaction_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id_z3,json_interaction_base,interaction_key)
   
    #member of cluster size information  : list of numbers for each cluster
    param_name = "member_cluster_size"
    if(param_name in parameters_list):
        base_interaction_path = base_interaction_json+"/pep-pep-contact/mol_in_cluster_info"
        json_interaction_base = "mol_belongs_cluster_of_size_info_"
        #interaction_key = "num_side_chain_protein_contacts"
        clusterwise_interaction_dict["member_cluster_size"],list_mol_time = calculate_member_cluster_size_for_list_of_mol_frame_range(mol_list_cluster,base_interaction_path,list_frame_id_z3,json_interaction_base)
    
    # =======================================================================================================================================================================================================================================================================================================================================================================================
   

    #if(list_mol_time1 == list_mol_time and list_mol_time1 == list_mol_time2):
       # print("ITS matches the order")
    
    #print(list_mol_time)
    segname_list_array = np.array(list_mol_time)
    segname_list = segname_list_array[:,0]

    # ==========================================================================================
    for idx in range(num_parameters):
        #print(idx,parameters_list[idx])
        parameters[idx]=clusterwise_interaction_dict[parameters_list[idx]]
    #print(parameters)
    
    xyz_file = peptide+"_"
    first_info = peptide+" "
    for label in parameters_list:
        label = label.replace("middle","mid")
        xyz_file = xyz_file+label+"_"
        first_info = first_info+label+" "
    xyz_file = xyz_file+"time.xyz"
    first_info = first_info+"time"
    gc.make_xyz_file_general(num_parameters,parameters,time,segname_list,xyz_file,first_info) 
    print("OUTPUT:",xyz_file)
    return xyz_file


##########################################################################################################################
#  It extract all the values of order parameter enlisted in parameters_list(KINDLY MENTAIN THE NAME CONVENSION example parameter names are given among 
#   line 311-319)
# Output : combines xyz file , outname of the file is in the order of the order parameter
##########################################################################################################################
peptide = "lll"

# number of peptides present in the system 
system_size = 600

if(peptide == 'fff'):
    basepath_interaction_json ="/home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/"+peptide+"/composite_OP/250A"
    fiberpath_seganme_json = "/home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/"+peptide+"/system/250A/segment_names_all_"+peptide+"_250A.json"
    num_frame = 400
else: 
    basepath_interaction_json ="/home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/"+peptide+"/composite_OP"
    fiberpath_seganme_json = "/home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/"+peptide+"/system/segment_names_all_"+peptide+".json"
    num_frame = 500


#parameters_list = ["pep_pep_contact","pep_wat_contact"]
#parameters_list = ["total_H_bond_contact","sidechain_pep_contact"]
#parameters_list = ["total_H_bond_contact","sidechain_pep_contact","pep_wat_contact"]
#parameters_list = ["sidechain_pep_contact","pep_wat_contact"]
#parameters_list = ["sidechain_pep_contact","pep_wat_contact"]
#parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact","pep_wat_contact"]
#parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact"]
parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact"]
#parameters_list = ["total_H_bond_contact","sc_wat_res1_contact","sc_wat_res2_contact","sc_wat_res3_contact"]


#--------------------------------------------------------------------------------------------------------------------------------------------------------------
xyz_file = composite_OP_extraction(peptide,basepath_interaction_json,fiberpath_seganme_json,num_frame,parameters_list)
#composite_OP_pep_electrostatic_H(peptide,base_interaction_json,fiber_seganme_json,frame_end)
#composite_OP_H_bond_classified(peptide,base_interaction_json,fiber_seganme_json,frame_end)
#composite_OP_H_bond_classified_contact_based_cluster_size(peptide,base_interaction_json,fiber_seganme_json,frame_end)


# MAKING the input file for onward calculation of TAMCIS
dict_TAMCIS = {}
dict_TAMCIS["peptide"] = peptide
dict_TAMCIS["basepath_interaction_json"] = basepath_interaction_json
dict_TAMCIS["fiberpath_seganme_json"] = fiberpath_seganme_json
dict_TAMCIS["num_frame"] = num_frame
dict_TAMCIS["parameters_list"] = parameters_list
dict_TAMCIS["xyz_raw_file"] = xyz_file
dict_TAMCIS["system_size"] = system_size

input_json = "input_TAMCIS.json"
with open(input_json,"w+")  as f:
    json.dump(dict_TAMCIS,f,indent=4)

print(input_json," IS CREATED !!! GREAT 1st step DONE")
    




        

