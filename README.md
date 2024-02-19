# TAMCIS
The Trajectory Analysis of Multidimensional Chemical Interaction Space (TAMCIS) is a systematic approach to deal with multidimeansional order parameter space
for the entire trajectory. It mainly tries to annote the trajeotory data in the chosen p-dimensional space ($\Gamma$ space) by filtering out the time homogeneous fluctuations.
The workflow for TAMCIS is given as follows:

1. Assumping that single molecule-wise chemical interactions are already calculated for entire trajectory length and stored in a organized format in a folder name "iii" for say tri-Isoleucine system
   Example Folder Structure :
   
2. Determination of order parameters in $\Gamma$ space and value of p.  This involved data extraction in xyz format :
 XYZ file : peptide_name pep_H_end_end_contact pep_H_end_mid_contact pep_H_mid_mid_contact sidechain_pep_contact time [ order parameter list with time]
```bash
# EDIT IN LINE 311-319
python3 step1_multi_D_OP_data_extraction_trajectory_xyz_general.py
```
Input : peptide name, base_interaction_json, fiber_seganme_json, frame_end
        parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact"]

 ---- It calls general_code.py
OUTPUT : input_TAMCIS.json , This file store every information so further TAMCIS code take this file as INPUT 

3. Look into time independent dataset of the all molecule for each order parameter seperately for taking decision about the manual bining step:
   Bar plot of the entire dataset. IT WILL TAKE INPUT "input_TAMCIS.json" automatically. JUST RUN the code
```bash
python3 step2_analysis_for_binning_parameter.py
```
4. Manually bin each parameter in a chemically meaningful way.
   INPUT : INSIDE CODE mention the binning for individual order parameter, as the INDEX OF THE ORDER PARAMETER STARTS FROM 1 and store in dictionary formate
    
```bash
# EDIT AT LINE : 10 
python3 step3_manual_binning_parameters_general.py
```  
5. Clustering on the binned data : 
```bash
# Copy the input_TAMCIS.json and XXX_BINNED.dat file formed after binning : to the location where clustering script is executed [High Memory Computation]
# maruti : node 37-43
module load codes/python-3.9.15
python3 step4_clustering_multi_D_binned_data_general.py
```
NOTE : Copy to your analysis location : output cluster dat file and ALSO the input_TAMCIS.json [this file have been updated]

6. Organize the clustering data  to json file format
```bash
python3 step5_organize_cluster_general.py
```
7. Extarct stray molecules information timewise
```bash
# LINE 52-53 provide the peptide name and num_frame
python3 step6_sizewise_cluster_gather.py
```
NOTE : This calculation is done on Argha's computer, example path : /media/argha/home_21/my_work/analysis/iii-600-mol-200A-box/pep-pep-contact/post_process/cluster_analysis
Copy : iii_all_cluster_of_size_1_1.json (example) to your working location.


8. This post processing step where stray molecules are removed from the clusters
```bash
# LINE 184 : Provide a shortcut tag for output file  
python3 step7_post_processing_stray_info_remove.py
```

9. General sanity check for any mismatch and just run this code 
```bash
python3 step8_check_postprocess.py
```
10. 
