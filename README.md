# TAMCIS
The Trajectory Analysis of Multidimensional Chemical Interaction Space (TAMCIS) is a systematic approach to deal with multidimeansional order parameter space
for the entire trajectory. It mainly tries to annote the trajeotory data in the chosen p-dimensional space ($\Gamma$ space) by filtering out the time homogeneous fluctuations.
The workflow for TAMCIS is given as follows:

1. Assumping that single molecule-wise chemical interactions are already calculated for entire trajectory length and stored in a organized format in a folder name "iii" for say tri-Isoleucine system
   Example Folder Structure :
   
2. Determination of order parameters in $\Gamma$ space and value of p.  This involved data extraction in xyz format :
 XYZ file : peptide_name pep_H_end_end_contact pep_H_end_mid_contact pep_H_mid_mid_contact sidechain_pep_contact time [ order parameter list with time]
```bash
# EDIT IN LINE 290-299
python3 multi_D_OP_data_extraction_trajectory_xyz_general_v2.py 
```
Input : peptide name, base_interaction_json, fiber_seganme_json, frame_end
        parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact"]

 ---- It calls general_code.py
OUTPUT : input_TAMCIS.json , This file store every information so further TAMCIS code take this file as INPUT 

3. Look into time independent dataset of the all molecule for each order parameter seperately for taking decision about the manual bining step:
   Bar plot of the entire dataset. IT WILL TAKE INPUT "input_TAMCIS.json" automatically. JUST RUN the code
   ```bash
   python3 analysis_for_binning_parameter_v1.py
   ```
4. Manually bin each parameter in a chemically meaningful way.
   INPUT : INSIDE CODE mention the binning for individual order parameter, as the INDEX OF THE ORDER PARAMETER STARTS FROM 1 and store in dictionary formate
    
   ```bash
   # EDIT AT LINE : 10 
   python3 manual_binning_parameters_general_v1.py
  
5. Clustering on the binned data :
```bash
python3 clustering_multi_D_binned_data_general_v3.py ```
```
5. Obtain organized json file
```bash
python3 DBSCAN_organize_cluster_general_v3.py
```
7. post_processing_stray_info_remove_v1.py
8. check_postprocess.py
9. 
