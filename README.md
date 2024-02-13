# TAMCIS
The Trajectory Analysis of Multidimensional Chemical Interaction Space (TAMCIS) is a systematic approach to deal with multidimeansional order parameter space
for the entire trajectory. It mainly tries to annote the trajeotory data in the chosen p-dimensional space ($/Gamma$ space) by filtering out the time homogeneous fluctuations.
The workflow for TAMCIS is given as follows:

1. Assumping that single molecule-wise chemical interactions are already calculated for entire trajectory length and stored in a organized format in a folder name "iii" for say tri-Isoleucine system
   Example Folder Structure :
   
2. Determination of order parameters in $\Gamma$ space and value of p.  This involved data extraction in xyz format :
 XYZ file : peptide_name pep_H_end_end_contact pep_H_end_mid_contact pep_H_mid_mid_contact sidechain_pep_contact time [ order parameter list with time]
```bash
python3 multi_D_OP_data_extraction_trajectory_xyz_general_v1.py 
```
Input : peptide name, base_interaction_json, fiber_seganme_json, frame_end
        parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact"]

 ---- It calls general_code.py

3. analysis_for_binning_parameter_v1.py 
4. manual_binning_parameters_general.py
5. DBSCAN_organize_cluster_general_v3.py
6. post_processing_stray_info_remove_v1.py
7. check_postprocess.py
8. 
