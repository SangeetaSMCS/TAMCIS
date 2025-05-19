# TAMCIS
The Trajectory Analysis of Multidimensional Chemical Interaction Space (TAMCIS) is a systematic approach to deal with multidimeansional order parameter space
for the entire trajectory. It mainly tries to annote the trajeotory data in the chosen p-dimensional space ($\Gamma$ space) by filtering out the time homogeneous fluctuations.
The workflow for TAMCIS is given as follows:



 Assumping that single molecule-wise chemical interactions are already calculated for entire trajectory length and stored in a organized format in a folder name "iii" for say tri-Isoleucine system
   Example Folder Structure :
   
1. Determination of order parameters in $\Gamma$ space and value of p.  This involved data extraction in xyz format :
   XYZ file : peptide_name pep_H_end_end_contact pep_H_end_mid_contact pep_H_mid_mid_contact sidechain_pep_contact time [ order parameter list with time]

  ```bash
  #INPUT : peptide name, base_interaction_json, fiber_seganme_json, frame_end
  #        parameters_list = ["pep_H_end_end_contact","pep_H_end_middle_contact","pep_H_middle_middle_contact","sidechain_pep_contact"]
  #---- It calls general_code.py
  #OUTPUT : input_TAMCIS.json , This xyz file store every information so further TAMCIS code take this file as INPUT eg. #"iii_pep_H_end_end_contact_pep_H_end_mid_contact_pep_H_mid_mid_contact_sidechain_pep_contact_time.xyz"
  # -----------------------
  # EDIT IN LINE 290-319
  #------------------------
  python3 step1_multi_D_OP_data_extraction_trajectory_xyz_general.py
  ```

2. Look into time independent dataset of the all molecule for each order parameter seperately for taking decision about the manual bining step:
   Bar plot of the entire dataset. IT WILL TAKE INPUT "input_TAMCIS.json" automatically. JUST RUN the code
```bash
 #INPUT : No manual input needed
 #OUPUT : png : "iii_analysis_4_bining_chbsc.png"
python3 step2_analysis_for_binning_parameter.py
```

3. Manually bin each parameter in a chemically meaningful way.
```bash
#INPUT : INSIDE CODE mention the binning for individual order parameter, as the INDEX OF THE ORDER PARAMETER STARTS FROM 1 and store in dictionary formate
#OUTPUT : updated input_TAMCIS.json and  the binned xyz file eg. "iii_pep_H_end_end_contact_pep_H_end_mid_contact_pep_H_mid_mid_contact_sidechain_pep_contact_time_BINNED.xyz"
#-------------------
# EDIT AT LINE : 10
#-------------------
python3 step3_manual_binning_parameters_general.py
```  
4. Clustering on the binned data :
```bash
#INPUT : input_TAMCIS.json and XXX_BINNED.xyz
#OUTPUT : updated input_TAMCIS.json and clustering dat file eg. "iii_cluster_data_ee_em_mm_sc_pep_min_sample_1_without_time_tgap_1_BINNED.dat"
# Copy the input_TAMCIS.json and XXX_BINNED.xyz file formed after binning : to the location where clustering script is executed [High Memory Computation]
# maruti : node 37-43
module load codes/python-3.9.15
python3 step4_clustering_multi_D_binned_data_general.py
```
NOTE : Copy to your analysis location : output cluster dat file and ALSO the input_TAMCIS.json [this file have been updated]

5. Organize the clustering data  to json file format
```bash
#INPUT : No manual input , it automatically takes input_TAMCIS.json
#OUTPUT : Sorted clusters json file eg. "iii_cluster_sorted_data_ee_em_mm_sc_pep_min_sample_1_without_time_tgap_1_BINNED.json" and updated input_TAMCIS.json
python3 step5_organize_cluster_general.py
```
6. Extarct stray molecules information timewise
Copy : iii_all_cluster_of_size_1_1.json (example) to your working location.
```bash
# LINE 52-53 provide the peptide name and num_frame
python3 step6_sizewise_cluster_gather.py
```

7. This post processing step where stray molecules are removed from the clusters
```bash
#INPUT : shortcut tag, alphanumeric name and stray molecule information eg. "iii_all_cluster_of_size_1_1.json"
#OUTPUT : clusters with alphanumeric name srored in json eg. "iii_cluster_alphanumeric_name_chbsp.json" and updated input_TAMCIS.json
#--------------------------------------------------
# LINE 184 : Provide a shortcut tag for output file
#--------------------------------------------------  
python3 step7_post_processing_stray_info_remove.py
```

8. General sanity check for any mismatch and just run this code 
```bash
INPUT : No manual input need
python3 step8_check_postprocess.py
```
9. Class-wise classification and sorting
```bash
#INPUT: order parameter name to understand the class and colors assigned to each class
#OUTPUT: json for classwise classification and respective gamma clusters in each class eg. "iii_sorted_dict_classes_gamma_cluster.json"
# EDIT AT LINE 110-117
python3 step9_sort_gamma_clusters_bin_and_classes.py
```
10. Classwise looking into the % population evalution of $\Gamma$ clusters 
```bash
#INPUT : Relative percentage population for heatmap dictionary for each class
# For time being it has been set so no need to change
#OUTPUT : $\leq 2^p-1$ Heatmap plots
#-----------------------------
# EDIT AT LINE 10
#-----------------------------
python3 step10_heatmap_time_all_possible_classwise_norm_sys_size_v1.py
```

11.
```bash
#INPUT : No manual input need
#OUTPUT : 
python3 analysis2_temporal_evol_all_classes_heatmap_n_lineplot.py
```
13. analysis_new1_heatmap_class_vs_clustersize_norm_csize_v2.py

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
11. Filtering those $\Gamma$ clusters based on  population over time and looking into the time evolution of those prominent $\Gamma$ clusters. For that we have put 2 filters : 1) population per timestep for $\gamma$ cluster atleast cross the  cutoff_num_mol_per_frame 2) $\Gamma$ cluster must present atleast  cutoff_num_occurences time out of the total length of the trajectory
```bash
#INPUT : percent_population ,cutoff_num_occurences, cutoff_num_mol_per_frame
#OUTPUT : Heatmap plot eg. "iii_CI_cluster_heatmap_time_evolution_chbsp_v2.png", update input_TAMCIS.json with import prominent classes
#------------------------
# EDIT AT LINE 9-11
#------------------------
python3 step11_heatmap_time_vs_class_clusters_filtered.py
```
12. Temporal every molecular hopping in prominent classes over time
```bash
#INPUT : no manual input
#OUTPUT : every molecular temporal hopping stored in csv file eg. "iii_temporal_hopping_chbsp.csv" and count for number of hopping exist for each prominent classes respect to $\Gamma$ clusters in json file eg.   
"iii_CI_cluster_hop_count_chbsp.json"
#ALERT : It is significant time consuming

python3 step12_temporal_molecular_hopping_imp_classes.py
```
13. heatmap plot every molecule in the system temporal hopping in classes
```bash
#INPUT:No manual input
#OUTPUT:png eg. "iii_temporal_hopping_for_fiber1_chbsp_v1.png"
python3 step13_all_molecular_hopping_plot_v1.py
```
14. heatmap plot selected molecule in the system temporal hopping in classes
```bash
#INPUT: json of segname list for selected molecules , selected indices if don't want to look at all segnames
#OUTPUT: png eg. "iii_temporal_hopping_for_fiber1_chbsp_v1.png"
#---------------------------
#EDIT AT LINE 39-43, 84-91
#--------------------------
python3 step14_selected_molecular_hopping_plot_v1.py
```
15. Assign VMD color ID for prominent class
```bash
#INPUT: Will ask for respective color VMD id through input with running the code
#OUTPUT: json file VMD_color_prominent_class.json

python3 step15_assign_VMD_color_for_pclasses.py
```
16. Move making : Extract tcl file of **selected** segnames (molecules participating in self assembly) for prominent classes for every timestep. Stored formate : e.g: frame_xxx/iii_category_segname_chbsp.tcl
```bash
#INPUT: json of segname list for selected molecules 
#OUTPUT: tcl files for every frame inside frame_id for respective class with color id set
#---------------------------
#EDIT AT LINE 56-60
#--------------------------
python3 step16_selected_molecular_hopping_for_movie.py
```
17. Render VMD picture for each frames with each classes 
```bash
#INPUT: set starting and ending frame number : LINE : 1 in step17_script_vmd_render_pic.sh
#INPUT: peptide tag LINE: 56 in step17_visualize_CI_clusters_movie.tcl
#OUTPUT : .tga files inside e.g. movie_iii
chmod 777 step17_script_vmd_render_pic.sh
./step17_script_vmd_render_pic.sh
```
 18. Look at VMD graphics at each frame for selected molecules
```bash
#INPUT: set peptide lll LINE:56
# set frame_id 499 LINE:66
vmd -e step18_visualize_CI_clusters_timepoint.tcl
```

