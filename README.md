# TAMCIS
The Trajectory Analysis of Multidimensional Chemical Interaction Space (TAMCIS) is a systematic approach to deal with multidimeansional order parameter space
for the entire trajectory. It mainly tries to annote the trajeotory data in the chosen p-dimensional space ($/Gamma$ space) by filtering out the time homogeneous fluctuations.
The workflow for TAMCIS is given as follows:

1. Determination of order parameters in $\Gamma$ space and value of p.  This involved data extraction in xyz format :
 XYZ file : peptide_name pep_H_end_end_contact pep_H_end_mid_contact pep_H_mid_mid_contact sidechain_pep_contact time [ order parameter list with time]
```bash
python3 multi_D_OP_data_extraction_trajectory_xyz_general_v1.py 
```
 ---- It calls general_code.py
