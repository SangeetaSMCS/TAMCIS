#source visual_input_render.tcl

proc disply_one_selection { mol_id sel representation_string color_string vmd_material } {
        mol selection $sel
        mol representation $representation_string
        mol color $color_string
        mol material $vmd_material
        mol addrep $mol_id
}



set vmd_material AOShiny

# VMD default colors. Integers next to the strings are VMD codes
set vmd_blue          0
set vmd_red           1
set vmd_gray          2
set vmd_orange        3
set vmd_yellow        4
set vmd_tan           5
set vmd_silver        6
set vmd_green         7
set vmd_white         8
set vmd_pink          9
set vmd_cyan          10
set vmd_purple        11
set vmd_lime          12
set vmd_mauve         13
set vmd_brown         14
set vmd_light_blue    15
set vmd_black         16

set newcartoon "NewCartoon 0.300000 40.000000 4.100000 0"
set newcartoon_thin "NewCartoon 0.29900000 40.000000 4.100000 0"
set newcartoon_thick "NewCartoon 0.3200000 40.000000 4.100000 0"
set newcartoon_very_thin "NewCartoon 0.1000000 40.000000 5.400000 0"
set newcartoon_thin_1 "NewCartoon 0.1300000 40.000000 5.400000 0"
set vdw "VDW 1.000000 30.000000"
set vdw_small "VDW 0.800000 30.000000"
set vdw_small_1 "VDW 0.700000 30.000000"
set ball_stick "CPK 1.000000 0.300000 30.000000 30.000000"
set cpk_stick "CPK 0.000000 0.300000 30.000000 30.000000"
#set licorice "Licorice 0.800000 30.000000 30.000000"
set licorice "Licorice 0.600000 30.000000 30.000000"
set cartoon "Cartoon 2.100000 40.000000 5.000000"
set cartoon_thin "Cartoon 2.000000 40.000000 4.9000000"
set cartoon_very_thin "Cartoon 1.000000 40.000000 2.5000000"
set lines "Lines 2.000000"
set surface "Surf 1.400000 0.000000"
set surface_1 "Surf 1.600000 0.000000"
set vmd_trace "Trace 0.300000 15.000000"
set quick_surf "QuickSurf 0.800000 0.500000 1.000000 1.000000"

##############USER INPUT#################################################################
set peptide lll

##########################################################################
#set base_path /home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/$peptide/pdb_file/500ns
#set psf_file $base_path/final_periodicity_corrected_fiber1_with_original_segments.psf
#set pdb_file $base_path/final_periodicity_corrected_fiber1_with_original_segments.pdb

#for {set frame_id 200}  {$frame_id < 300} {incr frame_id} {
source input_movie.tcl
#set frame_id 493
set tag chbsp
set base /home/sangeeta/sangeeta_adlab/clustering/final_self_assembled_peptides/$peptide
set psf_file $base/traj/trajectory.psf
set pdb_file $base/pdb_file/traj_pdb/frame_$frame_id.pdb
set viewpt_path $base/pdb_file/traj_pdb/final_viewpoint.tcl



set mol_id [mol load psf $psf_file pdb $pdb_file]

##########################################################################################



mol delrep 0 $mol_id


source number_prominent_class.tcl
for {set i 1} {$i <= $num_pclass} {incr i} {
    source frame_$frame_id/$peptide\_category_segname_$tag\_$i.tcl
    if {[llength $cluster_list_seg] > 0} {
        set sel0 "segname $cluster_list_seg and backbone"
        #set representation_string $lines
        set representation_string $licorice
        set color_string "ColorId $color_id"
        #disply_one_selection $mol_id $sel $representation_string $color_string Transparent
        disply_one_selection $mol_id $sel0 $representation_string $color_string $vmd_material
    }
}




####################################################################################
if 0 {
source $peptide\_core_intermediate_periphery.tcl
if {[llength $segment_list_core] > 0} {
    set sel7 "segname $segment_list_core"
    #set representation_string $lines
    set representation_string $quick_surf
    set color_string "ColorId 0"
    disply_one_selection $mol_id $sel7 $representation_string $color_string Transparent
    #disply_one_selection $mol_id $sel7 $representation_string $color_string $vmd_material
}

if {[llength $segment_list_intermediate] > 0} {
    set sel8 "segname $segment_list_intermediate"
    #set representation_string $lines
    set representation_string $quick_surf
    set color_string "ColorId 6"
    disply_one_selection $mol_id $sel8 $representation_string $color_string Transparent
    #disply_one_selection $mol_id $sel8 $representation_string $color_string $vmd_material
}

if {[llength $segment_list_periphery] > 0} {
    set sel9 "segname $segment_list_periphery"
    #set representation_string $lines
    set representation_string $quick_surf
    set color_string "ColorId 6"
    disply_one_selection $mol_id $sel9 $representation_string $color_string Transparent
    #disply_one_selection $mol_id $sel9 $representation_string $color_string $vmd_material
}
}

######################################################################################################################

display projection Orthographic
#display nearclip set 0.010000
#display cuedensity 0.300000
display ambientocclusion off
color Display Background white
axes location Off
color Labels Bonds black
label textthickness 1.515152

set tga_file movie_$peptide/frame_image_$frame_id.dat
source $viewpt_path
set viewpoints($mol_id) $viewpoint1
restore_viewpoint
render Tachyon $tga_file "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -o %s.tga
#render snapshot broken_mol_kmeans_v2/$sys\_$fiber\_broken_mol_kmeans.tga display %s

exit 
#}
