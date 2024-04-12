#!/bin/bash
for i in {0..200}; 
do
   echo "set frame_id $i" >> input_movie.tcl
   vmd -e step17_visualize_CI_clusters_movie.tcl
   rm movie_lll/*.dat
   rm input_movie.tcl
done 
