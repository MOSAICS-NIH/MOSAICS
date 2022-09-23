#!/bin/bash
core_name=$1
imin=$2
imax=$3
folder=$4
stride=$5 
dt=$6

i=0
for ((i=$imin; i<$imax; i+=$stride)); do
    in_name=$core_name'_'"$i"".dat"
    out_name=$core_name'_'"$i"".png"
    time=`echo "$i*$dt*0.000001" | bc`
    gnuplot -c heatmap_template_video.gnu [0:400]   [0:400] [0:1357]     $in_name   $out_name "Simulation time $time us"
    mv $out_name $folder
done



