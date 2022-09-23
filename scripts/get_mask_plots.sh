#!/bin/bash

core_name=$1
imin=$2
imax=$3
folder=$4

i=0
for ((i=$imin; i<$imax; i++)); do
    in_name=$core_name'_'"$i""_mask.dat"
    out_name=$core_name'_'"$i""_mask.png"
    echo $in_name $out_name

    dist=`echo "$i*0.1" | bc`
    gnuplot -c ../heatmap_template.gnu   [0:400]   [0:400] [:]     $in_name   $out_name "d = $dist +/- 0.5nm"
    mv $out_name $folder
done


