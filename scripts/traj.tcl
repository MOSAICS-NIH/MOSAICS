package require pbctools

set stride 1
set interval 16666
set num_sections 5
set ref_name  "clc_100_0.viparr_phi_psi.dms"
set traj_name "run.stk"

#set the working directory
cd /directory_with_traj/

#loop through chunks of trajectory
for {set i 0} {$i < $num_sections} {incr i} {

    #set start and end frames
    set start [expr $i * $interval]
    set end [expr ($i+1) * $interval - 1]

    #set the output file name
    set out_name "/name_of_output_directory/run"
    append out_name "_"
    append out_name $start
    append out_name "_"
    append out_name $end
    append out_name ".trr"

    #open reference molecule
    set trajMol [mol new $ref_name]

    #load trajectory
    mol addfile $traj_name first $start last $end step $stride waitfor all molid trajMol

    #write the output file
    animate write trr $out_name beg 1 end -1 waitfor all

    #delete the trajectory
    mol delete $trajMol
}

#quit vmd
quit

