#!/bin/sh


mpi=`grep "number of MPI processes =" $1 | tail -n 10 | awk '{print $6}'`
omp=`grep "number of threads       =" $1 | tail -n 10 | awk '{print $5}'`


adv=`grep "Advection       :" $1  | tail -n 10 | awk '{sum += $3; count += 1} END {print sum/count}'`
mac=`grep "MAC   Projection:" $1  | tail -n 10 | awk '{sum += $3; count += 1} END {print sum/count}'`
nodal=`grep "Nodal Projection:" $1  | tail -n 10 | awk '{sum += $3; count += 1} END {print sum/count}'`
react=`grep "Reactions       :" $1  | tail -n 10 | awk '{sum += $3; count += 1} END {print sum/count}'`
misc=`grep "Misc            :" $1  | tail -n 10 | awk '{sum += $3; count += 1} END {print sum/count}'`

tot=`grep "Time to advance timestep:" $1  | tail -n 10 | awk '{sum += $5; count += 1} END {print sum/count}'`

printf "# %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n" \
  "MPI" "threads" "advection" "MAC" "nodal" "reactions" "misc" "total"

printf " %10d  %10d  %10f  %10f  %10f  %10f  %10f  %10f\n" \
   $mpi   $omp  $adv  $mac  $nodal  $react  $misc  $tot


