#!/bin/bash

# first I'm going to copy everything over to $MEMBERWORK
MEMBERWORKDIR="$MEMBERWORK/ast106"
PROBLEMDIR="$MEMBERWORKDIR/MAESTROeX/Exec/SCIENCE/rotating_star"

cd $MEMBERWORKDIR
rm -rf MAESTROeX
rm -rf yt
#rm -rf Microphysics
#rm -rf amrex

cp -r ~/MAESTROeX .
cp -r ~/yt .
#cp -r ~/amrex .
#cp -r ~/Microphysics .

cd $PROBLEMDIR

# we're going to need to copy the actual helm_table.dat file here so that it can be accessed
HELM_TABLE=$(readlink -f helm_table.dat)
rm helm_table.dat
cp $HELM_TABLE helm_table.dat

# number of jobs to submit
njobs=10
submission_script="sub_rotstar.titan.pbs"

jobid=$(qsub $submission_script)

if [ $njobs -gt 1 ]; then 
    for i in $(seq 2 $njobs)
    do
	jobid=$(qsub -W depend=afterany:$jobid $submission_script)
    done
fi