#!/bin/bash

# number of jobs to submit
njobs=10
submission_script="sub_rotstar.titan.pbs"

jobid=$(qsub $submission_script)

if [ $njobs -gt 1 ]; then 
    for i in $(seq 2 $njobs)
    do
	jobid=$(qsub -W depend=afterok:$jobid $submission_script)
    done
fi