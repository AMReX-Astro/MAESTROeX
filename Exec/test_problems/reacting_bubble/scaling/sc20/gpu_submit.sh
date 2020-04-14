#!/bin/bash
#BSUB -P ast106
#BSUB -W 0:30
#BSUB -nnodes 1
#BSUB -alloc_flags smt1
#BSUB -J bubble
#BSUB -o bubble.%J
#BSUB -e bubble.%J

cd $LS_SUBCWD/../..

inputs_file=scaling/sc20/inputs_3d

n_mpi=6 # num nodes * 6 gpu per node
n_gpu=1
n_cores=1
n_rs_per_node=6

MAESTROeX_ex=./Maestro3d.pgi.MPI.CUDA.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu $MAESTROeX_ex $inputs_file 
