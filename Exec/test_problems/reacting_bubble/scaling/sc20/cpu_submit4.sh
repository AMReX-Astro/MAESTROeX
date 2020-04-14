#!/bin/bash
#BSUB -P ast106
#BSUB -W 0:10
#BSUB -nnodes 64
#BSUB -alloc_flags smt1
#BSUB -J cpu
#BSUB -o cpu.%J
#BSUB -e cpu.%J

cd $LS_SUBCWD/../..

inputs_file=scaling/sc20/inputs_3d_x4

n_mpi=448 # num nodes * 7 mpi per node 
n_threads=6
n_rs_per_node=7

MAESTROeX_ex=./Maestro3d.pgi.MPI.OMP.ex

export OMP_NUM_THREADS=6

jsrun -n $n_mpi -r $n_rs_per_node -c $n_threads -a 1 -d packed -b rs $MAESTROeX_ex $inputs_file
