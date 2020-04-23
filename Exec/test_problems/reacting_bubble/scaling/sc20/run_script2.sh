#!/bin/bash
#BSUB -P AST106
#BSUB -W 15
#BSUB -nnodes 8
#BSUB -alloc_flags smt1
#BSUB -J MAESTROeX
#BSUB -o MAESTROeX.%J
#BSUB -e MAESTROeX.%J

cd $LS_SUBCWD

inputs=inputs_3d_x2

n_mpi=48 # num nodes * 6 gpu per node
n_gpu=1
n_cores=1
n_rs_per_node=6

MAESTROeX_ex=./Maestro3d.pgi.TPROF.MPI.CUDA.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu $MAESTROeX_ex $inputs
