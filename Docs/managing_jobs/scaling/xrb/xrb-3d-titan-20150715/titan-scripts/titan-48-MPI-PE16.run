#!/bin/bash
#PBS -A ast106
#PBS -N xrb-48-MPI-PE16
#PBS -j oe
#PBS -l walltime=0:30:00,nodes=64
#PBS -q batch
#PBS -l gres=atlas1%atlas2

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

# Cray uses the term "processing elements" (PEs) to refer to the
# number of executables that are launched.  For our purposes, this is
# the number of MPI tasks.

# from https://www.olcf.ornl.gov/kb_articles/xk7-cpu-description/
#
# Each Titan compute node contains (1) AMD Opteron™ 6274 (Interlagos)
# CPU. Each CPU contains (2) die. Each die contains (4) “bulldozer”
# compute units and a shared L3 cache. Each compute unit contains (2)
# integer cores (and their L1 cache), a shared floating point
# scheduler, and shared L2 cache.
#
# To aid in task placement, each die is organized into a NUMA
# node. Each compute node contains (2) NUMA nodes. Each NUMA node
# contains a die’s L3 cache and its (4) compute units (8 cores).

# -n number of PEs (MPI tasks)
# -N number of PEs per node  (MPI tasks / node)
# -d number of CPUs per PE -- this should be the number of threads / MPI task
# -j number of CPUs to use per compute unit (since 2 CPUs share an FPU, we may want to reduce this)
# -S PEs to allocate per NUMA node (if we have 2 MPI tasks per node, we want -S 2 -- one on each NUMA node)
aprun -n 1024 -N 16 -d 1 -j 2 ./main.Linux.Cray.mpi.omp.exe inputs_3d_6.0cm.hi_dens
