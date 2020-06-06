Reacting bubble scaling test
============================

This directory contains a weak scaling setup for the reacting bubble problem.
This setup was designed for use on the OLCF Summit supercomputer.

To build and run, have the following modules loaded:

```
spectrum-mpi/10.3.1.2-20200121
cuda/10.1.243
pgi/19.10
```

`cd` up two directories, and `make -j16` will then build the executable `Maestro3d.pgi.TPROF.MPI.CUDA.ex`.
Then return to this directory and create symbolic links to `../../Maestro3d.pgi.TPROF.MPI.CUDA.ex`,
`../../helm_table.dat`, and `../../model.hse.cool.coulomb`. Then submit each of the five scaling runs
with `bsub run_script.sh`, `bsub run_script2.sh`, etc. Each script correspond to inputs files with
a different problem size. Note that you will need to ensure the project flag (`-P`) is valid for you
if running on Summit.

To generate a scaling plot based on the results in the `MAESTROeX.?????` stdout files from
the runs, simply run the `plot.py` script, which depends only on matplotlib and numpy. It is
recommended to have the `python/3.6.6-anaconda3-5.3.0` module loaded.
