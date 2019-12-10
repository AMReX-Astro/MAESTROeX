# Planar diagnostics test

Sets up a 2D isothermal atmosphere and records diagnostics for ten time-steps.

There is some stuff in the code about a symmetric perturbation that is related to a linear gravity wave problem from which this is derived. The amplitude is just set to zero so no perturbation is made to the base state.

A jupyter notebook is used to load the data and compare to the expected values. There are some analytically expected values to compare to, and also you can calculate e.g. total KE in python and compare to the version the code recorded. 
