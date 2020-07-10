# Rotating fully convective star diagnostics

Process a fully convective star problem to produce a set of diagnostic measures
including:

- Convection speed
- Ratio of mean rotation rates to surface rotation
- Meridional circulation velocities
- Br\"{u}nt-V\"{a}is\"{a}l\"{a} frequency

## Building & running

The n-dimensional diagnostic can be built by executing `make DIM=n`. This will
produce the executable `radial_nd.exe`. Note that `DIM=3` is default setting and
the only one that is currently supported. To run, the executable must be provided
with the name of the plotfile to be analyzed:
```
./radial_3d.exe infile=plotfile_name
```

Additional arguments are as follows:

- **model file diagnostics**: extra argument `modelfile=modelfile_name` need to be provided

