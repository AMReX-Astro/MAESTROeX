In an earlier version of xrb_layered, we directly converted the MESA profile to an eulerian grid using amrex initial models lagrangian routines. Instructions below, kept here for future reference. In the end, we decided to build the model from scratch, as explained in the parent directory readme.

## HSE model

To convert `profile37` (the pre-ignited data which we want to start the simulation from) into a amrex-astro usable format:

```
python convert_mesa.py profile37.data
```

creates `profile37.raw` (note: this script uses [py\_mesa\_reader](https://github.com/wmwolf/py_mesa_reader)). 

To convert the MESA lagrangian grid into eulerian (and readjust HSE if needed), we use the routines in [AMReX-ASTRO/initial\_models/lagrangian\_planar](https://github.com/AMReX-Astro/initial_models).  To get the correct list of species, the code there needs to be compiled with the CNO_extras network (`make NETWORK_DIR=CNO_extras`). The file `helm_table.dat` also needs to be copied over.

Values for gravity, domain dimensions, cutoffs, are set in `inputs_lagrangian_planar_xrb_layered`.

Run with:

```
{initial_models_directory}/lagrangian_planar/initialmodel1d.gnu.ex inputs_lagrangian_planar_xrb_layered
```

This creates the input file `xrb_layered_3cm.hse` used in the simulation.