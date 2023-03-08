## MESA

We first use MESA to create a model of a neutron star with mass $M=1.4M_\odot$ , radius $R=12$ km, and surface gravity $g=GM/R^2=1.29\times10^{14}$ g/cm $^2$ (MESA is Newtonian, and gravity has not been redshift-corrected). 

We then accrete a column $y\sim 2\times 10^4$ g/cm $^{2}$ of iron. This will be the substrate/buffer below the burning region.

We then accrete a "solar" mixture (70% hydrogen, 28% helium, 2% carbon-12) for a few days, until ignition in a helium layer. 

The last two profiles of the MESA run (`profile46.data` & `profile47.data`) are copied here. These steps are all done in [github.com/simonguichandut/mesa\_mixed\_burst](https://github.com/simonguichandut/mesa_mixed_burst). The relevant run is `run/make_ignition_model_for_hydro/`. 

## HSE model

To convert `profile46` (the pre-ignited data which we want to start the simulation from) into a amrex-astro usable format:

```
python convert_mesa.py profile46.data
```

creates `profile46.raw` (note: this script uses [py\_mesa\_reader](https://github.com/wmwolf/py_mesa_reader)). 

To convert the MESA lagrangian grid into eulerian (and readjust HSE if needed), we use the routines in [AMReX-ASTRO/initial\_models/lagrangian\_planar](https://github.com/AMReX-Astro/initial_models).  To get the correct list of species, the code there needs to be compiled with the CNO_extras network (`make NETWORK_DIR=CNO_extras`). The file `helm_table.dat` also needs to be copied over.

Values for gravity, domain dimensions, cutoffs, are set in `inputs_lagrangian_planar_xrb_layered`.

Run with:

```
{initial_models_directory}/lagrangian_planar/initialmodel1d.gnu.ex inputs_lagrangian_planar_xrb_layered
```

This creates the input file `xrb_layered_3cm.hse` used in the simulation.

## Figuring out the initial perturbation
The MESA profile we use is pre-ignition. We need to give it a "kick" (i.e. manually increase the temperature at the base) to get the burst to start.  We figure out how big this needs to be in `ignition_analysis.ipynb`.


