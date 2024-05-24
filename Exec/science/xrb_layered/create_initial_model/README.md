## MESA

We first use MESA to create a model of a neutron star with mass $M=1.4M_\odot$ , radius $R=12$ km, and surface gravity $g=GM/R^2=1.29\times10^{14}$ g/cm $^2$ (MESA is Newtonian, and gravity has not been redshift-corrected). 

We then accrete a column $y\sim 2\times 10^4$ g/cm $^{2}$ of iron. This will be the substrate/buffer below the burning region.

We then accrete a "solar" mixture (70% hydrogen, 28% helium, 2% carbon-12) for a few days, until ignition in a helium layer. 

The last two profiles of the MESA run (`profile37.data` & `profile38.data`) are copied in `mesa/`. These steps are all done in [github.com/simonguichandut/mesa\_mixed\_burst](https://github.com/simonguichandut/mesa_mixed_burst). The relevant run is `run/make_ignition_model_for_hydro/`. 

## Toy model

See the notebook `Initial_Model.ipynb` for a detailed explanation of how we use the MESA model (profile37) to construct a simplified toy model, with a temperature kick at the bottom of the fuel layer. This is done in `toy_atm/`, which is a modified version of https://github.com/AMReX-Astro/initial_models/tree/main/toy_atm .