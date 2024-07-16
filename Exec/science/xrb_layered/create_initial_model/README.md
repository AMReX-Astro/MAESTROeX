## MESA

We first use MESA to create a model of a neutron star with mass $M=1.4M_\odot$ , radius $R=12$ km, and surface gravity $g=GM/R^2=1.29\times10^{14}$ g/cm $^2$ (MESA is Newtonian, and gravity has not been redshift-corrected). Initially, the model is just an atmosphere of pure iron with a column depth of 5e9 g/cm2.

We then accrete a "solar" mixture (70% hydrogen, 28% helium, 2% carbon-12) for a few days, until ignition in a helium layer. 

The last two models of the MESA run (`profile4.data` & `profile5.data`) are copied in `mesa/`. These steps are all done in [github.com/simonguichandut/mesa\_mixed\_burst](https://github.com/simonguichandut/mesa_mixed_burst). The relevant run is `run/make_ignition_model_for_hydro/`. 

## Toy model

See the notebook `Initial_Model.ipynb` for a detailed explanation of how we use the pre-igntion MESA model (profile4) to construct a simplified toy model, with a temperature kick at the bottom of the fuel layer. This is done in `toy_atm/`, which is a modified version of https://github.com/AMReX-Astro/initial_models/tree/main/toy_atm .