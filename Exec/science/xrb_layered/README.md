# xrb_layered

Models the convection and burning in the envelope of a neutron star, where hydrogen
has previously been depleted by the hot CNO cycle, leaving a pure helium layer at 
large depths, where the burst will ignite.

This is the multidimensional version of a calculation previously done in MESA in:

  * `The imprint of convection on Type I X-ray bursts: Pauses in photospheric
     radius expansion lightcurves`, Guichandut, S.; Cumming, A, 2023. 
     https://arxiv.org/abs/2301.08769


## initial model

* `xrb_layered_3cm.hse`

An initial MESA model is placed into eulerian HSE using the lagrangian_planar routines
from AMREX/initial_models. More detail in create_initial_model/readme.
