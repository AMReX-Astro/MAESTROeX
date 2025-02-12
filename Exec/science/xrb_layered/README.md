# xrb_layered

Models the convection and burning in the envelope of a neutron star, where hydrogen
has previously been depleted by the hot CNO cycle, leaving a pure helium layer at 
large depths, where the burst will ignite.

Appears in:

  * `Hydrodynamical simulations of proton ingestion flashes in Type I X-ray Bursts`
     Guichandut, S.; Zingale, M.; Cumming, A, 2023. 
     https://arxiv.org/abs/2405.08952


This is the multidimensional version of a calculation previously done in MESA in:

  * `The imprint of convection on Type I X-ray bursts: Pauses in photospheric
     radius expansion lightcurves`, Guichandut, S.; Cumming, A, 2023. 
     https://arxiv.org/abs/2301.08769


## initial model

* `toy_atm_hot_3cm.hse`

Build using a MESA model as a guide. More details in `create_initial_model/`.
