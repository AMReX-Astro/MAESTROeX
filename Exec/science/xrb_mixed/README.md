# xrb_mixed

This setup models convection in the accreted layers (either pure He or
mixed H/He) on the surface of a neutron star, preceding an X-ray burst.

Science papers that used this setup include:

  * `Comparisons of Two- and Three-Dimensional Convection in Type I
    X-Ray Bursts`, Zingale, M.; Malone, C. M.; Nonaka, A.; Almgren,
    A. S.; Bell, J. B., 2015, The Astrophysical Journal, Volume 807,
    Issue 1, article id. 60

    https://doi.org/10.1088/0004-637X/807/1/60

  * `Multidimensional Modeling of Type I X-Ray Bursts. II.
    Two-dimensional Convection in a Mixed H/He Accretor`, Malone,
    C. M.; Zingale, M.; Nonaka, A.; Almgren, A. S.; Bell, J. B., 2014,
    The Astrophysical Journal, Volume 788, Issue 2, article id. 115

    https://doi.org/10.1088/0004-637X/788/2/115

  * `Multidimensional Modeling of Type I X-ray Bursts. I.
    Two-dimensional Convection Prior to the Outburst of a Pure 4He
    Accretor`, Malone, C. M.; Nonaka, A.; Almgren, A. S.; Bell, J. B.;
    Zingale, M., 2011, The Astrophysical Journal, Volume 728, Issue 2,
    article id. 118

    https://doi.org/10.1088/0004-637X/728/2/118


## initial model

* `inputs_2d.rprox`

  This initial model comes from the original MAESTRO.  It was generated
  with `initial_models/toy_atm` using the
  `_params.xrb_mixed.hi_dens.tall.CNO` inputs file.

  This is intended to be used with the rprox network.
