*****************
Reaction Networks
*****************

Introduction to MAESTROeX Networks
==================================

MAESTROeX models multiple species, described by the mass density of
the fluid, :math:`\rho`, and the mass fraction of the species, :math:`X_k \equiv
\rho_k/\rho`, where :math:`\rho_k` is the mass density of species :math:`k`. All
MAESTROeX problems, regardless of whether they model reactions, need a
network. In its most basic form, the network supplies the properties
of the species (atomic mass, atomic number) that are interpreted by
the equation of state to compute.

Additional networks will be made available in the Microphysics repo [1]_. These
will have interfaces for both MAESTROeX and CASTRO.

Notes of Specific Networks
==========================

general_null
------------

This is a ’null’ network – i.e. no burning, just define the
properties of the species for thermodynamics. The twist is that you
can create an inputs file to define what species you want to carry.
For example, the extern/networks/null/ network defines C12, O16, and
Mg24. To replicate this in general_null, we have the file
ignition.net with contents:

::

    # name       short name    aion     zion
     carbon-12      C12         12.0     6.0
     oxygen-16      O16         16.0     8.0
     magnesium-24   Mg24        24.0    12.0

To use this set of species in your problem, you would set:

::

    NETWORK_DIR := extern/networks/general_null
    GENERAL_NET_INPUTS := ignition.net

It is assumed that the \*.net files live in extern/networks/general_null/

Then at compile time, the network.f90 is created using these species and
compiled. (For the curious, the rule to build network.f90 lives in
extern/networks/general_null/GPackage.mak)

ignition_chamulak
-----------------

This network was introduced in our paper on convection in white dwarfs
as a model of Type Ia supernovae :cite:`wdconvect`. It models
carbon burning in a regime appropriate for a simmering white dwarf,
and captures the effects of a much larger network by setting the ash
state and energetics to the values suggested in :cite:`chamulak:2008`.

ignition_simple
---------------

This is the original network used in our white dwarf convection
studies :cite:`lowMach4`. It includes a single reaction,
:math:`^{12}\mathrm{C}({}^{12}\mathrm{C},\gamma){}^{24}\mathrm{Mg}`, using
the rate from Caughlin and Fowler :cite:`caughlan-fowler:1988`.

rprox
-----

This is a network introduced in a paper modeling mixed H/He X-ray
bursts :cite:`xrb2`. rprox that has 10 species approximates hit
CNO burning, triple-\ :math:`\alpha`, and rp-process breakout up through
:math:`^{56}\mathrm{Ni}`. Updated rates from ReacLib :cite:`ReacLib` are
used. The overall ideas in this network are based on Appendix C of
:cite:`wallacewoosley:1981`.

.. [1]
   https://github.com/AMReX-Astro/Microphysics
