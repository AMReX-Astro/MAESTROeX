***************************************************
Modifications for a Spherical Self-Gravitating Star
***************************************************

In papers II and III, we calculated the hydrostatic expansion of the base state
in plane-parallel geometry under the assumption that the weight of the
material above (or below) any given fluid parcel does not change
during hydrostatic expansion. This assumption holds when the
gravitational acceleration is independent of location. Here we discuss the
modifications to the algorithm in paper III required to treat a spherical
self-gravitating star.

One-dimensional Results
-----------------------

To test the spherical base state expansions, we inject heat at a
steady rate into a one-dimensional white dwarf model. This is similar
to the first test in paper II, except now in spherical coordinates.
As in that test, the compressible method with which we compare the low Mach number method
is the FLASH code's implementation of the
piecewise-parabolic method (PPM) in a one-dimensional spherical geometry.
The initial conditions for the white dwarf are those described in
Section 4.1 of paper III for the central region.

In the expansion of a plane-parallel atmosphere, heating at a
height :math:`r` above the base does not affect the pressure or density
below that height. By contrast, in a spherical symmetric
self-gravitating star, heating at a radius :math:`r` will lead to a pressure
and density decrease at the center in addition to the expansion of the
outer layers (see Schwarzchild & Harm, 1965, ApJ, 146, 855).

We apply a heating function of the form:

.. math::
    \Hext = H_0 \exp \left [-(r-r_0)^2 / W^2 \right ] ,

with :math:`r_0 = 4\times 10^7` cm, :math:`W = 10^7` cm, and :math:`H_0 = 1\times 10^{16}` 
erg g :math:`^{-1}` s :math:`^{-1}`. This is the same functional form as used
in the first test of paper II, but with a lower amplitude. Still, this
heating rate is far higher than what is expected during the convective
phase of Type Ia SNe. The heating term is added to the enthalpy
equation in the low Mach number equations in the same fashion as
described in paper II. In this test, we do not consider reactions.
Since this is a one-dimensional test all perturbational quantities,
as well as :math:`\Ubt`, are zero, so we are directly testing the computation of
:math:`w_0` as and the base state update as described in
the **Advect Base** procedure defined above.
Both the PPM and low Mach calculation use 768 zones in a domain :math:`5\times 10^8` cm high.

Figure [\[fig:spherical768\]](#fig:spherical768){reference-type="ref" reference="fig:spherical768"} shows the structure of the star after
heating for 10 s. The gray line is the initial star before any
heating.
We see that the compressible and low Mach number models
agree extremely well. Both capture the decrease in the density and
pressure at the center of the star and the considerable expansion in
radius. Only at the surface of the star do the temperatures differ slightly.
In all calculations, we set the minimum temperature to :math:`5\times 10^6` K. The PPM simulation required 13488 steps and the low Mach
(CFL :math:`=0.5`) calculation needed 203. Over the course of the
simulation, the Mach number of the flow remained less than :math:`0.35`, with the
maximum Mach number occurring at the surface of the star. This Mach
number pushes the limits of validity of the low Mach number model; a
smaller perturbation amplitude would result in a smaller Mach number.

Future improvements to the overall spherical base-state adjustment
algorithm will address the expansion in a simulation where the medium
outside the star is not brought down to arbitrarily low densities, but
instead a "cutoff density" is applied, as in the case of the
plane-parallel results presented in this paper. However, we expect
the changes to the overall method shown here to be small.

 Add a figure showing that we retain the correct solution
even when we place higher density material outside the star.

\clearpage
![image](\sphericalfigpath/spherical_adjust_768){width="5.0in"}
