***********************
Notes on Thermodynamics
***********************

Derivatives With Respect to Composition
=======================================

In the following we assume that the molar mass of species :math:`i` is given by its
atomic mass number, :math:`A_i = N_i + Z_i` where :math:`N_i` is the number of neutrons
and :math:`Z_i` is the number of protons in the isotope. This
is a slight approximation that ignores the mass difference between protons
and neutrons as well as some minor binding energy terms.

The number density [cm:math:`^-{3}`] of isotope :math:`i` is can be formed from the mass
density and the molar mass of that isotope as follows:

.. math::

   \label{eq:number density}
     n_i = \frac{\rho_i N_\text{A}}{A_i},

where :math:`N_\text{A}` is Avogadro’s number [# / mole]. The molar abundance,
:math:`Y_i`, is a measure of the number of moles of species :math:`i` per gram in the
system:

.. math::

   \label{eq:molar abundance}
     Y_i = \frac{n_i}{\rho N_\text{A}} = \frac{\rho_i}{\rho}\frac{1}{A_i}
     \equiv \frac{X_i}{A_i}

where we have defined the mass fraction, :math:`X_i = \frac{\rho_i}{\rho}`. Note

.. math::
   \sum_i X_i = 1.
   :label: eq:mass fraction sums to 1

We write the average molar mass and average proton number as:

.. math::
   \bar{A} = \frac{\sum_i A_i n_i}{\sum_i n_i} = \left(\sum_i X_i\right)
     \left(\sum_i \frac{X_i}{A_i}\right)^{-1}
   :label: eq:abar

.. math::
   \bar{Z} = \frac{\sum_i Z_i n_i}{\sum_i n_i} = \left(\sum_i Z_i
     \frac{X_i}{A_i}\right)\left(\sum_i \frac{X_i}{A_i}\right)^{-1}.
   :label: eq:zbar

Our algorithm requires terms involving the derivative thermodynamic variables
(:math:`p` or :math:`e`, e.g.) with respect to composition. Our EOS does not return such
derivatives but instead returns derivatives of these variables with respect
to :math:`\bar{A}` and :math:`\bar{Z}`. Using the chain rule, we have

.. math::

   \label{eq:p_xk}
     \frac{\ptl p}{\ptl X_k} = p_{X_k} =
     \frac{\ptl p}{\ptl \bar{A}}\frac{\ptl \bar{A}}{\ptl X_k} +
     \frac{\ptl p}{\ptl \bar{Z}}\frac{\ptl \bar{Z}}{\ptl X_k}.

From :eq:`eq:abar` and :eq:`eq:zbar` we have

.. math::
   \frac{\ptl \bar{A}}{\ptl X_k} = \left(\sum_i \frac{X_i}{A_i}\right)^{-1}
     - \frac{\bar{A}^2}{A_k} = -\frac{\bar{A}}{A_k}\left(\bar{A} - A_k\right)
   :label: eq:abar_X_k

.. math::
   \frac{\ptl \bar{Z}}{\ptl X_k} =
     \left(\frac{Z_k}{A_k}\right)\left(\sum_i \frac{X_i}{A_i}\right)^{-1}
     - \frac{\bar{Z}}{A_k}\left(\sum_i \frac{X_i}{A_i}\right)^{-1} =
     -\frac{\bar{A}}{A_k}\left(\bar{Z} - Z_k\right),
   :label: eq:zbar_X_k

where after differentiation we have used :eq:`eq:mass fraction sums to 1`
to write

.. math:: \left(\sum_i \frac{X_i}{A_i}\right)^{-1} = \bar{A}.

We therefore have

.. math::

   \label{eq:p_Xk_full}
     p_{X_k} = -\frac{\bar{A}}{A_k}\left(\bar{A} - A_k\right)
     \frac{\ptl p}{\ptl\bar{A}} - \frac{\bar{A}}{A_k}\left(\bar{Z} - Z_k\right)
     \frac{\ptl p}{\ptl\bar{Z}}.

Before it was brought to our attention by Frank Timmes, we were missing the
second term in :eq:`eq:abar_X_k`. The only place where such terms
appear in our algorithm is in a sum over all species, such as:

.. math::
   \sum_i p_{X_i}\dot{\omega}_i =
     -\bar{A}^2\frac{\ptl p}{\ptl\bar{A}}\sum_i \frac{\dot{\omega}_i}{A_i}
     +\bar{A}\frac{\ptl p}{\ptl\bar{A}}\sum_i \dot{\omega}_i
     -\bar{A}\bar{Z}\frac{\ptl p}{\ptl\bar{Z}}\sum_i \frac{\dot{\omega}_i}{A_i}
     +\bar{A}\frac{\ptl p}{\ptl \bar{Z}}\sum_i\frac{Z_i}{A_i}\dot{\omega}_i.
   :label: eq:sum over species

The second term in :eq:`eq:sum over species` is identically zero because

.. math:: \sum_k \dot{\omega}_k \equiv 0.

This second term arises from what was added to :eq:`eq:abar_X_k` by
Frank’s correction. Therefore, although important for individual derivatives
with respect to composition, this correction term has no effect on our
solution.

Convective stability criterion
==============================

Here we look at the criterion for convective stability in the case of
non-uniform chemical composition. This section follows Cox & Giuli
:raw-latex:`\cite{cg-ed2}` closely (see chapter 13).

Consider a fluid parcel that gets displaced upwards (against gravity) from
a radial location :math:`r` to :math:`r + \Delta r`.
The parcel is stable to convection if the displaced parcel’s density is
greater than
the surrounding fluid and gravity pushes the parcel back towards where it came
from. Then the criterion for stability should be

.. math::

   \begin{aligned}
    \rho_{parcel}(r+\Delta r) - \rho_{background}(r + \Delta r) &>& 0 \\
    \bigg[\rho_{parcel}(r) + \bigg(\frac{d\rho}{dr}\bigg)_{parcel}\Delta r\bigg] -
    \bigg[\rho_{background}(r) + \bigg(\frac{d\rho}{dr}\bigg)_{background}\Delta r\bigg] &>& 0 \end{aligned}

Since the parcel originates at r, :math:`\rho_{parcel}(r) = \rho_{background}(r)` and
so the stability criterion is

.. math::
   \bigg(\frac{d\rho}{dr}\bigg)_{parcel} > \bigg(\frac{d\rho}{dr}\bigg)_{background}
   :label: eqn:basicStability

Since the total pressure, :math:`P`, always increases inward in a star in hydrostatic
equilibrium, we can use :math:`P` instead of :math:`r` as the independent radial variable.
Then condition for stability can be written as

.. math:: \bigg( \frac{d \ln \rho}{d \ln P}\bigg )_{parcel} < \bigg(\frac{d \ln \rho}{d \ln P}\bigg)_{background}

Using the equation of state :math:`P = P( \rho, T, \bar{\mu})`, where
:math:`\bar{\mu}` is the average mass per molecule, we can write

.. math::
   d \ln P = \frac{\partial \ln P}{\partial \ln \rho} \bigg |_{T, \bar{\mu}}d \ln \rho + \frac{\partial \ln P}{\partial \ln T} \bigg |_{\rho, \bar{\mu}} d \ln T + \frac{\partial \ln P}{\partial \ln \bar{\mu}}\bigg |_{\rho, T} d \ln \bar{\mu}\
   :label: eqn:lnEOS

For convenience we introduce

.. math::

   \chi_{\rho} = \frac{\partial \ln P}{\partial \ln \rho}\bigg |_{T,\bar{\mu}} \qquad
     \chi_T = \frac{\partial \ln P}{\partial \ln T} \bigg |_{\rho,\bar{\mu}} \qquad
     \chi_{\bar{\mu}} = \frac{\partial \ln P}{\partial \ln \bar{\mu}} \bigg |_{\rho, T}

Then we can rearrange :eq:`eqn:lnEOS` to get

.. math::

   \frac{d \ln \rho}{\partial \ln P} = \frac{1}{\chi_\rho} -
     \frac{\chi_T}{\chi_\rho} \frac{d \ln T}{d \ln P}- \frac{\chi_{\bar{\mu}}}{\chi_\rho}
     \frac{d \ln \bar{\mu}}{d \ln P}

Then the general stability criterion is

.. math::
   \bigg ( \frac{1}{\chi_\rho} -
     \frac{\chi_T}{\chi_\rho} \frac{d \ln T}{d \ln P}- \frac{\chi_{\bar{\mu}}}{\chi_\rho}
     \frac{d \ln \bar{\mu}}{d \ln P} \bigg )_{parcel} <
     \bigg ( \frac{1}{\chi_\rho} -
     \frac{\chi_T}{\chi_\rho} \frac{d \ln T}{d \ln P}- \frac{\chi_{\bar{\mu}}}{\chi_\rho}
     \frac{d \ln \bar{\mu}}{d \ln P} \bigg )_{background}
   :label: eqn:genStability

Here’s where various assumptions/simplifications get used.

#. If no assumptions are made, you can’t get any further than
   :eq:`eqn:genStability`. Even in view of an infinitesimally small
   initial perturbation, you can’t, in general, assume the
   :math:`\chi`\ ’s in parcel are the same as the :math:`\chi`\ ’s in
   the background.  This applies in the case where nuclear reactions
   and/or ionization change the composition of the parcel. This case
   tends not to be of much interest for two reasons. Either
   composition effects get incorporated implicitly through assuming
   chemical equilibrium. Or both of these terms can be neglected in
   the rising parcel. This would be justified if the timescale for
   reactions is long compared to the convective timescale, and either
   the same is true for ionization or the fluid is fully ionized.

#. If we assume that :math:`\bar{\mu}` remains constant in the parcel, then
   :math:`\frac{d \ln \bar{\mu}}{d \ln P}` drops out for the parcel. In this case,
   we can assume, in view of the arbitrarily small initial perturbation of
   the parcel, that :math:`\chi_\rho` and :math:`\chi_T` to have the same values in the
   parcel as in the background. Then the stability criterion becomes

   .. math::

      \bigg (  \frac{d \ln T}{d \ln P} \bigg )_{parcel} >
        \bigg (  \frac{d \ln T}{d \ln P} + \frac{\chi_{\bar{\mu}}}{\chi_T}
        \frac{d \ln \bar{\mu}}{d \ln P} \bigg )_{background}
      \label{eqn:Ledoux}

   The Ledoux stability criterion is obtained by assuming that the parcel moves
   adiabatically.

#. If we assume
   that the background is in chemical equilibrium and the parcel achieves
   instantaneous chemical equilibrium, then :math:`\bar{\mu} = \bar{\mu}(\rho,T)` for
   the background and the parcel. (Note that we aren’t requiring constant
   composition in the parcel here.)
   The effect of variable composition are then absorbed into :math:`\chi_\rho` and
   :math:`\chi_T`. Again, we can take :math:`\chi_\rho` and :math:`\chi_T` to have the same values
   in the parcel as in the background. The criterion then is

   .. math::

      \bigg ( \frac{d \ln T}{d \ln P} \bigg )_{parcel} >
        \bigg ( \frac{d \ln T}{d \ln P} \bigg )_{background}
      \label{eqn:Schwarz}

   We obtain the Schwarzchild criterion for
   stability if we also assume the parcel moves adiabatically.

   The Scharwzchild criterion can be recast in terms of entropy if
   the EOS is taken as :math:`P(\rho, S)` instead of :math:`P(\rho, T)`. Then, in place
   of :eq:`eqn:lnEOS` we have

   .. math:: d \rho = \frac{\partial \rho}{\partial P} \bigg |_{S} d P + \frac{\partial \rho}{\partial S} \bigg |_{P} dS

   We can substitute this into :eq:`eqn:basicStability` for stability,
   and assuming the parcel moves adiabatically, we get

   .. math::

      \bigg ( \frac{\partial \rho}{\partial S} \bigg |_{P} \frac{dS}{dr}
        \bigg )_{background}< 0

   One of Maxwell’s relations is

   .. math:: \frac{\partial \rho^{-1}}{\partial S} \bigg |_{P} = \frac{\partial T}{\partial P} \bigg |_{S}

   All thermodynamically stable substances have temperatures that increase upon
   adiabatic compression, i.e. :math:`\frac{\partial T}{\partial P} \big |_{S} > 0`.
   So Maxwell’s relation implies that
   :math:`\frac{\partial \rho}{\partial S} \big |_{P} < 0`.
   The stability criterion then becomes

   .. math::
      \bigg ( \frac{d S}{d r} \bigg )_{background} > 0
      :label: eqn:stabilityEntr

Determining which stability criterion we want to enforce in creating the
initial model is complicated by the phenomenon of semiconvection, which
occurs when the Ledoux criterion is satisfied but the Schwarzchild is not,
i.e.

.. math::

   \bigg (  \frac{d \ln T}{d \ln P} \bigg )_{parcel} <
     \bigg (  \frac{d \ln T}{d \ln P} \bigg )_{background} <
     \bigg (  \frac{d \ln T}{d \ln P} \bigg )_{parcel} -
     \bigg ( \frac{\chi_{\bar{\mu}}}{\chi_T}
     \frac{d \ln \bar{\mu}}{d \ln P} \bigg )_{background}

(Note that :math:`\chi_{\bar{\mu}}` is negative, as pressure is inversely proportional
to mass per particle, and :math:`\frac{d \ln \bar{\mu}}{d \ln P}` is positive, since
nuclear reactions synthesize more massive particles in the center of the star.)
In this case, when a rising parcel eventually reaches neutral buoyancy, it will
have a temperature excess in comparison to it’s surroundings.
If the parcel can retain
it’s identity against diffusive mixing with the background long enough for
significant heat exchange to occur, then the parcel’s temperature will drop, it
will contract increasing it’s density, and the parcel will move inwards.
The time scale of semiconvection is much longer than the timescale of
traditional convection.

When we set up an initial model, we want to minimize any initial tendency
towards convective motions, as we want these to be driven by the heating due
to nuclear reactions,
not the initial configuration we supply. Thus I think we want to guard against
semiconvection as well as “traditional” convection by using the stability
criterion

.. math::

   \bigg ( \frac{d \ln T}{d \ln P} \bigg )_{parcel} =
   \frac{d \ln T}{d \ln P} \bigg |_{S,\bar{\mu}} >
     \bigg ( \frac{d \ln T}{d \ln P} \bigg )_{background}

Although this looks like the Schwarschild criterion (and, because I’m not
entirely sure on vocabulary, it might even be called the Schwarzchild
criterion), this does not simplify to :eq:`eqn:stabilityEntr`
because we need to keep the explicit :math:`\bar{\mu}` dependence in the EOS.

The question of whether we’re in chemical equilibrium or not might be a moot
point since our EOS (or any other part of the code) doesn’t enforce chemical
equilibrium. Thus, even
in the case of chemical equilbrium, we can’t in general
drop the explicit :math:`\bar{\mu}`
dependence from our equations. If we wanted to do that, then we would need
:math:`\bar{\mu}(\rho,T)` to be substituted for :math:`\bar{\mu}` inside the EOS.

.. _Sec:Adiabatic Excess:

Adiabatic Excess
================

The adiabatic excess, :math:`\Delta\nabla`, is a quantity used to determine
if a system is stable (:math:`\Delta\nabla < 0`) or unstable (:math:`\Delta\nabla
> 0`) to convection under the Schwarzschild criterion (i.e. neglecting
compositional gradients). Cox and Giuli (see chapter 9) define three
different “adiabatic exponents” that we will use:

.. math::

   \begin{aligned}
     \Gamma_1 &\equiv&   \left(\frac{d\ln p}{d\ln\rho}\right)_\text{ad} \\
     \frac{\Gamma_2}{\Gamma_2-1} &\equiv&
     \left(\frac{d\ln p}{d\ln T}\right)_\text{ad} \\
     \Gamma_3 - 1 &\equiv& \left(\frac{d\ln T}{d\ln\rho}\right)_\text{ad},\end{aligned}

where the subscript “ad” means along an adiabat. We can combine the
exponents to get the following relation

.. math::
   \Gamma_1 = \left(\frac{\Gamma_2}{\Gamma_2-1}\right)\left(\Gamma_3-1\right).
   :label: eq:Gamma relations


The adiabatic excess is defined as

.. math::

   \label{eq:adiabatic excess}
     \Delta\nabla = \nabla_\text{actual} - \nabla_\text{ad}

where

.. math::

   \label{eq:thermal gradient}
     \nabla \equiv \frac{d\ln T}{d\ln P}

is the thermal gradient. It is important to note that these thermal
gradients are only along the radial direction. The “actual”
gradient can be found from finite differencing the data whereas the
adiabatic term, :math:`\nabla_\text{ad} = \left(\Gamma_2-1\right) /
\Gamma_2`, will need to be calculated at each point using
thermodynamic relations. Our EOS only returns :math:`\Gamma_1` so we need
find another relation to use with :eq:`eq:Gamma relations` to solve
for the adiabatic excess.

The Schwarzschild criterion does not care about changes in composition
and we therefore write :math:`p = p(\rho,T)` and

.. math::
   d\ln p = \chi_\rho d\ln\rho + \chi_T d\ln T
   :label: eq:dp

where

.. math::

   \chi_\rho = \left(\frac{d\ln p}{d\ln\rho}\right)_T,\qquad
   \chi_T = \left(\frac{d\ln p}{d\ln T}\right)_\rho.

Dividing :eq:`eq:dp` by :math:`d\ln\rho` and taking this along an adiabat
we have

.. math::
   \left(\frac{d\ln p}{d\ln\rho}\right)_\text{ad} = \chi_\rho + \chi_T
     \left(\frac{d\ln T}{d\ln\rho}\right)_\text{ad}.
   :label: eq:dp2

Using the :math:`\Gamma`\ ’s, we have

.. math::
   \Gamma_1 = \chi_\rho + \chi_T\left(\Gamma_3-1\right).
   :label: eq:Gamma1 relation with Gamma2

Combining :eq:`eq:Gamma relations` and :eq:`eq:Gamma1 relation with Gamma2`
to eliminate :math:`\Gamma_3`, we have:

.. math::

   \label{eq:nabla_ad}
     \nabla_\text{ad} = \frac{\Gamma_1 - \chi_\rho}{\chi_T\Gamma_1}

which uses only terms which are easily returned from an EOS call.
