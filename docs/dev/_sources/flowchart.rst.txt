.. _ch:flowchart:

*******************
Governing Equations
*******************

The equation set and solution procedure used by MAESTROeX has changed
and improved over time.  In this chapter, we outline the model
equations and algorithmic options in the code.  The latest published
references for MAESTROeX are the multilevel paper :cite:`multilevel`
and the more recent :cite:`MAESTROeX`.  These two papers use the same
model equations, however the more recent paper, in addition to
retaining the original algorithmic capability of the previous code,
includes an option to use a new, simplified temporal integration
scheme.  We distinguish between the two temporal integration
strategies by referring to them as the "original temporal scheme" and
"new temporal scheme".  In this description, we make frequent
reference to papers I-IV and the multilevel paper (see §
:ref:`ch:intro`), which describe the developments of the original
temporal scheme.

Summary of the MAESTROeX Equation Set
=====================================

Here we summarize the equations solved by MAESTROeX. We refer the reader
to papers I through IV for the derivation
and motivation of the equation set.
We take the standard equations of reacting, compressible flow, and recast the
equation of state (EOS) as a divergence constraint on the velocity field.
The resulting model is a series of evolution equations for mass, energy,
and momentum, subject to an additional constraint on velocity, with a base state
density and pressure linked via hydrostatic equilibrium:

.. math::
   \frac{\partial \rho X_k}{\partial t} + \nabla \cdot (\rho \Ub X_k) = \rho \omegadot_k

.. math::
   \frac{\partial(\rho h)}{\partial t} =
      -\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt} + \rho\Hnuc + \rho\Hext,

.. math::
   \frac{\partial \Ub}{\partial t} + \Ub \cdot \nabla \Ub +
      \frac{\beta_0}{\rho} \nabla \left (\frac{p^\prime}{\beta_0} \right ) =
      -\frac{\rho^\prime}{\rho} |g| \er

.. math::
   \nabla \cdot (\beta_0 \Ub) =
      \beta_0 \left ( S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} \right )

.. math:: \nabla p_0 = -\rho_0 |g| \er

We discuss each of these equations in further detail below.


Lateral Average
---------------

A key concept in the MAESTROeX equation set and algorithm is
the lateral average.  The lateral average represents the average
value of a quantity at a given radius in spherical simulations
(or a given height in planar simulations).  We denote the
lateral average of a quantity with an overline, e.g.,
for any quantity :math:`\phi`, we denote
the average of :math:`\phi` over a layer at constant radius
as :math:`\overline{\phi}`.  For planar problems this routine is
a trivial average of all the values at a given height.
For spherical problems there is a
novel interpolation routine we use to average 3D data representing
a full spherical star into a 1D array representing the average.
Details can be found in :cite:`multilevel` and :cite:`MAESTROeX`.

For the velocity field, we can decompose the full velocity
field into a base state velocity and a local velocity,

.. math:: \Ub = w_0(r,t)\eb_r + \Ubt(\xb,t).

where :math:`r` is a 1D radial coordinate,
:math:`\xb` is a 3D Cartesian grid coordinate, and
:math:`\eb_r` is the unit vector in the outward radial direction.
Note that :math:`\overline{(\Ubt\cdot\eb_r)} = 0` and
:math:`w_0 = \overline{(\Ub\cdot\eb_r)}`.
In other words, the base state velocity can be thought of as the
lateral average of the outward radial velocity.
For the velocity decompsotion, we do not use the
same spatial averaging operators used
for all other variables; instead we derive an analytic expression
for the average expansion velocity and numerically integrate
this expression to obtain :math:`w_0`.

Mass
----

Conservation of mass gives the same continuity equation we have with
compressible flow:

.. math::

   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \Ub) = 0,
   \label{eq:flow:continuity}

where :math:`\rho` is the total mass density.
Additionally, we model the evolution of individual species that advect and react.
The creation and destruction
of the species is described by their create rate, :math:`\omegadot_k`, and the species
are defined by their mass fractions, :math:`X_k \equiv \rho_k / \rho`, giving

.. math::
   \frac{\partial \rho X_k}{\partial t} + \nabla \cdot (\rho \Ub X_k) = \rho \omegadot_k
   :label: eq:flow:rhoX

and

.. math:: \sum_k X_k = 1

In the original temporal scheme, we need to model the evolution of a base state density,
:math:`\rho_0`.  The governing equation can be obtained by laterally averaging the
full continuity equation, giving:

.. math::
   \frac{\partial\rho_0}{\partial t} = -\nabla\cdot(\rho_0 w_0 \eb_r),
   :label: eq:flow:base_density

Subtracting these two yields the evolution equation for the perturbational
density, :math:`\rho^\prime \equiv \rho - \rho_0`:

.. math::

   \frac{\partial\rho'}{\partial t} = -\Ub\cdot\nabla\rho' -
     \rho'\nabla\cdot\Ub - \nabla\cdot\left(\rho_0\Ubt\right) \label{eq:flow:rhoprime}

As first discussed in paper III and then refined in the multilevel paper, we capture
the changes that can occur due to significant convective
overturning by imposing the constraint that :math:`\overline{\rho'}=0`
for all time. This gives

.. math:: \frac{\partial\overline{\rho'}}{\partial t} = -\nabla\cdot(\etarho\eb_r).

where

.. math::
   \etarho = \overline{\left(\rho'\Ub\cdot\eb_r\right)}
   :label: eq:flow:etarho

In practice, we correct the drift by simply setting :math:`\rho_0 =
\overline{\rho}` after the advective update of :math:`\rho`. However we still need to
explicitly compute :math:`\etarho` since it appears in other equations.

Energy
------

We model the evolution of specific enthalpy, :math:`h`.
Strictly speaking this is not necessary to close the system,
but a user can enable the option to couple the energy with the
rest of the system by using the enthalpy to define the temperature.
The advantages of this coupling is an area of active research.
The evolution equation is

.. math::
   \frac{\partial(\rho h)}{\partial t} =
      -\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt} + \rho\Hnuc + \rho\Hext,
   :label: eq:flow:enthalpy

where :math:`p_0` is the 1D base state pressure, :math:`\Hnuc` and :math:`\Hext`
are energy sources due to reactions and user-defined external heating.

When we are using thermal diffusion, there will be an additional term in
the enthalpy equation (see § :ref:`sec:flow:diffusion`).

In the original temporal scheme, we utilized a base state enthlpy that effectively
represents the average over a layer; its evolution equation can be
found by laterally averaging :eq:`eq:flow:enthalpy`

.. math::
   \frac{\partial(\rho h)_0}{\partial t} = -\nabla\cdot\left[(\rho h)_0w_0\eb_r\right] +
     \psi + \overline{\rho \Hnuc} + \overline{\rho \Hext}.
   :label: eq:flow:enthalpy_base

We will often expand :math:`Dp_0/Dt` as

.. math:: \frac{Dp_0}{Dt} = \psi + (\Ubt \cdot \er) \frac{\partial p_0}{\partial r}

where we defined

.. math:: \psi \equiv \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial r}

In paper III, we showed that for a plane-parallel atmosphere with
constant gravity, :math:`\psi = \etarho g`

At times, we will define a temperature equation by writing :math:`h = h(T,p,X_k)`
and differentiating:

.. math::
   \frac{DT}{Dt} = \frac{1}{\rho c_p} \left\{ \left(1 - \rho h_p\right) \left
     [ \psi + (\Ubt \cdotb \er) \frac{\partial p_0}{\partial r} \right ]
    - \sum_k \rho \xi_k {\omegadot}_k
    + \rho \Hnuc + \rho \Hext \right \}   .
   :label: eq:flow:temp

Subtracting it from the full enthalpy equation gives:

.. math::
   \begin{align}
   \frac{\partial(\rho h)'}{\partial t} = &-\Ub\cdot\nabla(\rho h)' - (\rho h)'\nabla\cdot\Ub -
     \nabla\cdot\left[(\rho h)_0\Ubt\right] + \nonumber \\
   &\Ubt\cdot\nabla p_0
      + ( \rho\Hnuc - \overline{\rho \Hnuc}) + (\rho\Hext - \overline{\rho \Hext})
   \end{align}
   :label: eq:flow:rhohprime

Base State
----------

The stratified atmosphere is characterized by a one-dimensional
time-dependent base state, defined by a base state density, :math:`\rho_0`,
and a base state pressure, :math:`p_0`, in hydrostatic equilibrium:

.. math:: \nabla p_0 = -\rho_0 |g| \er

The gravitational acceleration, :math:`g` is either constant or a
point-mass with a :math:`1/r^2` dependence (see §
:ref:`sec:planarinvsqgravity`) for plane-parallel geometries, or a
monopole constructed by integrating the base state density for
spherical geometries.

For the time-dependence, we will define a base state velocity, :math:`w_0`,
which will adjust the base state from one hydrostatic equilibrium to
another in response to heating.

For convenience, we define a base state enthalphy, :math:`h_0`, as needed
by laterally averaging the full enthalpy, :math:`h`.

Base State Expansion
--------------------

In practice, we calculate :math:`w_0` by integrating
the one-dimensional divergence constraint. For a plane-parallel atmosphere, the
evolution is:

.. math::
   \frac{\partial w_0}{\partial r} = \Sbar - \frac{1}{\gammabar p_0} \etarho g
   :label: eq:flow:dw0dr_planar

Then we define

.. math::
   - \frac{\beta_0}{\rho_0} \frac{\partial (\pizero/\beta_0)}{\partial r} = \frac{\partial w_0}{\partial t} +
      w_0 \frac{\partial w_0}{\partial r} ,
   :label: eq:pizero

once :math:`w_0` at the old and new times is known, and the advective term is computed explicitly.
Then we can include this for completeness in the update for :math:`\ut.`


Momentum
--------

The compressible momentum equation (written in terms of velocity is):

.. math:: \rho \frac{\partial \Ub}{\partial t} + \rho \Ub \cdot \nabla \Ub + \nabla p = -\rho |g| \er

Subtracting off the equation of hydrostatic equilibrium, and defining the perturbational
pressure (sometimes called the dynamic pressure) as :math:`\pi \equiv p - p_0`,
and perturbational density as :math:`\rho' \equiv \rho - \rho_0`, we have:

.. math:: \rho \frac{\partial \Ub}{\partial t} + \rho \Ub \cdot \nabla \Ub + \nabla \pi = -\rho' |g| \er

or

.. math::

   \frac{\partial \Ub}{\partial t} + \Ub \cdot \nabla \Ub + \frac{1}{\rho} \nabla \pi =
      -\frac{\rho^\prime}{\rho} |g| \er

This is the form of the momentum equation that we solved in papers
I–IV and in the multilevel paper.

Several authors :cite:`KP:2012,VLBWZ:2013` explored the idea of energy
conservation in a low Mach number system and found that an additional
term (which can look like a buoyancy) is needed in the low Mach number
formulation, yielding:

.. math::
   \frac{\partial \Ub}{\partial t} + \Ub \cdot \nabla \Ub +
      \frac{\beta_0}{\rho} \nabla \left (\frac{p^\prime}{\beta_0} \right ) =
      -\frac{\rho^\prime}{\rho} |g| \er
   :label: eq:flow:newmomentum

This is the form that we enforce in MAESTROeX, and the choice is controlled
by ``use_alt_energy_fix``.

We decompose the full velocity field into a base state velocity,
:math:`w_0`, that governs the base state dynamics, and a local velocity,
:math:`\Ubt`, that governs the local dynamics, i.e.,

.. math:: \Ub = w_0(r,t)\eb_r + \Ubt(\xb,t).

with
:math:`\overline{(\Ubt\cdot\eb_r)} = 0` and
:math:`w_0 = \overline{(\Ub\cdot\eb_r)}`—the motivation for this splitting was given in paper II.
The velocity evolution equations are then

.. math::
   \frac{\partial w_0}{\partial t} = -w_0\frac{\partial w_0}{\partial
     r} - \frac{\beta_0}{\rho_0}\frac{\partial(\pi_0/\beta_0)}{\partial r}
   :label: eq:w0 evolution

.. math::
   \frac{\partial\Ubt}{\partial t} = -(\Ubt + w_0\er)\cdot\nabla\Ubt
     - \left(\Ubt\cdot\eb_r\right)\frac{\partial w_0}{\partial r}\eb_r -
   \frac{\beta_0}{\rho}\nabla\left(\frac{\pi}{\beta_0} \right) +
   \frac{\beta_0}{\rho_0}\frac{\partial(\pi_0/\beta_0)}{\partial r}\eb_r -
   \frac{\rho-\rho_0}{\rho}g\eb_r.
   :label: eq:flow:utildeupd

where :math:`\pi_0` is the base state component of the perturbational pressure.
By laterally averaging to :eq:`eq:U_divergence`,
we obtain a divergence constraint for :math:`w_0`:

.. math::
   \nabla\cdot(\beta_0 w_0 \eb_r) =
       \beta_0\left(\Sbar - \frac{1}{\gammabar p_0}
              \frac{\partial p_0}{\partial t}\right).
   :label: eq:w0 divergence

The divergence constraint for :math:`\Ubt` can be found by subtracting
:eq:`eq:w0 divergence` into :eq:`eq:U_divergence`, resulting in

.. math:: \nabla\cdot\left(\beta_0\Ubt\right) = \beta_0\left(S-\Sbar\right).\label{eq:utilde divergence}

Velocity Constraint
-------------------

The equation of state is cast into an elliptic constraint on the
velocity field by differentiating :math:`p_0(\rho, s, X_k)` along particle
paths, giving:

.. math::
   \nabla \cdot (\beta_0 \Ub) =
      \beta_0 \left ( S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} \right )
   :label: eq:U_divergence

where :math:`\beta_0` is a density-like variable that carries background
stratification, defined as

.. math:: \beta_0(r,t) = \rho_0(0,t)\exp\left(\int_0^r\frac{1}{\gammabar p_0}\frac{\partial p_0}{\partial r'}dr'\right),

and

.. math::
   S = -\sigma\sum_k\xi_k\omegadot_k + \frac{1}{\rho p_\rho}\sum_k p_{X_k}\omegadot_k + \sigma\Hnuc + \sigma\Hext + \frac{\sigma}{\rho} \nabla \cdot \kth \nabla T
   :label: eq:flow:S

where :math:`p_{X_k} \equiv \left. \partial p / \partial X_k
\right|_{\rho,T,X_{j,j\ne k}}`, :math:`\xi_k \equiv \left. \partial h /
\partial X_k \right |_{p,T,X_{j,j\ne k}},
p_\rho \equiv \left. \partial p/\partial \rho \right |_{T, X_k}`, and
:math:`\sigma \equiv p_T/(\rho c_p p_\rho)`, with
:math:`p_T \equiv \left. \partial p / \partial T \right|_{\rho, X_k}` and
:math:`c_p \equiv \left.  \partial h / \partial T
\right|_{p,X_k}` is the specific heat at constant pressure. The last
term is only present if we are using thermal diffusion (``use_thermal_diffusion = T``). In this term, :math:`\kth` is the thermal conductivity.

In this constraint, :math:`\gammabar` is the lateral average of
:math:`\Gamma_1 \equiv d\log p / d\log \rho |_s`. Using the lateral average
here makes it possible to cast the constraint as a
divergence. :cite:`KP:2012` discuss the general case where we want to
keep the local variations of :math:`\Gamma_1` (and we explored this in paper
III). We also look at this in § :ref:`sec_flow_gamma1vary`

Notation
========

Throughout the papers describing MAESTROeX, we’ve largely kept our
notation consistent. The table below defines the
frequently-used quantities and provides their units.

.. table:: Definition of symbols.

   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | symbol                | description                                                           | units                                |
   +=======================+=======================================================================+======================================+
   | :math:`c_p`           | specific heat at                                                      | erg g :math:`^{-1}` K :math:`^{-1}`  |
   |                       | constant pressure                                                     |                                      |
   |                       | (:math:`c_p \equiv \partial h / \partial T |_{p, X_k}`)               |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`f`             | volume discrepancy                                                    | –                                    |
   |                       | factor                                                                |                                      |
   |                       | (:math:`0 \le f \le 1`)                                               |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`g`             | gravitational                                                         | cm s  :math:`^{-2}`                  |
   |                       | acceleration                                                          |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`h`             | specific enthalpy                                                     | erg g :math:`^{-1}`                  |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\Hext`         | external heating                                                      | erg g :math:`^{-1}` s :math:`^{-1}`  |
   |                       | energy generation  rate                                               |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\Hnuc`         | nuclear energy                                                        | erg g :math:`^{-1}`                  |
   |                       | generation rate                                                       | s :math:`^{-1}`                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`h_p`           | :math:`h_p \equiv \partial h/\partial p |_{T,X_k}`                    | cm :math:`^{3}` g :math:`^{-1}`      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\kth`          | thermal conductivity                                                  | erg cm :math:`^{-1}`                 |
   |                       |                                                                       | s :math:`^{-1}` K :math:`^{-1}`      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`p_0`           | base state pressure                                                   | erg cm  :math:`^{-3}`                |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`p_T`           | :math:`p_T \equiv \partial p / \partial T |_{\rho,X_k}`               | erg cm :math:`^{-3}` K :math:`^{-1}` |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`p_{X_k}`       | :math:`p_{X_k}\equiv\partial p/\partial X_k|_{p,T,X_{j,j\ne k}}`      | erg cm :math:`^{-3}`                 |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`p_\rho`        | :math:`p_\rho \equiv \partial p/\partial \rho |_{T,X_k}`              | erg g :math:`^{-1}`                  |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`q_k`           | specific nuclear                                                      | erg g\ :math:`^{-1}`                 |
   |                       | binding energy                                                        |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`r`             | radial coordinate (direction of gravity)                              | cm                                   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`s`             | specific entropy                                                      | erg g :math:`^{-1}` K :math:`^{-1}`  |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`S`             | source term to the                                                    | s :math:`^{-1}`                      |
   |                       | divergence constraint                                                 |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`t`             | time                                                                  | s                                    |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`T`             | temperature                                                           | K                                    |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\Ub`           | total velocity                                                        | cm s :math:`^{-1}`                   |
   |                       | (:math:`\Ub = \Ubt + w_0 \eb_r`                                       |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\Ubt`          | local velocity                                                        | cm s :math:`^{-1}`                   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\uadv`         | advective velocity                                                    | cm s :math:`^{-1}`                   |
   |                       | (edge-centered)                                                       |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`w_0`           | base state expansion                                                  | cm s\ :math:`^{-1}`                  |
   |                       | velocity                                                              |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`X_k`           | mass fraction of the                                                  | –                                    |
   |                       | species                                                               |                                      |
   |                       | (:math:`\sum_k X_k = 1`)                                              |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\beta_0`       | coefficient to velocity in velocity constraint equation               | g cm :math:`^{-3}`                   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\Gamma_1`      | first adiabatic exponent                                              | –                                    |
   |                       | (:math:`\Gamma_1 \equiv d\log p/d\log \rho|_s`)                       |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\etarho`       | :math:`\etarho \equiv \overline{(\rho' \Ub \cdot \eb_r)}`             | g cm :math:`^{-2}` s :math:`^{-1}`   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\xi_k`         | :math:`\xi_k \equiv \partial h / \partial X_k |_{p,T,X_{j,j\ne k}}`   | erg g :math:`^{-1}`                  |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\pi`           | dynamic pressure                                                      | erg cm :math:`^{-3}`                 |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\pizero`       | base state dynamic pressure                                           | erg cm  :math:`^{-3}`                |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\rho`          | mass density                                                          | g cm :math:`^{-3}`                   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\rho_0`        | base state mass density                                               | g cm :math:`^{-3}`                   |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\rho'`         | perturbational density                                                | g cm :math:`^{-3}`                   |
   |                       | (:math:`\rho' = \rho - \rho_0`)                                       |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`(\rho h)_0`    | base state enthalpy density                                           | erg cm :math:`^{-3}`                 |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`(\rho h)'`     | perturbational enthalpy density                                       | erg cm :math:`^{-3}`                 |
   |                       | :math:`\left [(\rho h)' = \rho h - (\rho h)_0 \right ]`               |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\sigma`        | :math:`\sigma \equiv p_T/(\rho c_p p_\rho)`                           | erg :math:`^{-1}` g                  |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\psi`          | :math:`\psi \equiv D_0 p_0/Dt = \ptl p_0/\ptl t + w_0\ptl p_0/\ptl r` | erg cm :math:`^{-3}`                 |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+
   | :math:`\omegadot_k`   | creation rate for species :math:`k`                                   | s :math:`^{-1}`                      |
   |                       | (:math:`\omegadot_k \equiv DX_k/Dt`)                                  |                                      |
   +-----------------------+-----------------------------------------------------------------------+--------------------------------------+

.. [1]
   Here we see an unfortunate conflict
   of notation between the compressible hydro community and the
   incompressible community. In papers on compressible hydrodynamics,
   :math:`\Ub` will usually mean the vector of conserved quantities. In
   incompressible / low speed papers, :math:`\Ub` will mean the velocity vector.

.. _Sec:Time Advancement Algorithm:

Time Advancement Algorithm Ingredients
======================================

The full time advancement algorithm is detailed in :cite:`multilevel` for the original
algorithm and :cite:`MAESTROeX` for the new, simplified algorithm.  We do not repeat
that here.

The overall flow of the algorithm is depicted in the following flowcharts:

.. figure:: flowchart.png
   :align: center

   A flowchart of the algorithm. The thermodynamic state variables,
   base state variables, and local velocity are indicated in each
   step. Red text indicates that quantity was updated during that
   step. The predictor-corrector steps are outlined by the dotted
   box. The blue text indicates state variables that are the same in
   Step 6 as they are in Step 2, i.e., they are unchanged by the
   predictor steps. The diffusion steps (4a and 8a) are optional,
   depending on use_thermal_diffusion.

.. figure:: flowchart_4_8.png
   :align: center

   A flowchart for Steps 4 and 8. The thermodynamic state variables
   and base state variables are indicated in each step. Red text
   indicates that quantity was updated during that step. Note, for
   Step 4, the updated quantities should also have a :math:`⋆`
   superscript, e.g., Step 8I defines :math:`T^{(2)}` while Step 4I
   defines :math:`T^{(2),\star}`.



Here are some of the basic ingredients to the solver:


Definitions
-----------

Below we define operations that will be referenced in
§ :ref:`sec:flow:singlestep`.

**React State**\ :math:`[\rho^{\inp},(\rho h)^{\inp},X_k^{\inp},T^{\inp}, (\rho\Hext)^{\inp}, p_0^{\inp}] \rightarrow [\rho^{\outp}, (\rho h)^{\outp}, X_k^{\outp}, T^{\outp}, (\rho \omegadot_k)^{\outp}, (\rho\Hnuc)^{\outp}]`
evolves the species and enthalpy due to reactions through
:math:`\Delta t/2` according to:

.. math::

   \frac{dX_k}{dt} = \omegadot_k(\rho,X_k,T) ; \qquad
   \frac{dT}{dt}   = \frac{1}{c_p} \left ( -\sum_k \xi_k  \omegadot_k  + \Hnuc \right ).

Here the temperature equation comes from :eq:`eq:flow:temp` with :math:`Dp_0/Dt = 0` for
the burning part of the evolution.

Full details of the solution procedure can be found in Paper III. We
then define:

.. math::

   \begin{aligned}
   (\rho\omegadot_k)^{\outp} &=& \frac{\rho^{\outp} ( X_k^{\outp} - X_k^{\inp})}{\dt/2}, \\
   (\rho h)^{\outp} &=& (\rho h)^{\inp} + \frac{\dt}{2} (\rho\Hnuc)^{\outp} + \frac{\dt}{2} (\rho\Hext)^{\inp}.\end{aligned}

where the enthalpy update includes external heat sources :math:`(\rho\Hext)^{\inp}`.
As introduced in Paper IV, we update the temperature using :math:`T^{\outp} =
T(\rho^\outp,h^\outp,X_k^\outp)` for planar geometry or :math:`T^{\outp} =
T(\rho^\outp,p_0^\inp,X_k^\outp)` for spherical geometry, with this behavior
controlled by use_tfromp.
Note that the density remains unchanged within **React State**, i.e.,
:math:`\rho^{\outp} = \rho^{\inp}`.

The source code for this operation can be found in react_state.f90.

Advect Base Density
^^^^^^^^^^^^^^^^^^^

:math:`[\rho_0^\inp,w_0^\inp] \rightarrow [\rho_0^\outp, \rho_0^{\outp,\nph}]` is the process by which we
update the base state density through :math:`\dt` in time. We keep the
time-centered edge states, :math:`\rho_0^{\outp,\nph}`,
since they are used later in discretization of :math:`\etarho` for planar problems.

* planar:

  We discretize equation :eq:`eq:flow:base_density` to
  compute the new base state density,

  .. math:: \rho_{0,j}^{\outp} = \rho_{0,j}^{\inp} - \frac{\dt}{\dr} \left [ \left( \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j+\myhalf} - \left( \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j-\myhalf} \right ].

  We compute the time-centered edge states,
  :math:`{\rho_0}^{\outp,\nph}_{j\pm\myhalf}`, by discretizing an
  expanded form of :eq:`eq:flow:base_density`:

  .. math:: \frac{\partial \rho_0}{\partial t} + w_0 \frac{\partial \rho_0}{\partial r} = - \rho_0 \frac{\partial w_0}{\partial r},

  where the right hand side is used as the force term.

* spherical:

  The base state density update now includes the area factors in the
  divergences:

  .. math:: \rho_{0,j}^{\outp} = \rho_{0,j}^{\inp} - \frac{1}{r_j^2} \frac{\dt}{\dr} \left [ \left( r^2 \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j+\myhalf} - \left( r^2 \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j-\myhalf} \right].

  In order to compute the time-centered edge states, an additional
  geometric term is added to the forcing, due to the spherical
  discretization of :eq:`eq:flow:base_density`:

  .. math:: \frac{\partial \rho_0}{\partial t} + w_0 \frac{\partial \rho_0}{\partial r} = - \rho_0 \frac{\partial w_0}{\partial r} - \frac{2 \rho_0 w_0}{r}.



Enforce HSE
^^^^^^^^^^^

:math:`[p_0^{\inp},\rho_0^{\inp}] \rightarrow [p_0^{\outp}]` has replaced **Advect Base Pressure**
from Paper III as the process by which we update the base state
pressure. Rather than discretizing the evolution equation for
:math:`p_0`, we enforce hydrostatic equilibrium directly, which is numerically simpler
and analytically equivalent. We first set
:math:`p_{0,j=0}^{\outp} = p_{0,j=0}^{\inp}` and then update :math:`p_0^\outp` using:

.. math:: p_{0,j+1}^{\outp} = p_{0,j}^{\outp} + \Delta r g_{j+\myhalf}\frac{\left(\rho_{0,j+1}^{\inp}+\rho_{0,j}^{\inp}\right)}{2},

where :math:`g=g(\rho_0^\inp)`. As soon as :math:`\rho_{0,j}^\inp < \rho_{\rm cutoff}`, we set
:math:`p_{0,j+1}^\outp = p_{0,j}^\outp` for all remaining values of :math:`j`.
Then we compare :math:`p_{0,j_{\rm max}}^\outp` with :math:`p_{0,j_{\rm max}}^\inp` and offset
every element in :math:`p_0^\outp` so that :math:`p_{0,j_{\rm max}}^\outp = p_{0,j_{\rm max}}^\inp`.
We are effectively using the location where the :math:`\rho_0^\inp` drops below
:math:`\rho_{\rm cutoff}` as the starting point for integration.


Advect Base Enthalpy
^^^^^^^^^^^^^^^^^^^^
:math:`[(\rho h)_0^\inp,w_0^\inp,\psi^\inp] \rightarrow [(\rho h)_0^\outp]`
is the process by which we update the base state enthalpy through :math:`\dt` in time.

* planar:

  We discretize :eq:`eq:flow:enthalpy_base`, neglecting reaction source terms, to
  compute the new base state enthalpy,

  .. math:: (\rho h)_{0,j}^{\outp} = (\rho h)_{0,j}^{\inp} - \frac{\dt}{\Delta r} \left\{ \left[ (\rho h)_0^{\nph} w_0^{\inp}\right]_{j+\myhalf} - \left[ (\rho h)_0^{\nph} w_0^{\inp}\right]_{j-\myhalf} \right\} + \dt\psi_j^{\inp}.

  We compute the time-centered edge states, :math:`(\rho h)_0^{\nph}`, by discretizing
  an expanded form of :eq:`eq:flow:enthalpy_base`:

  .. math:: \frac{\partial (\rho h)_0}{\partial t} + w_0 \frac{\partial (\rho h)_0}{\partial r} = -(\rho h)_0 \frac{\partial w_0}{\partial r} + \psi.

* spherical:

  The base state enthalpy update now includes the area factors
  in the divergences:

  .. math::

       \begin{aligned}
       (\rho h)_{0,j}^{\outp} &= (\rho h)_{0,j}^{\inp} \nonumber \\
       & - \frac{1}{r_j^2} \frac{\dt}{\dr} \left \{ \left[ r^2 (\rho h)_0^{\nph} w_0^{\inp}\right]_{j+\myhalf} - \left[ r^2 (\rho h)_0^{\nph} w_0^{\inp}\right]_{j-\myhalf} \right\} +\dt\psi^{\inp,\nph}.\nonumber\\\end{aligned}

  In order to compute the time-centered edge states, an additional geometric
  term is added to the forcing, due to the spherical discretization of
  :eq:`eq:flow:enthalpy_base`:

  .. math:: \frac{\partial (\rho h)_0}{\partial t} + w_0 \frac{\partial (\rho h)_0}{\partial r} = -(\rho h)_0 \frac{\partial w_0}{\partial r} - \frac{2 (\rho h)_0 w_0}{r} + \psi.


Computing :math:`w_0`
^^^^^^^^^^^^^^^^^^^^^

Here we describe the process by which we compute :math:`w_0`. The arguments
are different for planar and spherical geometries.


* planar:

  :math:`[\Sbar^{\inp},\gammabar^{\inp}, p_0^{\inp},\psi^{\inp}]\rightarrow [w_0^{\outp}]`:

  In Paper III, we showed that :math:`\psi=\etarho g` for planar
  geometries, and derived derived :eq:`eq:flow:dw0dr_planar` as an
  alternate expression for eq:`eq:w0 divergence`. We discretize this
  as:

  .. math:: \frac{w_{0,j+\myhalf}^\outp-w_{0,j-\myhalf}^\outp}{\Delta r} = \left(\Sbar^{\inp} - \frac{1}{\gammabar^{\inp} p_0^{\inp}}\psi^{\inp}\right)_j,

  with :math:`w_{0,-\myhalf}=0`.


* spherical:

  :math:`[\Sbar^{\inp},\gammabar^{\inp},\rho_0^{\inp},p_0^{\inp},\etarho^{\inp}] \rightarrow[w_0^{\outp}]`:

  We begin with :eq:`eq:w0 divergence` written in spherical
  coordinates:

  .. math:: \frac{1}{r^2}\frac{\partial}{\partial r} \left (r^2 \beta_0 w_0 \right ) = \beta_0 \left ( \Sbar - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} \right ).

  We expand the spatial derivative and recall from Paper I that

  .. math:: \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial r} = \frac{1}{\beta_0} \frac{\partial \beta_0}{\partial r},

  giving:

  .. math::
     \frac{1}{r^2} \frac{\partial}{\partial r} \left (r^2 w_0 \right ) = \Sbar - \frac{1}{\gammabar p_0} \underbrace{\left( \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial r} \right)}_{\psi} .
    :label: eq:psi def

  We solve this equation for :math:`w_0` as described in Appendix B of the multilevel paper.


.. _sec:flow:singlestep:


Volume Discrepancy Changes
--------------------------

Chapter \ `[ch:volume] <#ch:volume>`__ describes the reasoning behind the volume discrepancy
term—a forcing term added to the constraint equation to bring us back to
the equation of state. This addition of this term (enabled with dpdt_factor= T) modifies our
equation set in the following way:

-  In **Step 2B**, to compute :math:`w_0`, we need to account for the volume discrepancy
   term by first defining :math:`p_{\rm EOS}^n = \overline{p(\rho,h,X_k)^n}`, and then using:

   .. math:: \frac{\partial w_0^{\nph,\star}}{\partial r} = \overline{S^{\nph,\star\star}} - \frac{1}{\overline{\Gamma_1^n}p_0^n}\psi^{\nmh} - \underbrace{\frac{f}{\overline{\Gamma_1^n}p_0^n}\left(\frac{p_0^n-\overline{p_{\rm EOS}^n}}{\Delta t}\right)}_{\delta\chi_{w_0}}.

-  In **Step 3**, the MAC projection should account for the volume discrepancy term:

   .. math::

      \nabla \cdot \left(\beta_0^n \uadvone\right) =
      \beta_0^n \left[ \left(S^{\nph,\star\star} - \overline{S^{\nph,\star\star}}\right)
      + \underbrace{\frac{f}{\gammabar^n p_0^n}
      \left(\frac{p_{\rm EOS}^n - \overline{p_{\rm EOS}^n}}{\Delta t^n}\right)}_{\delta\chi}\right].

-  In **Step 6B**, to compute :math:`w_0`, we need to account for the volume discrepancy
   term by first defining :math:`p_{\rm EOS}^{n+1,\star} = p(\rho,h,X_k)^{n+1,\star}`,
   :math:`\overline{\Gamma_1^{\nph,\star}} = (\overline{\Gamma_1^{n}}+\overline{\Gamma_1^{n+1,\star}})/2`,
   and :math:`p^{\nph,\star} = (p^{n}+p^{n+1,\star})/2`, and then using:

   .. math:: \frac{\partial w_0^{\nph}}{\partial r} = \overline{S^{\nph,\star}} - \frac{1}{\overline{\Gamma_1^{\nph,\star}}p_0^{\nph,\star}}\psi^{\nph,\star} - \frac{f}{\overline{\Gamma_1^{n+1,\star}}p_0^{n+1,\star}}\left(\frac{p_0^{n+1,\star}-\overline{p_{\rm EOS}^{n+1,\star}}}{\Delta t}\right) - \delta\chi_{w_0}

-  In **Step 7**, the MAC projection should account for the volume discrepancy term:

   .. math:: \nabla \cdot \left(\beta_0^{\nph,\star} \uadvtwo\right) = \beta_0^{\nph,\star}\left[\left(S^{\nph,\star} - \overline{S^{\nph,\star}}\right) + \frac{f}{\overline{\Gamma_1^{n+1,\star}} p_0^{n+1,\star}} \left(\frac{p_{\rm EOS}^{n+1,\star} - \overline{p_{\rm EOS}^{n+1,\star}}}{\Delta t^n}\right) + \delta\chi\right],

   where :math:`p(\rho,h,X_k)^{\nph,\star} = \left[p(\rho,h,X_k)^n +p(\rho,h,X_k)^{n+1,\star}\right]/2`.

-  In **Step 11**, the approximate projection should account for the volume
   discrepancy term:

   .. math::

      \nabla \cdot \left(\beta_0^{\nph} \Ubt^{n+1} \right)  = \beta_0^{\nph}\left\{  \left(S^{n+1} - \overline{S^{n+1}} \right)
      + \frac{f}{\overline{\Gamma_1^{n+1}} p_0^{n+1}}
      \left[\frac{p(\rho,h,X_k)^{n+1} - \overline{p(\rho,h,X_k)^{n+1}}}{\Delta t^n}\right]\right\}.


.. _sec_flow_gamma1vary:

:math:`\Gamma_1` Variation Changes
----------------------------------

The constraint we derive from requiring the pressure to be close to
the background hydrostatic pressure takes the form:

.. math:: \nablab \cdotb \Ub + \frac{1}{\Gamma_1 p_0} \frac{Dp_0}{Dt} = S  .

The default MAESTROeX algorithm replaces :math:`\Gamma_1` with :math:`\gammabar`,
allowing us to write this as a divergence constraint. In paper III,
we explored the effects of localized variations in :math:`\Gamma_1` by
writing :math:`\Gamma_1 = \gammabar + \delta \Gamma_1`. This gives us:

.. math::

   \nablab \cdotb \Ub + \frac{1}{(\gammabar + \delta \Gamma_1) \; p_0}
     \Ub \cdotb \nablab p_0 = S - \frac{1}{(\gammabar + \delta \Gamma_1) \; p_0}
     \frac{\partial p_0}{\partial t}  .

Assuming that :math:`\delta \Gamma_1 \ll \gammabar`, we then
have

.. math::

   \nablab \cdotb \Ub + \frac{1}{\gammabar p_0} \Ub \cdotb \nablab p_0
    \left [ 1 - \frac{\delta\Gamma_1}{\gammabar} + \frac{(\delta\Gamma_1)^2}{\gammabar^2} \right ]
   = S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t}
    \left [ 1 - \frac{\delta\Gamma_1}{\gammabar} + \frac{(\delta\Gamma_1)^2}{\gammabar^2} \right ] \\

Grouping by order of the correction, we have

.. math::
   \nablab \cdotb \Ub + \frac{1}{\gammabar p_0} \Ub \cdotb \nablab p_0
   = S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} +
   \underbrace{\frac{\delta \Gamma_1}{\gammabar^2 p_0}
               \left [\frac{\partial p_0}{\partial t} + \Ub \cdotb \nablab p_0\right ]}_{\mbox{first order corrections}}
    -
     \underbrace{\frac{(\delta\Gamma_1)^2}{\gammabar^3 p_0}
               \left [ \frac{\partial p_0}{\partial t} + \Ub \cdotb \nablab p_0 \right ]}_{\mbox{second order corrections}}
   :label: eq:gammafull

Keeping to First Order in :math:`\delta\Gamma_1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The base state evolution equation is the average of :eq:`eq:gammafull` over a layer

.. math::

   \nablab \cdot w_0 \er + \frac{1}{\gammabar p_0} w_0 \er \cdotb \nablab p_0 =
   \Sbar - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} + \overline{
     \left ( \frac{\delta \Gamma_1}{\gammabar^2 p_0} \Ubt \cdotb \nablab p_0
     \right ) }  .

where we see that the :math:`[\delta \Gamma_1/(\Gamma_1^2 p_0)] \partial p_0/\partial t` terms averages to zero, since the average of :math:`\delta\Gamma_1` term is zero.
Subtracting this from :eq:`eq:gammafull`, we have

.. math::

   \nablab \cdotb \Ubt + \frac{1}{\gammabar p_0} \Ubt \cdotb \nablab p_0 = S -
   \Sbar + \frac{\delta \Gamma_1}{\gammabar^2 p_0} \left (\psi + \Ubt \cdotb
   \nablab p_0 \right ) - \overline{ \left ( \frac{\delta \Gamma_1}{\gammabar^2 p_0} \Ubt
   \cdotb \nablab p_0 \right ) }  .

These can be written more compactly as:

.. math::

   \frac{\partial w_0}{\partial r} = \Sbar -\frac{1}{\gammabar p_0}\psi +
   \overline{ \left ( \frac{\delta \Gamma_1}{\gammabar^2 p_0} \Ubt \cdotb \nablab
     p_0 \right ) }  , \label{eq:base_w0_with_dgamma1}

for plane-parallel geometries (analogous to
:eq:`eq:flow:dw0dr_planar`, and

.. math::
   \nablab \cdotb (\beta_0 \Ubt) = \beta_0 \left [ S - \Sbar + \frac{\delta
       \Gamma_1}{\gammabar^2 p_0} \psi + \frac{\delta \Gamma_1}{\gammabar^2 p_0}
     \Ubt \cdotb \nablab p_0 - \overline{ \left ( \frac{\delta
         \Gamma_1}{\gammabar^2 p_0} \Ubt \cdotb \nablab p_0 \right ) } ~ \right ]
   :label: eq:constraint_with_delta_gamma

This constraint is not in a form that can be projected. To solve this
form, we need to use a lagged :math:`\Ubt` in the righthand side.

This change comes into MAESTROeX in a variety of steps, summarized here.
To enable this portion of the algorithm, set use_delta_gamma1_term = T.

-  In **Step 3**, we are doing the “predictor” portion of the
   MAESTROeX algorithm, getting the MAC velocity that satisfies the constraint,
   so we do not try to incorporate the :math:`\delta \Gamma_1` effect. We set
   all the :math:`\delta \Gamma_1` terms in
   :eq:`eq:constraint_with_delta_gamma` to zero.

-  In **Step 6**, we are computing the new time-centered source,
   :math:`S^{\nph,\star}` and the base state velocity, :math:`w_0^\nph`. Now we can
   incorporate the :math:`\delta \Gamma_1` effect. First we construct:

   .. math::

      \frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0 \approx
         \frac{\Gamma_1^{n+1,\star} - \overline{\Gamma_1^{n+1,\star}}}
               {{\overline{\Gamma_1^{n+1,\star}}}^2} \frac{1}{p_0^n} \Ubt^n \cdotb \nablab p_0^n

   Then we call **average** to construct the lateral average of this

   .. math::

      \overline{\frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0} =
         \mathbf{Avg} \left (\frac{\Gamma_1^{n+1,\star} - \overline{\Gamma_1^{n+1,\star}}}
               {{\overline{\Gamma_1^{n+1,\star}}}^2} \frac{1}{p_0^n} \Ubt^n \cdotb \nablab p_0^n \right )

   Since the average of this is needed in advancing :math:`w_0`, we modify :math:`\overline{S}`
   to include this average:

   .. math::

      \overline{S^{\nph,\star}} \leftarrow \overline{S^{\nph,\star}} +
         \overline{\frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0}

-  In **Step 7**, we now include the :math:`\delta \Gamma_1` term in the righthand
   side for the constraint by solving:

   .. math::

      \nabla \cdot \left(\beta_0^{\nph,\star} \uadvtwo\right) =
      \beta_0^{\nph,\star}\left( S^{\nph,\star} - \overline{S^{\nph,\star}} +
         \frac{\Gamma_1^{n+1,\star} - \overline{\Gamma_1^{n+1,\star}}}
               {{\overline{\Gamma_1^{n+1,\star}}}^2} \frac{1}{p_0^n}
               (\psi^{\nph,\star} + \Ubt^n \cdotb \nablab p_0^n)
      \right)

   We note that this includes the average of the correction term as shown
   in :eq:`eq:constraint_with_delta_gamma` because we modified
   :math:`\bar{S}` to include this already.

-  In **Step 10**, we do a construction much like that done in **Step 6**,
   but with the time-centerings of some of the quantities changed.
   First we construct:

   .. math::

      \frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0 \approx
         \frac{\Gamma_1^{n+1} - \overline{\Gamma_1^{n+1}}}
               {{\overline{\Gamma_1^{n+1}}}^2} \frac{1}{p_0^{n+1}} \Ubt^n \cdotb \nablab p_0^{n+1}

   Then we call **average** to construct the lateral average of this

   .. math::

      \overline{\frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0} =
         \mathbf{Avg} \left (\frac{\Gamma_1^{n+1} - \overline{\Gamma_1^{n+1}}}
               {{\overline{\Gamma_1^{n+1}}}^2} \frac{1}{p_0^{n+1}} \Ubt^n \cdotb \nablab p_0^{n+1} \right )

   Again we modify :math:`\overline{S}` to include this average:

   .. math::

      \overline{S^{n+1}} \leftarrow \overline{S^{n+1}} +
         \overline{\frac{\delta \Gamma_1}{\Gamma_1^2 p_0} \Ubt \cdotb \nablab p_0}

-  In **Step 11**, we modify the source of the constraint to include the
   :math:`\delta \Gamma_1` information. In particular, we solve:

   .. math::

      \nabla \cdot \left(\beta_0^{\nph} \Ubt^{n+1} \right)
      = \beta_0^{\nph} \left(S^{n+1} - \overline{S^{n+1}} +
         \frac{\Gamma_1^{n+1} - \overline{\Gamma_1^{n+1}}}
               {{\overline{\Gamma_1^{n+1}}}^2} \frac{1}{p_0^{n+1}}
               (\psi^{\nph} + \Ubt^n \cdotb \nablab p_0^{n+1})
      \right)

.. _sec:flow:diffusion:

Thermal Diffusion Changes
-------------------------

Thermal diffusion was introduced in the XRB :cite:`xrb`. This
introduces a new term to :math:`S` as well as the enthalpy equation.
Treating the enthalpy equation now requires a parabolic solve. We describe
that process here.

Immediately after **Step 4H**, diffuse the enthalpy through
a time interval of :math:`\dt`. First, define :math:`(\rho h)^{(1a),\star} = (\rho h)^{(2),\star}`.
We recompute :math:`(\rho h)^{(2),\star}` to account for thermal diffusion. Here we begin
with the enthalpy equation, but consider only the
diffusion terms,

.. math:: \frac{\partial (\rho h)}{\partial t} = \nabla\cdot\kth\nabla T.

We can recast this as an enthalpy-diffusion equation by applying the
chain-rule to :math:`h(p_0,T,X_k)`,

.. math:: \nabla h = h_p \nabla p_0 + c_p \nabla T + \sum_k \xi_k \nabla X_k  ,

giving

.. math::

   \frac{\partial (\rho h)}{\partial t}  =
    \nabla\cdot \frac{\kth}{c_p}\nabla h -
    \sum_k \nabla\cdot \frac{\xi_k \kth}{c_p}\nabla X_k -
    \nabla\cdot \frac{h_p \kth}{c_p}\nabla p_0.

Compute :math:`\kth^{(1)}, c_p^{(1)}`, and :math:`\xi_k^{(1)}` from :math:`\rho^{(1)}, T^{(1)}`, and :math:`X_k^{(1)}` as inputs to the equation of state. The update is given by

.. math::

   \begin{align}
   (\rho h)^{(2),\star} &= (\rho h)^{(1a),\star} + \frac{\dt}{2}\nabla\cdot\left(\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(2),\star} + \frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(2),\star} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n+1,\star} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),
   \end{align}

which is numerically implemented as a diffusion equation for :math:`h^{(2),\star}`,

.. math::

   \begin{align}
   \left(\rho^{(2),\star} - \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla\right)h^{(2),\star} &= (\rho h)^{(1a),\star} + \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\nonumber\\
   &- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(2),\star} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n+1,\star} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),
   \end{align}

Immediately after **Step 8H**, diffuse the enthalpy through a time interval of
:math:`\dt`. First, define :math:`(\rho h)^{(1a)} = (\rho h)^{(2)}`. We recompute :math:`(\rho h)^{(2)}` to
account for thermal diffusion. Compute :math:`\kth^{(2),\star}, c_p^{(2),\star}`, and
:math:`\xi_k^{(2),\star}`, from :math:`\rho^{(2),\star}, T^{(2),\star}`, and :math:`X_k^{(2),\star}` as inputs to
the equation of state. The update is given by

.. math::
   \begin{align}
   (\rho h)^{(2)} &= (\rho h)^{(1a)} + \frac{\dt}{2}\nabla\cdot\left(\frac{\kth^{(2),\star}}{c_p^{(2),\star}}\nabla h^{(2)} + \frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla X_k^{(2)} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla p_0^{n+1} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),
   \end{align}

which is numerically implemented as a diffusion equation for :math:`h^{(2)}`.

.. math::
   \begin{align}
   \left(\rho^{(2)} - \frac{\dt}{2}\nabla\cdot\frac{\kth^{(2),\star}}{c_p^{(2),\star}}\nabla\right)h^{(2)} &= (\rho h)^{(1a)} + \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\nonumber\\
   &- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla X_k^{(2)} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla p_0^{n+1} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),
   \end{align}
.. _sec:Initialization:

Initialization
==============

We start each calculation with user-specified initial values for
:math:`\rho`, :math:`X_k` and :math:`T,` as well as an initial background state. In
order for the low Mach number assumption to hold, the initial data
must be thermodynamically consistent with the initial background
state. In addition, the initial velocity field must satisfy an
initial approximation to the divergence constraint. We use an iterative
procedure to compute both an initial right-hand-side for the
constraint equation and an initial velocity field that satisfies
the constraint. The user specifies the number of iterations,
:math:`N_{\rm iters}^{S},` in this first step of the initialization procedure.

| The initial perturbational pressure also needs to be determined for
  use in **Steps 3**, **7**, and **11**.
  This is done through a second iterative procedure which follows the
  time advancement algorithm as described in **Steps 1-11** in
  §\ `2 <#Sec:Time Advancement Algorithm>`__.
  The user specifies the number of iterations,
  :math:`N_{\rm iters}^{\pi},` in this second step of the initialization procedure.
  The details for both iterations are given below.

First, we need to construct approximations to :math:`S^0, w_0^{-\myhalf}, \Delta t^0`,
and :math:`\Ub^0`. Start with initial data :math:`X_k^{\initp}, \rho^{\initp},` :math:`T^{\initp},`
an initial base state, and an initial guess for the velocity, :math:`\Ub^{\initp}`.
Set :math:`w_0^0 = 0` as an initial approximation. Use the equation of state to
determine :math:`(\rho h)^{\initp}`. Compute :math:`\beta_0^0` as a function of
the initial data. The next part of the initialization process
proceeds as follows.

#. *Initial Projection*: if ``do_initial_projection = T``, then we
   first project the velocity field with :math:`\rho = 1` and :math:`\beta_0^0`.
   The initial projection does not see reactions
   or external heating, and thus we set :math:`\dot\omega = \Hnuc = \Hext = 0` in :math:`S`.
   The reason for ignoring reactions and heating is that we need some kind of
   time scale over which to compute the effect of reactions, but we first
   need an estimate of the velocity
   field in order to derive the time step that will be used as a time scale.
   The elliptic equation we solve is

   .. math:: \nabla \cdot \beta_0^0 \nabla \phi = \underbrace{\beta_0^0(S - \Sbar)}_\mathrm{0~ except~ for~ diffusion} - \nabla \cdot (\beta_0^0 \Ub^{\initp})

   This :math:`\phi` is then used to correct the velocity field to obtain :math:`\Ub^{0,0}`.
   If do_initial_projection = F, set :math:`\Ub^{0,0} = \Ub^{\initp}`.

#. *“Divu” iterations*: Next we do ``init_divu_iter`` iterations
   to project the velocity field using a constraint that sees reactions
   and external heating.
   The initial timestep estimate is provided by ``firstdt`` and
   :math:`\Ub^{0,0}`, to allow us to compute the effect of reactions over :math:`\Delta t/2`.

   **Do** :math:`\nu = 1,...,N_{\rm iters}^{S}`.

   #. Estimate :math:`\Delta t^\nu` using :math:`\Ub^{0,\nu-1}` and :math:`w_0^{\nu-1}.`

   #. **React State**\ :math:`[ \rho^{\initp},(\rho h)^{\initp}, X_k^{\initp}, T^{\initp},
      (\rho^{\initp} \Hext), p_0^{\initp}] \rightarrow [\rho^{\outp}, (\rho h)^{\outp},
      X_k^{\outp}, T^{\outp}, (\rho \omegadot_k)^{0,\nu} ].`

   #. Compute :math:`S^{0,\nu}`
      using :math:`(\rho \omegadot_k)^{0,\nu}` and the initial data.

   #. Compute :math:`\overline{S^{0,\nu}} = {\mathrm{\bf Avg}} (S^{0,\nu}).`

   #. Compute :math:`w_0^{\nu}` as in **Step 2B** using :math:`\overline{S^{0,\nu}}` and :math:`\psi=0`.

   #. Project :math:`\Ub^{0,\nu-1}` using :math:`\beta_0^0` and
      :math:`(S^{0,\nu} - \overline{S^{0,\nu}})` as the source term.
      This yields :math:`\Ub^{0,\nu}.` In this projection, again the density is
      set to 1, and the elliptic equation we solve is:

      .. math:: \nabla \cdot \beta_0^0 \nabla \phi = \beta_0 (S - \Sbar)- \nabla \cdot (\beta_0^0 \Ub^{0,\nu-1})

   **End do.**

   Define :math:`S^0 = S^{0,N_{\rm iters}^S}`, :math:`w_0^{-\myhalf} = w_0^{N_{\rm iters}^S}`,
   :math:`\dt^0 = \Delta t^{N_{\rm iters}^S},` and :math:`\Ub^0 = \Ub^{0,N_{\rm iters}^S}.`

Next, we need to construct approximations to :math:`\etarho^{-\myhalf}, \psi^{-\myhalf}, S^1`,
and :math:`\pi^{-\myhalf}`. As initial approximations, set
:math:`\etarho^{-\myhalf}=0, \psi^{-\myhalf}=0, S^{1,0}=S^0`, and :math:`\pi^{-\myhalf}=0.`

#. *Pressure iterations*: Here we do ``init_iter`` iterations to get an
   approximation for the lagged pressure:

   **Do** :math:`\nu = 1,...,N_{\rm iters}^{\pi}`.

   #. Perform **Steps 1-11** as described above, using
      :math:`S^{\myhalf,\star\star} = (S^0 + S^{1,\nu-1})/2` in **Step 2** as described.
      The only other difference in the time advancement is that in **Step 11**
      we define :math:`{\bf V} = (\Ubt^{1,\star} - \Ubt^0)` and solve

      .. math:: L_\beta^\rho \phi = D \left ( \beta_0^{\myhalf} {\bf V} \right) - \beta_0^{\myhalf} \left[ \left(S^{1}-\overline{S^{1}}\right) - \left(S^{0}-\overline{S^{0}}\right) \right]  .

      (The motivation for this form of the projection in the initial pressure iterations
      is discussed in :cite:`almgren:bell:crutchfield`.)
      We discard the new velocity resulting from this, but keep the new
      value for :math:`\pi^{\myhalf} = \pi^{-\myhalf} + (1 / \dt) \; \phi.`
      These steps also yield new scalar data at time :math:`\dt,` which
      we discard, and new values for :math:`\etarho^{\myhalf}` (**Step 8C**),
      :math:`\psi^{\myhalf}` (**Step 8F**),
      :math:`S^{1,\nu}` (**Step 10A**), and :math:`\pi^{\myhalf}` (**Step 11**), which we keep.

   #. Set :math:`\pi^{-\myhalf} = \pi^{\myhalf}`, :math:`\etarho^{-\myhalf} = \etarho^{\myhalf}`,
      and :math:`\psi^{-\myhalf} = \psi^{\myhalf}`.

   **End do.**

   Finally, we define :math:`S^1 = S^{1,N_{\rm iters}^\pi}.`

The tolerances for these elliptic solves are described in § :ref:`sec:mg`.

Changes from Earlier Implementations
====================================

Changes Between Paper 3 and Paper 4
-----------------------------------

#. We defined the mapping of data between a 1D radial array and the 3D Cartesian
   grid for spherical problems (which we improve upon in the multilevel paper).

#. We update :math:`T` after the call to **React State**.

#. We have created a burning_cutoff_density, where the burning does
   not happen below this density. It is presently set to ``base_cutoff_density``.

#. Use corner coupling in advection.

#. We have an option, controlled by use_tfromp, to update temperature
   using :math:`T=T(\rho,X_k,p_0)` rather than :math:`T=T(\rho,h,X_k)`. The former completely
   decouples enthalpy from our system. For spherical problems, we use
   ``use_tfromp = TRUE``, for planar problems, we use use_tfromp = FALSE.

#. For spherical problems, we have changed the discretization of
   :math:`\Ubt\cdot\nabla p_0` in the enthalpy update to
   :math:`\nabla\cdot(\Ubt p_0) - p_0\nabla\cdot\Ubt`.

#. In paper III we discretized the enthalpy evolution equation in
   terms of :math:`T`. Since then we have discovered that
   discretizing the enthalpy evolution in perturbational form, :math:`(\rho h)'`,
   leads to better numerical properties. We use ``enthalpy_pred_type= 1``.
   This is more like paper II.

#. We have turned off the evolution of :math:`h` above the atmosphere and instead
   compute :math:`h` with the EOS using ``do_eos_h_above_cutoff = T``.

Changes Between Paper 4 and the Multilevel Paper
------------------------------------------------

See the multilevel paper for the latest.

Changes Between the Multilevel Paper and Paper 5
------------------------------------------------

#. Added rotation.

Changes Between Paper 5 and the XRB Paper
-----------------------------------------

#. We have added thermal diffusion, controlled by ``use_thermal_diffusion``,
   ``temp_diffusion_formulation``, and ``thermal_diffusion_type``.

#. We added the volume discrepancy term to the velocity constraint equation,
   controlled by the input parameter, ``dpdt_factor``.

#. For certain problems, we need to set ``do_eos_h_above_cutoff = F``
   to prevent large, unphysical velocities from appearing near the edge of the star.

Changes Since the XRB Paper
---------------------------

#. We switched to the new form of the momentum equation to
   :eq:`eq:flow:newmomentum` to conserve the low-Mach number form of
   energy.

#. We changed the form of the volume discrepancy term to get better
   agreement between the two temperatures.

Future Considerations
=====================

-  We are still exploring the effects of ``use_tfromp = F`` for spherical
   problems. We would eventually like to run in this mode, but :math:`T=T(\rho,X_k,p_0)`
   and :math:`T=T(\rho,h,X_k)` drift away from each other more than we would like. Our
   attempts at incorporating a ``dpdt_factor`` for spherical problems have not
   been successful.

.. _ch:methodology:

*********************
Numerical Methodology
*********************

We give an overview of the original MAESTRO algorithm :cite:`multilevel`,
as well as the recently introduced new temporal integration
scheme :cite:`MAESTROeX`.  These schemes share many of the same algorithmic
modules, however the new scheme is significantly simpler.

To summarize, in the original scheme we integrated the base state evolution
over the time step, which requires splitting the velocity equation into
more complicated average and perturbational equations.
In the new temporal integrator, we eliminated much of this complication
by introducing a predictor-corrector approach for the base state.
The key observation is that the base state density is simply the lateral
average of the full density, so we simply update the base state density
through averaging routines, rather than predicting the evolution with
a split velocity formulation, only to have to correct it with the
averaging operator in the end regardless.

Original Temporal Integrator
============================

New Temporal Integrator
=======================
