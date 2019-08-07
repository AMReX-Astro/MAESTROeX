.. _ch:flowchart:

*******************
MAESTROeX Flowchart
*******************

The equation set and solution procedure used by MAESTROeX has evolved
over time. In this chapter, we outline the algorithm currently
implemented in the code. The latest published reference for MAESTROeX
is the multilevel paper :raw-latex:`\cite{multilevel}`. In this description, we
make frequent reference to paper I :raw-latex:`\cite{lowMach}`,
paper II :raw-latex:`\cite{lowMach2}`, paper III :raw-latex:`\cite{lowMach3}`, and
paper IV :raw-latex:`\cite{lowMach4}`.

Summary of the MAESTROeX Equation Set
=====================================

Here we summarize the equations solved by MAESTROeX. We refer the reader
to papers I through IV and the multilevel paper for the derivation
and motivation of the equation set. (Note: this ‘traditional’ algorithm
uses Strang-splitting for the reactions. An alternate implementation, using
spectral deferred corrections is outlined in Chapter \ `[ch:sdc] <#ch:sdc>`__.)

Base State
----------

The stratified atmosphere is characterized by a one-dimensional
time-dependent base state, defined by a base state density, :math:`\rho_0`,
and a base state pressure, :math:`p_0`, in hydrostatic equilibrum:

.. math:: \nabla p_0 = -\rho_0 |g| \er

The gravitational acceleration, :math:`g` is either constant or a point-mass
with a :math:`1/r^2` dependence (see §\ `[sec:planarinvsqgravity] <#sec:planarinvsqgravity>`__) for plane-parallel geometries, or a monopole
constructed by integrating the base state density for spherical
geometries.

For the time-dependence, we will define a base state velocity, :math:`w_0`,
which will adjust the base state from one hydrostatic equilibrum to
another in response to heating.

For convenience, we define a base state enthalphy, :math:`h_0`, as needed
by laterally averaging the full enthalpy, :math:`h`.

Continuity
----------

Conservation of mass gives the same continuity equation we have with
compressible flow:

.. math::

   \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \Ub) = 0
   \label{eq:flow:continuity}

Additionally, we carry species around which can react. The creation and destruction
of the species is described by their create rate, :math:`\omegadot_k`, and the species
are defined by their mass fractions, :math:`X_k \equiv \rho_k / \rho`, giving

.. math::

   \frac{\partial \rho X_k}{\partial t} + \nabla \cdot (\rho \Ub X_k) = \rho \omegadot_k
   \label{eq:flow:rhoX}

and

.. math:: \sum_k X_k = 1

The base state density evolution equation can be defined by laterally averaging the
full continuity equation, giving:

.. math::

    \frac{\partial\rho_0}{\partial t} = -\nabla\cdot(\rho_0 w_0 \eb_r),
    \label{eq:flow:base_density}

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
   \label{eq:flow:etarho}

In practice, we correct the drift by simply setting :math:`\rho_0 =
\overline{\rho}` after the advective update of :math:`\rho`. However we still need to
explicitly compute :math:`\etarho` since it appears in other equations.

Constraint
----------

The equation of state is cast into an elliptic constraint on the
velocity field by differentiating :math:`p_0(\rho, s, X_k)` along particle
paths, giving:

.. math::

   \nabla \cdot (\beta_0 \Ub) =
      \beta_0 \left ( S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} \right )
   \label{eq:U divergence}

where :math:`\beta_0` is a density-like variable that carries background
stratification, defined as

.. math:: \beta_0(r,t) = \rho_0(0,t)\exp\left(\int_0^r\frac{1}{\gammabar p_0}\frac{\partial p_0}{\partial r'}dr'\right),

and

.. math::

   S = -\sigma\sum_k\xi_k\omegadot_k + \frac{1}{\rho p_\rho}\sum_k p_{X_k}\omegadot_k + \sigma\Hnuc + \sigma\Hext + \frac{\sigma}{\rho} \nabla \cdot \kth \nabla T
   \label{eq:flow:S}

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
divergence. :raw-latex:`\cite{KP:2012}` discuss the general case where we want to
keep the local variations of :math:`\Gamma_1` (and we explored this in paper
III). We also look at this in § \ `[sec:flow:gamma1vary] <#sec:flow:gamma1vary>`__.

Momentum
--------

The compressible momentum equation (written in terms of velocity is):

.. math:: \rho \frac{\partial \Ub}{\partial t} + \rho \Ub \cdot \nabla \Ub + \nabla p = -\rho |g| \er

Subtracting off the base state, and defining the perturbational
pressure (sometimes called the dynamic pressure) as :math:`\pi \equiv p - p_0`,
and perturbational density as :math:`\rho' \equiv \rho - \rho_0`, we have:

.. math:: \rho \frac{\partial \Ub}{\partial t} + \rho \Ub \cdot \nabla \Ub + \nabla \pi = -\rho' |g| \er

or

.. math::

   \frac{\partial \Ub}{\partial t} + \Ub \cdot \nabla \Ub + \frac{1}{\rho} \nabla \pi =
      -\frac{\rho^\prime}{\rho} |g| \er

This is the form of the momentum equation that we solved in papers
I–IV and in the multilevel paper.

Several authors :raw-latex:`\cite{KP:2012,VLBWZ:2013}` explored the idea of energy
conservation in a low Mach number system and found that an additional
term (which can look like a buoyancy) is needed in the low Mach number
formulation, yielding:

.. math::

   \frac{\partial \Ub}{\partial t} + \Ub \cdot \nabla \Ub +
      \frac{\beta_0}{\rho} \nabla \left (\frac{p^\prime}{\beta_0} \right ) =
      -\frac{\rho^\prime}{\rho} |g| \er
   \label{eq:flow:newmomentum}

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

   \begin{aligned}
   \frac{\partial w_0}{\partial t} &= -w_0\frac{\partial w_0}{\partial
     r} - \frac{\beta_0}{\rho_0}\frac{\partial(\pi_0/\beta_0)}{\partial r},\label{eq:w0
     evolution}\\
   %
   \frac{\partial\Ubt}{\partial t} &= -(\Ubt + w_0\er)\cdot\nabla\Ubt
     - \left(\Ubt\cdot\eb_r\right)\frac{\partial w_0}{\partial r}\eb_r -
   \frac{\beta_0}{\rho}\nabla\left(\frac{\pi}{\beta_0} \right) +
   \frac{\beta_0}{\rho_0}\frac{\partial(\pi_0/\beta_0)}{\partial r}\eb_r -
   \frac{\rho-\rho_0}{\rho}g\eb_r.\label{eq:flow:utildeupd}\end{aligned}

where :math:`\pi_0` is the base state component of the perturbational pressure.
By laterally averaging to equation (`[eq:U divergence] <#eq:U divergence>`__),
we obtain a divergence constraint for :math:`w_0`:

.. math::

   \nabla\cdot(\beta_0 w_0 \eb_r) =
       \beta_0\left(\Sbar - \frac{1}{\gammabar p_0}
              \frac{\partial p_0}{\partial t}\right).\label{eq:w0 divergence}

The divergence constraint for :math:`\Ubt` can be found by subtracting
(`[eq:w0 divergence] <#eq:w0 divergence>`__) into (`[eq:U divergence] <#eq:U divergence>`__), resulting in

.. math:: \nabla\cdot\left(\beta_0\Ubt\right) = \beta_0\left(S-\Sbar\right).\label{eq:utilde divergence}

Base State Expansion
--------------------

In practice, we calculate :math:`w_0` by integrating
the one-dimensional divergence constraint. For a plane-parallel atmosphere, the
evolution is:

.. math::

   \label{eq:flow:dw0dr_planar}
   \frac{\partial w_0}{\partial r} = \Sbar - \frac{1}{\gammabar p_0} \etarho g

Then we define

.. math::

   - \frac{\beta_0}{\rho_0} \frac{\partial (\pizero/\beta_0)}{\partial r} = \frac{\partial w_0}{\partial t} +
      w_0 \frac{\partial w_0}{\partial r} , \label{eq:pizero}

once :math:`w_0` at the old and new times is known, and the advective term is computed explicitly.
Then we can include this for completeness in the update for :math:`\ut.`

Energy
------

Finally, we add an equation for specific enthalpy evolution to our
system. Strictly speaking this is not necessary to close the system,
but it becomes convenient at times to define the temperature.

.. math::

   \frac{\partial(\rho h)}{\partial t} =
      -\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt} + \rho\Hnuc + \rho\Hext,\label{eq:flow:enthalpy}

We will often expand :math:`Dp_0/Dt` as

.. math:: \frac{Dp_0}{Dt} = \psi + (\Ubt \cdot \er) \frac{\partial p_0}{\partial r}

where we defined

.. math:: \psi \equiv \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial r}

When we are using thermal diffusion, there will be an additional term in
the enthalpy equation (see § \ `2.5 <#sec:flow:diffusion>`__).

In paper III, we showed that for a plane-parallel atmosphere with
constant gravity, :math:`\psi = \etarho g`

At times, we will define a temperature equation by writing :math:`h = h(T,p,X_k)`
and differentiating:

.. math::

   \label{eq:flow:temp}
   \frac{DT}{Dt} = \frac{1}{\rho c_p} \left\{ \left(1 - \rho h_p\right) \left
     [ \psi + (\Ubt \cdotb \er) \frac{\partial p_0}{\partial r} \right ]
    - \sum_k \rho \xi_k {\omegadot}_k
    + \rho \Hnuc + \rho \Hext \right \}   .

The base state evolution equations for density and enthalpy can be
found by averaging Eq \ `[eq:flow:enthalpy] <#eq:flow:enthalpy>`__
over a layer of constant radius, resulting in

.. math::

   \begin{aligned}
   \frac{\partial(\rho h)_0}{\partial t} &=& -\nabla\cdot\left[(\rho h)_0w_0\eb_r\right] +
     \psi + \overline{\rho \Hnuc} + \overline{\rho \Hext}. \label{eq:flow:enthalpy_base}\end{aligned}

Subtracting it from the full enthalpy equation gives:

.. math::

   \frac{\partial(\rho h)'}{\partial t} = -\Ub\cdot\nabla(\rho h)' - (\rho h)'\nabla\cdot\Ub -
     \nabla\cdot\left[(\rho h)_0\Ubt\right] + \Ubt\cdot\nabla p_0
      + ( \rho\Hnuc - \overline{\rho \Hnuc}) + (\rho\Hext - \overline{\rho \Hext})
   \label{eq:flow:rhohprime}

.. _Sec:Time Advancement Algorithm:

Time Advancement Algorithm
==========================

Here is the current description of the algorithm, based on the
description in the multilevel paper. The initialization has also been
included with more detail than was given in paper III. The main
driver for a single step of the algorithm is advance.f90—refer
to that code to see the sequence of functions called to implement each
step.

Definitions
-----------

Below we define operations that will be referenced in
§ \ `[sec:flow:singlestep] <#sec:flow:singlestep>`__.

**React State**\ :math:`[\rho^{\inp},(\rho h)^{\inp},X_k^{\inp},T^{\inp}, (\rho\Hext)^{\inp}, p_0^{\inp}] \rightarrow [\rho^{\outp}, (\rho h)^{\outp}, X_k^{\outp}, T^{\outp}, (\rho \omegadot_k)^{\outp}, (\rho\Hnuc)^{\outp}]`
evolves the species and enthalpy due to reactions through
:math:`\Delta t/2` according to:

.. math::

   \frac{dX_k}{dt} = \omegadot_k(\rho,X_k,T) ; \qquad
   \frac{dT}{dt}   = \frac{1}{c_p} \left ( -\sum_k \xi_k  \omegadot_k  + \Hnuc \right ).

Here the temperature equation comes from Eq. \ `[eq:flow:temp] <#eq:flow:temp>`__ with :math:`Dp_0/Dt = 0` for
the burning part of the evolution.

Full details of the
solution procedure can be found in Paper III. We then define:

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

**Advect Base Density**\ :math:`[\rho_0^\inp,w_0^\inp] \rightarrow [\rho_0^\outp, \rho_0^{\outp,\nph}]` is the process by which we
update the base state density through :math:`\dt` in time. We keep the
time-centered edge states, :math:`\rho_0^{\outp,\nph}`,
since they are used later in discretization of :math:`\etarho` for planar problems.

planar:
    We discretize equation (`[eq:flow:base_density] <#eq:flow:base_density>`__) to
    compute the new base state density,

    .. math:: \rho_{0,j}^{\outp} = \rho_{0,j}^{\inp} - \frac{\dt}{\dr} \left [ \left( \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j+\myhalf} - \left( \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j-\myhalf} \right ].

    We compute the time-centered edge states, :math:`{\rho_0}^{\outp,\nph}_{j\pm\myhalf}`,
    by discretizing an expanded form of equation (`[eq:flow:base_density] <#eq:flow:base_density>`__):

    .. math:: \frac{\partial \rho_0}{\partial t} + w_0 \frac{\partial \rho_0}{\partial r} = - \rho_0 \frac{\partial w_0}{\partial r},

    where the right hand side is used as the force term.

spherical:
    The base state density update now includes the area factors in the
    divergences:

    .. math:: \rho_{0,j}^{\outp} = \rho_{0,j}^{\inp} - \frac{1}{r_j^2} \frac{\dt}{\dr} \left [ \left( r^2 \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j+\myhalf} - \left( r^2 \rho_0^{\outp,\nph} w_0^{\inp}\right)_{j-\myhalf} \right].

    In order to compute the time-centered edge states, an additional geometric
    term is added to the forcing, due to the spherical discretization of
    (`[eq:flow:base_density] <#eq:flow:base_density>`__):

    .. math:: \frac{\partial \rho_0}{\partial t} + w_0 \frac{\partial \rho_0}{\partial r} = - \rho_0 \frac{\partial w_0}{\partial r} - \frac{2 \rho_0 w_0}{r}.

The source code for this operation can be found in advect_base.f90.

**Enforce HSE**\ :math:`[p_0^{\inp},\rho_0^{\inp}] \rightarrow [p_0^{\outp}]` has replaced **Advect Base Pressure**
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

The source code for this operation can be found in enforce_HSE.f90.

**Advect Base Enthalpy**\ :math:`[(\rho h)_0^\inp,w_0^\inp,\psi^\inp] \rightarrow [(\rho h)_0^\outp]`
is the process by which we update the base state enthalpy through :math:`\dt` in time.

planar:
    We discretize equation (`[eq:flow:enthalpy_base] <#eq:flow:enthalpy_base>`__), neglecting reaction
    source terms, to compute the new base state enthalpy,

    .. math:: (\rho h)_{0,j}^{\outp} = (\rho h)_{0,j}^{\inp} - \frac{\dt}{\Delta r} \left\{ \left[ (\rho h)_0^{\nph} w_0^{\inp}\right]_{j+\myhalf} - \left[ (\rho h)_0^{\nph} w_0^{\inp}\right]_{j-\myhalf} \right\} + \dt\psi_j^{\inp}.

    We compute the time-centered edge states, :math:`(\rho h)_0^{\nph}`, by discretizing
    an expanded form of equation (`[eq:flow:enthalpy_base] <#eq:flow:enthalpy_base>`__):

    .. math:: \frac{\partial (\rho h)_0}{\partial t} + w_0 \frac{\partial (\rho h)_0}{\partial r} = -(\rho h)_0 \frac{\partial w_0}{\partial r} + \psi.

spherical:
    The base state enthalpy update now includes the area factors
    in the divergences:

    .. math::

       \begin{aligned}
       (\rho h)_{0,j}^{\outp} &= (\rho h)_{0,j}^{\inp} \nonumber \\
       & - \frac{1}{r_j^2} \frac{\dt}{\dr} \left \{ \left[ r^2 (\rho h)_0^{\nph} w_0^{\inp}\right]_{j+\myhalf} - \left[ r^2 (\rho h)_0^{\nph} w_0^{\inp}\right]_{j-\myhalf} \right\} +\dt\psi^{\inp,\nph}.\nonumber\\\end{aligned}

    In order to compute the time-centered edge states, an additional geometric
    term is added to the forcing, due to the spherical discretization of
    (`[eq:flow:enthalpy_base] <#eq:flow:enthalpy_base>`__):

    .. math:: \frac{\partial (\rho h)_0}{\partial t} + w_0 \frac{\partial (\rho h)_0}{\partial r} = -(\rho h)_0 \frac{\partial w_0}{\partial r} - \frac{2 (\rho h)_0 w_0}{r} + \psi.

The source code for this operation can be found in advect_base.f90.

**Computing** :math:`w_0`\ [Sec:Computing w0]
Here we describe the process by which we compute :math:`w_0`. The arguments
are different for planar and spherical geometries.

**Compute** :math:`w_0` **Planar**
:math:`[\Sbar^{\inp},\gammabar^{\inp}, p_0^{\inp},\psi^{\inp}]\rightarrow [w_0^{\outp}]`:

In Paper III, we showed that :math:`\psi=\etarho g` for planar geometries,
and derived derived Eq. \ `[eq:flow:dw0dr_planar] <#eq:flow:dw0dr_planar>`__ as an alternate
expression for Eq. \ `[eq:w0 divergence] <#eq:w0 divergence>`__. We discretize this as:

.. math:: \frac{w_{0,j+\myhalf}^\outp-w_{0,j-\myhalf}^\outp}{\Delta r} = \left(\Sbar^{\inp} - \frac{1}{\gammabar^{\inp} p_0^{\inp}}\psi^{\inp}\right)_j,

with :math:`w_{0,-\myhalf}=0`.

**Compute** :math:`w_0` **Spherical**
:math:`[\Sbar^{\inp},\gammabar^{\inp},\rho_0^{\inp},p_0^{\inp},\etarho^{\inp}] \rightarrow[w_0^{\outp}]`:

We begin with equation (`[eq:w0 divergence] <#eq:w0 divergence>`__) written in spherical coordinates:

.. math:: \frac{1}{r^2}\frac{\partial}{\partial r} \left (r^2 \beta_0 w_0 \right ) = \beta_0 \left ( \Sbar - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} \right ).

We expand the spatial derivative and recall from Paper I that

.. math:: \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial r} = \frac{1}{\beta_0} \frac{\partial \beta_0}{\partial r},

giving:

.. math:: \frac{1}{r^2} \frac{\partial}{\partial r} \left (r^2 w_0 \right ) = \Sbar - \frac{1}{\gammabar p_0} \underbrace{\left( \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial r} \right)}_{\psi}.\label{eq:psi def}

We solve this equation for :math:`w_0` as described in Appendix B of the multilevel paper.

The source code for this operation can be found in make_w0.f90.

[sec:flow:singlestep] Single Step
---------------------------------

*Initialization*

This step remains unchanged from Paper III. See §\ `3 <#Sec:Initialization>`__
for details. The initialization step only occurs at the beginning of the simulation.
The initial values for :math:`\Ub^0, \rho^0, (\rho h)^0, X_k^0, T^0,
\rho_0^0, p_0^0`, and :math:`\overline{\Gamma_1^0}` are specified from the problem-dependent
initial conditions. The initial time step, :math:`\dt^0`, is computed as in
Paper III. Finally, initial
values for :math:`w_0^{-\myhalf}, \etarho^{-\myhalf}, \psi^{-\myhalf},
\pi^{-\myhalf}, S^0`, and :math:`S^1` come from a preliminary pass through
the algorithm.

*React the full state through the first time interval of* :math:`\dt / 2.`

Call **React State**\ :math:`[\rho^n, (\rho h)^n, X_k^n, T^n, (\rho\Hext)^n, p_0^n] \rightarrow [\rho^{(1)},(\rho h)^{(1)},X_k^{(1)},T^{(1)},(\rho \omegadot_k)^{(1)},(\rho \Hnuc)^{(1)}]`.

*Compute the provisional time-centered expansion,*
:math:`S^{\nph,\star\star}`, *provisional base state velocity,*
:math:`w_0^{\nph,\star}`, *and provisional base state velocity forcing.*

Compute :math:`S^{\nph,\star\star}`. We compute an estimate for the
time-centered expansion term in the velocity divergence constraint,
as given in Eq. \ `[eq:flow:S] <#eq:flow:S>`__. For the first time step (:math:`n=0`),
we set

.. math:: S^{n+\myhalf,\star\star} = \frac{S^0 + S^1}{2},

where :math:`S^1` is found during initialization. For other time steps
:math:`(n \ne 0)`, following :raw-latex:`\cite{SNe}`, we extrapolate
to the half-time using :math:`S` at the previous and current
time levels

.. math:: S^{\nph,\star\star} = S^n + \frac{\dt^n}{2} \frac{S^n - S^{n-1}}{\dt^{n-1}}.

Next, compute

.. math:: \overline{S^{\nph,\star\star}} = {\mathrm{\bf Avg}} \left(S^{\nph,\star\star}\right).

Compute :math:`w_0^{\nph,\star}`.

| For planar geometry, call
| **Compute** :math:`w_0` **Planar**\ :math:`[\overline{S^{\nph,\star\star}},\overline{\Gamma_1^n},p_0^n,\psi^{n-\myhalf}] \rightarrow [w_0^{\nph,\star}]`.

| For spherical geometry, call
| **Compute** :math:`w_0` **Spherical**\ :math:`[\overline{S^{\nph,\star\star}},\overline{\Gamma_1^n},\rho_0^n,p_0^n,\etarho^{n-\myhalf}] \rightarrow [w_0^{\nph,\star}]`.

Compute the provisional base state velocity forcing, using equation (38)
from paper III,

.. math:: -\frac{\beta_0}{\rho_0} \frac{\partial (\pi_0/\beta_0)}{\partial r} = \frac{\partial w_0}{\partial t} + w_0 \frac{\partial w_0}{\partial r},

with the following discretization:

.. math:: \left ( \frac{\beta_0}{\rho_0} \frac{\partial (\pi_0/\beta_0)}{\partial r} \right )^{n,\star} = -\frac{w_0^{\nph,\star} - w_0^\nmh}{(\dt^n+\dt^{n-1})/2} - w_0^{n,\star} \left(\frac{\partial w_0}{\partial r}\right)^{n,\star},

where :math:`w_0^{n,\star}` and :math:`(\partial w_0 / \partial r)^{n,\star}` are defined as

.. math::

   \begin{aligned}
   w_0^{n,\star} &=& \frac{\dt^{n} w_0^{\nmh} + \dt^{n-1} w_0^{\nph,\star}}{\dt^n+\dt^{n-1}}, \\
   \left(\frac{\partial w_0}{\partial r}\right)^{n,\star} &=& \frac{1}{\dt^n+\dt^{n-1}}\left [ \dt^{n} \left(\frac{\partial w_0 }{ \partial r}\right)^{\nmh} + \dt^{n-1} \left(\frac{\partial w_0 }{ \partial r}\right)^{\nph,\star} \right ].\nonumber \\\end{aligned}

If :math:`n=0`, we use :math:`\dt^{-1} = \dt^0`.

*Construct the provisional time-centered advective velocity on
edges,* :math:`\uadvone`.

The local velocity field is described by Eq. \ `[eq:flow:utildeupd] <#eq:flow:utildeupd>`__.
From this, we compute time-centered edge
velocities, :math:`\uadvonedag`, using
:math:`\Ub = \Ubt^n + w_0^{\nph,\star}`. The :math:`\dagger` superscript refers to the
fact that the predicted velocity field does not satisfy the divergence
constraint. We then construct :math:`\uadvone` from :math:`\uadvonedag`
using a MAC projection, as described in detail in Appendix B of Paper III.
We note that :math:`\uadvone` satisfies the discrete version of
:math:`\overline{(\uadvone\cdot\eb_r)}=0` as well as

.. math::

   \begin{aligned}
   \nabla \cdot \left(\beta_0^n \uadvone\right) &=& \beta_0^n \left(S^{\nph,\star\star} - \overline{S^{\nph,\star\star}}\right),\\
    \beta_0^n &=& \beta_0 \left(\rho_0^n, p_0^n, \overline{\Gamma_1^n}\right),\end{aligned}

where :math:`\beta_0` is computed as described in Appendix C of Paper III (although
note that an alternate procedure that uses linear reconstruction of :math:`g` in
constructing :math:`\beta_0` is enabled with ``use_linear_grav_in_beta``).

*Advect the base state and full state through a time interval of* :math:`\dt.`

Update :math:`\rho_0`, saving the time-centered density at radial edges by calling

**Advect Base Density**\ :math:`[\rho_0^{n},w_0^{\nph,\star}] \rightarrow [\rho_0^{(2a),\star}, \rho_0^{\nph,\star,\pred}]`.

Update :math:`(\rho X_k)` using a discretized version of the species
continuity equation, Eq. \ `[eq:flow:rhoX] <#eq:flow:rhoX>`__, where we omit the
reaction terms, which were already accounted for in **React State**:

.. math::

   \frac{\partial (\rho X_k)}{\partial t} + \nabla \cdot (\Ub \rho X_k) = 0  .
   \label{eq:species}

The update consists of two steps:

#. Compute the time-centered species edge states, :math:`(\rho X_k)^{\nph,\star,\pred}`,
   for the conservative update of :math:`(\rho X_k)^{(1)}`.

   | There are a variety of choices of quantities to predict to the
     edges (controlled by species_pred_type—see Chapter \ `[ch:pert] <#ch:pert>`__).
     By default, we use the equations

     .. math::

        \begin{aligned}
        \frac{\partial\rho'}{\partial t} &=& -\Ub\cdot\nabla\rho' -
             \rho'\nabla\cdot\Ub - \nabla\cdot\left(\rho_0\Ubt\right),
             \label{eq:Perturbational Density}  \\
        \frac{\partial X_k}{\partial t} &=& -\Ub\cdot\nabla X_k +
             \omegadot_k. \label{eq:Primitive Species}\end{aligned}

     to
     predict :math:`\rho^{'(1)} = \rho^{(1)} - \rho_0^n` and
     :math:`X_k^{(1)} = (\rho  X_k)^{(1)} / \rho^{(1)}` to time-centered edges using
     :math:`\Ub = \uadvone+w_0^{\nph,\star}\eb_r`, yielding :math:`\rho^{'\nph,\star,\pred}`
     and :math:`X_k^{\nph,\star,\pred}`.
     We convert the perturbational density to full state density using

     .. math:: \rho^{\nph,\star,\pred} = \rho^{'\nph,\star,\pred} + \frac{\rho_0^n + \rho_0^{(2a),\star}}{2},

     where the base state density terms are mapped to Cartesian edges.
     Then,
   | :math:`(\rho X_k)^{\nph,\star,\pred} = \rho^{\nph,\star,\pred} \, X_k^{\nph,\star,\pred}`.

#. Evolve :math:`(\rho X_k)^{(1)} \rightarrow (\rho X_k)^{(2),\star}` using

   .. math::

      \begin{aligned}
      (\rho X_k)^{(2),\star} &=& (\rho X_k)^{(1)} \nonumber \\
      && - \dt \left\{ \nabla \cdot \left[ \left(\uadvone+w_0^{\nph,\star} \eb_r\right) (\rho X_k)^{\nph,\star,\pred} \right] \right\},\nonumber \\\end{aligned}

   .. math::

      \rho^{(2),\star} = \sum_k (\rho X_k)^{(2),\star},
      \qquad
      X_k^{(2),\star} = (\rho X_k)^{(2),\star} / \rho^{(2),\star}.

Define a radial edge-centered :math:`\etarho^{\nph,\star}` (Eq. `[eq:flow:etarho] <#eq:flow:etarho>`__).

For planar geometry, since :math:`\etarho = \overline{\rho'(\Ub\cdot\eb_r)} = \overline{\rho(\Ub\cdot\eb_r)}-\overline{\rho_0(\Ub\cdot\eb_r}) = \overline{\rho(\Ub\cdot\eb_r)} - \rho_0w_0`,

.. math::

   \begin{aligned}
    \etarho^{\nph,\star} &=&  {\rm {\bf Avg}} \sum_k \left[ \left(\uadvone \cdot \eb_r + w_0^{\nph,\star}\right) (\rho X_k)^{\nph,\star,\pred} \right]\nonumber\\
   && - w_0^{\nph,\star} \rho_0^{\nph,\star,\pred},\end{aligned}

For spherical geometry, first construct
:math:`\etarho^{{\rm cart},\nph,\star} =
[\rho'(\Ub\cdot\eb_r)]^{\nph,\star}` on Cartesian cell centers using:

.. math::

   \begin{aligned}
   \etarho^{{\rm cart},\nph,\star} &=& \left[\left(\frac{\rho^{(1)}+\rho^{(2),\star}}{2}\right)-\left(\frac{\rho_0^n+\rho_0^{(2a),\star}}{2}\right)\right] \nonumber \\
   &&\cdot \left( \uadvone \cdot \eb_r  + w_0^{\nph,\star}\right).\end{aligned}

Then,

.. math:: \etarho^{\nph,\star} = {\rm {\bf Avg}}\left(\etarho^{{\rm cart},\nph,\star}\right).

This gives a radial cell-centered :math:`\etarho^{\nph,\star}`. To get
:math:`\etarho^{\nph,\star}` at radial edges, average the two neighboring
radial cell-centered values.

Correct :math:`\rho_0` by setting :math:`\rho_0^{n+1,\star} =` **Avg**\ :math:`(\rho^{(2),\star})`.

Update :math:`p_0` by calling
**Enforce HSE**\ :math:`[p_0^n,\rho_0^{n+1,\star}] \rightarrow [p_0^{n+1,\star}]`.

Compute :math:`\psi^{\nph,\star}`.

For planar geometry,

.. math::

   \psi_j^{\nph,\star} = \frac{1}{2} \left(\eta_{\rho,j-\myhalf}^{\nph,\star}
   + \eta_{\rho,j+\myhalf}^{\nph,\star}\right) g.

For spherical geometry, first compute:

.. math::

   \begin{aligned}
   \overline{\Gamma_1^{(1)}} &=& {\rm{\bf Avg}} \left[ \Gamma_1\left(\rho^{(1)}, p_0^{n}, X_k^{(1)}\right) \right]  , \\
   \overline{\Gamma_1^{(2),\star}} &=& {\rm{\bf Avg}} \left[ \Gamma_1\left(\rho^{(2),\star}, p_0^{n+1,\star}, X_k^{(2),\star}\right) \right].\end{aligned}

Then, define :math:`\psi^{\nph,\star}` using equation (`[eq:psi def] <#eq:psi def>`__)

.. math::

   \begin{aligned}
   \psi_j^{\nph,\star}
   &= \left(\frac{\overline{\Gamma_1^{(1)}}+\overline{\Gamma_1^{(2),\star}}}{2}\right)_j
   \left(\frac{p_0^n+p_0^{n+1,\star}}{2}\right)_j \nonumber \\
   & \left \{ \overline{S_j^{\nph,\star}} - \frac{1}{r_j^2} \left [ \left(r^2 w_0^{\nph,\star}\right)_{j+\myhalf} - \left(r^2 w_0^{\nph,\star}\right)_{j-\myhalf} \right ] \right \}.\nonumber \\\end{aligned}

| Update :math:`(\rho h)_0`. First, compute :math:`(\rho h)_0^n =` **Avg**\ :math:`[(\rho h)^{(1)}]`.
  Then, call
| **Advect Base Enthalpy**\ :math:`[(\rho h)_0^{n}, w_0^{\nph,\star}, \psi^{\nph,\star}] \rightarrow [(\rho h)_0^{n+1,\star}]`.

Update the enthalpy using a discretized version of the enthalpy
evolution equation (Eq. `[eq:flow:enthalpy] <#eq:flow:enthalpy>`__), again omitting the reaction and heating terms
since we already accounted for
them in **React State**. This equation takes the form:

.. math:: \frac{\partial (\rho h)}{\partial t}  = - \nabla \cdot (\Ub \rho h) + \psi + (\Ubt \cdot \eb_r) \frac{\partial p_0}{\partial r}.

For spherical geometry, we solve the
analytically equivalent form,

.. math:: \frac{\partial (\rho h)}{\partial t}  = - \nabla \cdot (\Ub \rho h) + \psi + \nabla \cdot (\Ubt p_0) - p_0 \nabla \cdot \Ubt,

which experience has shown to minimize the drift from thermodynamic
equilibrium. The update consists of two steps:

Compute the time-centered enthalpy edge state, :math:`(\rho h)^{\nph,\star,\pred},`
for the conservative update of :math:`(\rho h)^{(1)}`. There are a
variety of quantities that we can predict to the interfaces here
(controlled by enthalpy_pred_type—see
Chapter \ `[ch:pert] <#ch:pert>`__). For the default case, we use the
perturbational enthalpy equation, Eq. \ `[eq:flow:rhohprime] <#eq:flow:rhohprime>`__, neglecting reactions,

.. math::

   \frac{\partial(\rho h)'}{\partial t} = -\Ub\cdot\nabla(\rho h)' -
      (\rho h)'\nabla\cdot\Ub - \nabla\cdot\left[(\rho h)_0\Ubt\right] + \Ubt\cdot\nabla p_0
       \label{eq:Perturbational Enthalpy}.

to predict
:math:`(\rho h)' = (\rho h)^{(1)} - (\rho h)_0^n` to time-centered edges,
using :math:`\Ub = \uadvone+w_0^{\nph,\star} \eb_r`,
yielding :math:`(\rho h)^{'\nph,\star,\pred}`. We convert the perturbational
enthalpy to a full state enthalpy using

.. math:: (\rho h)^{\nph,\star,\pred} = (\rho h)^{'\nph,\star,\pred} + \frac{(\rho h)_0^n + (\rho h)_0^{n+1,\star}}{2}.

For planar geometry, we map :math:`(\rho h)_0` directly to Cartesian edges.
In spherical geometry, our experience has shown that a slightly different
approach leads to reduced discretization errors. We first map
:math:`h_0 \equiv (\rho h)_0/\rho_0` and :math:`\rho_0` to Cartesian edges separately,
and then multiply these terms to get :math:`(\rho h)_0`.

Evolve :math:`(\rho h)^{(1)} \rightarrow (\rho h)^{(2),\star}`.

For planar geometry,

.. math::

   \begin{aligned}
   (\rho h)^{(2),\star}
   &= (\rho h)^{(1)} \nonumber \\
   &- \dt \left\{ \nabla \cdot \left[ \left(\uadvone+w_0^{\nph,\star} \eb_r\right) (\rho h)^{\nph,\star,\pred} \right] \right\} \nonumber \\
   & + \dt \left(\uadvone \cdot \eb_r\right) \left(\frac{\partial p_0}{\partial r} \right)^{n} + \dt \psi^{\nph,\star},\end{aligned}

For spherical geometry,

.. math::

   \begin{aligned}
   (\rho h)^{(2),\star}
   &= (\rho h)^{(1)} \nonumber \\
   &- \dt \left\{ \nabla \cdot \left[ \left(\uadvone+w_0^{\nph,\star} \eb_r\right) (\rho h)^{\nph,\star,\pred} \right] \right\} \nonumber \\
   & + \dt \left \{ \nabla \cdot \left (\uadvone p_0^{n} \right ) - p_0^{n} \nabla \cdot \uadvone \right \} \nonumber \\
   &+ \dt \psi^{\nph,\star},\end{aligned}

Then, for each Cartesian cell where :math:`\rho^{(2),\star} < \rho_\mathrm{cutoff}`,
we recompute enthalpy using

.. math:: (\rho h)^{(2),\star} = \rho^{(2),\star}h\left(\rho^{(2),\star},p_0^{n+1,\star},X_k^{(2),\star}\right).

This behavior is controlled by ``do_eos_h_above_cutoff``.

Update the temperature using the equation of state:
:math:`T^{(2),\star} = T(\rho^{(2),\star}, h^{(2),\star}, X_k^{(2),\star})` (planar geometry) or
:math:`T^{(2),\star} = T(\rho^{(2),\star}, p_0^{n+1,\star}, X_k^{(2),\star})` (spherical geometry).

As before, this behavior is controlled by ``use_tfromp``.

*React the full state through a second time interval of* :math:`\dt / 2.`

| Call **React State**\ :math:`[ \rho^{(2),\star},(\rho h)^{(2),\star}, X_k^{(2),\star}, T^{(2),\star},(\rho\Hext)^{(2),\star}, p_0^{n+1,\star}]`
| :math:`\rightarrow [ \rho^{n+1,\star},(\rho h)^{n+1,\star}, X_k^{n+1,\star}, T^{n+1,\star}, (\rho \omegadot_k)^{(2),\star}, (\rho \Hnuc)^{(2),\star} ].`

*Compute the time-centered expansion,* :math:`S^{\nph,\star}`, *base state
velocity,* :math:`w_0^{\nph}`, *and base state velocity forcing.*

Compute :math:`S^{\nph,\star}`. First, compute :math:`S^{n+1,\star}` with

.. math::

   S^{n+1,\star} =  -\sigma  \sum_k  \xi_k  (\omegadot_k)^{(2),\star}  +
      \frac{1}{\rho^{n+1,\star} p_\rho} \sum_k p_{X_k}  ({\omegadot}_k)^{(2),\star} +
      \sigma \Hnuc^{(2),\star} + \sigma \Hext^{(2),\star},

where :math:`(\omegadot_k)^{(2),\star} = (\rho \omegadot_k)^{(2),\star} / \rho^{(2),\star}`
and the thermodynamic quantities are defined using
:math:`\rho^{n+1,\star}, X_k^{n+1,\star},` and :math:`T^{n+1,\star}` as inputs to
the equation of state. If we are using diffusion then we would include the diffusion
term in :math:`S`. Then, define

.. math::

   \overline{S^{\nph,\star}} = {\mathrm{\bf Avg}} (S^{\nph,\star}),
   \qquad
    S^{\nph.\star} = \frac{S^n + S^{n+1,\star}}{2},

Compute :math:`w_0^{\nph}`. First, define

.. math::

   \overline{\Gamma_1^{\nph,\star}} = \frac{\overline{\Gamma_1^n} + \overline{\Gamma_1^{n+1,\star}}}{2},
   \quad
   \rho_0^{\nph,\star} = \frac{\rho_0^{n} + \rho_0^{n+1,\star}}{2},
   \quad
   p_0^{\nph,\star} = \frac{p_0^{n} + p_0^{n+1,\star}}{2},

with

.. math:: \overline{\Gamma_1^{n+1,\star}} = {\rm{\bf Avg}} \left[ \Gamma_1\left(\rho^{n+1,\star}, p_0^{n+1,\star}, X_k^{n+1,\star}\right) \right].

| For planar geometry, call
| **Compute** :math:`w_0` **Planar**\ :math:`[\overline{S^{\nph,\star}},\overline{\Gamma_1^{\nph,\star}},p_0^{\nph,\star},\psi^{\nph,\star}]\rightarrow [w_0^{\nph}]`.

| For spherical geometry, call
| **Compute** :math:`w_0` **Spherical**\ :math:`[\overline{S^{\nph,\star}},\overline{\Gamma_1^{\nph,\star}},\rho_0^{\nph,\star},p_0^{\nph,\star},\etarho^{\nph,\star}]\rightarrow [w_0^{\nph}]`.

Compute the base state velocity forcing. Rearrange equation (`[eq:pizero] <#eq:pizero>`__),

.. math::

   \left ( \frac{\beta_0}{\rho_0} \frac{\partial (\pi_0/\beta_0)}{\partial r} \right )^n =
   -\frac{w_0^{\nph} - w_0^\nmh}{\myhalf(\dt^n+\dt^{n-1})}
   - w_0^n \left(\frac{\partial w_0}{\partial r}\right)^n,

where :math:`w_0^{n}` and :math:`(\partial w_0 / \partial r)^{n}` are defined as

.. math::

   \begin{aligned}
   w_0^n &=& \frac{\dt^{n} w_0^{\nmh} + \dt^{n-1} w_0^{\nph}}{\dt^n+\dt^{n-1}}, \\
   \left(\frac{\partial w_0}{\partial r}\right)^{n} &=& \frac{1}{\dt^n+\dt^{n-1} } \left [ \dt^{n} \left(\frac{\partial w_0 }{ \partial r}\right)^{\nmh} + \dt^{n-1} \left(\frac{\partial w_0 }{ \partial r}\right)^{\nph} \right ].\nonumber \\\end{aligned}

If :math:`n=0`, we use :math:`\dt^{-1} = \dt^0`.

*Construct the time-centered advective velocity on edges,* :math:`\uadvtwo`.

The procedure to construct :math:`\uadvtwodag` is identical to the procedure
for computing :math:`\uadvonedag` in **Step 3**, but uses
the updated values :math:`w_0^{\nph}` and :math:`\pi_0^n` rather than :math:`w_0^{\nph,\star}`
and :math:`\pi_0^{n,\star}`. We note that :math:`\uadvtwo` satisfies the discrete version of
:math:`\overline{(\uadvtwo\cdot\eb_r)}=0` as well as

.. math::

   \nabla \cdot \left(\beta_0^{\nph,\star} \uadvtwo\right) =
   \beta_0^{\nph,\star}\left(S^{\nph,\star} - \overline{S^{\nph,\star}}\right),

.. math::

   \beta_0^{\nph,\star} = \frac{ \beta_0^n +  \beta_0^{n+1,\star} }{2};
   \qquad
    \beta_0^{n+1,\star} = \beta_0 \left(\rho_0^{n+1,\star}, p_0^{n+1,\star}, \overline{\Gamma_1^{n+1,\star}}\right).

*Advect the base state and full state through a time interval of* :math:`\dt.`

Update :math:`\rho_0`, saving the time-centered density at radial edges by calling

**Advect Base Density**\ :math:`[\rho_0^{n},w_0^{\nph}] \rightarrow [\rho_0^{(2a)}, \rho_0^{\nph,\pred}]`.

Update :math:`(\rho X_k)`. This step is identical to **Step 4B** except we use
the updated values :math:`w_0^{\nph}, \uadvtwo`, and :math:`\rho_0^{(2a)}` rather than
:math:`w_0^{\nph,\star}, \uadvone`, and :math:`\rho_0^{(2a),\star}`. In particular:

#. Compute the time-centered species edge states, :math:`(\rho X_k)^{\nph,\pred}`,
   for the conservative update of :math:`(\rho X_k)^{(1)}`. We use equations
   (`[eq:Perturbational Density] <#eq:Perturbational Density>`__) and (`[eq:Primitive Species] <#eq:Primitive Species>`__) to
   predict :math:`\rho^{'(1)} = \rho^{(1)} - \rho_0^n` and
   :math:`X_k^{(1)} = (\rho  X_k)^{(1)} / \rho^{(1)}` to time-centered edges
   with :math:`\Ub = \uadvtwo+w_0^{\nph} \eb_r`,
   yielding :math:`\rho^{'\nph,\pred}` and :math:`X_k^{\nph,\pred}`.
   We convert the perturbational density to a full state density using

   .. math:: \rho^{\nph,\pred} = \rho^{'\nph,\pred} + \frac{\rho_0^n + \rho_0^{(2a)}}{2}.

   Then, :math:`(\rho X_k)^{\nph,\pred} = \rho^{\nph,\pred} \, X_k^{\nph,\pred}`.

#. Evolve :math:`(\rho X_k)^{(1)} \rightarrow (\rho X_k)^{(2)}` using

   .. math::

      (\rho X_k)^{(2)} = (\rho X_k)^{(1)}
      - \dt \left\{ \nabla \cdot \left[\left(\uadvtwo+w_0^{\nph} \eb_r\right)
      (\rho X_k)^{\nph,\pred} \right] \right\},

   .. math::

      \rho^{(2)} = \sum_k (\rho X_k)^{(2)},
      \qquad
      X_k^{(2)} = (\rho X_k)^{(2)} / \rho^{(2)}.

Define a radial edge-centered :math:`\etarho^{\nph}`.

For planar geometry,

.. math::

   \begin{aligned}
    \etarho^{\nph} &=& {\rm {\bf Avg}} \sum_k \left [\left(\uadvtwo \cdot \eb_r + w_0^{\nph}\right) (\rho X_k)^{\nph,\pred} \right] \nonumber \\
   &&- w_0^{\nph} \rho_0^{\nph,\pred},\end{aligned}

For spherical geometry, first construct
:math:`\etarho^{{\rm cart},\nph} = [\rho'(\Ub\cdot\eb_r)]^{\nph}` on Cartesian
cell centers using:

.. math:: \etarho^{{\rm cart},\nph} = \left[\left(\frac{\rho^{(1)}+\rho^{(2)}}{2}\right)-\left(\frac{\rho_0^n+\rho_0^{(2a)}}{2}\right)\right] \left(\uadvtwo \cdot \eb_r + w_0^{\nph}\right).

Then,

.. math:: \etarho^{\nph} = {\rm {\bf Avg}}\left(\etarho^{{\rm cart},\nph}\right).

This gives a radial cell-centered :math:`\etarho^{\nph}`. To get
:math:`\etarho^{\nph}` at radial edges, average the two neighboring
cell-centered values.

Correct :math:`\rho_0` by setting :math:`\rho_0^{n+1} =` **Avg**\ :math:`(\rho^{(2)})`.

Update :math:`p_0` by calling
**Enforce HSE**\ :math:`[p_0^n,\rho_0^{n+1}] \rightarrow [p_0^{n+1}]`.

Compute :math:`\psi^{\nph}`.

For planar geometry,

.. math::

   \psi_j^{\nph} = \frac{1}{2} \left(\eta_{\rho,j-\myhalf}^{\nph}
   + \eta_{\rho,j+\myhalf}^{\nph}\right) g.

For spherical geometry, first compute:

.. math::

   \overline{\Gamma_1^{(2)}} = {\rm{\bf Avg}} \left[ \Gamma_1\left(\rho^{(2)}, p_0^{n+1},
   X_k^{(2)}\right) \right].

Then, define :math:`\psi^{\nph}` using equation (`[eq:psi def] <#eq:psi def>`__):

.. math::

   \begin{aligned}
   \psi_j^{\nph}
   &= \left(\frac{\overline{\Gamma_1^{(1)}}+\overline{\Gamma_1^{(2)}}}{2}\right)_j \left(\frac{p_0^n+p_0^{n+1}}{2}\right)_j \nonumber \\
   & \left \{ \overline{S_j^{\nph}} - \frac{1}{r_j^2} \left [ \left(r^2 w_0^{\nph}\right)_{j+\myhalf} - \left(r^2 w_0^{\nph}\right)_{j-\myhalf} \right ] \right \}.\end{aligned}

Update :math:`(\rho h)_0` by calling
**Advect Base Enthalpy**\ :math:`[(\rho h)_0^n, w_0^{\nph}, \psi^{\nph}] \rightarrow [(\rho h)_0^{n+1}]`.

| Update the enthalpy. This step is identical to **Step 4H** except we use
  the updated values :math:`w_0^{\nph}, \uadvtwo, \rho_0^{n+1}, (\rho h)_0^{n+1}, p_0^{n+\myhalf}`,
  and :math:`\psi^{n+\myhalf}` rather than
| :math:`w_0^{\nph,\star}, \uadvone, \rho_0^{n+1,\star}, (\rho h)_0^{n+1,\star}, p_0^n`,
  and :math:`\psi^{n+\myhalf,\star}`. In particular:

Compute the time-centered enthalpy edge state, :math:`(\rho h)^{\nph,\pred},`
for the conservative update of :math:`(\rho h)^{(1)}`. We use equation
(`[eq:Perturbational Enthalpy] <#eq:Perturbational Enthalpy>`__) to predict
:math:`(\rho h)' = (\rho h)^{(1)} - (\rho h)_0^n` to time-centered edges
with :math:`\Ub = \uadvtwo+w_0^{\nph} \eb_r`,
yielding :math:`(\rho h)^{'\nph,\pred}`.
We convert the perturbational enthalpy to a full state enthalpy using

.. math:: (\rho h)^{\nph,\pred} = (\rho h)^{'\nph,\pred} + \frac{(\rho h)_0^n + (\rho h)_0^{n+1}}{2}.

Evolve :math:`(\rho h)^{(1)} \rightarrow (\rho h)^{(2)}`.

For planar geometry,

.. math::

   \begin{aligned}
   (\rho h)^{(2)}
   &=& (\rho h)^{(1)} - \dt \left\{ \nabla \cdot \left[ \left(\uadvtwo+w_0^{\nph} \eb_r\right)  (\rho h)^{\nph,\pred} \right] \right\} \nonumber \\
   && + \dt \left(\uadvtwo \cdot \eb_r\right) \left(\frac{\partial p_0}{\partial r} \right)^\nph + \dt \psi^{\nph},\end{aligned}

For spherical geometry,

.. math::

   \begin{aligned}
   (\rho h)^{(2)}
   &=& (\rho h)^{(1)} - \dt \left\{ \nabla \cdot \left[ \left(\uadvtwo+w_0^{\nph} \eb_r\right)  (\rho h)^{\nph,\pred} \right] \right\} \nonumber \\
   && + \dt \left[ \nabla \cdot \left (\uadvtwo p_0^{\nph} \right ) - p_0^{\nph} \nabla \cdot \uadvtwo \right] + \dt \psi^{\nph},\nonumber \\\end{aligned}

where :math:`p_0^\nph` is defined as :math:`p_0^\nph = (p_0^n+p_0^{n+1})/2`.

Then, for each Cartesian cell where :math:`\rho^{(2)} < \rho_\mathrm{cutoff}`, we recompute enthalpy using

.. math:: (\rho h)^{(2)} = \rho^{(2)}h\left(\rho^{(2)},p_0^{n+1},X_k^{(2)}\right).

Update the temperature using the equation of state:
:math:`T^{(2)} = T(\rho^{(2)}, h^{(2)}, X_k^{(2)})` (planar geometry) or
:math:`T^{(2)} = T(\rho^{(2)}, p_0^{n+1}, X_k^{(2)})` (spherical geometry).

Again, the actual inputs depend of ``use_tfromp``.

*React the full state through a second time interval of* :math:`\dt / 2.`

| Call **React State**\ :math:`[\rho^{(2)},(\rho h)^{(2)}, X_k^{(2)},T^{(2)}, (\rho\Hext)^{(2)}, p_0^{n+1}]`
| :math:`\rightarrow [\rho^{n+1}, (\rho h)^{n+1}, X_k^{n+1}, T^{n+1}, (\rho \omegadot_k)^{(2)}, (\rho \Hnuc)^{(2)} ].`

*Define the new time expansion,* :math:`S^{n+1}`, *and* :math:`\overline{\Gamma_1^{n+1}}`.

#. Define

   .. math::

      S^{n+1} =  -\sigma  \sum_k  \xi_k (\omegadot_k)^{(2)}  + \sigma \Hnuc^{(2)} +
        \frac{1}{\rho^{n+1} p_\rho} \sum_k p_{X_k}  ({\omegadot}_k)^{(2)}
         + \sigma \Hext^{(2)},

   where :math:`(\omegadot_k)^{(2)} = (\rho \omegadot_k)^{(2)} / \rho^{(2)}`
   and the thermodynamic quantities are defined using :math:`\rho^{n+1}`,
   :math:`X_k^{n+1}`, and :math:`T^{n+1}` as inputs to the equation of state.
   If we are doing thermal diffusion (use_thermal_diffusion= T)
   then we also include the diffusive term in :math:`S`.
   Then, compute

   .. math:: \overline{S^{n+1}} = {\mathrm{\bf Avg}} (S^{n+1}).

#. Define

   .. math::

      \overline{\Gamma_1^{n+1}} = {\rm{\bf Avg}}\left[\Gamma_1\left(\rho^{n+1}, p_0^{n+1},
      X_k^{n+1}\right) \right].

*Update the velocity*.

First, we compute the time-centered edge velocities, :math:`\Ubt^{\nph,\pred}`.
Then, we define

.. math:: \rho^\nph = \frac{\rho^n + \rho^{n+1}}{2}, \qquad \rho_0^\nph = \frac{\rho_0^n + \rho_0^{n+1}}{2}.

We update the velocity field :math:`\Ubt^n` to :math:`\Ubt^{n+1,\dagger}` by discretizing
equation (`[eq:flow:utildeupd] <#eq:flow:utildeupd>`__) as

.. math::

   \begin{aligned}
   \Ubt^{n+1,\dagger}
   &= \Ubt^n - \dt \left[\left(\uadvtwo+ w_0^{\nph} \eb_r\right) \cdot \nabla \Ubt^{\nph,\pred} \right] \nonumber \\
   &- \dt \left(\uadvtwo \cdot \eb_r\right)  \left(\frac{\partial w_0}{\partial r} \right)^\nph \eb_r \nonumber \\
   & + \dt \left[ - \frac{\beta_0^\nph}{\rho^\nph} \mathbf{G} \left ( \frac{\pi}{\beta_0}\right)^\nmh + \left(\frac{\beta_0}{\rho_0}\frac{\partial(\pi_0/\beta_0)}{\partial r}\right)^n \eb_r - \frac{\left(\rho^\nph-\rho_0^\nph\right)}{\rho^\nph} g^{\nph} \eb_r \right],\nonumber \\\end{aligned}

where :math:`\mathbf{G}` approximates a cell-centered gradient from nodal
data. Again, the :math:`\dagger` superscript refers
to the fact that the updated velocity does not satisfy the divergence
constraint.

Finally, we use an approximate nodal projection to define :math:`\Ubt^{n+1}`
from :math:`\Ubt^{n+1,\dagger},` such that :math:`\Ubt^{n+1}` approximately
satisfies

.. math::

   \nabla \cdot \left(\beta_0^{\nph} \Ubt^{n+1} \right)
   = \beta_0^{\nph} \left(S^{n+1} - \overline{S^{n+1}} \right),

where :math:`\beta_0^{\nph}` is defined as

.. math::

   \beta_0^{\nph} = \frac{\beta_0^n + \beta_0^{n+1}}{2}; \qquad
   \beta_0^{n+1} = \beta \left(\rho_0^{n+1}, p_0^{n+1}, \overline{\Gamma_1^{n+1}}, g^{n+1}\right).

As part of the projection we also define the new-time perturbational pressure,
:math:`\pi^\nph.` This projection necessarily differs from the MAC projection used in
**Step 3** and **Step 7** because the velocities in those steps are defined
on edges and :math:`\Ubt^{n+1}` is defined at cell centers, requiring different divergence
and gradient operators. Details of the approximate projection are given in Paper III.

*Compute a new* :math:`\dt.`

Compute :math:`\dt` for the next time step with the procedure described in
§3.4 of Paper III using :math:`w_0` as computed in **Step 6** and :math:`\Ubt^{n+1}`
as computed in **Step 11**.

This completes one step of the algorithm.

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

[sec:flow:gamma1vary] :math:`\Gamma_1` Variation Changes
--------------------------------------------------------

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
               \left [ \frac{\partial p_0}{\partial t} + \Ub \cdotb \nablab p_0 \right ]}_{\mbox{second order corrections}}  , \label{eq:gammafull}

Keeping to First Order in :math:`\delta\Gamma_1`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The base state evolution equation is the average of Eq. \ `[eq:gammafull] <#eq:gammafull>`__ over a layer

.. math::

   \nablab \cdot w_0 \er + \frac{1}{\gammabar p_0} w_0 \er \cdotb \nablab p_0 =
   \Sbar - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} + \overline{
     \left ( \frac{\delta \Gamma_1}{\gammabar^2 p_0} \Ubt \cdotb \nablab p_0
     \right ) }  .

where we see that the :math:`[\delta \Gamma_1/(\Gamma_1^2 p_0)] \partial p_0/\partial t` terms averages to zero, since the average of :math:`\delta\Gamma_1` term is zero.
Subtracting this from equation (`[eq:gammafull] <#eq:gammafull>`__), we have

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
Eq. \ `[eq:flow:dw0dr_planar] <#eq:flow:dw0dr_planar>`__), and

.. math::

   \nablab \cdotb (\beta_0 \Ubt) = \beta_0 \left [ S - \Sbar + \frac{\delta
       \Gamma_1}{\gammabar^2 p_0} \psi + \frac{\delta \Gamma_1}{\gammabar^2 p_0}
     \Ubt \cdotb \nablab p_0 - \overline{ \left ( \frac{\delta
         \Gamma_1}{\gammabar^2 p_0} \Ubt \cdotb \nablab p_0 \right ) } ~ \right ]
     ,  \label{eq:constraint_with_delta_gamma}

This constraint is not in a form that can be projected. To solve this
form, we need to use a lagged :math:`\Ubt` in the righthand side.

This change comes into MAESTROeX in a variety of steps, summarized here.
To enable this portion of the algorithm, set use_delta_gamma1_term = T.

-  In **Step 3**, we are doing the “predictor” portion of the
   MAESTROeX algorithm, getting the MAC velocity that satisfies the constraint,
   so we do not try to incorporate the :math:`\delta \Gamma_1` effect. We set
   all the :math:`\delta \Gamma_1` terms in
   Eq. \ `[eq:constraint_with_delta_gamma] <#eq:constraint_with_delta_gamma>`__ to zero.

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
   in Eq. \ `[eq:constraint_with_delta_gamma] <#eq:constraint_with_delta_gamma>`__ because we modified
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

Thermal diffusion was introduced in the XRB paper :raw-latex:`\cite{xrb}`. This
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

   \begin{aligned}
   (\rho h)^{(2),\star} &=& (\rho h)^{(1a),\star} + \frac{\dt}{2}\nabla\cdot\left(\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(2),\star} + \frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(2),\star} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n+1,\star} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),\end{aligned}

which is numerically implemented as a diffusion equation for :math:`h^{(2),\star}`,

.. math::

   \begin{aligned}
   \left(\rho^{(2),\star} - \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla\right)h^{(2),\star} &=& (\rho h)^{(1a),\star} + \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\nonumber\\
   &&- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(2),\star} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n+1,\star} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),\end{aligned}

Immediately after **Step 8H**, diffuse the enthalpy through a time interval of
:math:`\dt`. First, define :math:`(\rho h)^{(1a)} = (\rho h)^{(2)}`. We recompute :math:`(\rho h)^{(2)}` to
account for thermal diffusion. Compute :math:`\kth^{(2),\star}, c_p^{(2),\star}`, and
:math:`\xi_k^{(2),\star}`, from :math:`\rho^{(2),\star}, T^{(2),\star}`, and :math:`X_k^{(2),\star}` as inputs to
the equation of state. The update is given by

.. math::

   \begin{aligned}
   (\rho h)^{(2)} &=& (\rho h)^{(1a)} + \frac{\dt}{2}\nabla\cdot\left(\frac{\kth^{(2),\star}}{c_p^{(2),\star}}\nabla h^{(2)} + \frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla X_k^{(2)} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla p_0^{n+1} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),\end{aligned}

which is numerically implemented as a diffusion equation for :math:`h^{(2)}`.

.. math::

   \begin{aligned}
   \left(\rho^{(2)} - \frac{\dt}{2}\nabla\cdot\frac{\kth^{(2),\star}}{c_p^{(2),\star}}\nabla\right)h^{(2)} &=& (\rho h)^{(1a)} + \frac{\dt}{2}\nabla\cdot\frac{\kth^{(1)}}{c_p^{(1)}}\nabla h^{(1)}\nonumber\\
   &&- \frac{\dt}{2}\sum_k\nabla\cdot\left(\frac{\xi_k^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla X_k^{(2)} + \frac{\xi_k^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla X_k^{(1)}\right)\nonumber\\
   &&- \frac{\dt}{2}\nabla\cdot\left(\frac{h_p^{(2),\star}\kth^{(2),\star}}{c_p^{(2),\star}}\nabla p_0^{n+1} + \frac{h_p^{(1)}\kth^{(1)}}{c_p^{(1)}}\nabla p_0^{n}\right),\end{aligned}

.. raw:: latex

   \centering

.. figure:: \flowfigpath/flowchart
   :alt: [Fig:flowchart] A flowchart of the algorithm. The
   thermodynamic state variables, base state variables, and local velocity are
   indicated in each step. Red text indicates that quantity was
   updated during that step. The predictor-corrector steps are
   outlined by the dotted box. The blue text indicates state
   variables that are the same in **Step 6** as they are in
   **Step 2**, i.e., they are unchanged by the predictor steps.
   The diffusion steps (4a and 8a) are optional, depending on
   use_thermal_diffusion.

   [Fig:flowchart] A flowchart of the algorithm. The
   thermodynamic state variables, base state variables, and local velocity are
   indicated in each step. Red text indicates that quantity was
   updated during that step. The predictor-corrector steps are
   outlined by the dotted box. The blue text indicates state
   variables that are the same in **Step 6** as they are in
   **Step 2**, i.e., they are unchanged by the predictor steps.
   The diffusion steps (4a and 8a) are optional, depending on
   use_thermal_diffusion.

.. raw:: latex

   \centering

.. figure:: \flowfigpath/flowchart_4_8
   :alt: [Fig:flowchart48] A flowchart for **Steps 4** and **8**.
   The thermodynamic state variables and base state variables are
   indicated in each step. Red text indicates that quantity was
   updated during that step. Note, for **Step 4**, the updated
   quantities should also have a :math:`\star` superscript, e.g., **Step
   8I** defines :math:`T^{(2)}` while **Step 4I** defines :math:`T^{(2),\star}`
   .

   [Fig:flowchart48] A flowchart for **Steps 4** and **8**.
   The thermodynamic state variables and base state variables are
   indicated in each step. Red text indicates that quantity was
   updated during that step. Note, for **Step 4**, the updated
   quantities should also have a :math:`\star` superscript, e.g., **Step
   8I** defines :math:`T^{(2)}` while **Step 4I** defines :math:`T^{(2),\star}`
   .

.. _Sec:Initialization:

Initialization
==============

[sec:flow:initialization]

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

   #. Compute :math:`S^{0,\nu}` from equation (`[eq:defineS] <#eq:defineS>`__)
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
      is discussed in :raw-latex:`\cite{almgren:bell:crutchfield}`.)
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

The tolerances for these elliptic solves are described in § \ `[sec:mgtol] <#sec:mgtol>`__.

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

Changes Between the Multilevel Paper and Paper 5 :raw-latex:`\cite{wdconvect}`
------------------------------------------------------------------------------

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
   Eq. \ `[eq:flow:newmomentum] <#eq:flow:newmomentum>`__ to conserve the low-Mach number form of
   energy.

#. We changed the form of the volume discrepancy term to get better
   agreement between the two temperatures.

Future Considerations
=====================

-  Should we use a predictor-corrector for updating the full-state density?
   Specifically, after calling **Correct Base**, should we do a full-state density
   advance and **Correct Base** using the more accurate estimate of :math:`\rho_0^{n+1}`?

-  We are still exploring the effects of ``use_tfromp = F`` for spherical
   problems. We would eventually like to run in this mode, but :math:`T=T(\rho,X_k,p_0)`
   and :math:`T=T(\rho,h,X_k)` drift away from each other more than we would like. Our
   attempts at incorporating a ``dpdt_factor`` for spherical problems have not
   been successful.
