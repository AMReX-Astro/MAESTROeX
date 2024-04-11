.. _ch:intro:

*************************
Introduction to MAESTROeX
*************************

History of MAESTROeX
====================

MAESTROeX models the evolution of low Mach number astrophysical
flows that are in hydrostatic equilibrium.  The name MAESTROeX itself
derives from MAESTRO, the original (pure Fortran) low Mach number
stellar hydrodynamics code.  MAESTROeX began as a rewrite MAESTRO
in the C++/Fortran framework of AMReX, and has since then gained
additional algorithmic capabilities.  Henceforth, we will refer
only to MAESTROeX, although many of the ideas originated in
MAESTRO before the change to MAESTROeX.

The idea for MAESTROeX grew out of our success in applying low Mach
number combustion methods developed for terrestrial flames
:cite:`DayBell00` to small-scale (non-stratified) astrophysical
flames such as the Landau-Darrieus instability :cite:`SNld`,
Rayleigh-Taylor unstable flames :cite:`SNrt3d`, and flame-turbulence
interactions :cite:`SNturb`. Our original small-scale astrophysical
combustion algorithm is detailed in

-  *Adaptive Low Mach Number Simulations of Nuclear Flames,*
   J. B. Bell, M. S. Day, C. A. Rendleman, S. E. Woosley, & M. Zingale
   2004, JCP, 195, 2, 677 (henceforth BDRWZ) :cite:`SNe`

MAESTROeX was developed initially for modeling the convective phase in
a Chandrasekhar mass white dwarf preceding the ignition of a
Type Ia supernovae.  As
such, we needed to incorporate the compressibility effects due to
large-scale stratification in the star. The method closest in spirit
to MAESTROeX is the pseudo-incompressible method of Durran
:cite:`durran`, developed for terrestrial atmospheric flows (assuming
an ideal gas). Part of the complexity of the equations in MAESTROeX
stems from the need to describe a general equation of state.
Additionally, since reactions can significantly alter the hydrostatic
structure of a star, we incorporated extensions that capture the
expansion of the background state :cite:`almgren:2000`. The low Mach
number equations for stellar flows were developed in a series of
papers leading up to the first application to this problem:

-  *Low Mach Number Modeling of Type Ia
   Supernovae. I. Hydrodynamics,* A. S. Almgren, J. B. Bell,
   C. A. Rendleman, & M. Zingale 2006, ApJ, 637, 922 (henceforth
   paper I) :cite:`lowMach`

-  *Low Mach Number Modeling of Type Ia Supernovae. II. Energy
   Evolution,* A. S. Almgren, J. B. Bell, C. A. Rendleman, & M. Zingale
   2006, ApJ, 649, 927 (henceforth paper II) :cite:`lowMach2`

-  *Low Mach Number Modeling of Type Ia Supernovae. III. Reactions,*
   A. S. Almgren, J. B. Bell, A. Nonaka, & M. Zingale
   2008, ApJ, 684, 449 (henceforth paper III) :cite:`lowMach3`

-  *Low Mach Number Modeling of Type Ia Supernovae. IV. White Dwarf Convection,*
   M. Zingale, A. S. Almgren, J. B. Bell, A. Nonaka, & S. E. Woosley
   2009, ApJ, 704, 196 (henceforth paper IV) :cite:`lowMach4`

The adaptive mesh refinement version of the algorithm was presented in the
next papers in the series:

-  *MAESTRO: An Adaptive Low Mach Number Hydrodynamics Algorithm for Stellar
   Flows,* A. Nonaka, A. S. Almgren, J. B. Bell, M. J. Lijewski, C. M. Malone,
   & M. Zingale 2010, ApJS, 188, 358 (henceforth “the multilevel paper”) :cite:`multilevel`

The most recent developments for MAESTROeX are described in:

-  *MAESTROeX: A Massively Parallel Low Mach Number Astrophysical Solver,* D. Fan,
   A. Nonaka, A. S. Almgren, A. Harpole, & M. Zingale, 2019, submitted to ApJ
   :cite:`MAESTROeX`

We have many papers that describe applications of the method to Type
Ia supernovae, X-ray bursts, and stellar evolution. These are listed
on the main MAESTROeX website.  Some of these papers have appendices
that describe enhancements to the code—these are noted below.

-  *Multidimensional Modeling of Type I X-ray Bursts. I. Two-dimensional
   Convection Prior to the Outburst of a Pure 4He Accretor,*
   C. M. Malone, A. Nonaka, A. S. Almgren, J. B. Bell, & M. Zingale 2011,
   ApJ, 728, 118 (henceforth “the XRB paper”) :cite:`xrb`

   This introduces the thermal diffusion portion of the MAESTROeX algorithm.

-  *Comparisons of Two- and Three-Dimensional Convection in
   Type I X-ray Bursts* M. Zingale, C. M. Malone, A. Nonaka,
   A. S. Almgren, & J. B. Bell 2015, ApJ, 807, 60. :cite:`xrb3d`

   This has an appendix that describes the Godunov state construction in more
   detail than previous papers.

-  *Low Mach Number Modeling of Convection in Helium Shells on
   Sub-Chandrasekhar White Dwarfs II: Bulk Properties of Simple Models,*
   A. M. Jacobs, M. Zingale, A. Nonaka, A. S. Almgren, & J. B. Bell
   2016, ApJ, 827, 84. :cite:`subch2`

   This has an appendix that shows some test problems for the alternate energy
   formulation in MAESTROeX.

Brief Overview of Low Speed Approximations
==========================================

There are many low speed formulations of the equations of hydrodynamics
in use, each with their own applications. All of these methods share in
common a constraint equation on the velocity field that augments the
equations of motion.

Incompressible Hydrodynamics
----------------------------

The simplest low Mach number approximation is incompressible
hydrodynamics. This approximation is formally the zero
Mach number limit (:math:`M \rightarrow 0`)
of the Navier-Stokes equations. In incompressible hydrodynamics,
the velocity satisfies a constraint equation:

.. math:: \nabla \cdot \Ub = 0

which acts to instantaneously equilibrate the flow, thereby filtering
out soundwaves. The constraint equation implies that

.. math:: D\rho/Dt = 0

(through the continuity equation) which says that the density is
constant along particle paths. This means that there are no
compressibility effects modeled in this approximation.

Anelastic Hydrodynamics
-----------------------

In the anelastic approximation small amplitude thermodynamic
perturbations are carried with respect to a static hydrostatic
background (described by density :math:`\rho_0`). The density perturbation
is ignored in the continuity equation, resulting in a constraint
equation:

.. math:: \nabla \cdot (\rho_0 \Ub) = 0

This properly captures the compressibility effects due to the
stratification of the background. Because there is no evolution
equation for the perturbational density, approximations are made to
the buoyancy term in the momentum equation.

Low-Mach Number Combustion
--------------------------

In the low Mach number combustion model, the pressure is decomposed
into a dynamic, :math:`\pi`, and thermodynamic component, :math:`p_0`, the ratio
of which is :math:`O(M^2)`. The total pressure is replaced everywhere by the
thermodynamic pressure, except in the momentum equation. This
decouples the pressure and density and filters out the sound
waves. Large amplitude density and temperature fluctuations are
allowed. The only requirement is that the total pressure stay close to
the background pressure, which is assumed constant. This requirement
can be expressed as:

.. math:: p = p_0

and differentiating this along particle paths leads to a constraint on
the velocity field:

.. math:: \nabla \cdot \Ub = S

This looks like the constraint for incompressible hydrodynamics, but
now we have a source term, :math:`S`, representing the local compressibility
effects due to the energy generation and thermal diffusion. Since the
background pressure is taken to be constant, we cannot model flows
that cover a large fraction of a pressure scale height. However, this
method is ideal for exploring the physics of flames.

Pseudo-Incompressible Methods
-----------------------------

The pseudo-incompressible method incorporates both the local changes
to compressibility due to reaction/heat release, and the large-scale
changes due to the background stratification. This was originally
derived for an ideal gas equation of state for atmospherical flows.
Allowing the background pressure, :math:`p_0` to vary (e.g. in hydrostatic
equilibrium), differentiating pressure along particle paths gives:

.. math:: \nabla \cdot (p_0^{1/\gamma} \Ub) = H

where :math:`\gamma` is the ratio of specific heats and :math:`H` is the source.

MAESTROeX is based on this method, generalizing this constraint to an
arbitrary equation of state and allowing for the time-variation of the
base state.

Alternate Energy Formulation
----------------------------

Several authors :cite:`KP:2012,VLBWZ:2013` showed that with a slightly
different momentum equation, the low Mach number system can conserve
an energy (that is, a quantity that looks like the compressible
energy, but formed using the low Mach number quantities). This change
manifests itself as either a change to the buoyancy term or by
changing :math:`\nabla \pi` to :math:`\beta_0 \nabla (\pi/\beta_0)`. Furthermore,
:cite:`VLBWZ:2013` showed that the new formulation better captures the
vertical propagation of gravity waves. As of
Dec. 2013, this new formulation is the default in MAESTROeX.

Projection Methods 101
======================

Most astrophysical hydrodynamics codes
(e.g. CASTRO :cite:`castro` or FLASH :cite:`flash`) solve the
compressible Euler equations, which can be written in the form:

.. math:: \Ub_t + \nabla \cdot F(\Ub) = 0

where :math:`\Ub` is the vector of conserved quantities, :math:`\Ub = (\rho, \rho u,
\rho E)`, with :math:`\rho` the density, :math:`u` the velocity, and :math:`E` the total
energy per unit mass. This system of equations can be expressed
as a system of advection equations:

.. math:: {\bf q}_t + A({\bf q}) {\bf q}_x = 0

where :math:`{\bf q}` are called the primitive variables, and :math:`A` is the
Jacobian, :math:`A \equiv \partial F / \partial U`. The eigenvalues of the
matrix :math:`A` are the characteristic speeds—the speeds at which
information propagates. For the Euler equations, these are :math:`u` and :math:`u
\pm c`, where :math:`c` is the sound speed. Solution methods for the
compressible equations make use of this wave-nature to compute fluxes
at the interfaces of grid cells to update the state in time. An
excellent introduction to these methods is provided by LeVeque’s book
:cite:`leveque`. The timestep for these methods is limited by the time
it takes for the maximum characteristic speed to traverse one grid cell.
For very subsonic flows, this means that the timestep is dominated by
the propagation of soundwaves, which may not be important to the
overall dynamics of the flow.

In contrast, solving low Mach number systems (including the equations of
incompressible hydrodynamics) typically involves solving one or more
advection-like equations (representing, e.g. conservation of mass and
momentum) coupled with a divergence constraint on the velocity field.
For example, the equations of constant-density incompressible flow
are:

.. math::
   \Ub_t = -\Ub \cdot \nabla \Ub + \nabla p
   :label: incompressible_u

.. math::
   \nabla \cdot \Ub = 0


Here, :math:`\Ub` represents the velocity vector [1]_ and :math:`p` is
the dynamical pressure. The time-evolution equation for the velocity
(:eq:`incompressible_u`) can be solved
using techniques similar to those developed for compressible
hydrodynamics, updating the old velocity, :math:`\Ub^n`, to the new
time-level, :math:`\Ub^\star`.  Here the ‘:math:`^\star`’ indicates
that the updated velocity does not, in general, satisfy the divergence
constraint. A projection method will take this updated velocity and
force it to obey the constraint equation. The basic idea follows from
the fact that any vector field can be expressed as the sum of a
divergence-free quantity and the gradient of a scalar. For the
velocity, we can write:

.. math::
   \Ub^\star = \Ub^d + \nabla \phi
   :label: decomposition

where :math:`\Ub^d` is the divergence free portion of the velocity vector,
:math:`\Ub^\star`, and :math:`\phi` is a scalar. Taking the divergence of
:eq:`decomposition`, we have

.. math:: \nabla^2 \phi = \nabla \cdot \Ub^\star

(where we used :math:`\nabla \cdot \Ub^d = 0`).
With appropriate boundary conditions, this Poisson equation can be
solved for :math:`\phi`, and the final, divergence-free velocity can
be computed as

.. math:: \Ub^{n+1} = \Ub^\star - \nabla \phi

Because soundwaves are filtered, the timestep constraint now depends only
on :math:`|\Ub|`.

Extensions to variable-density incompressible
flows :cite:`bellMarcus:1992b` involve a slightly different
decomposition of the velocity field and, as a result, a slightly
different Poisson equation. There is also a variety of different ways
to express what is being projected :cite:`almgren:bell:crutchfield`,
and different discretizations of the divergence and gradient operators
lead to slightly different mathematical properties of the methods
(leading to “approximate
projections” :cite:`almgrenBellSzymczak:1996`). Finally, for
second-order methods, two projections are typically done per timestep.
The first (the ‘MAC’ projection :cite:`bellColellaHowell:1991`)
operates on the half-time, edge-centered advective velocities, making
sure that they satisfy the divergence constraint. These advective
velocities are used to construct the fluxes through the interfaces to
advance the solution to the new time. The second/final projection
operates on the cell-centered velocities at the new time, again
enforcing the divergence constraint. The MAESTROeX algorithm performs
both of these projections.

The MAESTROeX algorithm builds upon these ideas, using a different
velocity constraint equation that captures the compressibility
due to local sources and large-scale stratification.
