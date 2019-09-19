*****************************
Spectral Deferred Corrections
*****************************

SDC Overview
============

Spectral deferred corrections (SDC) is an iterative scheme used to integrate
the thermodynamic variables, :math:`(\rho X_k,\rho h)`, over a time step. It has
been shown to be more accurate and efficient than Strang splitting in a
terrestrial, non-stratified, low Mach number reacting flow solver :raw-latex:`\cite{Non11}`,
so we would like to develop an SDC version of MAESTRO.

MAESTRO integrates the following system of equations:

.. math::

   \begin{aligned}
   \frac{\partial\Ub}{\partial t} &=&
       -\Ub\cdot\nabla\Ub  - \frac{1}{\rho}\nabla\pi
       - \frac{\rho-\rho_0}{\rho} g\eb_r,\label{eq:momentum}\\
   \frac{\partial(\rho X_k)}{\partial t} &=&
       -\nabla\cdot(\rho X_k\Ub) + \rho\omegadot_k,\label{eq:species}\\
   \frac{\partial(\rho h)}{\partial t} &=&
       -\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt}
       + \rho\Hnuc + \nabla\cdot k_{\rm th}\nabla T,\label{eq:enthalpy}\end{aligned}

together with base state evolution equations and a constraint equation.
By default, MAESTRO advances the thermodynamic variables by coupling
the different physical processes together (advection, diffusion, reactions) using
Strang splitting. Specifically, we integrate the reaction terms
over half a time step ignoring contributions from advection and diffusion,
then from this intermediate state we integrate the advection and diffusion terms over
a full time step (while ignoring reactions), and finally
integrate the reactions over a half time step (ignoring advection and diffusion).
For problems where
the reactions and/or diffusion greatly alter the energy balance or composition
as compared to advection, this operator splitting approach can lead to large
splitting errors and highly inaccurate solutions.
This issue is can be particularly exasperating
for low Mach number methods that can take large advection-based time steps.

An alternate approach to advancing the thermodynamic variables is SDC.
SDC is an iterative approach to couple the various processes
together, with each process seeing an
approximation of the other processes as a source term. The SDC
algorithm converges to an integral representation of the solution in
time that couples all of the processes together in a self-consistent
fashion, see :raw-latex:`\cite{Non11}`.

As a first attempt, we will work on coupling advection and reactions only
via SDC, with a base state that is fixed in time, but not space.

Strang-Splitting Without Thermal Diffusion or Base State Evolution
==================================================================

In the Strang splitting version of MAESTRO, the reaction and advection
processes operate independent of one-another. The species and
enthalpy equations are integrated over :math:`\Delta t` using the following
sequence:

-  React for :math:`\Delta t/2`:

   In the Strang splitting version of MAESTRO, the reaction network solves just
   the reaction portion of the species evolution equations along with a
   temperature evolution equation. This is done using a standard stiff ODE solver
   package (like VODE) on the system:

   .. math::

      \begin{aligned}
      \frac{dX_k}{dt} &=& \omegadot_k(\rho,X_k,T), \\
      \frac{dT}{dt}   &=&
          \frac{1}{c_p} \left ( -\sum_k \xi_k  \omegadot_k  + \Hnuc \right ).\end{aligned}

   Here, :math:`T` is evolved solely to evaluate the reaction rates,
   :math:`\omegadot_k(\rho,X_k,T)`. Furthermore, we simplify the problem
   “freezing” the thermodynamics—i.e., :math:`c_p` and :math:`\xi_k` are evaluated at the
   start of the integration and held constant thereafter.
   The density remains constant during this step, i.e.,
   :math:`\rho^\mathrm{new} = \rho^\mathrm{old}`, and we
   update the enthalpy at the end of the integration as:

   .. math:: h^\mathrm{new} = h^\mathrm{old} + \frac{\Delta t}{2} \Hnuc.

-  Advect for :math:`\Delta t`:

   Beginning with the results from the previous reaction half time step, we integrate
   that state using the equations

   .. math::

      \begin{aligned}
      \frac{\partial(\rho X_k)}{\partial t} &=&
          -\nabla\cdot(\rho X_k\Ub), \\
      \frac{\partial(\rho h)}{\partial t} &=&
          -\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt}.\end{aligned}

   Note that no reaction terms appear here. Since the advection
   takes place using the state updated from the reaction step, the effect
   of the reactions is implicitly contained in the advective update.

-  React for :math:`\Delta t/2`:

   Finally, we react again, starting with the state left by the advection
   step.

Note that MAESTRO uses a predictor-corrector approach. After integrating :math:`(\rho X_k,\rho h)` over
the time step, we use this time-advanced state to get a better estimate of a time-centered :math:`\beta_0`
and :math:`S`. We re-compute the advective velocities using an updated divergence constraint and repeat
the thermodynamic variable advance.

SDC Without Thermal Diffusion or Base State Evolution
=====================================================

In the SDC version, the VODE integration at the end of an SDC
iteration is responsible for updating all the thermodynamic quantities
including both the advection (incorporated via a piecewise constant advective
flux divergence source term) and the reactions. This provides a much stronger coupling between
the physical processes. In particular, our system now looks like:

.. math::

   \begin{aligned}
   \frac{d(\rho X_k)}{dt} &=& \rho \omegadot_k(\rho,X_k,T) - \underbrace{\nabla\cdot(\rho X_k\Ub)}_{A_{\rho X_k}}\label{eq:sdc:rhoX} \\
   \frac{d(\rho h)}{dt}   &=& \rho \Hnuc - \underbrace{\nabla\cdot(\rho h\Ub) + \frac{Dp_0}{Dt}}_{A_{\rho h}} \label{eq:sdc:rhoh}\end{aligned}

Here, :math:`A_{\rho X_k}` and :math:`A_{\rho h}` are piecewise-constant (in time)
approximations to the change in :math:`{\rho X_k}` and :math:`{\rho h}` (respectively)
due to the advection. These are constructed by calling density_advance
and enthalpy_advance in MAESTRO and passed into the network solver
during the reaction step. A flowchart of the MAESTRO SDC algorithm is
shown in figure \ `[fig:sdc:flowchart] <#fig:sdc:flowchart>`__.

.. raw:: latex

   \centering

.. figure:: \sdcfigpath/flowchart_SDC
   :alt: [fig:sdc:flowchart] A flowchart of the MAESTRO SDC algorithm. The
   thermodynamic state variables and local velocity are
   indicated in each step. The base state is not shown as it is time-independent.
   Red text indicates that quantity was
   updated during that step. The predictor is
   outlined by the dotted box. The blue text indicates state
   variables that are the same in **Step 3** as they are in
   **Step 1**, i.e., they are unchanged by the predictor steps.
   The SDC loop is shown in the gray dotted box.

   [fig:sdc:flowchart] A flowchart of the MAESTRO SDC algorithm. The
   thermodynamic state variables and local velocity are
   indicated in each step. The base state is not shown as it is time-independent.
   Red text indicates that quantity was
   updated during that step. The predictor is
   outlined by the dotted box. The blue text indicates state
   variables that are the same in **Step 3** as they are in
   **Step 1**, i.e., they are unchanged by the predictor steps.
   The SDC loop is shown in the gray dotted box.

Advective Update
----------------

In the advective update, our goal is to compute :math:`A_{\rho X_k}` and
:math:`A_{\rho h}`. These terms approximate the following:

.. math::

   \begin{aligned}
   A_{\rho X_k} &=&  \left [- \nabla \cdot (\rho X_k \Ub) \right ]^{n+1/2} \\
   A_{\rho h}   &=&  \left [- \nabla \cdot (\rho h \Ub) + \frac{Dp_0}{Dt} \right ]^{n+1/2}\end{aligned}

The construction of the interface states used in the advective terms
uses either a time-lagged or iteratively-lagged approximation to the reaction
terms (:math:`I_{\rho X_k}` and :math:`I_{\rho h}`, see below) as a source term in the interface
prediction. This explicitly couples the reaction process to the
advection process.

Final Update
------------

The RHS routine that the ODE solver operates on will first construct
the density as:

.. math:: \rho = \sum_k (\rho X_k)

It will then derive the temperature from the equation of state. If we
are running with use_tfromp = T, then we do

.. math:: T = T(\rho, p_0, X_k)

otherwise, we do

.. math:: T = T(\rho, h, X_k)

Note that in contrast to the Strang splitting version, here we call the EOS
every time we enter the RHS routine, but here we call the EOS to compute temperature
rather than thermodynamic coefficients.

Finally we integrate the ODE system (Eqs. `[eq:sdc:rhoX] <#eq:sdc:rhoX>`__ and `[eq:sdc:rhoh] <#eq:sdc:rhoh>`__).
At the end of the integration, we define :math:`I_{\rho X_k}` and :math:`I_{\rho h}`. The actual
form of these depends on what quantities we predict to edges during
the construction of the advective fluxes.
Note that we only need :math:`I_{\rho X_k}` and :math:`I_{\rho h}` for the
prediction of the interface states, and not the VODE integration.
This is because all we need from the advection solver is the
approximation to :math:`A_{\rho X_k}` and :math:`A_{\rho h}` and not the final
updated state.

Species Source Terms.
---------------------

For the species prediction, the form of :math:`I` depends on
species_pred_type (see §\ `[sec:pred:density] <#sec:pred:density>`__).
We note that there is no :math:`I` term for :math:`\rho` or :math:`\rho'` prediction, since
the density evolution equation does not have a reaction source term.

-  species_pred_type = 1 (predict_rhoprime_and_X)
   or 3 (predict_rho_and_X)

   .. math::

      I_{X_k} = \frac{1}{\rho^{n+\myhalf}} \left [
            \frac{(\rho X_k)^\mathrm{new} -
                  (\rho X_k)^\mathrm{old}}{\Delta t} - A_{\rho X_k}  \right ].

   (Andy’s Idea) Define :math:`I_{X_k}` using

   .. math:: I_{X_k} = \frac{X_k^\mathrm{new} - X_k^\mathrm{old}}{\Delta t} - A_{X_k},

   where we first define a state that has only been updated with advection:

   .. math:: \frac{(\rho X_k)^{(1)} - (\rho X_k)^\mathrm{old}}{\Delta t} = A_{\rho X_k},

   and then define the species mass fractions,

   .. math::

      X_k^{(1)} = (\rho X_k)^{(1)} / \sum_k (\rho X_k)^{(1)}, \quad
      X_k^\mathrm{old} = (\rho X_k)^\mathrm{old} / \sum_k (\rho X_k)^\mathrm{old}, \quad
      X_k^\mathrm{new} = (\rho X_k)^\mathrm{new} / \sum_k (\rho X_k)^\mathrm{new},

   and finally define :math:`A_{X_k}` using

   .. math:: \frac{X^{(1)} - X^\mathrm{old}}{\Delta t}= A_{X_k}.

-  species_pred_type = 2 (predict_rhoX)

   .. math::

      I_{\rho X_k} = \frac{(\rho X_k)^\mathrm{new} - (\rho X_k)^\mathrm{old}}{\Delta t} - A_{\rho X_k}.
      \label{eq:sdc:Irhoo}

Enthalpy Source Terms.
----------------------

The appropriate constructions are:

-  enthalpy_pred_type = 0 (predict_rhoh)

   .. math:: I_{\rho h} = \frac{(\rho h)^{\rm new} - (\rho h)^{\rm old}}{\Delta t} - A_{\rho h}.

-  enthalpy_pred_type = 1 (predict_rhohprime, not implemented yet)

   (Andy’s Idea) Here we need an :math:`I_{\rho h}` term for the :math:`(\rho h)'` evolution
   equation (see Eq. \ `[rhohprime equation] <#rhohprime equation>`__). In this case we will use
   :math:`I_{(\rho h)'} = I_{\rho h}`. Since we are not evolving the base state, the PDE
   for :math:`(\rho h)_0` is simply :math:`\partial(\rho h)_0/\partial t = 0`, and thus the
   evolution equation for :math:`(\rho h)'` is the same as the evolution equation
   for :math:`\rho h`.

   In the future, when we enable base state evolution, the base state enthalpy
   evolution equation may need to know about the :math:`I_{\rho h}` source term.
   In particular, should :math:`(\rho h)_0` see a :math:`\overline{(\rho \Hnuc)}` term?
   what about an average thermal diffusion?

-  enthalpy_pred_type = 2 (predict_h )

   This is the most straightforward prediction type. The SDC solver
   integrates the equation for :math:`(\rho h)`:

   .. math:: \frac{\partial(\rho h)}{\partial t} = -\nabla\cdot(\rho h \Ub) + \frac{Dp_0}{Dt}  + \rho H_{\rm nuc}

   (shown here without diffusion or external heat sources). Expanding
   the time derivative and divergence, and using the continuity equation
   we see:

   .. math:: \frac{\partial h}{\partial t} = -\Ub \cdot \nabla h + \frac{1}{\rho} \frac{Dp_0}{Dt}  + \frac{1}{\rho} (\rho H_{\rm nuc}) \label{eq:sdc:h}

   Comparing these equations, we see that

   .. math::

      I_{h}  = \frac{1}{\rho^{n+\myhalf}} \left [
          \frac{(\rho h)^\mathrm{new} - (\rho h)^\mathrm{old}}{\Delta t} - A_{\rho h} \right ]

   (Andy’s Idea) Form :math:`I_h` in the same way we would form :math:`I_{X_k}` from above:

   .. math:: I_h = \frac{h^\mathrm{new} - h^\mathrm{old}}{\Delta t} - A_h,

   where we first define

   .. math:: \frac{(\rho h)^{(1)} - (\rho h)^\mathrm{old}}{\Delta t} = A_{\rho h},

   and then define :math:`h`,

   .. math:: h^{(1)} = (\rho h)^{(1)} / \sum_k(\rho X_k)^{(1)}, \quad h^\mathrm{old} = (\rho h)^\mathrm{old} / \sum_k(\rho X_k)^\mathrm{old}, \quad h^\mathrm{new} = (\rho h)^\mathrm{new} / \sum_k(\rho X_k)^\mathrm{new},

   and finally define :math:`A_h` using

   .. math:: I_h = \frac{h^{(1)} - h^\mathrm{old}}{\Delta t} = A_h.

-  enthalpy_pred_type = 3 (predict_T_then_rhoprime) or
   enthalpy_pred_type = 4 (predict_T_then_h )

   Both of these enthalpy_pred_types predict temperature. Expressing
   :math:`h = h(p_0,T,X_k)` and differentiating along particle paths:

   .. math::

      \begin{aligned}
      \frac{Dh}{Dt} &=& \left . \frac{\partial h}{\partial T} \right |_{p,X_k} \frac{DT}{Dt} +
                        \left . \frac{\partial h}{\partial p} \right |_{T,X_k} \frac{Dp_0}{Dt} +
                 \sum_k \left . \frac{\partial h}{\partial X_k} \right |_{p,T} \frac{DX_k}{Dt} \\
                    &=& c_p \frac{DT}{Dt} + h_p  \frac{Dp_0}{Dt} + \sum_k \xi_k \omegadot_k\end{aligned}

   where :math:`c_p`, :math:`h_p`, and :math:`\xi_k` are as defined in the table of symbols
   (Table `[table:sym] <#table:sym>`__), and we substitute :math:`DX_k/Dt = \omegadot_k` (from the species
   continuity equation, Eq. \ `[species equation] <#species equation>`__). Using Eq. \ `[eq:sdc:h] <#eq:sdc:h>`__, we have
   the familiar temperature evolution equation:

   .. math:: \rho c_p \frac{DT}{Dt} = \underbrace{(1 - \rho h_p) \frac{Dp_0}{Dt}}_{\begin{smallmatrix}\text{already~accounted~for} \\ \text{in~T~prediction}\end{smallmatrix}} - \sum_k \xi_k \rho \omegadot_k + \rho \Hnuc

   where the underbraced term is already present in mktempforce. Recognizing that
   Eq. \ `[eq:sdc:Irhoh] <#eq:sdc:Irhoh>`__ is the SDC approximation to :math:`(\rho \Hnuc)` and Eq. \ `[eq:sdc:Irhoo] <#eq:sdc:Irhoo>`__ is the
   SDC approximation to :math:`(\rho \omegadot_k)`, we can define

   .. math::

      I_T = \frac{1}{\rho^{n+\myhalf} c_p^{n+\myhalf}} \left \{
        \left [ \frac{(\rho h)^\mathrm{new} - (\rho h)^\mathrm{old}}{\Delta t} - A_{\rho h} \right ]
        - \sum_k \xi_k^{n+\myhalf} \left [      \frac{(\rho X_k)^\mathrm{new} -
                  (\rho X_k)^\mathrm{old}}{\Delta t} - A_{\rho X_k}  \right ] \right \}

   (Andy’s Idea) The idea is to advance the species and enthalpy with advection
   terms only, and compute the resulting temperature, :math:`T^{(1)}`. Compare that temperature
   with the final temperature computed by the SDC VODE call. The difference
   between these values is :math:`I_T`.

   .. math:: I_T = \frac{T^\mathrm{new} - T^\mathrm{old}}{\Delta t} - A_T,

   with :math:`A_T` given by

   .. math:: \frac{T^{(1)} - T^\mathrm{old}}{\Delta t} = A_T,

   and :math:`T^{(1)}` computed using the equation of state from :math:`\rho^{(1)}, X_k^{(1)}`,
   and :math:`h^{(1)}` (or :math:`p_0`, if use_tfromp = T).

Implementation
--------------

This is done in advance.f90 just after the call to react_state,
stored in the multifab called intra.
These terms are used as the source terms for the
advection step in the next SDC iteration.

Summary of Changes
------------------

The major changes from the non-SDC-enabled burners is the addition of
the advective terms to the system of ODEs, the fact that we integrate
:math:`(\rho X_k)` instead of just :math:`X_k`, integrate :math:`(\rho h)` instead
of :math:`T`, and the need to derive the
temperature from the input state for each RHS evaluation by VODE.

Note also that the SDC integration by VODE does not operate on
the velocities at all. That update is handled in the same fashion
as the Strang splitting version of the code.

The ignition_simple_SDC burner shows how to setup the system
for use_tfromp = T or F. Presently, this implementation
does not support evolve_base_state = T (in particular, we
need to evolve :math:`p_0` in the RHS routine).

Algorithm Flowchart - ADR with Fixed Base State
===============================================

We now include thermal diffusion and assume the base state is constant in time but not space:

.. math::

   \begin{aligned}
   \frac{\partial\Ub}{\partial t} =&
       -\Ub\cdot\nabla\Ub  - \frac{1}{\rho}\nabla\pi
       - \frac{\rho-\rho_0}{\rho} g\eb_r,\\
   \frac{\partial(\rho X_k)}{\partial t} =&
       -\nabla\cdot(\rho X_k\Ub) + \rho\omegadot_k,\label{eq:sdc species}\\
   \frac{\partial(\rho h)}{\partial t} =&
       -\nabla\cdot(\rho h\Ub) + \underbrace{\Ub\cdot\nabla p_0}_{Dp_0/Dt}
       + \rho\Hnuc
       + \underbrace{\nabla\cdot\frac{\kth}{c_p}\nabla h - \sum_k\nabla\cdot\frac{\xi_k k_{\rm th}}{c_p}\nabla X_k - \nabla\cdot\frac{h_p k_{\rm th}}{c_p}\nabla p_0}_{\nabla\cdot k_{\rm th}\nabla T}.\nonumber\\
   \label{eq:sdc enthalpy}\end{aligned}

| The time-advancement is divided into three major steps. The first step is the predictor, where we integrate the thermodynamic variables, :math:`(\rho,\rho X_k,\rho h)`, over the full time step. The second step is corrector, where we use the results from the predictor to perform a more accurate temporal integration of the thermodynamic variables. The third step is the velocity and dynamic pressure update.
| **Step 1:** (*Compute advection velocities*)
| Use :math:`\Ub^n` and a second-order Godunov method to compute time-centered edge velocities, :math:`\uadvsdcstar`, with time-lagged dynamic pressure and explicit buoyancy as forcing terms. The :math:`\star` superscript indicates that this field does not satisfy the divergence constraint. Compute :math:`S^{n+\myhalf,{\rm pred}}` by extrapolating in time,

  .. math:: S^{n+\myhalf,{\rm pred}} = S^n + \frac{\Delta t^n}{2}\frac{S^n - S^{n-1}}{\Delta t^{n-1}},

  and project :math:`\uadvsdcstar` to obtain :math:`\uadvsdcpred`, which satisfies

  .. math:: \nabla\cdot\left(\beta_0^n\uadvsdcpred\right) = S^{n+\myhalf,\rm{pred}}.

| **Step 2:** (*Predictor*)
| In this step, we integrate :math:`(\rho, \rho X_k, \rho h)` over the full time step. The quantities :math:`(S, \beta_0, k_{\rm th}, c_p, \xi_k, h_p)^n` are computed from the the thermodynamic variables at :math:`t^n`. This step is divided into several sub-steps:
| **Step 2A:** (*Compute advective flux divergences*)
| Use :math:`\uadvsdcpred` and a second-order Godunov integrator to compute time-centered edge states, :math:`(\rho X_k, \rho h)^{n+\myhalf,(0)}`, with time-lagged reactions (:math:`I^{\rm lagged} = I^{(j_{\rm max})}` from the previous time step), explicit diffusion, and time-centered thermodynamic pressure as source terms. Define the advective flux divergences as

  .. math::

     \begin{aligned}
     A_{\rho X_k}^{(0)} &=& -\nabla\cdot\left[\left(\rho X_k\right)^{n+\myhalf,{(0)}}\uadvsdcpred\right],\\
     A_{\rho h}^{(0)} &=& -\nabla\cdot\left[\left(\rho h\right)^{n+\myhalf,(0)}\uadvsdcpred\right] + \uadvsdcpred\cdot\nabla p_0.\end{aligned}

  Next, use these fluxes to compute the time-advanced density,

  .. math:: \frac{\rho^{n+1} - \rho^n}{\Delta t} = \sum_k A_{\rho X_k}^{(0)}.

  Then, compute preliminary, time-advanced species using

  .. math:: \frac{\rho^{n+1}\widehat{X}_k^{n+1,(0)} - (\rho X_k)^n}{\Delta t} = A_{\rho X_k}^{(0)} + I_{\rho X_k}^{\rm lagged}.\label{eq:sdc species 2}

| **Step 2B:** (*Compute diffusive flux divergence*)
| Solve a Crank-Nicolson-type diffusion equation for :math:`\widehat{h}^{n+1,(0)}`, using transport coefficients evaluated at :math:`t^n` everywhere,

  .. math::

     \begin{aligned}
     \frac{\rho^{n+1}\widehat{h}^{n+1,(0)} - (\rho h)^n}{\Delta t} =& A_{\rho h}^{(0)} + I_{\rho h}^{\rm lagged}\nonumber\\
     & + \half\left(\nabla\cdot\frac{\kth^n}{c_p^n}\nabla h^n + \nabla\cdot\frac{\kth^n}{c_p^n}\nabla \widehat{h}^{n+1,(0)}\right)\nonumber\\
     & - \half\left(\sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla X_k^n + \sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla\widehat{X}_k^{n+1,(0)}\right)\nonumber\\
     & - \half\left(\nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0 + \nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0\right),\label{eq:sdc enthalpy 2}\end{aligned}

  which is equivalent to

  .. math::

     \begin{aligned}
     \left(\rho^{n+1} - \frac{\Delta t}{2}\nabla\cdot\frac{k_{\rm th}^n}{c_p^n}\nabla\right)\widehat{h}^{n+1,(0)} =& (\rho h)^n + \Delta t\Bigg[A_{\rho h}^{(0)} + I_{\rho h}^{\rm lagged} + \left(\half\nabla\cdot\frac{\kth^n}{c_p^n}\nabla h^n\right)\nonumber\\
     &\hspace{0.85in} - \half\left(\sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla X_k^n + \sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla\widehat{X}_k^{n+1,(0)}\right)\nonumber\\
     &\hspace{0.85in} - \half\left(\nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0 + \nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0\right)\Bigg].\end{aligned}

| **Step 2C:** (*Advance thermodynamic variables*)
| Define :math:`Q_{\rho X_k}^{(0)}` as the right hand side of (`[eq:sdc species 2] <#eq:sdc species 2>`__) without the :math:`I_{\rho X_k}^{\rm lagged}` term, and define :math:`Q_{\rho h}^{(0)}` as the right hand side of (`[eq:sdc enthalpy 2] <#eq:sdc enthalpy 2>`__) without the :math:`I_{\rho h}^{\rm lagged}` term. Use VODE to integerate (`[eq:sdc species] <#eq:sdc species>`__) and (`[eq:sdc enthalpy] <#eq:sdc enthalpy>`__) over :math:`\Delta t` to advance :math:`(\rho X_k, \rho h)^n` to :math:`(\rho X_k, \rho h)^{n+1,(0)}` using the piecewise-constant advection and diffusion source terms:

  .. math::

     \begin{aligned}
     \frac{\partial(\rho X_k)}{\partial t} &=& Q_{\rho X_k}^{(0)} + \rho\dot\omega_k\\
     \frac{\partial(\rho h)}{\partial t} &=& Q_{\rho h}^{(0)} + \rho\Hnuc.\end{aligned}

  At this point we can define :math:`I_{\rho X_k}^{(0)}` and :math:`I_{\rho h}^{(0)}`, or whatever term we need depending on our species and enthalpy edge state prediction types, for use in the corrector step. In our first implementation, we are predicting :math:`\rho X_k` and :math:`\rho h`, in which case we define:

  .. math::

     \begin{aligned}
     I_{\rho X_k}^{(0)} &=& \frac{(\rho X_k)^{n+1,(0)} - (\rho X_k)^n}{\Delta t} - Q_{\rho X_k}^{(0)}\\
     I_{\rho h}^{(0)} &=& \frac{(\rho h)^{n+1,(0)} - (\rho h)^n}{\Delta t} - Q_{\rho h}^{(0)}.\end{aligned}

| **Step 3:** (*Update advection velocities*)
| First, compute :math:`S^{n+\myhalf}` and :math:`\beta_0^{n+\myhalf}` using

  .. math:: S^{n+\myhalf} = \frac{S^n + S^{n+1,(0)}}{2}, \qquad \beta_0^{n+\myhalf} = \frac{\beta_0^n + \beta_0^{n+1,(0)}}{2}.

  Then, project :math:`\uadvsdcstar` to obtain :math:`\uadvsdc`, which satisfies

  .. math:: \nabla\cdot\left(\beta_0^{n+\myhalf}\uadvsdc\right) = S^{n+\myhalf}.

| **Step 4:** (*Corrector Loop*)
| We loop over this step from :math:`j=1,j_{\rm max}` times. In the corrector, we use the time-advanced state from the predictor to perform a more accurate integration of the thermodynamic variables. The quantities :math:`(S, \beta_0, k_{\rm th}, c_p, \xi_k, h_p)^{n+1,(j-1)}` are computed from :math:`(\rho,\rho X_k,\rho h)^{n+1,(j-1)}`. This step is divided into several sub-steps:
| **Step 4A:** (*Compute advective flux divergences*)
| Use :math:`\uadvsdc` and a second-order Godunov integrator to compute time-centered edge states, :math:`(\rho X_k, \rho h)^{n+\myhalf}`, with iteratively-lagged reactions (:math:`I^{(j-1)}`), explicit diffusion, and time-centered thermodynamic pressure as source terms. Define the advective flux divergences as

  .. math::

     \begin{aligned}
     A_{\rho X_k}^{(j)} &=& -\nabla\cdot\left[\left(\rho X_k\right)^{n+\myhalf,(j)}\uadvsdc\right],\\
     A_{\rho h}^{(j)} &=& -\nabla\cdot\left[\left(\rho h\right)^{n+\myhalf,(j)}\uadvsdc\right] + \uadvsdc\cdot\nabla p_0.\end{aligned}

  Then, compute preliminary, time-advanced species using

  .. math:: \frac{\rho^{n+1}\widehat{X}_k^{n+1,(j)} - (\rho X_k)^n}{\Delta t} = A_{\rho X_k}^{(j)} + I_{\rho X_k}^{(j-1)}.\label{eq:sdc species 3}

| **Step 4B:** (*Compute diffusive flux divergence*)
| Solve a backward-Euler-type correction equation for :math:`\widehat{h}^{n+1,(j)}`,

  .. math::

     \begin{aligned}
     \frac{\rho^{n+1}\widehat{h}^{n+1,(j)} - (\rho h)^n}{\Delta t} =& A_{\rho h}^{(j)} + I_{\rho h}^{(j-1)}\nonumber\\
     & + \nabla\cdot\frac{\kth^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla\widehat{h}^{n+1,(j)} + \half\left(\nabla\cdot\frac{\kth^n}{c_p^n}\nabla h^n - \nabla\cdot\frac{\kth^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla h^{n+1,(j-1)}\right)\nonumber\\
     & - \half\left(\sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla X_k^n + \sum_k\nabla\cdot\frac{\xi_k^{n+1,(j-1)} k_{\rm th}^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla\widehat{X}_k^{n+1,(j)}\right)\nonumber\\
     & - \half\left(\nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0 + \nabla\cdot\frac{h_p^{n+1,(j-1)}k_{\rm th}^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla p_0\right),\label{eq:sdc enthalpy 3}\end{aligned}

  which is equivalent to

  .. math::

     \begin{aligned}
     \left(\rho^{n+1} - \Delta t\nabla\cdot\frac{\kth^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla\right)\widehat{h}^{n+1,(j)} =& (\rho h)^n + \Delta t\Bigg[A_{\rho h}^{(j)} + I_{\rho h}^{(j-1)} \nonumber\\
     & + \half\left(\nabla\cdot\frac{\kth^n}{c_p^n}\nabla h^n - \nabla\cdot\frac{\kth^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla h^{n+1,(j-1)}\right)\nonumber\\
     & - \half\left(\sum_k\nabla\cdot\frac{\xi_k^n k_{\rm th}^n}{c_p^n}\nabla X_k^n + \sum_k\nabla\cdot\frac{\xi_k^{n+1,(j-1)} k_{\rm th}^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla\widehat{X}_k^{n+1,(j)}\right)\nonumber\\
     & - \half\left(\nabla\cdot\frac{h_p^n k_{\rm th}^n}{c_p^n}\nabla p_0 + \nabla\cdot\frac{h_p^{n+1,(j-1)}k_{\rm th}^{n+1,(j-1)}}{c_p^{n+1,(j-1)}}\nabla p_0\right)\Bigg].\end{aligned}

| **Step 4C:** (*Advance thermodynamic variables*)
| Define :math:`Q_{\rho X_k}^{(j)}` as the right hand side of (`[eq:sdc species 3] <#eq:sdc species 3>`__) without the :math:`I_{\rho X_k}^{(j-1)}` term, and define :math:`Q_{\rho h}^{(j)}` as the right hand side of (`[eq:sdc enthalpy 3] <#eq:sdc enthalpy 3>`__) without the :math:`I_{\rho h}^{(j-1)}` term. Use VODE to integerate (`[eq:sdc species] <#eq:sdc species>`__) and (`[eq:sdc enthalpy] <#eq:sdc enthalpy>`__) over :math:`\Delta t` to advance :math:`(\rho X_k, \rho h)^n` to :math:`(\rho X_k, \rho h)^{n+1,(j)}` using the piecewise-constant advection and diffusion source terms:

  .. math::

     \begin{aligned}
     \frac{\partial(\rho X_k)}{\partial t} &=& Q_{\rho X_k}^{(j)} + \rho\dot\omega_k\\
     \frac{\partial(\rho h)}{\partial t} &=& Q_{\rho h}^{(j)} + \rho\Hnuc.\end{aligned}

  At this point we can define :math:`I_{\rho X_k}^{(j)}`, :math:`I_{\rho h}^{(j)}`, and any other :math:`I` terms we need depending on
  our species and enthalpy edge state prediction types, for use in the predictor in the next time step. In our first implementation, we are predicting :math:`\rho X_k` and :math:`\rho h`, in which case we define:

  .. math::

     \begin{aligned}
     I_{\rho X_k}^{(j)} &=& \frac{(\rho X_k)^{n+1,(j)} - (\rho X_k)^n}{\Delta t} - Q_{\rho X_k}^{(j)}\\
     I_{\rho h}^{(j)} &=& \frac{(\rho h)^{n+1,(j)} - (\rho h)^n}{\Delta t} - Q_{\rho h}^{(j)}.\end{aligned}

| **Step 5:** (*Advance velocity and dynamic pressure*)
| Similar to the original MAESTRO algorithm, more to come.
