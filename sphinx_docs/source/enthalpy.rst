*****************
Notes on Enthalpy
*****************

Evolution Equations
===================

The compressible and low Mach number formulations of the governing
equations both share the unapproximated continuity and momentum equations
shown here (where :math:`p = p_0 + \pi`).

.. math::

   \begin{aligned}
   \frac{\partial(\rho X_k)}{\partial t} &=& -\nabla\cdot(\rho X_k\Ub) +
   \rho\omegadot_k,\label{enth:eq:species}\\
   \frac{\partial\Ub}{\partial t} &=& -\Ub\cdot\nabla\Ub  -
     \frac{1}{\rho}\nabla\pi -
     \frac{\rho-\rho_0}{\rho} g\eb_r,\label{eq:momentum}\end{aligned}

In the compressible formulation we complete the system with an energy
equation as well as an equation of state; in the low Mach number
formulation we can derive a constraint on the velocity by setting
:math:`p_{EOS} = p(\rho,h,X_k) = p_0` and differentiating the equation of
state along particle paths. In this case adding the energy equation
would over-constrain the system, but its solution can help us in
providing a diagnostic capability for the solution. We can write the
energy equation in terms of internal energy, :math:`e`, or enthalpy, :math:`h`;
these two equations are analytically equivalent.

.. math::

   \begin{aligned}
   \frac{\partial(\rho h)}{\partial t} + \nabla\cdot(\rho h\Ub)
    &=& \frac{Dp}{Dt} + \rho\Hnuc \nonumber \\
   %
   \frac{\partial(\rho e)}{\partial t} + \nabla\cdot(\rho e \Ub)
    &=& - p\nabla\cdot \Ub + \rho\Hnuc \nonumber \end{aligned}

In low Mach number combustion, :math:`Dp/Dt = 0,` which makes the enthalpy
equation preferable to the internal energy equation because we don’t
need to evaluate :math:`p \nabla \cdot \Ub.`

Derivation of Velocity Constraint
=================================

Differentiating the equation of state, written in the form, :math:`p =
p(\rho,T,X_k),` along particle paths, we can write

.. math:: \frac{D p}{Dt}  = p_\rho \frac{D \rho}{Dt} + p_T \frac{D T}{Dt} + \sum_k p_{X_k} \frac{D X_k}{Dt}

Then, by rearranging the terms, we get

.. math::

   \frac{D \rho}{Dt}  = -\rho \nabla \cdot \Ub =
       \frac{1}{p_\rho}
       \left( \frac{D p}{Dt} - p_T \frac{D T}{Dt}
                             - \sum_k p_{X_k} {\omegadot}_k \right)  ,

with :math:`p_\rho = \left.\partial p/\partial \rho\right|_{X_k,T}`,
:math:`p_{X_k} = \left.\partial p/\partial X_k \right|_{T,\rho,(X_j,j\ne k)}`,
and :math:`p_T = \left.\partial p/\partial T\right|_{\rho,X_k}`.

Using the Enthalpy Equation
---------------------------

Now writing :math:`h = h(p,T,X_k)` and expanding :math:`Dh/Dt` and using the
enthalpy evolution equation:

.. math::

   \rho \frac{D h}{Dt}  = \rho \left( h_p \frac{D p}{Dt} + c_p \frac{D T}{Dt} + \sum_k h_{X_k} \frac{D X_k}{Dt} \right)
                        = \frac{D p}{D t} + \rho \Hnuc

we can express :math:`DT/Dt` in terms of

.. math::

   \frac{DT}{Dt} = \frac{1}{\rho c_p} \left( (1 - \rho h_p) \frac{D p}{D t}
   - \sum_k \rho \xi_k \omegadot_k + \rho \Hnuc \right)  , \label{eq:dTdt}

where :math:`c_p = \left.\partial h/\partial T\right|_{p,X_k}` is the
specific heat at constant pressure,
:math:`\xi_k = \left.\partial h/\partial X_k \right|_{p,T}`,
and :math:`h_p = \left.\partial h/\partial p\right|_{T,X_k}.`

We could then write,

.. math::

   \begin{aligned}
   \nabla \cdot \Ub &=& \frac{1}{\rho p_\rho} \left(
   - \frac{D p}{D t} + \frac{p_T}{\rho c_p}
     \left( (1 - \rho h_p) \frac{D p}{D t} - \rho \sum_k \xi_k \omegadot_k + \rho \Hnuc \right)
   + \sum_k p_{X_k} \omegadot_k \right)  \\
                    &=& \frac{1}{\rho p_\rho}
     \left( \frac{p_T}{\rho c_p}(1  - \rho h_p) - 1 \right) \frac{D p}{D t}
    + \frac{1}{\rho p_\rho} \left(
     \frac{p_T}{\rho c_p} (\rho \Hnuc - \rho \sum_k \xi_k   \omegadot_k)
                                  + \sum_k p_{X_k} \omegadot_k \right)  .\end{aligned}

When we derived this expression we explicitly retained the dependence
of :math:`h` on :math:`p`, as shown by the presence of the :math:`h_p` term.

Then, replacing :math:`p` by :math:`p_0(r)`, :math:`Dp/Dt` becomes :math:`\Ub \cdot
\nabla p_0`, and the divergence constraint can be written

.. math::

   \nabla \cdot \Ub + \alpha \Ub \cdot \nabla p_0 =
   \frac{1}{\rho p_\rho} \left(
      \frac{p_T}{\rho c_p} \left(
     - \sum_k\rho  \xi_k \omegadot_k + \rho \Hnuc \right)
    + \sum_k p_{X_k} \omegadot_k \right)  \equiv \tilde{S}  \label{eq:full_divu_constraint} ,

where we define

.. math::

   \alpha(\rho,T) \equiv - \left( \frac{(1 - \rho h_p )p_T - \rho c_p}{\rho^2
     c_p p_\rho} \right)  . \label{eq:alphadef}

Using the Energy Equation
-------------------------

Now writing :math:`e = e(p,T,X_k)` and expanding :math:`De/Dt` and using the
energy evolution equation:

.. math::

   \rho \frac{D e}{Dt}  = \rho \left( e_p \frac{D p}{Dt} + e_T \frac{D T}{Dt} + \sum_k e_{X_k} \frac{D X_k}{Dt} \right)
                        = -p \nabla \cdot \Ub + \rho \Hnuc

we can express :math:`DT/Dt` in terms of

.. math::

   \begin{aligned}
   \frac{DT}{Dt} &=& \frac{1}{\rho e_T} \left(
                     -p \nabla \cdot \Ub + \rho \Hnuc
                    - \rho e_p \frac{D p}{Dt}
                    - \rho \sum_k e_{X_k} \omegadot_k   \right)\end{aligned}

where :math:`e_T = \left.\partial e/\partial T\right|_{p,X_k}`
and :math:`e_p = \left.\partial e/\partial p\right|_{T,X_k}.`
Then

.. math::

   -\rho \nabla \cdot \Ub =
       \frac{1}{p_\rho} \left( \frac{D p}{Dt}
       - \frac{p_T}{\rho e_T} \left(
                    -p \nabla \cdot \Ub + \rho \Hnuc
                    - \rho e_p \frac{D p}{Dt}
                    - \rho \sum_k e_{X_k} \omegadot_k   \right)
       - \sum_k p_{X_k} {\omegadot}_k \right)  ,

which leads to

.. math::

   \left(-\rho -  \frac{p p_T}{\rho e_T p_\rho} \right) \nabla \cdot \Ub =
       \frac{1}{p_\rho} \left( (1 + \frac{e_p p_T}{e_T}) \frac{D p}{Dt}
                              - \frac{p_T}{\rho e_T} \left(
                                           \rho \Hnuc
                                           - \rho \sum_k e_{X_k} \omegadot_k   \right)
       - \sum_k p_{X_k} {\omegadot}_k \right)  ,

Note that we can replace :math:`p` by :math:`p_0` in the coefficient on the
l.h.s. as well as on the r.h.s.

Comparison of Constraints
=========================

If we set :math:`\omegadot_k = \Hnuc = 0` for simplicity, then the
constraint as derived using :math:`h` can be written

.. math::

   \nabla \cdot \Ub
   + \left( \frac{(1 - \rho h_p )p_T - \rho c_p}{\rho^2
     c_p p_\rho} \right) \frac{D p_0}{D t} = 0

and the constraint derived using :math:`e` can be written

.. math::

   \nabla \cdot \Ub
   + \left( \frac{\rho e_T + \rho e_p p_T}{\rho^2 e_T p_\rho + p p_T} \right) \frac{D p_0}{D t} = 0

We note that if we evaluate both constraints for :math:`p = \rho R T,` with
:math:`h = c_p T,` :math:`e = c_v T,` :math:`c_p = c_v + R` and :math:`\gamma = c_p / c_v,`
then both constraints reduce to

.. math:: \nabla \cdot \Ub + \frac{1}{\gamma p} \frac{D p_0}{D t} = 0

Enthalpy vs Energy Equation
===========================

The full enthalpy equation, with no approximations, appears as:

.. math::

   \frac{\partial(\rho h)}{\partial t} = -\nabla\cdot(\rho h\Ub) +
     \frac{Dp}{Dt} + \rho\Hnuc \label{eq:enthalpy}

Here, :math:`h = e + p/\rho` is the specific enthalpy, with :math:`e` the specific
internal energy. In the low Mach number formulation, we replace :math:`p`
with :math:`p_0` in the :math:`Dp/Dt` term, however, the definition of enthalpy
implicitly contains a pressure. When calling the equation of state,
we take :math:`h` and :math:`\rho` as inputs. The equation of state is expressed
in terms of :math:`T` and :math:`\rho`, so it iterates until it finds the :math:`h` that
we desire. This :math:`h` will be of the form :math:`h = e + p_\mathrm{EOS}/\rho`,
where :math:`p_\mathrm{EOS}` is the pressure returned from the EOS. Note that
:math:`p_\mathrm{EOS}` may not be equal to :math:`p_0`—this may be what
causes us to drive off of the constraint.

The mismatch between the pressure implicit in the definition of :math:`h`
and :math:`p_0` can be seen by substituting :math:`h = e + p/\rho` into the
enthalpy equation, where we replace :math:`p` with :math:`p_0` in the :math:`Dp/Dt` term:

.. math::

   \begin{aligned}
   \frac{\partial(\rho h)}{\partial t} &=& -\nabla\cdot(\rho h\Ub) +
     \frac{Dp_0}{Dt} + \rho\Hnuc \nonumber \\
   %
   \frac{\partial(\rho e)}{\partial t} + \frac{\partial p}{\partial t} &=&
    -\nabla\cdot(\rho e\Ub) -\nabla\cdot(p\Ub) + \frac{Dp_0}{Dt} + \rho\Hnuc \nonumber \\
   %
   \frac{\partial(\rho e)}{\partial t} &=&
    -\nabla\cdot(\rho e\Ub) - p\nabla\cdot\Ub + \rho\Hnuc +
     \left \{ \frac{Dp_0}{Dt} - \frac{Dp}{Dt} \right \} \nonumber \end{aligned}

However, if we solve the evolution equation for :math:`e` we would
substitute :math:`p_0` for :math:`p` in this equation as well. Thus, we can pose
the situation as the following. If we solve the evolution equation
for :math:`h` then we effectively are solving

.. math::

   \frac{\partial(\rho e)}{\partial t} +
     \nabla\cdot(\rho e\Ub) = -p \; \nabla\cdot\Ub + \rho\Hnuc +
     \left \{ \frac{Dp_0}{Dt} - \frac{Dp}{Dt} \right \} \nonumber

but if we solve the evolution equation for :math:`e` we are effectively solving

.. math::

   \frac{\partial(\rho e)}{\partial t} +
     \nabla\cdot(\rho e\Ub) = p_0 \nabla\cdot\Ub + \rho\Hnuc \nonumber

The second equation subtracted from the first gives:

.. math::
   \frac{D (p_0 - p)}{Dt} - (p_0 - p)  \nabla \cdot \Ub = 0,
   :label: eq:difference between h and e equations


but this equation is only true, in general, if :math:`p=p_0`.

Suppose we solve the current enthalpy equation, but when we call the
EOS, we subtract :math:`p_0` from :math:`\rho h` and then call the EOS with :math:`e`
instead of :math:`h`. This is equivalent to:

.. math::

   \begin{aligned}
   \frac{\partial(\rho h)}{\partial t} &=& -\nabla\cdot(\rho h\Ub) +
     \frac{Dp_0}{Dt} + \rho\Hnuc \nonumber \\
   %
   \frac{\partial(\rho e)}{\partial t} + \frac{\partial p_0}{\partial t} &=&
    -\nabla\cdot(\rho e\Ub) -\nabla\cdot(p_0 \Ub) + \frac{Dp_0}{Dt} + \rho\Hnuc \nonumber \\
   %
   \frac{\partial(\rho e)}{\partial t} &=&
    -\nabla\cdot(\rho e\Ub) - p_0 \nabla\cdot\Ub + \rho\Hnuc  \nonumber\end{aligned}

which is identical to solving the energy equation with :math:`p\to p_0`.
This option is enabled in MAESTROeX via
use_eos_e_instead_of_h = T.

Constant :math:`\gamma` Gas
---------------------------

Going back to the constant :math:`\gamma`, ideal gas EOS, we can rewrite the
enthalpy equation as a pressure evolution equation

.. math::
   \begin{aligned}
     \frac{\partial\rho h}{\partial t} + \nabla\cdot\left(\rho h\Ub\right) &=& \frac{Dp_0}{Dt} + \rho H {} \nonumber\\
     \frac{\gamma}{\gamma-1}\frac{\partial p}{\partial t} + \frac{\gamma}{\gamma-1}\nabla\cdot\left(p\Ub\right) &=& \frac{Dp_0}{Dt} + \rho H {} \nonumber\\
     \frac{Dp}{Dt} &=& -p\nabla\cdot\Ub + \frac{\gamma}{\gamma-1}\left(\frac{Dp_0}{Dt}+\rho H\right)
   \end{aligned}
   :label: eq:H:pressure evolution

Similarly, we can derive a pressure evolution equation from the energy equation

.. math::
   \frac{Dp}{Dt} = -p\nabla\cdot\Ub - \left(\gamma-1\right)p_0\nabla\cdot\Ub + \rho H
   :label: eq:e:pressure evolution

Now, if we further make the assumption that :math:`p_0` is constant,
:math:`Dp_0/Dt = 0`, then the divergence constraint for such a gas reduces
to

.. math::
   \nabla\cdot\Ub = \frac{\gamma-1}{\gamma p_0}\rho H .
   :label: eq:div constraint for constant gamma

Plugging this back into either of :eq:`eq:H:pressure evolution` or
:eq:`eq:e:pressure evolution` gives

.. math::
   \frac{Dp}{Dt} = \frac{\gamma-1}{\gamma}\left(1-\frac{p}{p_0}\right)\rho H.
   :label: eq:pressure evolution constant gamma

If :math:`p_0` is assumed constant and using :eq:`div constraint for
constant gamma`, the difference between the enthalpy equation and the
energy equation, :eq:`eq:difference between h and e equations`, can be
rewritten as

.. math:: -\frac{Dp}{Dt} - \frac{\gamma}{\gamma-1}\left(1-\frac{p}{p_0}\right)\rho H = 0,

where the equality holds from :eq:`eq:pressure evolution constant
gamma`. In other words, for the constant :math:`\gamma` gas we have
:math:`p=p_0` as expected.

Outstanding Questions
=====================

#. Why do we want to start with enthalpy instead of internal energy?

   We believe that the original desire stems from our experience with
   smallscale combustion. There, stratification is not important and
   :math:`Dp_0/Dt = 0`, so the enthalpy equation becomes a conservation
   equation for :math:`(\rho h)`.

#. Should we call the EOS with :math:`h` as is, or call the EOS with :math:`e =
   h - p_0 / \rho`?

#. When we stay on the constraint, i.e. :math:`p_{EOS} = p_0`, then the
   equations for :math:`e` and for :math:`h` are identical. However, once we are
   off the constraint, do the terms in the current evolution equation
   for :math:`h` serve to drive us back to the constraint? Recall our
   current "volume discrepancy factor" acts as a source term which
   modifies the divergence constraint, which effectively modifies both
   :math:`\rho` and :math:`T` (or :math:`e` or :math:`h`). The term in the enthalpy equation
   only modifies :math:`\rho.` Is this relevant and/or useful?? Recall that
   the current "volume discrepancy factor" takes the form of adding to
   the r.h.s. of the constraint:

   .. math::

      \nabla \cdot(\beta_0\Ub) = \beta_0\left(S-\frac{1}{\overline{\Gamma_1} p_0}
             \frac{\partial p_0}{\partial t} - \frac{f}{\overline{\Gamma_1} p_0}
             \frac{p_0-p_\text{EOS}}{\dt}\right)

#. Suppose we corrected the :math:`h` (or :math:`e` equation) by using the full
   :math:`p` instead of :math:`p_0`? Would this be more or less consistent (one
   could imagine doing this as a correction after solving for :math:`\pi`
   earlier in the timestep).

#. In computing the thermodynamic coefficients in :math:`S` for the
   projection, don’t we need these to be in terms of :math:`p_0` instead of
   :math:`p(h,\rho)`?
