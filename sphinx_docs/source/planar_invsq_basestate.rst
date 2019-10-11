.. _sec:planarinvsqgravity:

**********************************************************
Modifications for a :math:`1/r^2` Plane-Parallel Basestate
**********************************************************

In the plane parallel assumption, we model a layer of thickness
:math:`\Delta R`` a distance :math:`R_\mathrm{base}`` from the center of a star,
under the assumption the :math:`\Delta R \ll R_\mathrm{base}``. In this
assumption, we can neglect the curvature of the atmosphere. Here, we
extend that basic assumption to allow for the gravitational acceleration
to still fall off as :math:`1/r^2`` as we move outward in the envelope. We assume
that the mass of the envelope is insignificant, and that only the mass of
the underlying star contributes to the gravitational acceleration.

Constraint Equation
-------------------

We begin with the :math:`w_0` constraint equation (derived elsewhere), including the volume-discrepancy term:

.. math::
    \nabla \cdot  ({w_0 \er}) = \bar{S} - \frac{1}{\gammabar p_0} \psi - \frac{f}{\gammabar p_0} \frac{p_0 - \overline{{p_0}_\mathrm{EOS}}}{\dt}

where

.. math::
    \psi = \frac{D_0 p_0}{D t} = \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial t}.

In Cartesian geometry, the divergence on the left hand side is simply :math:`\partialw_0 / \partial r`.
Let us now define :math:`\ow` and :math:`\dw` such that

.. math::
    w_0 = \ow + \dw

We take :math:`\ow` to satisfy

.. math::
    \frac{\partial \ow}{\partial r} = \bar{S} - \frac{f}{\gammabar p_0} \frac{p_0 - \overline{{p_0}_\mathrm{EOS}}}{\dt} ,

leaving the equation

.. math::
    \frac{\partial \dw}{\partial r} = - \frac{1}{\gammabar p_0} \left ({ \frac{\partial p_0}{\partial t} + w_0 \frac{\partial p_0}{\partial r}} \right ).

If we multiply by :math:`\gammabar p_0`, differentiate by r, then switch the order of :math:`\partial t`
and :math:`\partial r` where they appear in the same term, we get

.. math::
    \frac{\partial}{\partial r} \left[ \gammabar p_0 \frac{\partial \delta w_0}{\partial r} \right] =
   - \frac{\partial}{\partial t} \frac{\partial p_0}{\partial r} -  \frac{\partial}{\partial r} \left (w_0 \frac{\partial p_0}{\partial r} \right) .

We then substitute in the equation of hydrostatic equilibrium:

.. math::
    \frac{\partial p_0}{\partial r} = -\rho_0 g \quad\quad \mbox{where} \quad\quad
   g = \frac{G m}{r^2}.

We will assume that the mass of the atmosphere is negligible relative to the
mass of the core, which is outside of the simulation domain. This simplifies
the equation by allowing us to assume that :math:`m`, the enclosed mass, is constant.
So we now have

.. math::
    \frac{\partial}{\partial r} \left[ \gammabar p_0 \frac{\partial}{\partial r}{({\delta w_0})} \right]
   = \frac{\partial}{\partial t}{({\rho_0 g}}) + \frac{\partial}{\partial r}{({w_0 \rho_0 g})}
   = \rho_0 \left ({\frac{\partial}{\partial t}{g} + w_0 \frac{\partial}{\partial r}{g}} \right )
         + g \left ({\frac{\partial}{\partial t}{\rho_0} + \frac{\partial}{\partial r}{({w_0 \rho_0})}} \right ).

We now recall Equation 29 from Paper III:

.. math::
    \frac{\partial}{\partial t}{\rho_0} = - \nabla \cdot ({w_0 \rho_0 \er})
                  - \nabla \cdot ({\etarho \er}),

which is, in Cartesian geometry,

.. math::
    \frac{\partial \rho_0}{\partial t} = - \frac{\partial}{\partial r}{({w_0 \rho_0})}
                  - \frac{\partial \etarho}{\partial r}.

Substituting this expression yields

.. math::
    \frac{\partial}{\partial r}{} \left[ \gammabar p_0 \frac{\partial}{\partial r}{({\delta w_0})} \right]
   = \rho_0 \frac{D_0 g}{D t} + g \left ({- \frac{\partial}{\partial r}{({w_0 \rho_0})} - \frac{\partial \etarho}{\partial r} + \frac{\partial}{\partial r}{({w_0 \rho_0})}} \right)
   = \rho_0 \frac{D_0 g}{D t} - g \frac{\partial \etarho}{\partial r}.

We then differentiate the gravitational acceleration:

.. math::
    \begin{aligned}
    \frac{D_0 g}{D t}
   & = & \frac{D_0}{Dt} \left ({\frac{G m}{r^2}} \right ) \nonumber \\
   & = & G m \left ({\frac{\partial}{\partial t}{({r^{-2}})} + w_0 \frac{\partial}{\partial r}{({r^{-2}})}} \right ) \nonumber \\
   & = & - \frac{2 G m w_0}{r^3} \nonumber \\
   & = & - \frac{2 w_0 g}{r}.\end{aligned}

Substituting in this expression gives our final result:

.. math::
    \frac{\partial}{\partial r}{} \left[ \gammabar p_0 \frac{\partial}{\partial r}{({\dw})} \right]
   = - \frac{2 w_0 \rho_0 g}{r} - g \frac{\partial \etarho}{\partial r}

Uniform :math:`\dr` Discretization
----------------------------------

Collecting all of the :math:`\dw` terms on the left side, our constraint equation
appears as:

.. math::
    \frac{\partial}{\partial r} \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]
 + \frac{2 \rho_0 g \dw}{r} = -\frac{2 \rho_0 g \ow}{r} - g \frac{\partial \etarho}{\partial r}

On a uniform mesh (constant :math:`\dr`, we would discretize this as:

.. math::
    \begin{aligned}
    \frac{1}{\dr} \left \{ \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]_j
                     - \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]_{j-1}
              \right \}
            &+& \left [ \frac{2 \rho_0 g}{r} \dw \right ]_{j-\half} \nonumber \\
            = &-& \left [ \frac{2 \rho_0 g}{r} \ow \right ]_{j-\half}
              - \frac{g_{j-\half}}{\dr} \left [ {\etarho}_j - {\etarho}_{j-1} \right ]\end{aligned}

Expanding the :math:`\partial \dw / \partial r` terms, we have:

.. math::
    \begin{aligned}
    \frac{1}{(\dr)^2} \left \{ \left [ (\gammabar p_0)_j
                               \left ({\dw}_{j+\half} - {\dw}_{j-\half} \right ) \right ]
                      \right . &-& \left . \left [ (\gammabar p_0)_{j-1}
                               \left ({\dw}_{j-\half} - {\dw}_{j-\thalf} \right ) \right ]
              \right \}
            + \left [ \frac{2 \rho_0 g}{r} \dw \right ]_{j-\half} \nonumber \\
            = &-& \left [ \frac{2 \rho_0 g}{r} \ow \right ]_{j-\half}
              - \frac{g_{j-\half}}{\dr} \left [ {\etarho}_j - {\etarho}_{j-1} \right ]\end{aligned}

As with the spherical case (multilevel paper, appendix B), we write this in the form:

.. math::
    A_j (\dw)_{j-\thalf} + B_j (\dw)_{j-\myhalf} + C_j (\dw)_{j+\myhalf} = F_j,

then:

.. math::
    \begin{aligned}
    A_j &=& \frac{1}{\dr^2} \left( {\gammabar p_0}\right)_{j-1}, \\
    B_j &=& -\frac{1}{\dr^2} \left[ \left( {\gammabar p_0}\right)_{j}  + \left( {\gammabar p_0}\right)_{j-1} \right] +  \frac{2}{r_{j-\myhalf}} \left (\rho_0 g \right )_{j-\half}  , \\
    C_j &=& \frac{1}{\dr^2} \left( {\gammabar p_0}\right)_{j}  , \\
    F_j &=&  -\frac{2}{r_{j-\myhalf}} (\rho_0 g)_{j-\half}   (\ow)_{j-\half} - \frac{g_{j-\half}}{\dr} \left[ \left( \etarho \right)_{j} - \left( \etarho \right)_{j-1} \right] \end{aligned}

Non-Uniform :math:`\dr` Discretization
--------------------------------------

\centering
![image](\planeinvsqfigpath/grid2){width="4in"}

Consider the above non-uniform grid spacing,
where :math:`\dr_c = 2 \dr_f`. Here, the discretization of the Laplacian-like term is more complex.
We want to compute

.. math::
    \frac{\partial}{\partial r} \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]_{j-\half}

This is to be centered at :math:`j-\half`, which we accomplish by averaging the two fine grids and then
differencing:

.. math::
    \frac{\partial}{\partial r} \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]_{j-\half} =
    \frac{1}{\dr_c} \left \{ \frac{1}{2} \left [
      \left ( \gammabar p_0 \frac{\partial \dw}{\partial r} \right )_{j+1} +
      \left ( \gammabar p_0 \frac{\partial \dw}{\partial r} \right )_{j} \right ]
   - \left ( \gammabar p_0 \frac{\partial \dw}{\partial r} \right )_{j-1}
   \right \}

Expanding the :math:`\partial \dw / \partial r` terms results in a equation depending on :math:`\dw` at
4 different edge locations---this no longer fits into the tri-diagonal format used in the
uniform grid case. In detail, it becomes:

.. math::
    \begin{aligned}
    \frac{\partial}{\partial r} \left [ \gammabar p_0 \frac{\partial \dw}{\partial r} \right ]_{j-\half} &=&
    \frac{1}{\dr_c} \left \{  \frac{1}{2} \left [
      \left ( \gammabar p_0 \right )_{j+1} \frac{(\dw)_{j+\thalf} - (\dw)_{j+\half}}{\dr_f} +
      \left ( \gammabar p_0 \right )_{j}   \frac{(\dw)_{j+\half} - (\dw)_{j-\half}}{\dr_f}  \right ] \right .  \nonumber \\
      && \qquad \left .  - \left ( \gammabar p_0 \right )_{j-1} \frac{(\dw)_{j-\half} - (\dw)_{j-\thalf}}{\dr_c}
      \right \}\end{aligned}

which has terms proportional to :math:`(\dw)_{j-\thalf}`, :math:`(\dw)_{j-\half}`, :math:`(\dw)_{j+\half}`, and :math:`(\dw)_{j+\thalf}`

Boundary Conditions
-------------------

Together with the assumption that the mass of the envelope does not
contribute to the gravitational acceleration, we assume that as we move
a fluid element in the atmosphere, it does not drive a velocity at the very
base of the layer. Therefore, we take :math:`w_0(r_\mathrm{base}) = 0`.
