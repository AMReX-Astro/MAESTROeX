************************
Godunov Interface States
************************

These are working notes for the Godunov step in MAESTROeX and VARDEN.

MAESTROeX Notation
==================

-  For 2D, :math:`\Ub = (u,w)` and :math:`\Ubt = (\ut,\wt)`.
   Note that :math:`u = \ut`. We will use the shorthand :math:`\ib = (x,r)`.

-  For 3D plane parallel, :math:`\Ub = (u,v,w)`
   and :math:`\Ubt = (\ut,\vt,\wt)`. Note that :math:`u = \ut` and :math:`v = \vt`.
   We will use the shorthand :math:`\ib = (x,y,r)`.

-  For 3D spherical, :math:`\Ub = (u,v,w)`
   and :math:`\Ubt = (\ut,\vt,\wt)`. We will use the shorthand
   :math:`\ib = (x,y,z)`.

Computing :math:`\Ub` From :math:`\Ubt`
---------------------------------------

For plane-parallel problems, in order to compute :math:`w` from
:math:`\wt`, we use simple averaging

.. math:: w_{\ib}^n = \wt_{\ib}^n + \frac{w_{0,r-\half} + w_{0,r+\half}}{2}.

For spherial problems, in order to compute :math:`\Ub` from :math:`\Ubt`,
we first map :math:`w_0` onto :math:`w_0^{\mac}` using put_w0_on_edges,
where :math:`w_0^{\mac}` only contains normal velocities at each face.
Then we construct :math:`\Ub` by using

.. math:: u_{\ib} = \ut_{\ib} + \frac{w_{0,\ib+\half\eb_x}^{\mac} + w_{0,\ib-\half\eb_x}^{\mac}}{2},

.. math:: v_{\ib} = \vt_{\ib} + \frac{w_{0,\ib+\half\eb_y}^{\mac} + w_{0,\ib-\half\eb_y}^{\mac}}{2},

.. math:: w_{\ib} = \wt_{\ib} + \frac{w_{0,\ib+\half\eb_z}^{\mac} + w_{0,\ib-\half\eb_z}^{\mac}}{2}.

To compute full edge-state velocities, simply add :math:`w_0`
(for plane-parallel) or w0mac to the perturbational
velocity directly since only edge-based quantities are involved.

Computing :math:`\partial w_0/\partial r`
-----------------------------------------

For plane-parallel problems, the spatial derivatives of :math:`w_0`
are given by the two-point centered difference:

.. math:: \left(\frac{\partial w_0}{\partial r}\right)_{\ib} = \frac{w_{0,r+\half}-w_{0,r-\half}}{h}.

For spherical problems, we compute the radial bin centered gradient using

.. math:: \left(\frac{\partial w_0}{\partial r}\right)_{r} = \frac{w_{0,r+\half}-w_{0,r-\half}}{\Delta r}.

Then we put :math:`\partial w_0/\partial r` onto a Cartesian grid
using put_1d_array_on_cart_3d_sphr.


Computing :math:`\Ubt^{\trans}` in MAESTROeX
============================================

| In advance_premac, we call mkutrans, to compute
  :math:`\Ubt^{\trans}`. We will only compute the normal
  component of velocity at each face.
  These transverse velocities do not contain :math:`w_0`, so immediately
  following the call to mkutrans, we call addw0 to compute
  :math:`\Ub^{\trans}` from :math:`\Ubt^{\trans}`.
| The evolution equation for the perturbational velocity is:

  .. math:: \frac{\partial\Ubt}{\partial t} = -\Ub\cdot\nabla\Ubt \underbrace{- (\Ubt\cdot\eb_r)\frac{\partial w_0}{\partial r}\eb_r - \frac{1}{\rho}\nabla\pi + \frac{1}{\rho_0}\frac{\partial\pi_0}{\partial r}\eb_r - \frac{(\rho-\rho_0)}{\rho}g\eb_r}_{\hbox{forcing terms}}.\label{Perturbational Velocity Equation}

  We extrapolate each velocity component to edge-centered, time-centered locations. For example,

  .. math::

     \begin{aligned}
     \ut_{R,\ib-\half\eb_x} &=& \ut_{\ib}^n + \frac{h}{2}\frac{\partial\ut_{\ib}^n}{\partial x} + \frac{\dt}{2}\frac{\partial\ut_{\ib}^n}{\partial t} \nonumber \\
     &=& \ut_{\ib}^n + \frac{h}{2}\frac{\partial\ut_{\ib}^n}{\partial x} + \frac{\dt}{2}
     \left(-\ut_{\ib}^n\frac{\partial\ut_{\ib}^n}{\partial x} - \wt_{\ib}^n\frac{\partial\ut_{\ib}^n}{\partial r} + \text{forcing terms}\right)\end{aligned}

  We are going to use a 1D Taylor series extrapolation in space and time.
  By 1D, we mean that we omit any spatial derivatives that are not in the
  direction of the extrapolation. We also omit the underbraced forcing terms.
  We also use characteristic tracing.

  .. math:: \ut_{R,\ib-\half\eb_x} = \ut_{\ib}^n + \left[\frac{1}{2} - \frac{\dt}{2h}\min(0,\ut_{\ib}^n)\right]\partial\ut_{\ib}^n


2D Cartesian Case
-----------------

We predict :math:`\ut` to x-faces using a 1D extrapolation:

.. math::

   \begin{aligned}
   \ut_{L,\ib-\half\eb_x} &=& \ut_{\ib-\eb_x}^n + \left[\half - \frac{\dt}{2h}\max(0,u_{\ib-\eb_x}^n)\right]\Delta_x \ut_{\ib-\eb_x}^n,\\
   \ut_{R,\ib-\half\eb_x} &=& \ut_{\ib}^n - \left[\half + \frac{\dt}{2h}\min(0,u_{\ib}^n)\right]\Delta_x \ut_{\ib}^n.\end{aligned}

We pick the final trans states using a Riemann solver:

.. math::

   \ut^{\trans}_{\ib-\half\eb_x} =
   \begin{cases}
   0, & \left(\ut_{L,\ib-\half\eb_x} \le 0 ~~ {\rm AND} ~~ \ut_{R,\ib-\half\eb_x} \ge 0\right) ~~ {\rm OR} ~~ \left|\ut_{L,\ib-\half\eb_x} + \ut_{R,\ib-\half\eb_x}\right| < \epsilon, \\
   \ut_{L,\ib-\half\eb_x}, & \ut_{L,\ib-\half\eb_x} + \ut_{R,\ib-\half\eb_x} > 0, \\
   \ut_{R,\ib-\half\eb_x}, & \ut_{L,\ib-\half\eb_x} + \ut_{R,\ib-\half\eb_x} < 0, \\
   \end{cases}

We predict :math:`\wt` to r-faces using a 1D extrapolation:

.. math::

   \begin{aligned}
   \wt_{L,\ib-\half\eb_r} &=& \wt_{\ib-\eb_r}^n + \left[\half - \frac{\dt}{2h}\max(0,w_{\ib-\eb_r}^n)\right]\Delta_r \wt_{\ib-\eb_r}^n,\\
   \wt_{R,\ib-\half\eb_r} &=& \wt_{\ib}^n - \left[\half + \frac{\dt}{2h}\min(0,w_{\ib}^n)\right]\Delta_r \wt_{\ib}^n.\end{aligned}

We pick the final :math:`\trans` states using a Riemann solver, noting
that we upwind based on the full velocity.

.. math::

   \wt^{\trans}_{\ib-\half\eb_r} =
   \begin{cases}
   0, & \left(w_{L,\ib-\half\eb_r} \le 0 ~~ {\rm AND} ~~ w_{R,\ib-\half\eb_r} \ge 0\right) ~~ {\rm OR} ~~ \left|w_{L,\ib-\half\eb_r} + w_{R,\ib-\half\eb_r}\right| < \epsilon, \\
   \wt_{L,\ib-\half\eb_r}, & w_{L,\ib-\half\eb_r} + w_{R,\ib-\half\eb_r} > 0, \\
   \wt_{R,\ib-\half\eb_r}, & w_{L,\ib-\half\eb_r} + w_{R,\ib-\half\eb_r} < 0, \\
   \end{cases}


.. _d-cartesian-case-1:

3D Cartesian Case
-----------------

We use the exact same procedure in 2D and 3D to compute :math:`\ut^{\trans}` and
:math:`\wt^{\trans}`. The procedure for computing :math:`\vt^{\trans}` is analogous to
computing :math:`\ut^{\trans}`. We predict :math:`\vt` to y-faces using the
1D extrapolation:

.. math::

   \begin{aligned}
   \vt_{L,\ib-\half\eb_y} &=& \vt_{\ib-\eb_y}^n + \left[\half - \frac{\dt}{2h}\max(0,v_{\ib-\eb_y}^n)\right]\Delta_y \vt_{\ib-\eb_y}^n, \\
   \vt_{R,\ib-\half\eb_y} &=& \vt_{\ib}^n - \left[\half + \frac{\dt}{2h}\min(0,v_{\ib}^n)\right]\Delta_y \vt_{\ib}^n,\end{aligned}

.. math::

   \vt^{\trans}_{\ib-\half\eb_y} =
   \begin{cases}
   0, & \left(v_{L,\ib-\half\eb_y} \le 0 ~~ {\rm AND} ~~ v_{R,\ib-\half\eb_y} \ge 0\right) ~~ {\rm OR} ~~ \left|v_{L,\ib-\half\eb_y} + v_{R,\ib-\half\eb_y}\right| < \epsilon, \\
   \vt_{L,\ib-\half\eb_y}, & v_{L,\ib-\half\eb_y} + v_{R,\ib-\half\eb_y} > 0, \\
   \vt_{R,\ib-\half\eb_y}, & v_{L,\ib-\half\eb_y} + v_{R,\ib-\half\eb_y} < 0. \\
   \end{cases}


3D Spherical Case
-----------------

We predict the normal components of velocity to the normal faces
using a 1D extrapolation. The equations for all three directions
are identical to those given in the 2D and 3D plane-parallel
sections. As in the plane-parallel case, make sure
that the advection velocities, as well as
the upwind velocity, is done with the full velocity, not the
perturbational velocity.


Computing :math:`\Ubt^{\mac,*}` in MAESTROeX
============================================

| In advance_premac, we call velpred to compute
  :math:`\Ubt^{\mac,*}`. We will only compute the normal component of
  velocity at each face.
| For reference, here is the perturbational velocity equation from before:

  .. math:: \frac{\partial\Ubt}{\partial t} = -\Ub\cdot\nabla\Ubt \underbrace{- (\Ubt\cdot\eb_r)\frac{\partial w_0}{\partial r}\eb_r \underbrace{- \frac{1}{\rho}\nabla\pi + \frac{1}{\rho_0}\frac{\partial\pi_0}{\partial r}\eb_r - \frac{(\rho-\rho_0)}{\rho}g\eb_r}_{\hbox{terms included in $\fb_{\Ubt}$}}}_{\hbox{forcing terms}}.

  Note that the :math:`\partial w_0/\partial r` term is treated like a forcing
  term, but it is not actually part of :math:`\fb_{\Ubt}`. We make use of the 1D
  extrapolations used to compute :math:`\Ubt^{\trans}`
  (:math:`\ut_{L/R,\ib-\half\eb_x}`, :math:`\vt_{L/R,\ib-\half\eb_y}`,
  and :math:`\wt_{L/R,\ib-\half\eb_r}`), as well as the “:math:`\trans`” states
  (:math:`\ut_{\ib-\half\eb_x}^{\trans}`, :math:`\vt_{\ib-\half\eb_y}^{\trans}`,
  and :math:`\wt_{\ib-\half\eb_r}^{\trans}`)


.. _d-cartesian-case-2:

2D Cartesian Case
-----------------

#. Predict :math:`\ut` to r-faces using a 1D extrapolation.

#. Predict :math:`\ut` to x-faces using a full-dimensional extrapolation.

#. Predict :math:`\wt` to x-faces using a 1D extrapolation.

#. Predict :math:`\wt` to r-faces using a full-dimensional extrapolation.

Predict :math:`\ut` to r-faces using a 1D extrapolation:

.. math::

   \begin{aligned}
   \ut_{L,\ib-\half\eb_r} &=& \ut_{\ib-\eb_r}^n + \left[\half - \frac{\dt}{2h}\max(0,w_{\ib-\eb_r}^n)\right]\Delta_r \ut_{\ib-\eb_r}^n, \\
   \ut_{R,\ib-\half\eb_r} &=& \ut_{\ib} - \left[\half + \frac{\dt}{2h}\min(0,w_{\ib}^n)\right]\Delta_r \ut_{\ib}^n.\end{aligned}

Upwind based on :math:`w^{\trans}`:

.. math::

   \ut_{\ib-\half\eb_r} =
   \begin{cases}
   \half\left(\ut_{L,\ib-\half\eb_r} + \ut_{R,\ib-\half\eb_r}\right), & \left|w^{\trans}_{\ib-\half\eb_r}\right| < \epsilon \\
   \ut_{L,\ib-\half\eb_r}, & w^{\trans}_{\ib-\half\eb_r} > 0, \\
   \ut_{R,\ib-\half\eb_r}, & w^{\trans}_{\ib-\half\eb_r} < 0. \\
   \end{cases}

Predict :math:`\ut` to x-faces using a full-dimensional extrapolation,

.. math::

   \begin{aligned}
   \ut_{L,\ib-\half\eb_x}^{\mac,*} &=& \ut_{L,\ib-\half\eb_x} - \frac{\dt}{4h}\left(w_{\ib-\eb_x+\half\eb_r}^{\trans}+w_{\ib-\eb_x-\half\eb_r}^{\trans}\right)\left(\ut_{\ib-\eb_x+\half\eb_r} - \ut_{\ib-\eb_x-\half\eb_r}\right) + \frac{\dt}{2}f_{\ut,\ib-\eb_x}, \nonumber \\
   && \\
   \ut_{R,\ib-\half\eb_x}^{\mac,*} &=& \ut_{R,\ib-\half\eb_x} - \frac{\dt}{4h}\left(w_{\ib+\half\eb_r}^{\trans}+w_{\ib-\half\eb_r}^{\trans}\right)\left(\ut_{\ib+\half\eb_r} - \ut_{\ib-\half\eb_r}\right) + \frac{\dt}{2}f_{\ut,\ib}.\end{aligned}

Solve a Riemann problem:

.. math::

   \ut_{\ib-\half\eb_x}^{\mac,*} =
   \begin{cases}
   0, & \left(u_{L,\ib-\half\eb_x}^{\mac,*} \le 0 ~~ {\rm AND} ~~ u_{R,\ib-\half\eb_x}^{\mac,*} \ge 0\right) ~~ {\rm OR} ~~ \left|u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*}\right| < \epsilon, \\
   \ut_{L,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} > 0, \\
   \ut_{R,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} < 0.
   \end{cases}

Predict :math:`\wt` to x-faces using a 1D extrapolation:

.. math::

   \begin{aligned}
   \wt_{L,\ib-\half\eb_x} &=& \wt_{\ib-\eb_x}^n + \left[\half - \frac{\dt}{2h}\max(0,u_{\ib-\eb_x}^n)\right]\Delta_x \wt_{\ib-\eb_x}^n, \\
   \wt_{R,\ib-\half\eb_x} &=& \wt_{\ib} - \left[\half + \frac{\dt}{2h}\min(0,u_{\ib}^n)\right]\Delta_x \wt_{\ib}^n.\end{aligned}

Upwind based on :math:`u^{\trans}`:

.. math::

   \wt_{\ib-\half\eb_x} =
   \begin{cases}
   \half\left(\wt_{L,\ib-\half\eb_x} + \wt_{R,\ib-\half\eb_x}\right), & \left|u^{\trans}_{\ib-\half\eb_x}\right| < \epsilon \\
   \wt_{L,\ib-\half\eb_x}, & u^{\trans}_{\ib-\half\eb_x} > 0, \\
   \wt_{R,\ib-\half\eb_x}, & u^{\trans}_{\ib-\half\eb_x} < 0. \\
   \end{cases}

Predict :math:`\wt` to r-faces using a full-dimensional extrapolation:

.. math::

   \begin{aligned}
   \wt_{L,\ib-\half\eb_r}^{\mac,*} = \wt_{L,\ib-\half\eb_r} &-& \frac{\dt}{4h}\left(u_{\ib-\eb_r+\half\eb_x}^{\trans}+u_{\ib-\eb_r-\half\eb_x}^{\trans}\right)\left(\wt_{\ib-\eb_r+\half\eb_x} - \wt_{\ib-\eb_r-\half\eb_x}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(\wt_{\ib-\half\eb_r}^{\trans}+\wt_{\ib-\frac{3}{2}\eb_r}^{\trans}\right)\left(w_{0,\ib-\half\eb_r} - w_{0,\ib-\frac{3}{2}\eb_r}\right) + \frac{\dt}{2}f_{\wt,\ib-\eb_r}, \nonumber \\
   && \\
   \wt_{R,\ib-\half\eb_r}^{\mac,*} = \wt_{R,\ib-\half\eb_r} &-& \frac{\dt}{4h}\left(u_{\ib+\half\eb_x}^{\trans}+u_{\ib-\half\eb_x}^{\trans}\right)\left(\wt_{\ib+\half\eb_x} - \wt_{\ib-\half\eb_x}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(\wt_{\ib+\half\eb_r}^{\trans}+\wt_{\ib-\half\eb_r}^{\trans}\right)\left(w_{0,\ib+\half\eb_r} - w_{0,\ib-\half\eb_r}\right) + \frac{\dt}{2}f_{\wt,\ib}.\end{aligned}

Solve a Riemann problem:

.. math::

   \wt_{\ib-\half\eb_r}^{\mac,*} =
   \begin{cases}
   0, & \left(w_L^{\mac,*} \le 0 ~~ {\rm AND} ~~ w_R^{\mac,*} \ge 0\right) ~~ {\rm OR} ~~ \left|w_L^{\mac,*} + w_R^{\mac,*}\right| < \epsilon, \\
   \wt_{L,\ib-\half\eb_r}^{\mac,*}, & w_L^{\mac,*} + w_R^{\mac,*} > 0, \\
   \wt_{R,\ib-\half\eb_r}^{\mac,*}, & w_L^{\mac,*} + w_R^{\mac,*} < 0.
   \end{cases}


.. _d-cartesian-case-3:

3D Cartesian Case
-----------------

This algorithm is more complicated than the 2D case since we include
the effects of corner coupling.

#. Predict :math:`\ut` to y-faces using a 1D extrapolation.

#. Predict :math:`\ut` to r-faces using a 1D extrapolation.

#. Predict :math:`\vt` to x-faces using a 1D extrapolation.

#. Predict :math:`\vt` to r-faces using a 1D extrapolation.

#. Predict :math:`\wt` to x-faces using a 1D extrapolation.

#. Predict :math:`\wt` to y-faces using a 1D extrapolation.

#. Update prediction of :math:`\ut` to y-faces by accounting for :math:`r`-derivatives.

#. Update prediction of :math:`\ut` to r-faces by accounting for :math:`y`-derivatives.

#. Update prediction of :math:`\vt` to x-faces by accounting for :math:`r`-derivatives.

#. Update prediction of :math:`\vt` to r-faces by accounting for :math:`x`-derivatives.

#. Update prediction of :math:`\wt` to x-faces by accounting for :math:`y`-derivatives.

#. Update prediction of :math:`\wt` to y-faces by accounting for :math:`x`-derivatives.

#. Predict :math:`\ut` to x-faces using a full-dimensional extrapolation.

#. Predict :math:`\vt` to y-faces using a full-dimensional extrapolation.

#. Predict :math:`\wt` to r-faces using a full-dimensional extrapolation.

* Predict :math:`\ut` to y-faces using a 1D extrapolation.

  .. math::

     \begin{aligned}
     \ut_{L,\ib-\half\eb_y} &=& \ut_{\ib-\eb_y}^n + \left[\half - \frac{\dt}{2h}\max(0,v_{\ib-\eb_y}^n)\right]\Delta_y \ut_{\ib-\eb_y}^n, \\
     \ut_{R,\ib-\half\eb_y} &=& \ut_{\ib} - \left[\half + \frac{\dt}{2h}\min(0,v_{\ib}^n)\right]\Delta_y \ut_{\ib}^n.\end{aligned}

  Upwind based on :math:`v^{\trans}`:

  .. math::

     \ut_{\ib-\half\eb_y} =
     \begin{cases}
     \half\left(\ut_{L,\ib-\half\eb_y} + \ut_{R,\ib-\half\eb_y}\right), & \left|v^{\trans}_{\ib-\half\eb_y}\right| < \epsilon \\
     \ut_{L,\ib-\half\eb_y}, & v^{\trans}_{\ib-\half\eb_y} > 0, \\
     \ut_{R,\ib-\half\eb_y}, & v^{\trans}_{\ib-\half\eb_y} < 0. \\
     \end{cases}

* Predict :math:`\ut` to r-faces using a 1D extrapolation.

* Predict :math:`\vt` to x-faces using a 1D extrapolation.

* Predict :math:`\vt` to r-faces using a 1D extrapolation.

* Predict :math:`\wt` to x-faces using a 1D extrapolation.

* Predict :math:`\wt` to y-faces using a 1D extrapolation.

* Update prediction of :math:`\ut` to y-faces by accounting for :math:`r`-derivatives.
  The notation :math:`\ut_{\ib-\half\eb_y}^{y|r}` means state :math:`\ut_{\ib-\half\eb_y}` that has been updated to account for transverse derives in the r-direction.

  .. math::

     \begin{aligned}
     \ut_{L,\ib-\half\eb_y}^{y|r} &=& \ut_{L,\ib-\half\eb_y} - \frac{\dt}{6h}\left(w_{\ib-\eb_y+\half\eb_r}^{\trans}+w_{\ib-\eb_y-\half\eb_r}^{\trans}\right)\left(\ut_{\ib-\eb_y+\half\eb_r}-\ut_{\ib-\eb_y-\half\eb_r}\right), \\
     \ut_{R,\ib-\half\eb_y}^{y|r} &=& \ut_{R,\ib-\half\eb_y} - \frac{\dt}{6h}\left(w_{\ib+\half\eb_r}^{\trans}+w_{\ib-\half\eb_r}^{\trans}\right)\left(\ut_{\ib+\half\eb_r}-\ut_{\ib-\half\eb_r}\right).\end{aligned}

  Upwind based on :math:`v^{\trans}`:

  .. math::

     \ut_{\ib-\half\eb_y}^{y|r} =
     \begin{cases}
     \half\left(\ut_{L,\ib-\half\eb_y}^{y|r} + \ut_{R,\ib-\half\eb_y}^{y|r}\right), & \left|v_{\ib-\half\eb_y}^{\trans}\right| < \epsilon \\
     \ut_{L,\ib-\half\eb_y}^{y|r}, & v_{\ib-\half\eb_y}^{\trans} > 0, \\
     \ut_{R,\ib-\half\eb_y}^{y|r}, & v_{\ib-\half\eb_y}^{\trans} < 0.
     \end{cases}

* Update prediction of :math:`\ut` to r-faces by accounting for :math:`y`-derivatives.
* Update prediction of :math:`\vt` to x-faces by accounting for :math:`r`-derivatives.
* Update prediction of :math:`\vt` to r-faces by accounting for :math:`x`-derivatives.
* Update prediction of :math:`\wt` to x-faces by accounting for :math:`y`-derivatives.
* Update prediction of :math:`\wt` to y-faces by accounting for :math:`x`-derivatives.
* Predict :math:`\ut` to x-faces using a full-dimensional extrapolation.

  .. math::

     \begin{aligned}
     \ut_{L,\ib-\half\eb_x}^{\mac,*} = \ut_{L,\ib-\half\eb_x} &-& \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{\trans}+v_{\ib-\eb_x-\half\eb_y}^{\trans}\right)\left(\ut_{\ib-\eb_x+\half\eb_y}^{y|r}-\ut_{\ib-\eb_x-\half\eb_y}^{y|r}\right) \nonumber \\
     &-& \frac{\dt}{4h}\left(w_{\ib-\eb_x+\half\eb_r}^{\trans}+w_{\ib-\eb_x-\half\eb_r}^{\trans}\right)\left(\ut_{\ib-\eb_x+\half\eb_r}^{r|y}-\ut_{\ib-\eb_x-\half\eb_r}^{r|y}\right) + \frac{\dt}{2}f_{u,\ib-\eb_x}, \nonumber \\
     && \\
     \ut_{R,\ib-\half\eb_x}^{\mac,*} = \ut_{R,\ib-\half\eb_x} &-& \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{\trans}+v_{\ib-\half\eb_y}^{\trans}\right)\left(\ut_{\ib+\half\eb_y}^{y|r}-\ut_{\ib-\half\eb_y}^{y|r}\right) \nonumber \\
     &-& \frac{\dt}{4h}\left(w_{\ib+\half\eb_r}^{\trans}+w_{\ib-\half\eb_r}^{\trans}\right)\left(\ut_{\ib+\half\eb_r}^{r|y}-\ut_{\ib-\half\eb_r}^{r|y}\right) + \frac{\dt}{2}f_{u,\ib}.\end{aligned}

  Solve a Riemann problem:

  .. math::

     \ut_{\ib-\half\eb_x}^{\mac,*} =
     \begin{cases}
     0, & \left(u_{L,\ib-\half\eb_x}^{\mac,*} \le 0 ~~ {\rm AND} ~~ u_{R,\ib-\half\eb_x}^{\mac,*} \ge 0\right) ~~ {\rm OR} ~~ \left|u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*}\right| < \epsilon, \\
     \ut_{L,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} > 0, \\
     \ut_{R,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} < 0.
     \end{cases}

* Predict :math:`\vt` to y-faces using a full-dimensional extrapolation.
* Predict :math:`\wt` to r-faces using a full-dimensional extrapolation.
  In this step, make sure you account for the :math:`\partial w_0/\partial r`
  term before solving the Riemann problem:

  .. math::

     \begin{aligned}
     \wt_{L,\ib-\half\eb_r}^{\mac,*} &=& \wt_{L,\ib-\half\eb_r}^{\mac,*} -
     \frac{\dt}{4h}\left(\wt^{\trans}_{\ib+\half\eb_r} + \wt^{\trans}_{\ib-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\
     \wt_{R,\ib-\half\eb_r}^{\mac,*} &=& \wt_{R,\ib-\half\eb_r}^{\mac,*} -
     \frac{\dt}{4h}\left(\wt^{\trans}_{\ib-\half\eb_r} + \wt^{\trans}_{\ib-\frac{3}{2}\eb_r}\right)\left(w_{0,\ib-\half\eb_r}-w_{0,\ib-\frac{3}{2}\eb_r}\right)\end{aligned}


.. _d-spherical-case-1:

3D Spherical Case
-----------------

The spherical case is the same as the plane-parallel 3D Cartesian
case, except the :math:`\partial w_0/\partial r` term enters
in the full dimensional extrapolation for each direction.
As in the plane-parallel case, make sure to upwind using the
full velocity.


.. _Scalar Edge State Prediction in MAESTROeX:

Computing :math:`\rho^{'\edge}, X_k^{\edge},(\rho h)^{'\edge}`, and :math:`\Ubt^{\edge}` in MAESTROeX
=====================================================================================================

We call make_edge_scal to compute :math:`\rho^{'\edge}, X_k^{\edge},
(\rho h)^{'\edge}`, and :math:`\Ubt^{\edge}` at each edge.
The procedure is the same for each quantity, so we shall simply denote
the scalar as :math:`s`. We always need to compute :math:`\rho'` and :math:`X_k` to faces,
and the choice of energy prediction is as follows:

-  For enthalpy_pred_type = 1, we predict :math:`(\rho h)'` to faces.

-  For enthalpy_pred_type = 2, we predict :math:`h` to faces.

-  For enthalpy_pred_type = 3 and 4, we predict :math:`T` to faces.

-  For enthalpy_pred_type = 5, we predict :math:`h'` to faces.

-  For enthalpy_pred_type = 6, we predict :math:`T'` to faces.

We are using enthalpy_pred_type = 1 for now. The equations
of motion are:

.. math::

   \begin{aligned}
   \frac{\partial \rho'}{\partial t} &=& -\Ub\cdot\nabla\rho' \underbrace{- \rho'\nabla\cdot\Ub - \nabla\cdot\left(\rho_0\Ubt\right)}_{f_{\rho'}}, \\
   \frac{\partial X_k}{\partial t} &=& -\Ub\cdot\nabla X_k ~~~ \text{(no forcing)}, \\
   \frac{\partial(\rho h)'}{\partial t} &=& -\Ub\cdot\nabla(\rho h)' \underbrace{- (\rho h)'\nabla\cdot\Ub - \nabla\cdot\left[(\rho h)_0\Ubt\right] + \left(\Ubt\cdot\eb_r\right)\frac{\partial p_0}{\partial r} + \nabla\cdot\kth\nabla T}_{f_{(\rho h)'}}, \nonumber \\
   && \\
   \frac{\partial\Ubt}{\partial t} &=& -\Ub\cdot\nabla\Ubt \underbrace{- \left(\Ubt\cdot\eb_r\right)\frac{\partial w_0}{\partial r}\eb_r \underbrace{- \frac{1}{\rho}\nabla\pi + \frac{1}{\rho_0}\frac{\partial\pi_0}{\partial r}\eb_r - \frac{(\rho-\rho_0)}{\rho}g\eb_r}_{\hbox{terms included in $\fb_{\Ubt}$}}}_{\hbox{forcing terms}}.\end{aligned}


.. _d-cartesian-case-4:

2D Cartesian Case
-----------------

#. Predict :math:`s` to r-faces using a 1D extrapolation.

#. Predict :math:`s` to x-faces using a full-dimensional extrapolation.

#. Predict :math:`s` to x-faces using a 1D extrapolation.

#. Predict :math:`s` to r-faces using a full-dimensional extrapolation.

* Predict :math:`s` to r-faces using a 1D extrapolation:

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_r} &=& s_{\ib-\eb_r}^n + \left(\half - \frac{\dt}{2h}w_{\ib-\half\eb_r}^{\mac}\right)\Delta_r s_{\ib-\eb_r}^n, \\
     s_{R,\ib-\half\eb_r} &=& s_{\ib} - \left(\half + \frac{\dt}{2h}w_{\ib-\half\eb_r}^{\mac}\right)\Delta_r s_{\ib}^n.\end{aligned}

  Upwind based on :math:`w^{\mac}`:

  .. math::

     s_{\ib-\half\eb_r} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_r} + s_{R,\ib-\half\eb_r}\right), & \left|w^{\mac}_{\ib-\half\eb_r}\right| < \epsilon \\
     s_{L,\ib-\half\eb_r}, & w^{\mac}_{\ib-\half\eb_r} > 0, \\
     s_{R,\ib-\half\eb_r}, & w^{\mac}_{\ib-\half\eb_r} < 0. \\
     \end{cases}

  Predict :math:`s` to x-faces using a full-dimensional extrapolation. First, the normal derivative and forcing terms:

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} &=& s_{\ib-\eb_x}^n + \left(\half - \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib-\eb_x}^n + \frac{\dt}{2}f_{\ib-\eb_x}^n \\
     s_{R,\ib-\half\eb_x}^{\edge} &=& s_{\ib}^n - \left(\half + \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib}^n + \frac{\dt}{2}f_{\ib}^n \end{aligned}

  Account for the transverse terms:

  **if** is_conservative **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} &=& s_{L,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{2h}\left[\left(w^{\mac}s\right)_{\ib-\eb_x+\half\eb_r} - \left(w^{\mac}s\right)_{\ib-\eb_x-\half\eb_r}\right] - \frac{\dt}{2h}s_{\ib-\eb_x}^{n}\left(u_{\ib-\half\eb_x}^{\mac}-u_{\ib-\frac{3}{2}\eb_x}^{\mac}\right)\nonumber \\
     &&\\
     s_{R,\ib-\half\eb_x}^{\edge} &=& s_{R,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{2h}\left[\left(w^{\mac}s\right)_{\ib+\half\eb_r} - \left(w^{\mac}s\right)_{\ib-\half\eb_r}\right] - \frac{\dt}{2h}s_{\ib}^{n}\left(u_{\ib+\half\eb_x}^{\mac}-u_{\ib-\half\eb_x}^{\mac}\right)\end{aligned}

  **else**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} &=& s_{L,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(w^{\mac}_{\ib-\eb_x+\half\eb_r} + w^{\mac}_{\ib-\eb_x-\half\eb_r}\right)\left(s_{\ib-\eb_x+\half\eb_r} - s_{\ib-\eb_x-\half\eb_r}\right)\\
     s_{R,\ib-\half\eb_x}^{\edge} &=& s_{R,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(w^{\mac}_{\ib+\half\eb_r} + w^{\mac}_{\ib-\half\eb_r}\right)\left(s_{\ib+\half\eb_r} - s_{\ib-\half\eb_r}\right)\end{aligned}

  **end if**

* Account for the :math:`\partial w_0/\partial r` term:

  **if** is_vel **and** comp = 2 **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} &=& s_{L,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib-\eb_x+\half\eb_r} + \wt^{\mac}_{\ib-\eb_x-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\
     s_{R,\ib-\half\eb_x}^{\edge} &=& s_{R,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib+\half\eb_r} + \wt^{\mac}_{\ib-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\\end{aligned}

  **end if**

* Upwind based on :math:`u^{\mac}`.

  .. math::

     s_{\ib-\half\eb_x}^{\edge} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_x}^{\edge} + s_{R,\ib-\half\eb_x}^{\edge}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
     s_{L,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
     s_{R,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} < 0.
     \end{cases}

  Predict :math:`s` to x-faces using a 1D extrapolation:

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x} &=& s_{\ib-\eb_x}^n + \left(\half - \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib-\eb_x}^n, \\
     s_{R,\ib-\half\eb_x} &=& s_{\ib} - \left(\half + \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib}^n.\end{aligned}

  Upwind based on :math:`u^{\mac}`:

  .. math::

     s_{\ib-\half\eb_x} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_x} + s_{R,\ib-\half\eb_x}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
     s_{L,\ib-\half\eb_x}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
     s_{R,\ib-\half\eb_x}, & u^{\mac}_{\ib-\half\eb_x} < 0. \\
     \end{cases}

  Predict :math:`s` to r-faces using a full-dimensional extrapolation. First, the normal derivative and forcing terms:

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_r}^{\edge} &=&  s_{\ib-\eb_r}^n + \left(\half - \frac{\dt}{2h}w_{\ib-\half\eb_r}^{\mac}\right)\Delta_r s_{\ib-\eb_r}^n + \frac{\dt}{2}f_{\ib-\eb_r}^n \\
     s_{R,\ib-\half\eb_r}^{\edge} &=&  s_{\ib}^n - \left(\half + \frac{\dt}{2h}w_{\ib-\half\eb_r}^{\mac}\right)\Delta_r s_{\ib}^n + \frac{\dt}{2}f_{\ib}^n \end{aligned}

  Account for the transverse terms:
  **if** is_conservative **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_r}^{\edge} &=& s_{L,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{2h}\left[\left(u^{\mac}s\right)_{\ib-\eb_r+\half\eb_x} - \left(u^{\mac}s\right)_{\ib-\eb_r-\half\eb_x}\right] - \frac{\dt}{2h}s_{\ib-\eb_r}^{n}\left(w_{\ib-\half\eb_r}^{\mac}-w_{\ib-\frac{3}{2}\eb_r}^{\mac}\right)\nonumber\\
     && \\
     s_{R,\ib-\half\eb_r}^{\edge} &=& s_{R,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{2h}\left[\left(u^{\mac}s\right)_{\ib+\half\eb_x} - \left(u^{\mac}s\right)_{\ib-\half\eb_x}\right] - \frac{\dt}{2h}s_{\ib}^{n}\left(w_{\ib+\half\eb_r}^{\mac}-w_{\ib-\half\eb_r}^{\mac}\right)\end{aligned}

  **else**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_r}^{\edge} &=& s_{L,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{4h}\left(u^{\mac}_{\ib-\eb_r+\half\eb_x} + u^{\mac}_{\ib-\eb_r-\half\eb_x}\right)\left(s_{\ib-\eb_r+\half\eb_x} - s_{\ib-\eb_r-\half\eb_x}\right)\\
     s_{R,\ib-\half\eb_r}^{\edge} &=& s_{R,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{4h}\left(u^{\mac}_{\ib+\half\eb_x} + u^{\mac}_{\ib-\half\eb_x}\right)\left(s_{\ib+\half\eb_x} - s_{\ib-\half\eb_x}\right)\end{aligned}

  **end if**

* Account for the :math:`\partial w_0/\partial r` term:
  **if** is_vel **and** comp = 2 **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_r}^{\edge} &=& s_{L,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib-\half\eb_r} + \wt^{\mac}_{\ib-\frac{3}{2}\eb_r}\right)\left(w_{0,\ib-\half\eb_r}-w_{0,\ib-\frac{3}{2}\eb_r}\right) \\
     s_{R,\ib-\half\eb_r}^{\edge} &=& s_{R,\ib-\half\eb_r}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib+\half\eb_r} + \wt^{\mac}_{\ib-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\\end{aligned}

  **end if**

* Upwind based on :math:`w^{\mac}`:

  .. math::

     s_{\ib-\half\eb_r} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_r} + s_{R,\ib-\half\eb_r}\right), & \left|w^{\mac}_{\ib-\half\eb_r}\right| < \epsilon \\
     u_{L,\ib-\half\eb_r}, & w^{\mac}_{\ib-\half\eb_r} > 0, \\
     u_{R,\ib-\half\eb_r}, & w^{\mac}_{\ib-\half\eb_r} < 0. \\
     \end{cases}


.. _d-cartesian-case-5:

3D Cartesian Case
-----------------

This algorithm is more complicated than the 2D case since we include
the effects of corner coupling.

#. Predict :math:`s` to x-faces using a 1D extrapolation.

#. Predict :math:`s` to y-faces using a 1D extrapolation.

#. Predict :math:`s` to r-faces using a 1D extrapolation.

#. Update prediction of :math:`s` to x-faces by accounting for y-derivatives.

#. Update prediction of :math:`s` to x-faces by accounting for r-derivatives.

#. Update prediction of :math:`s` to y-faces by accounting for x-derivatives.

#. Update prediction of :math:`s` to y-faces by accounting for r-derivatives.

#. Update prediction of :math:`s` to r-faces by accounting for x-derivatives.

#. Update prediction of :math:`s` to r-faces by accounting for y-derivatives.

#. Predict :math:`s` to x-faces using a full-dimensional extrapolation.

#. Predict :math:`s` to y-faces using a full-dimensional extrapolation.

#. Predict :math:`s` to r-faces using a full-dimensional extrapolation.

* Predict :math:`s` to x-faces using a 1D extrapolation.

  .. math::
     s_{L,\ib-\half\eb_x} &=& s_{\ib-\eb_x}^n + \left(\half - \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib-\eb_x}^n,
    :label: 3D predict s to left

  .. math::
     s_{R,\ib-\half\eb_x} &=& s_{\ib} - \left(\half + \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib}^n.
     :label: 3D predict s to right

  Upwind based on :math:`u^{\mac}`:

  .. math::

     s_{\ib-\half\eb_x} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_x} + s_{R,\ib-\half\eb_x}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
     s_{L,\ib-\half\eb_x}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
     s_{R,\ib-\half\eb_x}, & u^{\mac}_{\ib-\half\eb_x} < 0. \\
     \end{cases}

* Predict :math:`s` to y-faces using a 1D extrapolation.

* Predict :math:`s` to r-faces using a 1D extrapolation.

* Update prediction of :math:`s` to x-faces by accounting for y-derivatives.
  The notation :math:`s_{\ib-\half\eb_x}^{x|y}` means “state :math:`s_{\ib-\half\eb_x}`
  that has been updated to account for the transverse derivatives in
  the :math:`y`-direction”.

  **if** is_conservative **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{x|y} &=& s_{L,\ib-\half\eb_x} - \frac{\dt}{3h}\left[(sv^{\mac})_{\ib-\eb_x+\half\eb_y}-(sv^{\mac})_{\ib-\eb_x-\half\eb_y}\right], \\
     s_{R,\ib-\half\eb_x}^{x|y} &=& s_{R,\ib-\half\eb_x} - \frac{\dt}{3h}\left[(sv^{\mac})_{\ib+\half\eb_y}-(sv^{\mac})_{\ib-\half\eb_y}\right].\end{aligned}

  **else**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{x|y} &=& s_{L,\ib-\half\eb_x} - \frac{\dt}{6h}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac} + v_{\ib-\eb_x-\half\eb_y}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_y} - s_{\ib-\eb_x-\half\eb_y}\right), \\
     s_{R,\ib-\half\eb_x}^{x|y} &=& s_{R,\ib-\half\eb_x} - \frac{\dt}{6h}\left(v_{\ib+\half\eb_y}^{\mac} + v_{\ib-\half\eb_y}^{\mac}\right)\left(s_{\ib+\half\eb_y} - s_{\ib-\half\eb_y}\right).\end{aligned}

  **end if**

* Upwind based on :math:`u^{\mac}`:

  .. math::

     s_{\ib-\half\eb_x}^{x|y} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_x}^{x|y} + s_{R,\ib-\half\eb_x}^{x|y}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
     s_{L,\ib-\half\eb_x}^{x|y}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
     s_{R,\ib-\half\eb_x}^{x|y}, & u^{\mac}_{\ib-\half\eb_x} < 0.
     \end{cases}

* Update prediction of :math:`s` to x-faces by accounting for r-derivatives.

* Update prediction of :math:`s` to y-faces by accounting for x-derivatives.

* Update prediction of :math:`s` to y-faces by accounting for r-derivatives.

* Update prediction of :math:`s` to r-faces by accounting for x-derivatives.

* Update prediction of :math:`s` to r-faces by accounting for y-derivatives.

* Predict :math:`s` to x-faces using a full-dimensional extrapolation.

  **if** is_conservative **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} = s_{L,\ib-\half\eb_x} &-& \frac{\dt}{2h}\left[(s^{y|r}v^{\mac})_{\ib-\eb_x+\half\eb_y}-({s^{y|r}v^{\mac})_{\ib-\eb_x-\half\eb_y}}\right] \nonumber \\
     &-& \frac{\dt}{2h}\left[(s^{r|y}w^{\mac})_{\ib-\eb_x+\half\eb_r}-({s^{r|y}w^{\mac})_{\ib-\eb_x-\half\eb_r}}\right] \nonumber \\
     &-& \frac{\dt}{2h}s_{\ib-\eb_x}\left(u_{\ib-\half\eb_x}^{\mac}-u_{\ib-\frac{3}{2}\eb_x}^{\mac}\right) + \frac{\dt}{2}f_{\ib-\eb_x}, \\
     s_{R,\ib-\half\eb_x}^{\edge} = s_{R,\ib-\half\eb_x} &-& \frac{\dt}{2h}\left[(s^{y|r}v^{\mac})_{\ib+\half\eb_y}-({s^{y|r}v^{\mac})_{\ib-\half\eb_y}}\right] \nonumber \\
     &-& \frac{\dt}{2h}\left[(s^{r|y}w^{\mac})_{\ib+\half\eb_r}-({s^{r|y}w^{\mac})_{\ib-\half\eb_r}}\right] \nonumber \\
     &-& \frac{\dt}{2h}s_{\ib}\left(u_{\ib+\half\eb_x}^{\mac}-u_{\ib-\half\eb_x}^{\mac}\right) + \frac{\dt}{2}f_{\ib}.\end{aligned}

  **else**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} = s_{L,\ib-\half\eb_x} &-& \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac}+v_{\ib-\eb_x-\half\eb_y}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_y}^{y|r}-s_{\ib-\eb_x-\half\eb_y}^{y|r}\right) \nonumber \\
     &-& \frac{\dt}{4h}\left(w_{\ib-\eb_x+\half\eb_r}^{\mac}+w_{\ib-\eb_x-\half\eb_r}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_r}^{r|y}-s_{\ib-\eb_x-\half\eb_r}^{r|y}\right) + \frac{\dt}{2}f_{\ib-\eb_x}, \nonumber \\
     && \\
     s_{R,\ib-\half\eb_x}^{\edge} = s_{R,\ib-\half\eb_x} &-& \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{\mac}+v_{\ib-\half\eb_y}^{\mac}\right)\left(s_{\ib+\half\eb_y}^{y|r}-s_{\ib-\half\eb_y}^{y|r}\right) \nonumber \\
     &-& \frac{\dt}{4h}\left(w_{\ib+\half\eb_r}^{\mac}+w_{\ib-\half\eb_r}^{\mac}\right)\left(s_{\ib+\half\eb_r}^{r|y}-s_{\ib-\half\eb_r}^{r|y}\right) + \frac{\dt}{2}f_{\ib}.\end{aligned}

  **end if**

* Account for the :math:`\partial w_0/\partial r` term:

  **if** is_vel **and** comp = 2 **then**

  .. math::

     \begin{aligned}
     s_{L,\ib-\half\eb_x}^{\edge} &=& s_{L,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib-\eb_x+\half\eb_r} + \wt^{\mac}_{\ib-\eb_x-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\
     s_{R,\ib-\half\eb_x}^{\edge} &=& s_{R,\ib-\half\eb_x}^{\edge} -
     \frac{\dt}{4h}\left(\wt^{\mac}_{\ib+\half\eb_r} + \wt^{\mac}_{\ib-\half\eb_r}\right)\left(w_{0,\ib+\half\eb_r}-w_{0,\ib-\half\eb_r}\right) \\\end{aligned}

  **end if**

* Upwind based on :math:`u^{\mac}`:

  .. math::

     s_{\ib-\half\eb_x}^{\edge} =
     \begin{cases}
     \half\left(s_{L,\ib-\half\eb_x}^{\edge} + s_{R,\ib-\half\eb_x}^{\edge}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
     s_{L,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
     s_{R,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} < 0.
     \end{cases}

* Predict :math:`s` to y-faces using a full-dimensional extrapolation.
* Predict :math:`s` to r-faces using a full-dimensional extrapolation.


.. _d-spherical-case-2:

3D Spherical Case
-----------------

The spherical case is the same as the plane-parallel 3D Cartesian
case, except the :math:`\partial w_0/\partial r` term enters in the full
dimensional extrapolation for each direction when predicting velocity
to faces. As in the plane-parallel case, make sure upwind based on
the full velocity.

Computing :math:`\Ub^{\mac,*}` in VARDEN
========================================

.. _d-cartesian-case-6:

2D Cartesian Case
-----------------

* We do a 1D Taylor series extrapolation to get both components of velocity at the x-face:

  .. math::
     u_{L,\ib-\half\eb_x}^{1D} = u_{\ib-\eb_x} + \left[\half - \frac{\dt}{2h}{\rm max}(0,u_{\ib-\eb_x})\right]\Delta_xu_{\ib-\eb_x},
     :label:  varden U_L^1D

  .. math::
     u_{R,\ib-\half\eb_x}^{1D} = u_{\ib} + \left[\half - \frac{\dt}{2h}{\rm min}(0,u_{\ib})\right]\Delta_xu_{\ib} .

  .. math::
     v_{L,\ib-\half\eb_x}^{1D} = v_{\ib-\eb_x} + \left[\half - \frac{\dt}{2h}{\rm max}(0,v_{\ib-\eb_x})\right]\Delta_xv_{\ib-\eb_x},

  .. math::
     v_{R,\ib-\half\eb_x}^{1D} &=& v_{\ib} + \left[\half - \frac{\dt}{2h}{\rm min}(0,v_{\ib})\right]\Delta_xv_{\ib}.

  We obtain the normal velocity using the Riemann problem:

  .. math::

     u_{\ib-\half\eb_x}^{1D} =
     \begin{cases}
     0, & \left(u_{L,\ib-\half\eb_x}^{1D} \le 0 ~~ {\rm AND} ~~ u_{R,\ib-\half\eb_x}^{1D} \ge 0\right) ~~ {\rm OR} ~~ \left|u_{L,\ib-\half\eb_x}^{1D} + u_{R,\ib-\half\eb_x}^{1D}\right| < \epsilon, \\
     u_{L,\ib-\half\eb_x}^{1D}, & u_{L,\ib-\half\eb_x}^{1D} + u_{R,\ib-\half\eb_x}^{1D} > 0, \\
     u_{R,\ib-\half\eb_x}^{1D}, & u_{L,\ib-\half\eb_x}^{1D} + u_{R,\ib-\half\eb_x}^{1D} < 0.
     \end{cases}

  We obtain the transverse velocity by upwinding based on
  :math:`u_{\ib-\half\eb_x}^{1D}`:

  .. math::
     v_{\ib-\half\eb_x}^{1D} =
     \begin{cases}
     \half\left(v_{L,\ib-\half\eb_x}^{1D} + v_{R,\ib-\half\eb_x}^{1D}\right), & \left|u_{\ib-\half\eb_x}^{1D}\right| < \epsilon \\
     v_{L,\ib-\half\eb_x}^{1D}, & u_{\ib-\half\eb_x}^{1D} > 0, \\
     v_{R,\ib-\half\eb_x}^{1D}, & u_{\ib-\half\eb_x}^{1D} < 0.
     \end{cases}
     :label: Transverse Velocity Riemann Problem

* We perform analogous operations to compute both components of velocity
  at the y-faces, :math:`\Ub_{\ib-\half\eb_y}^{1D}`.

* Now we do a full-dimensional extrapolation to get the MAC velocity at
  the x-faces (note that we only compute the normal components):

  .. math::

     \begin{aligned}
     u_{L,\ib-\half\eb_x}^{\mac,*} &=& u_{L,\ib-\half\eb_x}^{1D} - \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{1D}+v_{\ib-\eb_x-\half\eb_y}^{1D}\right)\left(u_{\ib-\eb_x+\half\eb_y}^{1D} - u_{\ib-\eb_x-\half\eb_y}^{1D}\right) + \frac{\dt}{2}f_{u,\ib-\eb_x}, \\
     u_{R,\ib-\half\eb_x}^{\mac,*} &=& u_{R,\ib-\half\eb_x}^{1D} - \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{1D}+v_{\ib-\half\eb_y}^{1D}\right)\left(u_{\ib+\half\eb_y}^{1D} - u_{\ib-\half\eb_y}^{1D}\right) + \frac{\dt}{2}f_{u,\ib}.\end{aligned}

  Then we solve a Riemann problem:

  .. math::
     u_{\ib-\half\eb_x}^{\mac,*} =
     \begin{cases}
     0, & \left(u_{L,\ib-\half\eb_x}^{\mac,*} \le 0 ~~ {\rm AND} ~~ u_{R,\ib-\half\eb_x}^{\mac,*} \ge 0\right) ~~ {\rm OR} ~~ \left|u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*}\right| < \epsilon, \\
     u_{L,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} > 0, \\
     u_{R,\ib-\half\eb_x}^{\mac,*}, & u_{L,\ib-\half\eb_x}^{\mac,*} + u_{R,\ib-\half\eb_x}^{\mac,*} < 0.
     \end{cases}
     :label: umac Riemann Problem

* We perform analogous operations to compute the normal velocity at the
  y-faces, :math:`v^{\mac,*}_{\ib-\half\eb_y}`.

.. _d-cartesian-case-7:

3D Cartesian Case
-----------------

This is more complicated than the 2D case because we include corner
coupling. We compute :math:`\Ub_{\ib-\half\eb_x}^{1D},
\Ub_{\ib-\half\eb_y}^{1D}`, and :math:`\Ub_{\ib-\half\eb_z}^{1D}` in
an analogous manner as :eq:`varden U_L^1D]`-:eq:`Transverse Velocity
Riemann Problem`. Then we compute an intermediate state,
:math:`u_{\ib-\half\eb_y}^{y|z}`, which is described as “state
:math:`u_{\ib-\half\eb_y}^{1D}` that has been updated to account for
the transverse derivatives in the z direction”, using:

.. math::

   \begin{aligned}
   u_{L,\ib-\half\eb_y}^{y|z} &=& u_{L,\ib-\half\eb_y}^{1D} - \frac{\dt}{6h}\left(w_{\ib-\eb_y+\half\eb_z}^{1D}+w_{\ib-\eb_y-\half\eb_z}^{1D}\right)\left(u_{\ib-\eb_y+\half\eb_z}^{1D}-u_{\ib-\eb_y-\half\eb_z}^{1D}\right), \\
   u_{R,\ib-\half\eb_y}^{y|z} &=& u_{R,\ib-\half\eb_y}^{1D} - \frac{\dt}{6h}\left(w_{\ib+\half\eb_z}^{1D}+w_{\ib-\half\eb_z}^{1D}\right)\left(u_{\ib+\half\eb_z}^{1D}-u_{\ib-\half\eb_z}^{1D}\right).\end{aligned}

Then upwind based on :math:`v_{\ib-\half\eb_y}^{1D}`:

.. math::

   u_{\ib-\half\eb_y}^{y|z} =
   \begin{cases}
   \half\left(u_{L,\ib-\half\eb_y}^{y|z} + u_{R,\ib-\half\eb_y}^{y|z}\right), & \left|v_{\ib-\half\eb_y}^{1D}\right| < \epsilon \\
   u_{L,\ib-\half\eb_y}^{y|z}, & v_{\ib-\half\eb_y}^{1D} > 0, \\
   u_{R,\ib-\half\eb_y}^{y|z}, & v_{\ib-\half\eb_y}^{1D} < 0.
   \end{cases}

We use an analogous procedure to compute five more intermediate states,
:math:`u_{\ib-\half\eb_z}^{z|y}, v_{\ib-\half\eb_x}^{x|z},
v_{\ib-\half\eb_z}^{z|x}, w_{\ib-\half\eb_x}^{x|y}`, and
:math:`w_{\ib-\half\eb_y}^{y|x}`. Then we do a full-dimensional
extrapolation to get the MAC velocities at normal faces:

.. math::

   \begin{aligned}
   u_{L,\ib-\half\eb_x}^{\mac,*} = u_{L,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{1D}+v_{\ib-\eb_x-\half\eb_y}^{1D}\right)\left(u_{\ib-\eb_x+\half\eb_y}^{y|z}-u_{\ib-\eb_x-\half\eb_y}^{y|z}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(w_{\ib-\eb_x+\half\eb_z}^{1D}+w_{\ib-\eb_x-\half\eb_z}^{1D}\right)\left(u_{\ib-\eb_x+\half\eb_z}^{z|y}-u_{\ib-\eb_x-\half\eb_z}^{z|y}\right) + \frac{\dt}{2}f_{u,\ib-\eb_x}, \\
   u_{R,\ib-\half\eb_x}^{\mac,*} = u_{R,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{1D}+v_{\ib-\half\eb_y}^{1D}\right)\left(u_{\ib+\half\eb_y}^{y|z}-u_{\ib-\half\eb_y}^{y|z}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(w_{\ib+\half\eb_z}^{1D}+w_{\ib-\half\eb_z}^{1D}\right)\left(u_{\ib+\half\eb_z}^{z|y}-u_{\ib-\half\eb_z}^{z|y}\right) + \frac{\dt}{2}f_{u,\ib}.\end{aligned}

Then we use the Riemann solver given above for the 2D case (:eq:`umac Riemann Problem`) to compute
:math:`u_{\ib-\half\eb_x}^{\mac,*}`. We use an analogous procedure to
obtain :math:`v_{\ib-\half\eb_y}^{\mac,*}` and
:math:`w_{\ib-\half\eb_z}^{\mac,*}`.


Computing :math:`\Ub^{\edge}` and :math:`\rho^{\edge}` in VARDEN
================================================================

To compute :math:`\Ub^{\edge}`, VARDEN uses the exact same algorithm
as the :math:`s^{\edge}` case in MAESTROeX. The algorithm for
:math:`\rho^{\edge}` in VARDEN is slightly different than the
:math:`s^{\edge}` case in MAESTROeX since it uses a “conservative”
formulation. Here, :math:`s` is used in place of either :math:`\rho, u, v`, or
:math:`w` (in 3D).

.. _d-cartesian-case-8:

2D Cartesian Case
-----------------

The 1D extrapolation is:

.. math::
   s_{L,\ib-\half\eb_x}^{1D} = s_{\ib-\eb_x}^n + \left(\half - \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib-\eb_x}^n,
   :label: varden s_L^1D

.. math::
   s_{R,\ib-\half\eb_x}^{1D} = s_{\ib} - \left(\half + \frac{\dt}{2h}u_{\ib-\half\eb_x}^{\mac}\right)\Delta_x s_{\ib}^n.
   :label: varden s_R^1D

Then we upwind based on :math:`u^{\mac}`:

.. math::

   s_{\ib-\half\eb_x}^{1D} =
   \begin{cases}
   \half\left(s_{L,\ib-\half\eb_x}^{1D} + s_{R,\ib-\half\eb_x}^{1D}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
   s_{L,\ib-\half\eb_x}^{1D}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
   s_{R,\ib-\half\eb_x}^{1D}, & u^{\mac}_{\ib-\half\eb_x} < 0. \\
   \end{cases}

We use an analogous procedure to obtain :math:`s_{\ib-\half\eb_y}^{1D}`.
Now we do a full-dimensional extrapolation of :math:`s` to each face. The
extrapolation of a “non-conserved” :math:`s` to x-faces is:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{\edge} &=& s_{L,\ib-\half\eb_x}^{1D} - \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac}+v_{\ib-\eb_x-\half\eb_y}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_y}^{1D} - s_{\ib-\eb_x-\half\eb_y}^{1D}\right) + \frac{\dt}{2}f_{s,\ib-\eb_x}, \\
   s_{R,\ib-\half\eb_x}^{\edge} &=& s_{R,\ib-\half\eb_x}^{1D} - \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{\mac}+v_{\ib-\half\eb_y}^{\mac}\right)\left(s_{\ib+\half\eb_y}^{1D} - s_{\ib-\half\eb_y}^{1D}\right) + \frac{\dt}{2}f_{s,\ib}.\end{aligned}

The extrapolation of a “conserved” :math:`s` to x-faces is:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{\edge} = s_{L,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{2h}\left[(s^{1D} v^{\mac})_{\ib-\eb_x+\half\eb_y} - (s^{1D} v^{\mac})_{\ib-\eb_x-\half\eb_y}\right] \nonumber \\
   &-& \frac{\dt}{2}s_{\ib-\eb_x}(\nabla\cdot\Ub^{\mac})_{\ib-\eb_x} + \frac{\dt}{2h}s_{\ib-\eb_x}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac} - v_{\ib-\eb_x-\half\eb_y}^{\mac}\right) + \frac{\dt}{2}f_{s,\ib-\eb_x}, \\
   s_{R,\ib-\half\eb_x}^{\edge} = s_{R,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{2h}\left[(s^{1D} v^{\mac})_{\ib+\half\eb_y} - (s^{1D} v^{\mac})_{\ib-\half\eb_y}\right] \nonumber \\
   &-& \frac{\dt}{2}s_{\ib}(\nabla\cdot\Ub^{\mac})_{\ib} + \frac{\dt}{2h}s_{\ib}\left(v_{\ib+\half\eb_y}^{\mac} - v_{\ib-\half\eb_y}^{\mac}\right) + \frac{\dt}{2}f_{s,\ib}.\end{aligned}

Then we upwind based on :math:`u^{\mac}`.

.. math::
   s_{\ib-\half\eb_x}^{\edge} =
   \begin{cases}
   \half\left(s_{L,\ib-\half\eb_x}^{\edge} + s_{R,\ib-\half\eb_x}^{\edge}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
   s_{L,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
   s_{R,\ib-\half\eb_x}^{\edge}, & u^{\mac}_{\ib-\half\eb_x} < 0.
   \end{cases}
   :label: varden s^edge upwind

We use an analogous procedure to compute :math:`s_{\ib-\half\eb_y}^{\edge}`.

.. _d-cartesian-case-9:

3D Cartesian Case
-----------------

This is more complicated than the 2D case because we include corner
coupling. We first compute :math:`s_{\ib-\half\eb_x}^{1D}`,
:math:`s_{\ib-\half\eb_y}^{1D}`, and :math:`s_{\ib-\half\eb_z}^{1D}` in an
analogous manner to :eq:`varden s_L^1D` and :eq:`varden s_R^1D`. Then we compute six intermediate states,
:math:`s_{\ib-\half\eb_x}^{x|y}, s_{\ib-\half\eb_x}^{x|z},
s_{\ib-\half\eb_y}^{y|x}, s_{\ib-\half\eb_y}^{y|z},
s_{\ib-\half\eb_z}^{z|x}`, and :math:`s_{\ib-\half\eb_z}^{z|y}`. For the
“non-conservative case”, we use, for example:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{x|y} &=& s_{L,\ib-\half\eb_x}^{1D} - \frac{\dt}{6h}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac} + v_{\ib-\eb_x-\half\eb_y}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_y}^{1D} - s_{\ib-\eb_x-\half\eb_y}^{1D}\right), \\
   s_{R,\ib-\half\eb_x}^{x|y} &=& s_{R,\ib-\half\eb_x}^{1D} - \frac{\dt}{6h}\left(v_{\ib+\half\eb_y}^{\mac} + v_{\ib-\half\eb_y}^{\mac}\right)\left(s_{\ib+\half\eb_y}^{1D} - s_{\ib-\half\eb_y}^{1D}\right).\end{aligned}

For the “conservative” case, we use, for example:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{x|y} &=& s_{L,\ib-\half\eb_x}^{1D} - \frac{\dt}{3h}\left[(sv^{\mac})_{\ib-\eb_x+\half\eb_y}-(sv^{\mac})_{\ib-\eb_x-\half\eb_y}\right], \\
   s_{R,\ib-\half\eb_x}^{x|y} &=& s_{R,\ib-\half\eb_x}^{1D} - \frac{\dt}{3h}\left[(sv^{\mac})_{\ib+\half\eb_y}-(sv^{\mac})_{\ib-\half\eb_y}\right].\end{aligned}

Then we upwind based on :math:`u^{\mac}`:

.. math::

   s_{\ib-\half\eb_x}^{x|y} =
   \begin{cases}
   \half\left(s_{L,\ib-\half\eb_x}^{x|y} + s_{R,\ib-\half\eb_x}^{x|y}\right), & \left|u^{\mac}_{\ib-\half\eb_x}\right| < \epsilon \\
   s_{L,\ib-\half\eb_x}^{x|y}, & u^{\mac}_{\ib-\half\eb_x} > 0, \\
   s_{R,\ib-\half\eb_x}^{x|y}, & u^{\mac}_{\ib-\half\eb_x} < 0.
   \end{cases}

We use an analogous procedure to compute the other five intermediate
states. Now we do a full-dimensional extrapolation of :math:`s` to each
face. The extrapolation of a “non-conserved” :math:`s` to x-faces is:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{\edge} = s_{L,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{4h}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac}+v_{\ib-\eb_x-\half\eb_y}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_y}^{y|z}-s_{\ib-\eb_x-\half\eb_y}^{y|z}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(w_{\ib-\eb_x+\half\eb_z}^{\mac}+w_{\ib-\eb_x-\half\eb_z}^{\mac}\right)\left(s_{\ib-\eb_x+\half\eb_z}^{z|y}-s_{\ib-\eb_x-\half\eb_z}^{z|y}\right) \nonumber \\
   &+& \frac{\dt}{2}f_{s,\ib-\eb_x}, \\
   s_{R,\ib-\half\eb_x}^{\edge} = s_{R,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{4h}\left(v_{\ib+\half\eb_y}^{\mac}+v_{\ib-\half\eb_y}^{\mac}\right)\left(s_{\ib+\half\eb_y}^{y|z}-s_{\ib-\half\eb_y}^{y|z}\right) \nonumber \\
   &-& \frac{\dt}{4h}\left(w_{\ib+\half\eb_z}^{\mac}+w_{\ib-\half\eb_z}^{\mac}\right)\left(s_{\ib+\half\eb_z}^{z|y}-s_{\ib-\half\eb_z}^{z|y}\right) \nonumber \\
   &+& \frac{\dt}{2}f_{s,\ib}.\end{aligned}

The extrapolation of a “conserved” :math:`s` to x-faces is:

.. math::

   \begin{aligned}
   s_{L,\ib-\half\eb_x}^{\edge} = s_{L,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{2h}\left[(s^{y|z}v^{\mac})_{\ib-\eb_x+\half\eb_y}-({s^{y|z}v^{\mac})_{\ib-\eb_x-\half\eb_y}}\right] \nonumber \\
   &-& \frac{\dt}{2h}\left[(s^{z|y}w^{\mac})_{\ib-\eb_x+\half\eb_z}-({s^{z|y}w^{\mac})_{\ib-\eb_x-\half\eb_z}}\right] \nonumber \\
   &-& \frac{\dt}{2}s_{\ib-\eb_x}(\nabla\cdot\Ub^{\mac})_{\ib-\eb_x} \nonumber \\
   &+& \frac{\dt}{2h}s_{\ib-\eb_x}\left(v_{\ib-\eb_x+\half\eb_y}^{\mac}-v_{\ib-\eb_x-\half\eb_y}^{\mac}+w_{\ib-\eb_x+\half\eb_z}^{\mac}-w_{\ib-\eb_x-\half\eb_z}^{\mac}\right) \nonumber \\
   &+& \frac{\dt}{2}f_{s,\ib-\eb_x}, \\
   s_{R,\ib-\half\eb_x}^{\edge} = s_{R,\ib-\half\eb_x}^{1D} &-& \frac{\dt}{2h}\left[(s^{y|z}v^{\mac})_{\ib+\half\eb_y}-({s^{y|z}v^{\mac})_{\ib-\half\eb_y}}\right] \nonumber \\
   &-& \frac{\dt}{2h}\left[(s^{z|y}w^{\mac})_{\ib+\half\eb_z}-({s^{z|y}w^{\mac})_{\ib-\half\eb_z}}\right] \nonumber \\
   &-& \frac{\dt}{2}s_{\ib}(\nabla\cdot\Ub^{\mac})_{\ib} \nonumber \\
   &+& \frac{\dt}{2h}s_{\ib}\left(v_{\ib+\half\eb_y}^{\mac}-v_{\ib-\half\eb_y}^{\mac}+w_{\ib+\half\eb_z}^{\mac}-w_{\ib-\half\eb_z}^{\mac}\right) \nonumber \\
   &+& \frac{\dt}{2}f_{s,\ib}.\end{aligned}

Then we upwind based on :math:`u^{\mac}`, as in :eq:`varden s^edge upwind`.
We use an analogous procedure to compute both
:math:`s_{\ib-\half\eb_y}^{\edge}` and :math:`s_{\ib-\half\eb_z}`.


ESTATE_FPU in GODUNOV_2D/3D.f
=============================

* First, the normal predictor.

  .. math::

     \begin{aligned}
     s_L^x &=& s_{\ib-\eb_x} + \left(\half - \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib-\eb_x} + \underbrace{\frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x}}_{\text{IF USE\_MINION}} \\
     s_R^x &=& s_{\ib} - \left(\half + \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \underbrace{\frac{\dt}{2}\text{TFORCES}_{\ib}}_{\text{IF USE\_MINION}}\end{aligned}

  **If** USE_MINION **and** ICONSERVE **then:**

  .. math::

     \begin{aligned}
     s_L^x &=& s_L^x - \frac{\dt}{2}s_{\ib-\eb_x}\text{DIVU}_{\ib-\eb_x} \\
     s_R^x &=& s_R^x - \frac{\dt}{2}s_{\ib}\text{DIVU}_{\ib}\end{aligned}

  Apply boundary conditions on :math:`s_L^x` and :math:`s_R^x`. Then,

  .. math::
     \text{s}_{\ib-\half\eb_x}^x =
     \begin{cases}
     s_L^x, & \text{UEDGE}_{\ib-\half\eb_x} > 0, \\
     s_R^x, & \text{else}. \\
     \end{cases}
     :label: ESTATE_FPU Upwind

* Then, if :math:`|\text{UEDGE}_{\ib-\half\eb_x}| \le \epsilon`, we set :math:`s_{\ib-\half\eb_x}^x = (s_L^x+s_R^x)/2`. The procedure to obtain :math:`s_{\ib-\half\eb_y}^y` is analogous.

* Now, the transverse terms.

  **If** ICONSERVE **then:**

  .. math::

     \begin{aligned}
     \text{sedge}_L^x &=& s_{\ib-\eb_x} + \left(\half - \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib-\eb_x} + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{\text{VEDGE}_{\ib-\eb_x+\half\eb_y}s_{\ib-\eb_x+\half\eb_y}^y - \text{VEDGE}_{\ib-\eb_x-\half\eb_y}s_{\ib-\eb_x-\half\eb_y}^y}{h_y}\right.\nonumber\\
     && ~~~~~~~~~~ \left. - \frac{s_{\ib-\eb_x}(\text{VEDGE}_{\ib-\eb_x+\half\eb_y}-\text{VEDGE}_{\ib-\eb_x-\half\eb_y})}{h_y}+s_{\ib-\eb_x}\text{DIVU}_{\ib-\eb_x}\right]\\
     \text{sedge}_R^x &=& s_{\ib} - \left(\half + \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{\text{VEDGE}_{\ib+\half\eb_y}s_{\ib+\half\eb_y}^y - \text{VEDGE}_{\ib-\half\eb_y}s_{\ib-\half\eb_y}^y}{h_y}\right.\nonumber\\
     && ~~~~~~~~~~ \left. - \frac{s_{\ib}(\text{VEDGE}_{\ib+\half\eb_y}-\text{VEDGE}_{\ib-\half\eb_y})}{h_y}+s_{\ib}\text{DIVU}_{\ib}\right]\end{aligned}

* Now, define :math:`\text{VBAR}_{\ib} = (\text{VEDGE}_{\ib+\half\eb_y}+\text{VEDGE}_{\ib-\half\eb_y})/2`.

  **If** NOT ICONSERVE **and** :math:`\text{VEDGE}_{\ib+\half\eb_y}\cdot\text{VEDGE}_{\ib-\half\eb_y} < 0` **and** :math:`\text{VBAR}_{\ib} < 0` **then:**

  .. math::
     \begin{align}
     \text{sedge}_L^x = s_{\ib-\eb_x} &+ \left(\half - \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \nonumber\\
     & - \frac{\dt}{2}\left[\frac{\text{VBAR}_{\ib-\eb_x}(s_{\ib-\eb_x+\eb_y}-s_{\ib-\eb_x})}{h_y}\right]
     \end{align}
     :label: transverse upwinding 1

  .. math::
     \begin{align}
     \text{sedge}_R^x = s_{\ib} &- \left(\half + \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib} \nonumber\\
     & - \frac{\dt}{2}\left[\frac{\text{VBAR}_{\ib}(s_{\ib+\eb_y}-s_{\ib})}{h_y}\right]
     \end{align}

  **Else If** NOT ICONSERVE **and** :math:`\text{VEDGE}_{\ib+\half\eb_y}\cdot\text{VEDGE}_{\ib-\half\eb_y} < 0` **and** :math:`\text{VBAR}_{\ib} \ge 0` **then:**

  .. math::

     \begin{aligned}
     \text{sedge}_L^x = s_{\ib-\eb_x} &+& \left(\half - \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib-\eb_x} + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{\text{VBAR}_{\ib-\eb_x}(s_{\ib-\eb_x}-s_{\ib-\eb_x-\eb_y})}{h_y}\right] \\
     \text{sedge}_R^x = s_{\ib} &-& \left(\half + \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{\text{VBAR}_{\ib}(s_{\ib}-s_{\ib-\eb_y})}{h_y}\right]\end{aligned}

  **Else If** NOT ICONSERVE **and** :math:`\text{VEDGE}_{\ib+\half\eb_y}\cdot\text{VEDGE}_{\ib-\half\eb_y} \ge 0` **then:**

  .. math::
     \begin{align}
     \text{sedge}_L^x &= s_{\ib-\eb_x} + \left(\half - \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib-\eb_x} + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{(\text{VEDGE}_{\ib-\eb_x+\half\eb_y}+\text{VEDGE}_{\ib-\eb_x-\half\eb_y})(s_{\ib-\eb_x+\half\eb_y}-s_{\ib-\eb_x-\half\eb_y})}{2h_y}\right] \\
     \text{sedge}_R^x &=& s_{\ib} - \left(\half + \frac{\dt}{h_x}\text{UEDGE}_{\ib-\half\eb_x}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib} \nonumber\\
     && - \frac{\dt}{2}\left[\frac{(\text{VEDGE}_{\ib+\half\eb_y}+\text{VEDGE}_{\ib-\half\eb_y})(s_{\ib+\half\eb_y}-s_{\ib-\half\eb_y})}{2h_y}\right]
     \end{align}
     :label: transverse upwinding 6

*  Finally, upwind analogous to :eq:`ESTATE_FPU Upwind` to get :math:`\text{sedge}_{\ib-\half\eb_x}`.


ESTATE in GODUNOV_2D/3D.f
=========================

First, the normal predictor.

.. math::

   \begin{aligned}
   s_L^x &=& s_{\ib-\eb_x} + \left(\half - \frac{\dt}{h_x}u_{\ib-\eb_x}\right)\Delta^x s_{\ib-\eb_x} \\
   s_R^x &=& s_{\ib} - \left(\half + \frac{\dt}{h_x}u_{\ib}\right)\Delta^x s_{\ib}\end{aligned}

**If** USE_MINION **then:**

.. math::

   \begin{aligned}
   s_L^x &=& s_L^x + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \\
   s_R^x &=& s_R^x + \frac{\dt}{2}\text{TFORCES}_{\ib}\end{aligned}

Apply boundary conditions on :math:`s_L^x` and :math:`s_R^x`. Then,

.. math::
   \text{s}_{\ib-\half\eb_x}^x =
   \begin{cases}
   s_L^x, & \text{UAD}_{\ib-\half\eb_x} > 0, \\
   s_R^x, & \text{else}. \\
   \end{cases}
   :label: ESTATE Upwind

Then, if :math:`|\text{UAD}_{\ib-\half\eb_x}| \le \epsilon`, we set :math:`s_{\ib-\half\eb_x}^x = (s_L^x+s_R^x)/2`.

.. math::

   \begin{aligned}
   \text{sedge}_L^x = s_{\ib-\eb_x} &+& \left(\half - \frac{\dt}{h_x}u_{\ib-\eb_x}\right)\Delta^x s_{\ib-\eb_x} + \frac{\dt}{2}\text{TFORCES}_{\ib-\eb_x} \nonumber\\
   && - \frac{\dt}{2}\left[\frac{(\text{VAD}_{\ib-\eb_x+\half\eb_y}+\text{VAD}_{\ib-\eb_x-\half\eb_y})(s_{\ib-\eb_x+\half\eb_y}-s_{\ib-\eb_x-\half\eb_y})}{2h_y}\right] \\
   \text{sedge}_R^x = s_{\ib} &-& \left(\half + \frac{\dt}{h_x}u_{\ib}\right)\Delta^x s_{\ib} + \frac{\dt}{2}\text{TFORCES}_{\ib} \nonumber\\
   && - \frac{\dt}{2}\left[\frac{(\text{VAD}_{\ib+\half\eb_y}+\text{VAD}_{\ib-\half\eb_y})(s_{\ib+\half\eb_y}-s_{\ib-\half\eb_y})}{2h_y}\right]\end{aligned}

Note that the 2D and 3D algorithms are different - in 3D the transverse
terms use upwinding analogous to :eq:`transverse upwinding 1`-:eq:`transverse upwinding 6`, using UAD
instead of UEDGE. Finally, upwind analogous to :eq:`ESTATE Upwind`
to get :math:`\text{sedge}_{\ib-\half\eb_x}`, but use UEDGE instead of UAD.


Piecewise Parabolic Method (PPM)
================================

Consider a scalar, :math:`s`, which we wish to predict to
time-centered edges.  The PPM method is an improvement over the
piecewise-linear method.  Using our notation, we modify equations
(:eq:`3D predict s to left` and :eq:`3D predict s to right`) in Section `4 <#Scalar Edge
State Prediction in MAESTROeX>`__ to obtain better estimates for the
time-centered 1D edge states, :math:`s_{L/R,\ib-\myhalf\eb_x}`,
etc.. Once these states are obtained, we continue with the
full-dimensional extrapolations as described before.

The PPM method is described in a series of papers:

-  Colella and Woodward 1984 - describes the basic method.

-  Miller and Colella 2002 - describes how to apply PPM to a multidimensional
   unsplit Godunov method and generalizes the characteristic tracing for more complicated
   systems. Note that we only need to upwind based on the fluid velocity, so we don’t
   need to use fancy characteristic tracing.

-  Colella and Sekora 2008 - describes new fancy quadratic limiters. There are
   several errors in the printed text, so we have implemented a corrected version from
   Phil Colella.

Here are the steps for the :math:`x`-direction. For simplicity, we replace the vector index notation with a simple scalar notation (:math:`\ib+\eb_x \rightarrow i+1`, etc.).

-  **Step 1:** Compute :math:`s_{i,+}` and :math:`s_{i,-}`, which are spatial interpolations of
   :math:`s` to the hi and lo faces of cell :math:`i`, respectively. See Sections
   `9.1 <#Sec:ColellaWoodward>`__ and `9.2 <#Sec:ColellaSekora>`__ for the two options.

-  **Step 2:** Construct a quadratic profile within each cell.

   .. math:: s_i^I(x) = s_{i,-} + \xi\left[s_{i,+} - s_{i,-} + s_{6,i}(1-\xi)\right],\label{Quadratic Interp}

   .. math:: s_{6,i}= 6s_{i} - 3\left(s_{i,-}+s_{i,+}\right),

   .. math:: \xi = \frac{x - (i-\myhalf)h}{h}, ~ 0 \le \xi \le 1.

-  **Step 3:** Integrate quadratic profiles to get the average value swept over the face
   over time.
   Define the following integrals, where :math:`\sigma = |u|\Delta t/h`:

   .. math::

        \begin{aligned}
        \mathcal{I}_{i,+}(\sigma) &=& \frac{1}{\sigma h}\int_{(i+\myhalf)h-\sigma h}^{(i+\myhalf)h}s_i^I(x)dx \\
        \mathcal{I}_{i,-}(\sigma) &=& \frac{1}{\sigma h}\int_{(i-\myhalf)h}^{(i-\myhalf)h+\sigma h}s_i^I(x)dx\end{aligned}

   Plugging in (`[Quadratic Interp] <#Quadratic Interp>`__) gives:

   .. math::

        \begin{aligned}
        \mathcal{I}_{i,+}(\sigma) &=& s_{j,+} - \frac{\sigma}{2}\left[s_{j,+}-s_{j,-}-\left(1-\frac{2}{3}\sigma\right)s_{6,i}\right], \\
        \mathcal{I}_{i,-}(\sigma) &=& s_{j,-} + \frac{\sigma}{2}\left[s_{j,+}-s_{j,-}+\left(1-\frac{2}{3}\sigma\right)s_{6,i}\right].\end{aligned}

-  **Step 4:** Obtain 1D edge states.
   Perform a 1D extrapolation, without source terms, to get
   left and right edge states. Add the source terms later if desired/necessary.

   .. math::

        \begin{aligned}
        s_{L,i-\myhalf} &=&
        \begin{cases}
        \mathcal{I}_{i-1,+}(\sigma), & u_{i-1} ~ \text{or} ~ u_{i-\myhalf}^{\mac} > 0 \\
        s_{i-1}, & \text{else}.
        \end{cases}\\
        s_{R,i-\myhalf} &=&
        \begin{cases}
        \mathcal{I}_{i,-}(\sigma), & u_{i} ~ \text{or} ~ u_{i-\myhalf}^{\mac} < 0 \\
        s_{i}, & \text{else}.
        \end{cases}\end{aligned}

.. _Sec:ColellaWoodward:

Colella and Woodward Based Approach
-----------------------------------

Spatially interpolate :math:`s` to edges.
Use a 4th-order interpolation in space with van Leer limiting to obtain edge values:

.. math::
   s_{i+\myhalf} = \frac{1}{2}\left(s_{i} + s_{i+1}\right) - \frac{1}{6}\left(\delta s_{i+1}^{vL} - \delta s_{i}^{vL}\right),
   :label: eq:CW Edge

.. math::

   \delta s_i =
   \frac{1}{2}\left(s_{i+1}-s_{i-1}\right),

.. math::

   \delta s_i^{vL} =
   \begin{cases}
   \text{sign}(\delta s_i)\min\left(|\delta s_i|, ~ 2|s_{i+1}-s_{i}|, ~ 2|s_i-s_{i-1}|\right), & {\rm if} ~ (s_{i+1}-s_i)(s_i-s_{i-1}) > 0,\\
   0, & {\rm otherwise}.
   \end{cases}

A more compact way of writing this is

.. math:: s = \text{sign}(\delta s_i),

.. math:: \delta s_i^{vL} = s\max\left\{\min\left[s\delta s_i, 2s(s_{i+1}-s_i),2s(s_i-s_{i-1})\right],0\right\}

Without the limiters, :eq:`eq:CW Edge` is the familiar 4th-order spatial interpolation formula:

.. math:: s_{i+\myhalf} = \frac{7}{12}\left(s_{i+1}+s_i\right) - \frac{1}{12}\left(s_{i+2}+s_{i-1}\right).

Next, we must ensure that :math:`s_{i+\myhalf}` lies between the adjacent
cell-centered values:

.. math:: \min\left(s_{i},s_{i+1}\right) \le s_{i+\myhalf} \le \max\left(s_{i},s_{i+1}\right).

In anticipation of further limiting, we set double-valued face-centered values:

.. math:: s_{i,+} = s_{i+1,-} = s_{i+\myhalf}.

Modify :math:`s_{i,\pm}` using a quadratic limiter.
First, we test whether
:math:`s_i` is a local extreumum with the condition:

.. math:: \left(s_{i,+}-s_i\right)\left(s_i-s_{i,-}\right) \le 0,

If this condition is true, we constrain :math:`s_{i,\pm}` by setting
:math:`s_{i,+} = s_{i,-} = s_i`. If not, we then apply a second test to determine
whether :math:`s_i` is sufficiently close to :math:`s_{i,\pm}` so that a quadratic
interpolate would contain a local extremum. We define
:math:`\alpha_{i,\pm} = s_{i,\pm} - s_i`. If one of :math:`|\alpha_{i,\pm}| \ge 2|\alpha_{i,\mp}|`
holds, then for that choice of :math:`\pm = +,-` we set:

.. math:: s_{i,\pm} = 3s_i - 2s_{i,\mp}.

.. _Sec:ColellaSekora:

Colella and Sekora Based Approach
---------------------------------

* Spatially interpolate :math:`s` to edges.

  Use a 4th-order interpolation in space to obtain edge values:

  .. math:: s_{i+\myhalf} = \frac{7}{12}\left(s_{i+1}+s_i\right) - \frac{1}{12}\left(s_{i+2}+s_{i-1}\right).

  Then, if :math:`(s_{i+\myhalf}-s_i)(s_{i+1}-s_{i+\myhalf}) < 0`, we limit :math:`s_{i+\myhalf}` using
  a nonlinear combination of approximations to the second derivative.
  First, define:

  .. math::

     \begin{aligned}
     (D^2s)_{i+\myhalf} &=& 3\left(s_{i}-2s_{i+\myhalf}+s_{i+1}\right) \\
     (D^2s)_{i+\myhalf,L} &=& s_{i-1}-2s_{i}+s_{i+1} \\
     (D^2s)_{i+\myhalf,R} &=& s_{i}-2s_{i+1}+s_{i+2}\end{aligned}

  Then, define

  .. math:: s = \text{sign}\left[(D^2s)_{i+\myhalf}\right],

  .. math:: (D^2s)_{i+\myhalf,\text{lim}} = s\max\left\{\min\left[Cs(D^2s)_{i+\myhalf,L},Cs(D^2s)_{i+\myhalf,R},s(D^2s)_{i+\myhalf}\right],0\right\},

  where :math:`C=1.25` was used in Colella and Sekora. Then,

  .. math:: s_{i+\myhalf} = \frac{1}{2}\left(s_{i}+s_{i+1}\right) - \frac{1}{6}(D^2s)_{i+\myhalf,\text{lim}}.

  Now we implement Phil’s new version of the algorithm to eliminate sensitivity to roundoff.
  First we need to detect whether a particular cell corresponds to an “extremum”. There
  are two tests. For the first test, define

  .. math:: \alpha_{i,\pm} = s_{i\pm\myhalf} - s_i.

  If :math:`\alpha_{i,+}\alpha_{i,-} \ge 0`, then we are at an extremum. We apply the second
  test if either :math:`|\alpha_{i,\pm}| > 2|\alpha_{i,\mp}|`. Then, we define:

  .. math::

     \begin{aligned}
     (Ds)_{i,{\rm face},-} &=& s_{i-\myhalf} - s_{i-\sfrac{3}{2}} \\
     (Ds)_{i,{\rm face},+} &=& s_{i+\sfrac{3}{2}} - s_{i+\myhalf}\end{aligned}

  .. math:: (Ds)_{i,{\rm face,min}} = \min\left[\left|(Ds)_{i,{\rm face},-}\right|,\left|(Ds)_{i,{\rm face},+}\right|\right].

  .. math::

     \begin{aligned}
     (Ds)_{i,{\rm cc},-} &=& s_{i} - s_{i-1} \\
     (Ds)_{i,{\rm cc},+} &=& s_{i+1} - s_{i}\end{aligned}

  .. math:: (Ds)_{i,{\rm cc,min}} = \min\left[\left|(Ds)_{i,{\rm cc},-}\right|,\left|(Ds)_{i,{\rm cc},+}\right|\right].

  If :math:`(Ds)_{i,{\rm face,min}} \ge (Ds)_{i,{\rm cc,min}}`, set
  :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm face},\pm}`. Otherwise, set
  :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm cc},\pm}`. Finally, we are at an extreumum if
  :math:`(Ds)_{i,+}(Ds)_{i,-} \le 0`.

* Now that we have finished the extremum tests, if we are at an extremum,
  we scale :math:`\alpha_{i,\pm}`. First, we define

  .. math::

     \begin{aligned}
     (D^2s)_{i} &=& 6(\alpha_{i,+}+\alpha_{i,-}) \\
     (D^2s)_{i,L} &=& s_{i-2}-2s_{i-1}+s_{i} \\
     (D^2s)_{i,R} &=& s_{i}-2s_{i+1}+s_{i+2} \\
     (D^2s)_{i,C} &=& s_{i-1}-2s_{i}+s_{i+1}\end{aligned}

  Then, define

  .. math:: s = \text{sign}\left[(D^2s)_{i}\right],

  .. math:: (D^2s)_{i,\text{lim}} = \max\left\{\min\left[s(D^2s)_{i},Cs(D^2s)_{i,L},Cs(D^2s)_{i,R},Cs(D^2s)_{i,C}\right],0\right\}.

  Then,

  .. math:: \alpha_{i,\pm} = \frac{\alpha_{i,\pm}(D^2s)_{i,\text{lim}}}{\max\left[\left|(D^2s)_{i}\right|,1\times 10^{-10}\right]}

  Otherwise, if we are not at an extremum and :math:`|\alpha_{i,\pm}| > 2|\alpha_{i,\mp}|`,
  then define

  .. math:: s = \text{sign}(\alpha_{i,\mp})

  .. math:: \delta\mathcal{I}_{\text{ext}} = \frac{-\alpha_{i,\pm}^2}{4\left(\alpha_{j,+}+\alpha_{j,-}\right)},

  .. math:: \delta s = s_{i\mp 1} - s_i,

  If :math:`s\delta\mathcal{I}_{\text{ext}} \ge s\delta s`, then we perform the following test.
  If :math:`s\delta s - \alpha_{i,\mp} \ge 1\times 10^{-10}`, then

  .. math:: \alpha_{i,\pm} =  -2\delta s - 2s\left[(\delta s)^2 - \delta s \alpha_{i,\mp}\right]^{\myhalf}

  otherwise,

  .. math:: \alpha_{i,\pm} =  -2\alpha_{i,\mp}

  Finally, :math:`s_{i,\pm} = s_i + \alpha_{i,\pm}`.
