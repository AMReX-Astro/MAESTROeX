We carry around three different 1D radial quantities: :math:`\etarhoec`
(edge-centered), :math:`\etarhocc` (cell-centered), and :math:`\divetarho`
(cell-centered). These notes discuss when each of these is used, and
how they are computed, in both plane-parallel and spherical.

The Mixing Term, :math:`\etarho`
================================

The base state evolves in response to heating and mixing in the star.
The density evolution is governed by

.. math::
   \frac{\partial \rho_0}{\partial t} = -
    \nablab \cdotb \left( \rho_0 w_0 \er \right)
   - \nablab \cdotb \left( \etarho \er \right)  ,
   :label: eq:rho0upd_new

with

.. math::
   \etarho(r) = \overline{\left(\rhop \Ubt \cdot \er \right)} = \frac{1}{A(\Omega_H)}
    \int_{\Omega_H}  (\rhop \Ubt \cdot \er ) \; dA  ,
   :label: eq:eta

designed to keep the average value of the full density, :math:`\rho`, over a
layer of constant radius in the star equal to :math:`\rho_0`. To complete
the update of the base state, we need evolution equations for the
pressure, :math:`p_0`, and velocity, :math:`w_0`. For spherical geometry, the
derivation of :math:`w_0` constraint equation is shown in
the multilevel paper, resulting in the following system

.. math::

   \begin{aligned}
   w_0 &= \ow + \delta w_0 \\
   \frac{1}{r^2} \frac{\partial}{\partial r} \left (r^2 \ow \right ) &= \Sbar \\
   \frac{\partial}{\partial r} \left[ \frac{\gammabar p_0}{r^2} \frac{\partial}{\partial r} (r^2 \dw) \right] &= - \frac{g}{r^2} \frac{\partial (r^2 \etarho)}{\partial r} - \frac{4 (\ow + \dw) \rho_0 g}{r}
   - 4 \pi G \rho_0 \etarho \label{eq:dw0constraint}\end{aligned}

In paper III, we introduced a mixing term, :math:`\etarho`, to the density
evolution equation (:eq:`eq:rho0upd_new`), with the objective of
keeping the base state density equal to the average of the density
over a layer. For a spherical base state, it is best to define the
average in terms of spherical coordinates,

.. math:: \overline{q} = \frac{1}{4\pi} \int_{\Omega_H} q(r,\theta,\phi) \; d\Omega

where :math:`\int_{\Omega_H} d\Omega = 4\pi` represents the integral over
the spherical :math:`\theta` and :math:`\phi` angles at constant radius.

Recall from Paper III that if initially :math:`\overline{\rho'} = 0`,
there is no guarantee that :math:`\overline{\rho'} = 0` will hold at later
time. To see this, start with the equation for the perturbational density

.. math::
   \frac{\partial\rho'}{\partial t} = -\Ub\cdot\nabla\rho' -
     \rho'\nabla\cdot\Ub - \nabla\cdot\left(\rho_0\Ubt\right),
   :label: eq:Perturbational Density

written in a slightly different form:

.. math:: \frac{\partial\rho'}{\partial t} + \nabla\cdot\left(\rho'\Ub\right) = -\nabla\cdot\left(\rho_0\Ubt\right).

We integrate this over a spherical shell of thickness :math:`2h` at radius :math:`r_0`, i.e.,
:math:`\Omega_H \times (r_0-h, r_0+h)`, and normalize by the integration
volume, which is :math:`\sim 4\pi r_0^2  2h` for small :math:`h`, to obtain:

.. math::

   \begin{aligned}
   \frac{1}{4\pi r_0^2 2h}\int_{r_0-h}^{r_0+h}r^2 dr\int_{\Omega_H}\left[\frac{\partial\rho'}{\partial t} + \nabla\cdot\left(\rho'\Ub\right)\right]d\Omega
   =& - \frac{1}{4\pi r_0^2 2h}\int_{r_0-h}^{r_0+h}r^2 dr\int_{\Omega_H}\nabla\cdot\left(\rho_0\Ubt\right)d\Omega \nonumber \\
   =& \left. -\frac{1}{4\pi r_0^2 2h}\int_{\Omega_H}\left[\rho_0\left(\Ubt\cdot\eb_r\right)\right]r^2 d\Omega\right|_{r_0-h}^{r_0+h} \nonumber \\
   =& 0,\end{aligned}

where we have used the divergence theorem in spherical coordinates to transform
the volume integral on the right hand side into an area integral over :math:`\Omega_H`.
We see that the right hand side disappears since
:math:`\int_{\Omega_H}\rho_0(\Ubt\cdot\eb_r)d\Omega=0`, which follows from the definition
of :math:`\Ubt`. Now, expanding the remaining terms and taking the limit as
:math:`h\rightarrow 0,` we can write

.. math::

   \begin{aligned}
   0 =& \lim_{h\rightarrow 0} \frac{1}{4\pi r_0^2 2h} \int_{r_0-h}^{r_0+h} r^2  dr \int_{\Omega_H} \left[\frac{\partial \rho'}{\partial t} + \nabla\cdot(\rho'\Ub)\right] d\Omega \nonumber \\
   =& \frac{\partial}{\partial t} \left( \lim_{h\rightarrow 0} \frac{1}{4\pi r_0^2} \frac{1}{2h} \int_{r_0-h}^{r_0+h} r^2  dr \int_{\Omega_H}  \rho' d\Omega \right) + \lim_{h\rightarrow 0} \left[\frac{1}{4\pi r_0^2  2h} \int_{r_0-h}^{r_0+h} r^2  dr \int_{\Omega_H}  \nabla \cdot ( \rho' \Ub )  d\Omega\right] \nonumber \\
   =& \frac{\partial}{\partial t} \left(\frac{1}{4\pi} \int_{\Omega_H}  \rho'  d\Omega \right) + \lim_{h\rightarrow 0} \left\{\left.\frac{1}{r_0^2 2h}\left[\frac{1}{4\pi}\int_{\Omega_H} \rho' (\Ub \cdot \eb_r)  d\Omega  \right]  r^2 \right |_{r_0-h}^{r_0+h} \right\} \nonumber  \\
   =&  \frac{\partial}{\partial t} \overline{\rho'} + \lim_{h\rightarrow 0} \left\{ \left .  \frac{1}{r_0^2  2h} \left[ \overline{\rho' (\Ub \cdot \eb_r)} \right]  r^2  \right |_{r_0-h}^{r_0+h} \right\} \nonumber  \\
   =&  \frac{\partial}{\partial t} \overline{\rho'} + \lim_{h\rightarrow 0} \frac{1}{r_0^2 2h} \int_{r_0-h}^{r_0+h} \nabla \cdot \left[\overline{\rho' (\Ub \cdot \eb_r)} \eb_r \right]  r^2  dr \nonumber \\
   =&  \frac{\partial}{\partial t} \overline{\rho'} + \nabla \cdot \left[ \overline{\rho' (\Ub \cdot \eb_r)} \eb_r \right]\end{aligned}

again using the divergence theorem, extracting the time derivative from the spatial integral,
and switching the order of operations as appropriate.

In short,

.. math::
   \frac{\partial}{\partial t} \overline{\rho'} = - \nabla\cdot\left[\overline{\rho'\left(\Ub\cdot\eb_r\right)}\eb_r\right] = -\nabla\cdot\left(\etarho\eb_r\right),
   :label: eq:rhopbar

and thus, :math:`\etarho = \overline{\left(\rho'\Ub\cdot\eb_r\right)}`.

We need both :math:`\etarho` alone and its divergence for the
various terms in the construction of :math:`w_0` and the correction to
:math:`\rho_0`. The quantity :math:`\etarho` is edge-centered on our grid, and
for Cartesian geometries, we constructed it by averaging the
appropriate fluxes through the grid boundaries. For a spherical base
state, this does not work, since the spherical shells do not line up
with the Cartesian grid boundaries.

Therefore, we take a different approach. We compute :math:`\etarho` by
constructing the quantity :math:`\rhop \Ubt \cdot \er` in each cell, and then
use our average routine to construct a 1-d, cell-centered
:math:`\eta_{\rho,r}` (this is essentially numerically solving the integral
in :eq:`eq:eta`. The edge-centered values of :math:`\etarho`,
:math:`\eta_{\rho,r+1/2}` are then constructed by simple
averaging:

.. math:: \eta_{\rho,r+1/2} = \frac{\eta_{\rho,r} + \eta_{\rho,r+1}}{2}  .

Instead of differencing :math:`\eta_{\rho,r+1/2}` to construct the
divergence, we instead use :eq:`eq:rhopbar` directly, by writing:

.. math::

   \left [ \nabla \cdot (\etarho \er ) \right ]^{n+1/2}
   = - \frac{\overline{\rhop^{n+1}} - \overline{\rhop^n}}{\Delta t}
   = - \frac{\overline{\rhop^{n+1}}}{\Delta t}  ,

where we have made use of the fact that :math:`\overline{\rhop^n} = 0` by construction.

:math:`\eta` Flow Chart
=======================

#. Enter ``advance_timestep`` with :math:`[\etarhoec, \etarhocc]^{n-\myhalf}`.

#. Call ``make_w0``. The spherical version uses uses :math:`\etarho^{{\rm ec},n-\myhalf}` and :math:`\etarho^{{\rm cc},n-\myhalf}`.

#. Call ``density_advance``. The plane-parallel version computes :math:`\etarho^{{\rm flux},n+\myhalf,*}`.

#. Call ``make_etarho`` to compute :math:`[\etarhoec, \etarhocc]^{n+\myhalf,*}`. The plane-parallel version uses :math:`\etarho^{{\rm flux},n+\myhalf,*}`.

#. Call ``make_psi``. The plane-parallel version uses :math:`\etarho^{{\rm cc},n+\myhalf,*}`.

#. Call ``make_w0``. The spherical version uses uses :math:`\etarho^{{\rm ec},n+\myhalf,*}` and :math:`\etarho^{{\rm cc},n+\myhalf,*}`.

#. Call ``density_advance``. The plane-parallel version computes :math:`\etarho^{{\rm flux},n+\myhalf}`.

#. Call ``make_etarho`` to compute :math:`[\etarhoec, \etarhocc]^{n+\myhalf}`. The plane-parallel version uses :math:`\etarho^{{\rm flux},n+\myhalf}`.

#. Call ``make_psi``. The plane-parallel version uses :math:`\etarho^{{\rm cc},n+\myhalf}`.

Computing :math:`\etarhoec` and :math:`\etarhocc`
=================================================

This is done in ``make_eta.f90``.

Plane-Parallel
--------------

We first compute a radial edge-centered multifab, :math:`\etarho^{\rm flux}`, using

.. math:: \eta_{\rho,\ib+\myhalf\eb_r}^{\rm flux} = \left[\left(\Ubt_{\ib+\myhalf\eb_r}^{n+\myhalf}\cdot\eb_r\right) + w_{0,r+\myhalf}^{n+\myhalf}\right] \rho_{\ib+\myhalf\eb_r}^{n+\myhalf} - w_{0,r+\myhalf}^{n+\myhalf}\rho_{0,r+\myhalf}^{n+\myhalf, {\rm pred}}

:math:`\etarhoec` is the edge-centered “average” value of :math:`\eta_{\rho}^{\rm flux}`,

.. math:: \eta_{\rho,r+\myhalf}^{\rm ec} = \overline{\eta_{\rho,\ib+\myhalf\eb_r}^{\rm flux}}

:math:`\etarhocc` is a cell-centered average of :math:`\etarhoec`,

.. math:: \eta_{\rho,r}^{\rm cc} = \frac{\eta_{\rho,r+\myhalf}^{\rm ec} + \eta_{\rho,r-\myhalf}^{\rm ec}}{2}.

.. _Sec:eta Spherical:

Spherical
---------

First, construct :math:`\eta_{\rho}^{\rm cart} =
[\rho'(\Ubt\cdot\eb_r)]^{n+\myhalf}` using:

.. math:: \left[\frac{\rho^n+\rho^{n+1}}{2}-\left(\frac{\rho_0^n+\rho_0^{n+1}}{2}\right)^{\rm cart}\right] \sum_d\left(\frac{\Ubt_{\ib+\myhalf\eb_d}^{n+\myhalf}\cdot\eb_d+\Ubt_{\ib-\myhalf\eb_d}^{n+\myhalf}\cdot\eb_d}{2}\right)n_d.

Then, :math:`\etarhocc` is the cell-centered average of :math:`\eta_{\rho}^{\rm cart}`,

.. math:: \etarhocc = \overline{\eta_{\rho}^{\rm cart}}.

On interior faces, :math:`\etarhoec` is the average of :math:`\etarhocc`,

.. math:: \eta_{\rho,r-\myhalf}^{\rm ec} = \frac{\eta_{\rho,r-1}^{\rm cc} + \eta_{\rho,r}^{\rm cc}}{2}.

At the upper and lower boundaries, we use

.. math::

   \begin{aligned}
   \eta_{\rho,-\myhalf}^{\rm ec} &=& 0, \\
   \eta_{\rho,{\rm nr}-\myhalf}^{\rm ec} &=& \eta_{\rho,{\rm nr}-1}^{\rm cc}.\end{aligned}

Using :math:`\etarhoec`
=======================

.. _plane-parallel-1:

Plane-Parallel
--------------

NOT USED.

Spherical
---------

In ``make_w0``, :math:`\etarhoec` is used in the construction of the RHS
for the :math:`\delta w_0` equation.

Using :math:`\etarhocc`
=======================

.. _plane-parallel-2:

Plane-Parallel
--------------

In ``make_psi``, :math:`\psi = \etarhocc g`.

.. _spherical-1:

Spherical
---------

In ``make_w0``, :math:`\etarhocc` is used in the construction of the RHS
for the :math:`\delta w_0` equation.
