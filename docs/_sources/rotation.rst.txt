*********************
Rotation in MAESTROeX
*********************

Rotation
========

To handle rotation in MAESTROeX we move to the co-rotating reference frame. Time
derivatives of a vector in the inertial frame are related to those in the
co-rotating frame by

.. math::
   \left(\frac{D}{Dt}\right)_\text{i} =
   \left[\left(\frac{D}{Dt}\right)_\text{rot} + \Omegab\times\ \right],
   :label: eq:derivative relations

where i(rot) refers to the inertial(co-rotating) frame and :math:`\Omegab` is
the angular velocity vector. Using :eq:`eq:derivative relations` and
assuming :math:`\Omegab` is constant, we have

.. math::
   \frac{Dv_\text{i}}{Dt} = \frac{Dv_\text{rot}}{Dt} +
   2\Omegab \times \mathbf{v_\text{rot}} +
   \Omegab\times\left(\Omegab\times r_\text{rot}\right).
   :label: eq:rotational velocity relation

Plugging this into the momentum equation and making use of the continuity
equation we have a momentum equation in the co-rotating frame:

.. math::
   \frac{\partial(\rho\mathbf{U})}{\partial t} +
   \nabla\cdot\left(\mathbf{U}(\rho\mathbf{U}) + p \right) =
   \underbrace{-2\rho\Omegab\times\mathbf{U}}_{\text{Coriolis}} -
   \underbrace{\rho\Omegab\times\left(\Omegab\times
     \rb\right)}_{\text{Centrifugal}} -
   \rho \mathbf{g}.
   :label: eq:momentum equation with rotation

The Coriolis and Centrifugal terms are force terms that will be added to
right hand side of the equations in mk_vel_force. Note that the
Centrifugal term can be rewritten as a gradient of a potential and absorbed
into either the pressure or gravitational term to create an effective pressure
or effective gravitational potential. Furthermore, because it can be written
as a gradient, this term would drop out completely in the projection if we were
doing incompressible flow.
In what follows we will include the Centrifugal term explicitly.

.. _Sec:Using Spherical Geometry:

Using Spherical Geometry
------------------------

In spherical geometry implementing rotation is straightforward as we don’t have
to worry about special boundary conditions. We assume a geometry as shown in
:numref:`Fig:rotation in spherical` where the primed coordinate system is
the simulation domain coordinate system, the unprimed system is the primed
system translated by the vector :math:`\rb_\text{c}` and is added here for
clarity. The :math:`\rb_\text{c}` vector is given by center in the
geometry module. In spherical, the only runtime parameter of importance
is rotational_frequency in units of Hz. This specifies the angular
velocity vector which is assumed to be along the **k** direction:

.. math::

   \Omegab \equiv \Omega {\bf k} = 2\pi *
   \left(\text{rotational_frequency}\right) {\bf k}.

The direction of :math:`\rb` is given as the normal vector which is
passed into mk_vel_force; in particular

.. math::

   \cos\theta \equiv \frac{\rb\cdot {\bf k}}{r} =
   {\bf normal}\cdot {\bf k} = \text{normal}(3).

The magnitude of :math:`\rb` is calculated based on the current zone’s location with
respect to the center.
Using this notation we can write the Centrifugal term as

.. math::
   \begin{align}
   \Omegab\times\left(\Omegab\times\rb\right) &=
   \left(\Omegab\cdot\rb\right)\Omegab - \left(\Omegab\cdot\Omegab\right)\rb\\
   &= \Omega^2 r *\left[\text{normal}(3)\right]{\bf k} -
   \Omega^2 r *{\bf normal} = \left(
   \begin{array}{c}
   \Omega^2r*\left[\text{normal}(1)\right]\\
   \Omega^2r*\left[\text{normal}(2)\right]\\
   0 \end{array}\right).
   \end{align}

The Coriolis term is a straightforward cross-product:

.. math::
   \begin{align}
   \Omegab \times \Ub &= \left|
   \begin{array}{ccc}
     {\bf{i}}&{\bf{j}}&{\bf{k}}\\
     0 & 0 & \Omega\\
     u & v & w
   \end{array}\right|\\
   &= \left(
   \begin{array}{c}
   -\Omega v\\ \Omega u \\ 0
   \end{array}
   \right).
   \end{align}

.. _Fig:rotation in spherical:
.. figure:: rotation_spherical.png
   :align: center

   Geometry of rotation when spherical_in :math:`=1`. We assume the
   star to be rotating about the :math:`z` axis with rotational frequency :math:`\Omega`.


.. _Sec:Using Plane-Parallel Geometry:

Using Plane-Parallel Geometry
-----------------------------
