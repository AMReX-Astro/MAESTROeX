*************************
Volume Discrepancy Factor
*************************

The volume discrepancy term is used in the constraint equation to
force the system back to the equation of state. We write our velocity
constraint equation as

.. math::
   \nablab \cdotb (\beta_0 \Ub)  = \beta_0 \left(S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} - \frac{f}{\gammabar p_0} \frac{p_0 - p_\mathrm{EOS}}{\Delta t} \right )  .
   :label: eq:fullconstraint

Here, :math:`f` is the volume discrepancy factor and ranges from 0 to 1, and
:math:`p_\mathrm{EOS}` is the thermodynamic pressure as returned by the EOS,
using the full state as input.
In practice we evaluate this as

.. math:: p_\mathrm{EOS} = p(\rho,h,X_k)

The idea behind this forcing term is that if :math:`p_\mathrm{EOS} > p_0` then
:math:`\nablab \cdotb (\beta_0 \Ub) > 0`, and the cell will evacuate. This
has the effect of returning us to :math:`p_\mathrm{EOS} = p_0`.

In MAESTROeX, we decomponse the velocity into a base state component
and a local component. The base state constraint equation is simply
the horizontal average of the full constraint. Starting with
:math:`\Ub = \Ubt + w_0 \er` in :eq:`eq:fullconstraint`, we have

.. math::
   \nablab \cdotb (\beta_0 w_0 \er) + \nablab \cdotb (\beta_0 \Ubt)  = \beta_0 \left(S - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} - \frac{f}{\gammabar p_0} \frac{p_0 - p_\mathrm{EOS}}{\Delta t} \right )  .
   :label: eq:fullconstraint2

Averaging this over a layer, we note that :math:`\overline{\nablab \cdotb (\beta_0 \Ubt)} = 0`,
yielding

.. math::
   \nablab \cdotb (\beta_0 w_0 \er)  = \beta_0 \left(\overline{S} - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} - \frac{f}{\gammabar p_0} \frac{p_0 - \overline{p_\mathrm{EOS}}}{\Delta t} \right )
   :label: eq:w0constraint_vd

and

.. math::
   \nablab \cdotb (\beta_0 \Ubt)  = \beta_0 \left(S - \overline{S} + \frac{f}{\gammabar p_0} \frac{p_\mathrm{EOS} - \overline{p_\mathrm{EOS}}}{\Delta t} \right )  .
   :label: eq:Utconstraint_vd

In solving the :math:`w_0` evolution :eq:`eq:w0constraint_vd`, we expand the divergence, giving

.. math::

   \nablab \cdotb (w_0 \er)  = \overline{S} - \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial t} -
   w_0 \er  \cdotb \frac{1}{\beta_0} \nablab \beta_0  - \frac{f}{\gammabar p_0} \frac{p_0 - \overline{p_\mathrm{EOS}}}{\Delta t}  .

Recalling that

.. math::

   \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial r} =
    \frac{1}{\beta_0} \frac{\partial \beta_0}{\partial r}

(see Paper I), and that :math:`\psi \equiv D_0 p_0 / Dt \equiv \partial p_0 / \partial t +
w_0 \partial p_0 / \partial r`, we have

.. math::

   \nablab \cdotb (w_0 \er)  = \overline{S} - \frac{1}{\gammabar p_0} \psi -
    \frac{f}{\gammabar p_0} \frac{p_0 - \overline{p_\mathrm{EOS}}}{\Delta t}  .

This is the form that is solved in ``make_w0``.
