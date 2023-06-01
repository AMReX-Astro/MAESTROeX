************************
Notes on :math:`\beta_0`
************************

The goal of :math:`\beta_0` is to capture the expansion of a displaced fluid
element due to the stratification of the atmosphere. MAESTROeX computes
:math:`\beta_0` as:

.. math::

   \label{eq:beta_0}
   \beta_0(r,t) = \rho_0 \exp\left (  \int_0^r  \frac{1}{\gammabar p_0} \frac{\partial p_0}{\partial r^\prime} dr^\prime \right )

Constant Composition
====================

Consider an isentropically stratified atmosphere, with a constant
composition as a function of :math:`r`. If you displace a parcel of fluid
upwards, it will expand adabatically and continue to rise until its
density matches the background density. Even if :math:`\gammabar` is not
constant in :math:`r`, following from the definition of :math:`\beta_0`,

.. math:: \frac{1}{\beta_0} \frac{d\beta_0}{dr} = \frac{1}{\gammabar p_0} \frac{dp_0}{dr}

and the definition of :math:`\Gamma_1`,

.. math:: \Gamma_1 = \left . \frac{d \log p}{d \log \rho} \right |_s

So, at constant entropy, from the definition of :math:`\Gamma_1`, it must hold
that

.. math:: \frac{1}{\rho} \frac{d \rho}{dr} = \frac{1}{\Gamma_1 p} \frac{d p}{dz}  .

Comparing to the definition of :math:`\beta_0` then

.. math:: \frac{1}{\beta_0} \frac{\beta_0}{dr} =\frac{1}{\gammabar p_0}\frac{dp_0}{dr} = \frac{1}{\rho_0} \frac{d\rho_0}{dr}   .

Therefore, :math:`\beta_0 = \rho_0`.

This means that if we have a constant composition and an
isentropically stratified atmosphere, as we displace a fluid element,
it will always remain neutrally buoyant.

Composition Gradient
====================

If there is a change in composition with :math:`r`, the situation is more
complicated. Consider again an isentropically stratified atmosphere,
now with a composition gradient. If you displace a parcel of fluid
upwards, it will rise. If there are no processes that change the
composition (e.g. reactions), then the composition in the fluid
element will remain fixed. As it rises, it will the ambient medium
will have a different composition that it has. In this case, what is
the path to equilibrium?

.. _Sec:On the Affect of Chemical Potential:

On the Effect of Chemical Potential
===================================

In MAESTRO, we do things in an operator split fashion — the hydro is
de-coupled from the burning. This means that during the hydro parts
of the algorithm (where :math:`\beta_0` is used), the system is fixed in
chemical equilibrium. For completeness, however, here we describe the
effects of the species’ chemical potentials, which were neglected in
the original derivation of :math:`\beta_0`. Note that similar terms appear
in the calculation of things such as specific heats,
which *may* be important in the burning step — there appears to
be very little about this in the literature, but everyone seems to
assume it makes little difference.

.. _Sec:Derivation of alpha:

Derivation of :math:`\alpha`
----------------------------

In paper I, :math:`\alpha` is defined as

.. math::
   \alpha\equiv -\left(
   \frac{(1-\rho h_p)p_T-\rho c_p}{\rho^2c_pp_\rho}\right)
   :label: eq:alpha

where

.. math::

   h_p \equiv \left(\frac{\partial h}{\partial p}\right)_{T,X}, \quad
   c_p \equiv \left(\frac{\partial h}{\partial T}\right)_{p,X}, \quad
   p_T \equiv \left(\frac{\partial p}{\partial T}\right)_{\rho,X}, \quad
   p_\rho \equiv \left(\frac{\partial p}{\partial \rho}\right)_{T,X}

where the subscript :math:`X` means holding all :math:`X_i` constant. In the
absence of reactions, the :math:`X` subscript can be dropped from all
derivatives and with the use of the equation of state :math:`p=p(\rho,T)`,
:math:`\alpha` can be written as :math:`\alpha=\alpha(\rho,T)`. Such a system
without reactions and in thermal equilibrium could be either a pure
system of one species, or a system of many species in chemical (and
therefore *thermodynamic*) equilibrium. Cox and Giuli (hereafter
CG) call the former type of system a “simple system” and therefore
the latter a “non-simple system” in chemical equilibrium. The
analysis in paper I that reduced :eq:`eq:alpha` to

.. math::
   \alpha = \frac{1}{\Gamma_1p_0}
   :label: eq:alpha_simp_no_rxn

used CG’s discussion of the various adiabatic :math:`\Gamma`\ ’s. However,
their discussion only pertains to “simple systems” or “non-simple
systems” in chemical equilibrium. In general, nuclear reactions will
be important and therefore this analysis needs to be reformed.

Even in the presence of reactions, :eq:`eq:alpha` can be rewritten
as was done in paper I:

.. math::
   \alpha = -\frac{1}{p\chi_\rho c_p}\left[\left(\frac{1}{\rho\chi_\rho}
   - \frac{\rho e_\rho}{p\chi_\rho}\right)\frac{p\chi_T}{T} - c_p\right],
   :label: eq:alpha2

where

.. math::

   \begin{aligned}
   \chi_{\rho} &\equiv \left(\frac{\partial\ln p}{\partial\ln\rho}
   \right)_{T,X} \\
   \chi_{T} &\equiv \left(\frac{\partial\ln p}{\partial\ln T}
   \right)_{\rho,X}.\end{aligned}

Following the results of paper I, we want to find a relation
between :math:`p\chi_\rho` and :math:`\Gamma_1`.

For an equation of state :math:`p=p(\rho,T,X)` we have

.. math::

   d\ln p = \left(\frac{\partial\ln p}{\partial\ln\rho}\right)_{T,X}d\ln\rho +
   \left(\frac{\partial\ln p}{\partial\ln T}\right)_{\rho,X}d\ln T +
   \sum_i\left(\frac{\partial\ln p}{\partial\ln X_i}\right)_{\rho,T,(X_j,j
   \neq i)} d\ln X_i.

We define another logarithmic derivative

.. math::

   \begin{aligned}
   \chi_{X_{i}} &\equiv \left(\frac{\partial\ln p}{\partial\ln X_i}
   \right)_{\rho,T,(X_j,j\neq i)}\end{aligned}

and therefore

.. math::

   d\ln p = \chi_\rho \ d\ln\rho + \chi_T \ d\ln T + \sum_i \chi_{X_i}\
   d\ln X_i.

From here we get the general statement

.. math::

   \frac{\partial\ln p}{\partial \ln \rho} = \chi_\rho +
   \chi_T\frac{\partial \ln T}{\partial\ln \rho} +
   \sum_i\chi_{X_i}\frac{\partial\ln X_i}{\partial\ln \rho}

which must hold for an adiabatic process as well, and therefore we have

.. math::
   \Gamma_1 = \chi_\rho + \chi_T\left(\Gamma_3-1\right)
     + \sum_i\chi_{X_i}\Gamma_{4,i}
   :label: eq:gamma1

where we use CG’s definition of :math:`\Gamma_1` and :math:`\Gamma_3` and introduce a
fourth gamma function:

.. math::

   \Gamma_1 \equiv \left(
   \frac{\partial \ln p}{\partial \ln \rho}\right)_{\text{AD}},\quad
   \Gamma_3-1\equiv \left(
   \frac{\partial \ln T}{\partial \ln \rho}\right)_{\text{AD}},\quad
   \Gamma_{4,i}\equiv \left(
   \frac{\partial\ln X_i}{\partial\ln\rho}\right)_{\text{AD}},

where the subscript AD means along an adiabat. We now derive an expression
for :math:`\Gamma_3`.

The first law of thermodynamics can be written as

.. math:: dQ = dE + pdV - \sum_i\mu_idN_i

where :math:`\mu_i=\left(
\frac{\partial E}{\partial N_i}\right)_{\text{AD},\rho,(N_j,j\neq i)}` is
the chemical potential; or per unit mass we have

.. math::

   \begin{aligned}
     dq &= de - \frac{p}{\rho^2}d\rho - \sum_i\mu_id
     \left(\frac{n_i}{\rho}\right)\\
     &= de - \frac{p}{\rho^2}d\rho - \sum_i
     \left(
     \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}dX_i\end{aligned}

where we have used :math:`X_i \equiv \rho_i/\rho = A_in_i/\rho N_\text{A}`
and the chemical potential has been replaced with
:math:`\mu_i = \frac{A_i}{N_\text{A}}\left(\frac{\partial e}{\partial X_i}  \right)_{\rho,\text{AD},(X_j,j\neq i)}`.
Using this and expressing the specific internal energy as :math:`e=e(\rho,T,X)`
we then have

.. math::

   dq = c_vdT +
   \left[\left(\frac{\partial e}{\partial \rho}\right)_{T,X} -\frac{p}{\rho^2}
     \right]d\rho +
   \sum_i\left[
     \left(\frac{\partial e}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)} -
     \left(\frac{\partial e}{\partial X_i}
     \right)_{\rho,\text{AD},(X_j,j\neq i)}\right]dX_i

and

.. math::
   \begin{aligned}
   \left(\frac{d\ln T}{d\ln\rho}\right)_\text{AD} \equiv \Gamma_3-1
   &= \frac{1}{c_vT}\left[
   \frac{p}{\rho} - \left(\frac{\partial e}{\partial\ln\rho}\right)_{T,X} +
   \right.{}\nonumber\\
   &\qquad\qquad  \left.\sum_i \left[
       \left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       -
       \left(\frac{\partial e}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)}
       \right]X_i\Gamma_{4,i}\right]
   \end{aligned}
   :label: eq:gamma3_first

Now we need to evaluate :math:`\left(\partial e/\partial \ln\rho\right)_{T,X}`.
Again using the first law and the fact that :math:`ds=dq/T` is an exact
differential (i.e. mixed derivatives are equal) we have

.. math::

   \begin{aligned}
   \label{eq:dedlnrho}
     \left(
     \frac{\partial}{\partial\rho}\left[\frac{c_v}{T}\right]\right)_{T,X} &=
     \left(\frac{\partial}{\partial T}\left[\frac{1}{T}
       \left(\frac{\partial e}{\partial\rho}\right)_{T,X} - \frac{p}{T\rho^2}
       \right]\right)_{\rho,X}{}\nonumber\\
     \frac{1}{T}\left(\frac{\partial}{\partial\rho}\left(
     \frac{\partial e}{\partial T}\right)_{\rho,X}\right)_{T,X} &=
     -\frac{1}{T^2}\left(\frac{\partial e}{\partial\rho}\right)_{T,X} +
     \frac{1}{T}\left(\frac{\partial}{\partial T}\left(
     \frac{\partial e}{\partial\rho}\right)_{T,X}\right)_{\rho,X}
     +\frac{p}{T^2\rho^2} -
     \frac{1}{T\rho^2}\left(\frac{\partial p}{\partial T}\right)_{\rho,X}
     {}\nonumber\\
     \therefore\quad \left(\frac{\partial e}{\partial\ln \rho}\right)_{T,X} &=
     \frac{p}{\rho}\left(1-\chi_T\right),\end{aligned}

exactly the same result if we were to exclude species information.
Similarly, we can find an expression for the derivative of energy with
respect to composition

.. math::

   \begin{aligned}
     \left(\frac{\partial}{\partial X_i}\left[
       \frac{c_v}{T}\right]\right)_{\rho,T,(X_j,j\neq i)} &=
     \left(\frac{\partial}{\partial T}\left[
       \frac{1}{T}\left(\frac{\partial e}{\partial X_i}
         \right)_{\rho,T,(X_j,j\neq i)} - \frac{1}{T}\left(
         \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right]\right)_{\rho,X}\\
     \frac{1}{T}\left(\frac{\partial }{\partial X_i}\left(
     \frac{\partial e}{\partial T}\right)_{\rho,X}\right)_{\rho,T,(X_j,j\neq i)}
     &= \frac{1}{T^2}\left[\left(\frac{\partial e}{\partial X_i}
       \right)_{\rho,\text{AD},(X_j,j\neq i)} -
       \left(\frac{\partial e}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)}
       \right] + \\
     &\ \ \ \ \ \frac{1}{T}\left[
       \left(\frac{\partial}{\partial T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)}
       \right)_{\rho,X} -
       \left(\frac{\partial }{\partial T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}\right]\\
     \therefore\quad
     \left(\frac{\partial e}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)} &=
     \left(\frac{\partial e}{\partial X_i}
     \right)_{\rho,\text{AD},(X_j,j\neq i)} - \left(
     \frac{\partial}{\partial\ln T}\left(
     \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
     \right)_{\rho,X}.\end{aligned}

Plugging these back into :eq:`eq:gamma3_first` we have

.. math::
   \Gamma_3-1 = \frac{1}{c_vT}\left[\frac{p}{\rho}\chi_T +\sum_i
       \left(\frac{\partial}{\partial\ln T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}X_i
       \Gamma_{4,i}\right],
   :label: eq:gamma3_second

or

.. math::
   c_v = \frac{1}{T(\Gamma_3-1)}\left[\frac{p}{\rho}\chi_T +\sum_i
       \left(\frac{\partial}{\partial\ln T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}X_i
       \Gamma_{4,i}\right].
   :label: eq:cv

We can obtain an expression for the specific heat at constant pressure
from the enthalpy

.. math::

   \begin{aligned}
     c_p \equiv \left(\frac{\partial h}{\partial T}\right)_{p,X} &=
     \left(\frac{\partial e}{\partial T}\right)_{p,X} - \frac{p}{\rho^2}
     \left(\frac{\partial \rho}{\partial T}\right)_{p,X}\\
     &= \left(\frac{\partial e}{\partial T}\right)_{p,X} + \frac{p}{\rho^2}
     \left(\frac{\partial p}{\partial T}\right)_{\rho,X}
     \left(\frac{\partial \rho}{\partial p}\right)_{T,X}\\
     &=\left(\frac{\partial e}{\partial T}\right)_{p,X} + \frac{p}{\rho T}
     \frac{\chi_t}{\chi_\rho}.\end{aligned}

The first term on the rhs can be obtained from writing :math:`e=e(p,T,X)` and
:math:`p=p(\rho,T,X)`:

.. math::

   \begin{aligned}
     de &= \left(\frac{\partial e}{\partial p}\right)_{T,X}dp
     + \left(\frac{\partial e}{\partial T}\right)_{p,X}dT +
     \sum_i \left(\frac{\partial e}{\partial X_i}\right)_{p,T,(X_j,j\neq i)}
     dX_i\\
     dp &= \left(\frac{\partial p}{\partial \rho}\right)_{T,X}d\rho +
     \left(\frac{\partial p}{\partial T}\right)_{\rho,X}dT + \sum_i
     \left(\frac{\partial p}{\partial X_i}\right)_{\rho,T,(X_j,j\neq i)}dX_i\\
     \therefore \ \left(\frac{\partial e}{\partial T}\right)_{\rho,X} &=
     \left(\frac{\partial e}{\partial p}\right)_{T,X}
     \left(\frac{\partial p}{\partial T}\right)_{\rho,X} +
     \left(\frac{\partial e}{\partial T}\right)_{p,X}\\
     \Rightarrow \ \left(\frac{\partial e}{\partial T}\right)_{p,X} &= c_v -
     \left(\frac{\partial e}{\partial \rho}\right)_{T,X}
     \left(\frac{\partial \rho}{\partial p}\right)_{T,X}
     \left(\frac{\partial p}{\partial T}\right)_{\rho,X}\\
     &= c_v - \frac{p\chi_T}{\rho T\chi_\rho}\left(1-\chi_T\right)\end{aligned}

and

.. math:: c_p = \frac{p}{\rho T}\frac{\chi_T^2}{\chi_\rho} + c_v

Dividing this by :eq:`eq:cv` and using the relation between the
:math:`\Gamma` ’s, :eq:`eq:gamma1`, we then have

.. math::

   \begin{aligned}
   \label{eq:pchirho}
     \gamma \equiv \frac{c_p}{c_v} &= 1 + \frac{p(\Gamma_3-1)}{\rho }
     \frac{\chi_T^2}{\chi_\rho}\left[\frac{p}{\rho}\chi_T +\sum_i
       \left(\frac{\partial}{\partial\ln T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}X_i
       \Gamma_{4,i}\right]^{-1}{}\nonumber\\
     &= 1 + \frac{p\chi_T\left(\Gamma_1 - \chi_\rho -
       \sum_i \chi_{X_i}\Gamma_{4,i}\right)}{p\chi_\rho\chi_T + \rho
       \chi_\rho\sum_i \left(
       \frac{\partial}{\partial \ln T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}X_i\Gamma_{4,i}}{}\nonumber\\
     &= \frac{p\chi_T\Gamma_1 + \sum_i \left[\rho\chi_\rho\left(
       \frac{\partial}{\partial \ln T}\left(\frac{\partial e}{\partial X_i}
       \right)_{\rho,\text{AD},(X_j,j\neq i)}\right)_{\rho,X}X_i - p\chi_T
       \chi_{X_i}\right]\Gamma_{4,i}}{p\chi_\rho\chi_T + \rho
       \chi_\rho\sum_i \left(
       \frac{\partial}{\partial \ln T}\left(
       \frac{\partial e}{\partial X_i}\right)_{\rho,\text{AD},(X_j,j\neq i)}
       \right)_{\rho,X}X_i\Gamma_{4,i}}{}\nonumber\\
     \Rightarrow p\chi_\rho &= \frac{1}{\chi_T\gamma}\left[p\chi_T\Gamma_1 +
       \sum_i \left[\rho\chi_\rho\left(1-\gamma\right)\left(
         \frac{\partial}{\partial \ln T}\left(\frac{\partial e}{\partial X_i}
         \right)_{\rho,\text{AD},(X_j,j\neq i)}\right)_{\rho,X}X_i - p\chi_T
         \chi_{X_i}\right]\Gamma_{4,i}\right].\end{aligned}

Plugging `[eq:pchirho] <#eq:pchirho>`__ into :eq:`eq:alpha2` and rewriting the
partial derivative of :math:`e` with the help of `[eq:dedlnrho] <#eq:dedlnrho>`__ we have

.. math::

   \begin{aligned}
   \alpha &= -\frac{1}{p\chi_\rho c_p}\left[\left(\frac{1}{\rho\chi_\rho}
     - \frac{\rho e_\rho}{p\chi_\rho}\right)\frac{p\chi_T}{T} - c_p\right] \\
   &=\frac{\gamma}{c_p}\frac{c_p\chi_T + \left(\rho
     \left(\frac{\partial e}{\partial\ln\rho}\right)_{T,X}-p\right)
     \frac{\chi_T^2}{T\rho\chi_\rho}}
   {p\chi_T\Gamma_1 +
       \sum_i \left[\rho\chi_\rho\left(1-\gamma\right)\left(
         \frac{\partial}{\partial \ln T}\left(\frac{\partial e}{\partial X_i}
         \right)_{\rho,\text{AD},(X_j,j\neq i)}\right)_{\rho,X}X_i - p\chi_T
         \chi_{X_i}\right]\Gamma_{4,i}}\\
   &=\frac{\gamma}{\Gamma_1 p c_p}\left[\frac{c_p - \frac{p\chi_T^2}
       {T\rho\chi_\rho}}
     {1 + \sum_i \left[\frac{\rho\chi_\rho}{p\chi_T}
         \left(1-\gamma\right)\left(
         \frac{\partial}{\partial \ln T}\left(\frac{\partial e}{\partial X_i}
         \right)_{\rho,\text{AD},(X_j,j\neq i)}\right)_{\rho,X}X_i -
         \chi_{X_i}\right]\frac{\Gamma_{4,i}}{\Gamma_1}}\right]\\
   &=\left(\frac{1}{\Gamma_1p}\right)
   \left[1 + \sum_i \left[\frac{\rho\chi_\rho}{p\chi_T}
       \left(1-\gamma\right)\left(
       \frac{\partial}{\partial \ln T}\left(\frac{\partial e}{\partial X_i}
       \right)_{\rho,\text{AD},(X_j,j\neq i)}\right)_{\rho,X}X_i -
       \chi_{X_i}\right]\frac{\Gamma_{4,i}}{\Gamma_1}\right]^{-1}\\\end{aligned}

.. math::

   \boxed{
     \alpha = \frac{1}{\Gamma_1p}\left[1 + \sum_i\left[\frac{\rho^2p_\rho}
         {pp_T}(1-\gamma)\frac{N_\text{A}}{A_i}
         \left(\frac{\partial\mu_i}{\partial T}\right)_{\rho,X}X_i -
         \chi_{X_i}\right]\frac{\Gamma_{4,i}}{\Gamma_1}\right]^{-1}
   }

.. _Recalling Derivation of beta0:

Recalling Derivation of :math:`\beta_0`
---------------------------------------

Recall from paper I that :math:`\beta_0` was derived from the equation

.. math:: \nabla\cdot\mathbf{U} + \alpha\mathbf{U}\cdot\nabla p_0 = \tilde{S}

in such a fashion that we ended up with an equation of the form

.. math::

   \label{eq:beta constraint}
   \nabla\cdot\left(\beta_0(r)\mathbf{U}\right) = \beta_0\tilde{S}.

The derivation in Appendix B of paper I for a :math:`\beta_0` that
satisfies `[eq:beta constraint] <#eq:beta constraint>`__ automatically assumed :math:`\alpha
= \left(\Gamma_{1_0}p_0\right)^{-1}`. This would have to be modified
with the above derivation of :math:`\alpha` to be correct in a non-operator
split fashion.
