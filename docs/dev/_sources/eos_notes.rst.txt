***********************
Equation of State Notes
***********************

EOS Calls
=========

Initialization
--------------

advance_timestep
----------------

Step 1.
    *Define the average expansion at time* :math:`t^\nph` *and the new* :math:`w_0.` *
    if ``dpdt_factor``* :math:`>` 0 *then*

    -  In ``makePfromRhoH``, we compute a thermodynamic :math:`p^n` for the volume discrepancy
       term using :math:`(\rho,h,X)^n`.

    end if

Step 2.
    *Construct the provisional edge-based advective velocity,* :math:`\uadvone.`

Step 3.
    *React the full state and base state through the first time interval
    of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` burner, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^n`.

    -  In ``react_state``, we compute :math:`T^{(1)}` using :math:`(\rho,h,X)^{(1)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(1)}` and :math:`p_0^n` (if ``use_tfromp = T``)

    if ``evolve_base_state = T`` then

    -  In ``make_gamma``, we compute :math:`\Gamma_1` using :math:`(\rho,X)^{(1)}` and :math:`p_0^n`.

    end if

Step 4.
    | *Advect the base state, then the full state, through a time interval
      of* :math:`\dt.`
    | if ``use_thermal_diffusion = T`` then

    -  In advance before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{(1)}`.

    -  In ``enthalpy_advance`` :math:`\rightarrow` update_scal, we compute :math:`h` above
       the ``base_cutoff_density`` using :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}`.

    -  In ``thermal_conduct``, we compute :math:`T^{(2),*}` using :math:`(\rho,h,X)^{(2),*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}` (if ``use_tfromp = T``).

    else

    -  In ``enthalpy_advance`` :math:`\rightarrow` ``update_scal``, we compute :math:`h` above
       the ``base_cutoff_density`` using :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}`.

    -  In ``enthalpy_advance``, we compute :math:`T^{(2),*}` using :math:`(\rho,h,X)^{(2),*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}` (if ``use_tfromp = T``).

    end if

Step 5.
    *React the full state through a second time interval of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` burner, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^{(2),*}`.

    -  In ``react_state``, we compute :math:`T^{n+1,*}` using :math:`(\rho,h,X)^{n+1,*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{n+1,*}` and :math:`p_0^{n+1,*}` (if use_tfromp = T).

    if ``evolve_base_state`` then

    -  In ``make_gamma``, we compute :math:`\Gamma_1` using :math:`(\rho,X)^{n+1,*}` and :math:`p_0^{n+1,*}`.

    end if

Step 6.
    | *Define a new average expansion rate at time* :math:`t^\nph.`
    | if ``use_thermal_diffusion`` then

    -  In advance before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{n+1,*}`.

    end if

    -  In ``make_S``, we compute thermodynamic variables using :math:`(\rho,T,X)^{n+1,*}`.

    if ``dpdt_factor`` :math:`>` 0 then

    -  In ``makePfromRhoH``, we compute a thermodynamic :math:`p^{n+1,*}` for the volume
       discrepancy term using :math:`(\rho,h,X)^{n+1,*}`.

    end if

Step 7.
    *Construct the final edge-based advective velocity,* :math:`\uadvtwo.`

Step 8.
    | *Advect the base state, then the full state, through a time interval
      of* :math:`\dt.`
    | if ``use_thermal_diffusion = T`` then

    -  In ``enthalpy_advance`` :math:`\rightarrow` ``update_scal``, we compute :math:`h` above
       the ``base_cutoff_density`` using :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}`.

    -  In ``advance`` before ``thermal_conduct``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{(2),*}`.

    -  In ``thermal_conduct``, we compute :math:`T^{(2)}` using :math:`(\rho,h,X)^{(2)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

    else

    -  In ``enthalpy_advance`` :math:`\rightarrow` ``update_scal``, we compute :math:`h` above
       the ``base_cutoff_density`` using :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}`.

    -  In ``enthalpy_advance``, we compute :math:`T^{(2)}` using :math:`(\rho,h,X)^{(2)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

    end if

Step 9.
    *React the full state and base state through a second time interval
    of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` burner, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^{(2)}`.

    -  In ``react_state``, we compute :math:`T^{n+1}` using :math:`(\rho,h,X)^{n+1}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{n+1}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

    if ``evolve_base_state = T`` then

    -  In ``make_gamma``, we compute :math:`\Gamma_1` using :math:`(\rho,X)^{n+1}` and :math:`p_0^{n+1}`.

    end if

Step 10.
    | *Compute* :math:`S^{n+1}` *for the final projection.*
    | if ``make_explicit_thermal`` then

    -  In ``advance`` before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{n+1}`.

    end if

    -  In ``make_S``, we compute thermodynamic variables using :math:`(\rho,T,X)^{n+1}`.

Step 11.
    | *Update the velocity.*
    | if ``dpdt_factor`` :math:`>` 0 then

    -  In ``makePfromRhoH``, we compute a thermodynamic :math:`p^{n+1}` for the volume
       discrepancy term using :math:`(\rho,h,X)^{n+1}`.

    end if

Step 12.
    *Compute a new* :math:`\dt.`

make_plotfile
-------------

Temperature Usage
=================

.. _advance_timestep-1:

advance_timestep
----------------

Step 1.
    *Define the average expansion at time* :math:`t^\nph` *and the new* :math:`w_0.`

Step 2.
    *Construct the provisional edge-based advective velocity,* :math:`\uadvone.`

Step 3.
    *React the full state and base state through the first time interval
    of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` ``burner``, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^n`.

    -  In ``react_state``, we compute :math:`T^{(1)}` using :math:`(\rho,h,X)^{(1)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(1)}` and :math:`p_0^n` (if ``use_tfromp = T``).

Step 4.
    | *Advect the base state, then the full state, through a time interval
      of* :math:`\dt.`
    | if ``use_thermal_diffusion = T`` then

    -  In ``advance`` before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{(1)}`.

    -  In ``thermal_conduct``, we compute :math:`T^{(2),*}` using :math:`(\rho,h,X)^{(2),*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}` (if ``use_tfromp = T``).

    else

    -  In ``enthalpy_advance``, we compute :math:`T^{(2),*}` using :math:`(\rho,h,X)^{(2),*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2),*}` and :math:`p_0^{n+1,*}` (if ``use_tfromp = T``).

    end if

Step 5.
    *React the full state through a second time interval of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` ``burner``, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^{(2),*}`.

    -  In ``react_state``, we compute :math:`T^{n+1,*}` using :math:`(\rho,h,X)^{n+1,*}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{n+1,*}` and :math:`p_0^{n+1,*}` (if ``use_tfromp = T``).

Step 6.
    | *Define a new average expansion rate at time* :math:`t^\nph.`
    | if ``use_thermal_diffusion = T`` then

    -  In advance before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{n+1,*}`.

    end if

    -  In ``make_S``, we compute thermodynamic variables using :math:`(\rho,T,X)^{n+1,*}`.

Step 7.
    *Construct the final edge-based advective velocity,* :math:`\uadvtwo.`

Step 8.
    | *Advect the base state, then the full state, through a time interval
      of* :math:`\dt.`
    | if ``use_thermal_diffusion = T`` then

    -  In ``advance`` before ``thermal_conduct``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{(2),*}`.

    -  In ``thermal_conduct``, we compute :math:`T^{(2)}` using :math:`(\rho,h,X)^{(2)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

    else

    -  In ``enthalpy_advance``, we compute :math:`T^{(2)}` using :math:`(\rho,h,X)^{(2)}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{(2)}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

    end if

Step 9.
    *React the full state and base state through a second time interval
    of* :math:`\dt/2.`

    -  In ``react_state`` :math:`\rightarrow` ``burner``, we compute :math:`c_p` and :math:`\xi_k`
       for inputs to VODE using :math:`(\rho,T,X)^{(2)}`.

    -  In ``react_state``, we compute :math:`T^{n+1}` using :math:`(\rho,h,X)^{n+1}`
       (if ``use_tfromp = F``) or :math:`(\rho,X)^{n+1}` and :math:`p_0^{n+1}` (if ``use_tfromp = T``).

Step 10.
    | *Compute* :math:`S^{n+1}` *for the final projection.*
    | if ``make_explicit_thermal`` then

    -  In ``advance`` before ``make_explicit_thermal``, we compute the coefficients for
       :math:`\nabla\cdot(\kth/c_p)\nabla h + \cdots` using :math:`(\rho,T,X)^{n+1}`.

    end if

    -  In ``make_S``, we compute thermodynamic variables using :math:`(\rho,T,X)^{n+1}`.

Step 11.
    *Update the velocity.*

Step 12.
    *Compute a new* :math:`\dt.`
