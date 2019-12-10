module model_util_module

  use network, only : nspec
  use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
  use eos_type_module
  use eos_module
  use make_grav_module, only : grav_zone
  use probin_module, only: gamma1

  implicit none

contains

  pure function set_species(y) result (xn)

    double precision, intent(in) :: y
    double precision :: xn(nspec)

    xn(:) = ZERO
    xn(1) = ONE - fv(y)
    xn(2) = fv(y)

  end function set_species

  pure function fv(y) result (f_v)

    double precision, intent(in) :: y
    double precision :: f_v

    if (y < 1.9375d0 * 4.d8) then
       f_v = ZERO

    else if (y > 2.0625d0 * 4.d8) then
       f_v = ONE

    else
       f_v = HALF * (ONE + sin(8.d0 * M_PI * (y/4.d8 - 2.d0)))
    endif

  end function fv

  function dUdy(y, U) result (dU)

    double precision, intent(in) :: y, U(2)
    double precision :: dU(2), gamma, gamma0

    type(eos_t) :: eos_state

    eos_state % rho = U(1)
    eos_state % p = U(2)
    eos_state % xn = set_species(y)

    call eos(eos_input_rp, eos_state)

    gamma0 = eos_state % gam1
    gamma = gamma0 + fv(y) * (gamma1 - gamma0)

    dU(2) = exp(U(1)) * grav_zone(y) / exp(U(2))
    dU(1) = dU(2) / gamma

  end function dUdy

end module model_util_module
