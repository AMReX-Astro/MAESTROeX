module problem_functions_module

  use amrex_constants_module
  use probin_module, only: min_time_step, react_its
  ! use meth_params_module, only: prob_lo

  implicit none

  private

contains

  subroutine get_min_timestep(dt) bind(C, name="get_min_timestep")

    double precision, intent(inout) :: dt

    dt = min_time_step

  end subroutine get_min_timestep

  subroutine get_react_its(its) bind(C, name="get_react_its")

    integer, intent(inout) :: its

    its = react_its

  end subroutine get_react_its

end module problem_functions_module
