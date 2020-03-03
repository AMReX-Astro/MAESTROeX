module state_indices_module
   use network, only : nspec, naux
   implicit none

   integer, parameter :: Nscal = 4 + nspec

   ! scalars
   integer, parameter :: rho_comp = 1
   integer, parameter :: rhoh_comp = 2
   integer, parameter :: temp_comp = 3
   integer, parameter :: pi_comp = 4
   integer, parameter :: spec_comp = 5

end module state_indices_module
