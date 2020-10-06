module state_indices_module
   use network, only : nspec, naux
   implicit none

   integer, parameter :: Nscal = 4 + nspec + naux

   ! scalars
   integer, parameter :: rho_comp = 1
   integer, parameter :: rhoh_comp = 2
   integer, parameter :: temp_comp = 3
   integer, parameter :: pi_comp = 4
   integer, parameter :: spec_comp = 5
   integer, parameter :: aux_comp = 5 + nspec

end module state_indices_module
