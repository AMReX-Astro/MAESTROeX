
#include <Maestro.H>

using namespace amrex;

void
Maestro::MakeVelForce (Vector<MultiFab>& vel_force,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& uedge,
                       const Vector<MultiFab>& rho,
                       const Vector<Real>& rho0,
                       const Vector<Real>& grav_cell,
                       const Vector<Real>& w0_force,
                       int do_add_utilde_force)
{
    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        MultiFab& vel_force_mf = vel_force[lev];
        const MultiFab& gpi_mf = gpi[lev];
        const MultiFab& rho_mf = rho[lev];
        const MultiFab& uedge_mf = uedge[lev][0];
        const MultiFab& vedge_mf = uedge[lev][1];
#if (AMREX_SPACEDIM == 3)
        const MultiFab& wedge_mf = uedge[lev][2];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(vel_force_mf); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // We will also pass "validBox", which specifies the "valid" region.
            make_vel_force(lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                           BL_TO_FORTRAN_FAB(vel_force_mf[mfi]),
                           BL_TO_FORTRAN_FAB(gpi_mf[mfi]),
                           BL_TO_FORTRAN_N_3D(rho_mf[mfi],Rho),
                           BL_TO_FORTRAN_3D(uedge_mf[mfi]),
                           BL_TO_FORTRAN_3D(vedge_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                           BL_TO_FORTRAN_3D(wedge_mf[mfi]),
#endif
                           w0.dataPtr(),
                           w0_force.dataPtr(),
                           rho0.dataPtr(),
                           grav_cell.dataPtr(),
                           do_add_utilde_force);
        }
    }



    for (int lev=finest_level-1; lev>=0; --lev)
    {
        AverageDownTo(lev,vel_force,0,AMREX_SPACEDIM); // average lev+1 down to lev
    }

}
