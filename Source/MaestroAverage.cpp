
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// Given a multifab of data (phi), average down to a base state quantity, phibar.
// If we are in plane-parallel, the averaging is at constant height.
// If we are spherical, then the averaging is done at constant radius.

void Maestro::Average (const Vector<MultiFab>& phi,
                       RealVector& phibar,
                       int comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Average()",Average);

    for (int i=0; i<phibar.size(); ++i) {
        phibar[i] = 0.0;
    }

    if (spherical == 0) {

        // planar case

        // phibar is dimensioned to "max_radial_level" so we must mimic that for phisum
        // so we can simply swap this result with phibar
        RealVector phisum((max_radial_level+1)*nr_fine,0.0);

        // this stores how many cells there are laterally at each level
        Vector<int> ncell(max_radial_level+1);

        // loop is over the existing levels (up to finest_level)
        for (int lev=0; lev<=finest_level; ++lev) {

            // Get the index space of the domain
            const Box domainBox = geom[lev].Domain();

            // compute number of cells at any given height for each level
            if (AMREX_SPACEDIM==2) {
                ncell[lev] = domainBox.bigEnd(0)+1;
            }
            else if (AMREX_SPACEDIM==3) {
                ncell[lev] = (domainBox.bigEnd(0)+1)*(domainBox.bigEnd(1)+1);
            }

            // get references to the MultiFabs at level lev
            const MultiFab& phi_mf = phi[lev];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
            for ( MFIter mfi(phi_mf); mfi.isValid(); ++mfi )
            {

                // Get the index space of the valid region
                const Box& validBox = mfi.validbox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                average(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                        BL_TO_FORTRAN_N_3D(phi_mf[mfi],comp),
                        phisum.dataPtr());
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(max_radial_level+1)*nr_fine);

        // divide phisum by ncell so it stores "phibar"
        divide_phisum_by_ncell(phisum.dataPtr(),ncell.dataPtr());

        // swap pointers so phibar contains the computed average
        std::swap(phisum,phibar);

    }
    else if (spherical == 1 && use_exact_base_state) {
        // spherical case with uneven base state spacing

        // phibar is dimensioned to "max_radial_level" so we must mimic that for phisum
        // so we can simply swap this result with phibar
        RealVector phisum((max_radial_level+1)*nr_fine,0.0);

        // this stores how many cells there are at each level
        Vector<int> ncell((max_radial_level+1)*nr_fine,0);

        // loop is over the existing levels (up to finest_level)
        for (int lev=0; lev<=finest_level; ++lev) {

            // get references to the MultiFabs at level lev
            const MultiFab& phi_mf = phi[lev];
            const MultiFab& cc_to_r = cell_cc_to_r[lev];

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
            for ( MFIter mfi(phi_mf); mfi.isValid(); ++mfi )
            {

                // Get the index space of the valid region
                const Box& validBox = mfi.validbox();

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                average_sphr_irreg(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                   BL_TO_FORTRAN_N_3D(phi_mf[mfi],comp),
                                   phisum.dataPtr(), ncell.dataPtr(),
                                   BL_TO_FORTRAN_3D(cc_to_r[mfi]));
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(max_radial_level+1)*nr_fine);
        ParallelDescriptor::ReduceIntSum(ncell.dataPtr(),(max_radial_level+1)*nr_fine);

        // divide phisum by ncell so it stores "phibar"
        divide_phisum_by_ncell_irreg(phisum.dataPtr(),ncell.dataPtr());

        // swap pointers so phibar contains the computed average
        std::swap(phisum,phibar);
    }
    else {
        // spherical case with even base state spacing

        // For spherical, we construct a 1D array at each level, phisum, that has space
        // allocated for every possible radius that a cell-center at each level can
        // map into.  The radial locations have been precomputed and stored in radii.
        Vector<Real> phisum((finest_level+1)*(nr_irreg+2),0.0);
        Vector<Real>  radii((finest_level+1)*(nr_irreg+3));
        Vector<int> ncell((finest_level+1)*(nr_irreg+2),0);

        // radii contains every possible distance that a cell-center at the finest
        // level can map into
        for (int lev=0; lev<=finest_level; ++lev) {

            // Get the index space of the domain
            const Real* dx = geom[lev].CellSize();

            compute_radii_sphr(&lev, radii.dataPtr(), &finest_level, dx);
        }

        // loop is over the existing levels (up to finest_level)
        for (int lev=finest_level; lev>=0; --lev) {

            // Get the grid size of the domain
            const Real* dx = geom[lev].CellSize();

            // get references to the MultiFabs at level lev
            const MultiFab& phi_mf = phi[lev];

            // create mask assuming refinement ratio = 2
            int finelev = lev+1;
            if (lev == finest_level) finelev = finest_level;

            const BoxArray& fba = phi[finelev].boxArray();
            const iMultiFab& mask = makeFineMask(phi_mf, fba, IntVect(2));


            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
            for ( MFIter mfi(phi_mf); mfi.isValid(); ++mfi )
            {

                // Get the index space of the valid region
                const Box& validBox = mfi.validbox();

                int use_mask = !(lev==finest_level);

                // call fortran subroutine
                // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
                // lo/hi coordinates (including ghost cells), and/or the # of components
                // We will also pass "validBox", which specifies the "valid" region.
                // we include the mask so we don't double count; i.e., we only consider
                // cells that are not covered by finer cells when constructing the sum
                sum_phi_3d_sphr(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                                BL_TO_FORTRAN_N_3D(phi_mf[mfi],comp),
                                phisum.dataPtr(), radii.dataPtr(), &finest_level,
                                dx, ncell.dataPtr(),
                                BL_TO_FORTRAN_3D(mask[mfi]), &use_mask);
            }
        }

        // reduction over boxes to get sum
        ParallelDescriptor::ReduceRealSum(phisum.dataPtr(),(finest_level+1)*(nr_irreg+2));
        ParallelDescriptor::ReduceIntSum(ncell.dataPtr(),(finest_level+1)*(nr_irreg+2));

        average_sphr(phisum.dataPtr(),phibar.dataPtr(),ncell.dataPtr(),radii.dataPtr(),&finest_level);

    }


}
