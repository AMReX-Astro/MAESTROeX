
#include <Maestro.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write plotfile to disk
void
Maestro::DiagFile (const int step,
                   const Real t_in,
                   const Vector<Real>& rho0_in,
                   const Vector<Real>& p0_in,
                   const Vector<MultiFab>& u_in,
                   const Vector<MultiFab>& s_in,
                   int& index)
{
    if (spherical == 0) {
	Print() << "ERROR: WriteDiagFile() not written for non-spherical geometry" << std::endl;
        return;
    }

    // timer for profiling
    BL_PROFILE_VAR("Maestro::DiagFile()",DiagFile);

    // w0mac will contain an edge-centered w0 on a Cartesian grid,
    // for use in computing divergences.
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

    // w0r_cart is w0 but onto a Cartesian grid in cell-centered as
    // a scalar.  Since w0 is the radial expansion velocity, w0r_cart
    // is the radial w0 in a zone
    Vector<MultiFab> w0r_cart(finest_level+1);

    if (spherical == 1) {

        for (int lev=0; lev<=finest_level; ++lev) {
            AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                         w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                         w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                w0mac[lev][idim].setVal(0.);
            }

            w0r_cart[lev].define(grids[lev], dmap[lev], 1, 0);
            w0r_cart[lev].setVal(0.);
        }

        // put w0 on Cartesian edges as a vector
        MakeW0mac(w0mac);

        // put w0 in Cartesian cell-centers as a scalar (the radial
        // expansion velocity)
        Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0);

    } else {
        Abort("ERROR: WriteDiagFile() not supported for non-spherical geometry");
    }

    // initialize diagnosis variables
    Real T_max=0.0, T_center=0.0;
    int ncenter=0;
    Vector<Real> coord_Tmax(AMREX_SPACEDIM,0.0);
    Vector<Real> vel_Tmax(AMREX_SPACEDIM,0.0);
    Real Mach_max=0.0;
    Real Rloc_Tmax, vr_Tmax;

    for (int lev=0; lev<=finest_level; ++lev) {

        // diagnosis variables at each level
        Real T_max_local=0.0;
        Real T_center_level=0.0;
        int ncenter_level=0;
        Vector<Real> coord_Tmax_local(AMREX_SPACEDIM,0.0);
        Vector<Real> vel_Tmax_local(AMREX_SPACEDIM,0.0);
        Real Mach_max_level=0.0;

        // get references to the MultiFabs at level lev
        const MultiFab& sin_mf = s_in[lev];
        const MultiFab& uin_mf = u_in[lev];
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];
        const MultiFab& w0rcart_mf = w0r_cart[lev];
        const MultiFab& normal_mf = normal[lev];
        const Real* dx = geom[lev].CellSize();

        // create mask assuming refinement ratio = 2
        int finelev = lev+1;
        if (lev == finest_level) finelev = finest_level;

        const BoxArray& fba = s_in[finelev].boxArray();
        const iMultiFab& mask = makeFineMask(sin_mf, fba, IntVect(2));


        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(sin_mf, true); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            int use_mask = !(lev==finest_level);
            
            // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data,
            // lo/hi coordinates (including ghost cells), and/or the # of components
            // we include the mask so we don't double count; i.e., we only consider
            // cells that are not covered by finer cells
            diag_sphr(&lev, ARLIM_3D(tileBox.loVect()), ARLIM_3D(tileBox.hiVect()),
                      BL_TO_FORTRAN_FAB(sin_mf[mfi]),
                      rho0_in.dataPtr(), p0_in.dataPtr(),
                      BL_TO_FORTRAN_3D(uin_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0rcart_mf[mfi]),
                      dx,
                      BL_TO_FORTRAN_3D(normal_mf[mfi]),
                      &T_max_local, coord_Tmax_local.dataPtr(), vel_Tmax_local.dataPtr(),
                      &ncenter_level, &T_center_level, &Mach_max_level,
                      BL_TO_FORTRAN_3D(mask[mfi]), &use_mask);
        }

        // sum ncenter and T_center over all processors
        ParallelDescriptor::ReduceRealSum(T_center_level);
        ParallelDescriptor::ReduceIntSum(ncenter_level);

        // find the largest Mach number over all processors
        ParallelDescriptor::ReduceRealMax(Mach_max_level);

        // for T_max, we want to know where the hot spot is, so we do a
        // gather on the temperature and find the index corresponding to
        // the maxiumum.  We then pack the coordinates and velocities
        // into a local array and gather that to the I/O processor and
        // pick the values corresponding to the maximum.
        int nprocs = ParallelDescriptor::NProcs();
        int ioproc = ParallelDescriptor::IOProcessorNumber();

        Vector<Real> T_max_data(nprocs);

        if (nprocs == 1) {
            T_max_data[0] = T_max_local;
        } else {
            ParallelDescriptor::Gather(&T_max_local, 1, &T_max_data[0], 1, ioproc);
        }

        // determine index of max global T
        int index_max = 0;
        Real T_max_level = 0.0;
        for (int ip=0; ip<nprocs; ++ip) {
            if (T_max_data[ip] > T_max_level) {
                T_max_level = T_max_data[ip];
                index_max = ip;
            }
        }

        // T_max_coords will contain both the coordinate information and
        // the velocity information, so there are 2*dm values on each
        // proc
        Vector<Real> T_max_coords(2*AMREX_SPACEDIM);
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            T_max_coords[i] = coord_Tmax_local[i];
            T_max_coords[i+3] = vel_Tmax_local[i];
        }

        Vector<Real> T_max_coords_level(2*AMREX_SPACEDIM*nprocs);

        if (nprocs == 1) {
            for (int i=0; i<2*AMREX_SPACEDIM; ++i) {
                T_max_coords_level[i] = T_max_coords[i];
            }
        } else {
            ParallelDescriptor::Gather(&T_max_coords[0], 6, &T_max_coords_level[0], 6, ioproc);
        }

        // initialize global variables
        Vector<Real> coord_Tmax_level(AMREX_SPACEDIM);
        Vector<Real> vel_Tmax_level(AMREX_SPACEDIM);

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            coord_Tmax_level[i] = T_max_coords_level[2*AMREX_SPACEDIM*index_max+i];
            vel_Tmax_level[i] = T_max_coords_level[2*AMREX_SPACEDIM*index_max+i+3];
        }

        //
        // reduce the current level's data with the global data
        //
        if (ParallelDescriptor::IOProcessor()) {

            // if T_max_level is the new max, then copy the location as well
            if ( T_max_level > T_max ) {
                T_max = T_max_level;

                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    coord_Tmax[i] = coord_Tmax_level[i];
                    vel_Tmax[i] = vel_Tmax_level[i];
                }

                // compute center of domain
                Vector<Real> center(3, 0.0);
                const Real* probLo = geom[0].ProbLo();
                const Real* probHi = geom[0].ProbHi();

                for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                    center[idim] = 0.5*(*(probLo+idim) + *(probHi+idim));
                }

                // compute the radius of the bubble from the center
                Rloc_Tmax = sqrt( (coord_Tmax[0] - center[0])*(coord_Tmax[0] - center[0]) +
                                  (coord_Tmax[1] - center[1])*(coord_Tmax[1] - center[1]) +
                                  (coord_Tmax[2] - center[2])*(coord_Tmax[2] - center[2]) );

                // use the coordinates of the hot spot and the velocity components
                // to compute the radial velocity at the hotspot
                vr_Tmax = ((coord_Tmax[0] - center[0])/Rloc_Tmax)*vel_Tmax[0] +
                          ((coord_Tmax[1] - center[1])/Rloc_Tmax)*vel_Tmax[1] +
                          ((coord_Tmax[2] - center[2])/Rloc_Tmax)*vel_Tmax[2];
            }

            T_center = T_center + T_center_level;

            ncenter = ncenter + ncenter_level;

            if ( Mach_max_level > Mach_max ) {
                Mach_max = Mach_max_level;
            }
        }
    }

    // normalize
    if (ParallelDescriptor::IOProcessor()) {

        // for a full star ncenter should be 8 -- there are only 8 zones
        // that have a vertex at the center of the star.  For an octant,
        // ncenter should be 1
        if ( !((ncenter == 8 && !octant) ||
               (ncenter == 1 && octant)) ) {
            Abort("ERROR: ncenter invalid in Diag()");
        } else {
            T_center = T_center/ncenter;
        }
    }

    // write out diagnosis data if at initialization
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename = "diag_temp.out";
        std::ofstream diagfile;

        if (step == 0) {
            // create file after initialization
            diagfile.open(diagfilename, std::ofstream::out |
                          std::ofstream::trunc | std::ofstream::binary);

            // write variable names
            diagfile << std::setw(20) << std::left << "time";
            diagfile << std::setw(20) << std::left << "max{T}";
            diagfile << std::setw(20) << std::left << "x(max{T})";
            diagfile << std::setw(20) << std::left << "y(max{T})";
            diagfile << std::setw(20) << std::left << "z(max{T})";
            diagfile << std::setw(20) << std::left << "vx(max{T})";
            diagfile << std::setw(20) << std::left << "vy(max{T})";
            diagfile << std::setw(20) << std::left << "vz(max{T})";
            diagfile << std::setw(20) << std::left << "R(max{T})";
            diagfile << std::setw(20) << std::left << "vr(max{T})";
            diagfile << std::setw(20) << std::left << "T_center";
	    diagfile << std::setw(20) << std::left << "max{Mach}" << std::endl;

            diagfile.precision(10);
            diagfile << std::scientific;
            diagfile << std::setw(20) << std::left << t_in;
            diagfile << std::setw(20) << std::left << T_max;
            diagfile << std::setw(20) << std::left << coord_Tmax[0];
            diagfile << std::setw(20) << std::left << coord_Tmax[1];
            diagfile << std::setw(20) << std::left << coord_Tmax[2];
            diagfile << std::setw(20) << std::left << vel_Tmax[0];
            diagfile << std::setw(20) << std::left << vel_Tmax[1];
            diagfile << std::setw(20) << std::left << vel_Tmax[2];
            diagfile << std::setw(20) << std::left << Rloc_Tmax;
            diagfile << std::setw(20) << std::left << vr_Tmax;
            diagfile << std::setw(20) << std::left << T_center;
            diagfile << std::setw(20) << std::left << Mach_max << std::endl;

        } else {

            // store variable values in data array to be written later
            diagfile_data[index*12  ] = t_in;
            diagfile_data[index*12+1] = T_max;
            diagfile_data[index*12+2] = coord_Tmax[0];
            diagfile_data[index*12+3] = coord_Tmax[1];
            diagfile_data[index*12+4] = coord_Tmax[2];
            diagfile_data[index*12+5] = vel_Tmax[0];
            diagfile_data[index*12+6] = vel_Tmax[1];
            diagfile_data[index*12+7] = vel_Tmax[2];
            diagfile_data[index*12+8] = Rloc_Tmax;
            diagfile_data[index*12+9] = vr_Tmax;
            diagfile_data[index*12+10] = T_center;
	    diagfile_data[index*12+11] = Mach_max;
            index += 1;
        }
    }
}


// put together a vector of multifabs for writing
void
Maestro::WriteDiagFile (int& index)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WriteDiagFile()",WriteDiagFile);

    // time
    // T_max
    // coord_Tmax (3)
    // vel_Tmax (3)
    // Rloc_Tmax
    // vr_Tmax
    // T_center
    // Mach_max

    // write out diagnosis data
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename = "diag_temp.out";
        std::ofstream diagfile(diagfilename, std::ofstream::out |
                               std::ofstream::app | std::ofstream::binary);

        diagfile.precision(10);
        diagfile << std::scientific;
        for (int ii=0; ii<index; ++ii) {
            for (int icomp=0; icomp<12; ++icomp) {
                diagfile << std::setw(20) << std::left << diagfile_data[ii*12+icomp];
            }
            diagfile << std::endl;
        }

        // close file
        diagfile.close();

        // reset buffer array
        index = 0;
    }

}
