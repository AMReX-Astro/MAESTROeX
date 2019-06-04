
#include <Maestro.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write diagnostics files to disk
// We hold many timesteps-worth of diagnostic information in a buffer
// and output to the files only when flush_diag() is called.  This
// gives better performance on large machines with slow filesystems.
void
Maestro::DiagFile (const int step,
                   const Real t_in,
                   const RealVector& rho0_in,
                   const RealVector& p0_in,
                   const Vector<MultiFab>& u_in,
                   const Vector<MultiFab>& s_in,
                   int& index)
{
    if (spherical == 0) {
        if (ParallelDescriptor::IOProcessor()) {
            Warning("WARNING: WriteDiagFile() not written for non-spherical geometry");
        }
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

    // rho_Hnuc and rho_Hext are used to determine energy generation
    Vector<MultiFab> stemp             (finest_level+1);
    Vector<MultiFab> rho_Hext          (finest_level+1);
    Vector<MultiFab> rho_omegadot      (finest_level+1);
    Vector<MultiFab> rho_Hnuc          (finest_level+1);

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
        Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0,1);

        // compute rho_Hext and rho_Hnuc
        for (int lev=0; lev<=finest_level; ++lev) {
            stemp             [lev].define(grids[lev], dmap[lev],   Nscal, 0);
            rho_Hext          [lev].define(grids[lev], dmap[lev],       1, 0);
            rho_omegadot      [lev].define(grids[lev], dmap[lev], NumSpec, 0);
            rho_Hnuc          [lev].define(grids[lev], dmap[lev],       1, 0);
        }

        if (dt < small_dt) {
            React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, small_dt);
        } else {
            React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, dt*0.5);
        }

    } else {
        Abort("ERROR: WriteDiagFile() not supported for non-spherical geometry");
    }

    // initialize diagnosis variables
    // diag_temp.out
    Real T_max=0.0, T_center=0.0;
    int ncenter=0;
    Vector<Real> coord_Tmax(AMREX_SPACEDIM,0.0);
    Vector<Real> vel_Tmax(AMREX_SPACEDIM,0.0);
    Real Rloc_Tmax, vr_Tmax;

    // diag_vel.out
    Real U_max=0.0, Mach_max=0.0;
    Real kin_ener=0.0, int_ener=0.0;
    Vector<Real> vel_center(AMREX_SPACEDIM,0.0);

    // diag_enuc.out
    Real enuc_max=0.0;
    Vector<Real> coord_enucmax(AMREX_SPACEDIM,0.0);
    Vector<Real> vel_enucmax(AMREX_SPACEDIM,0.0);
    Real Rloc_enucmax, vr_enucmax;
    Real nuc_ener=0.0;


    for (int lev=0; lev<=finest_level; ++lev) {

        // diagnosis variables at each level
        // diag_temp.out
        Real T_max_local=0.0;
        Real T_center_level=0.0;
        int ncenter_level=0;
        Vector<Real> coord_Tmax_local(AMREX_SPACEDIM,0.0);
        Vector<Real> vel_Tmax_local(AMREX_SPACEDIM,0.0);

        // diag_vel.out
        Real U_max_level=0.0;
        Real Mach_max_level=0.0;
        Real kin_ener_level=0.0;
        Real int_ener_level=0.0;
        Vector<Real> vel_center_level(AMREX_SPACEDIM,0.0);

        // diag_enuc.out
        Real enuc_max_local=0.0;
        Vector<Real> coord_enucmax_local(AMREX_SPACEDIM,0.0);
        Vector<Real> vel_enucmax_local(AMREX_SPACEDIM,0.0);
        Real nuc_ener_level=0.0;


        // get references to the MultiFabs at level lev
        const MultiFab& sin_mf = s_in[lev];
        const MultiFab& uin_mf = u_in[lev];
        const MultiFab& w0macx_mf = w0mac[lev][0];
        const MultiFab& w0macy_mf = w0mac[lev][1];
        const MultiFab& w0macz_mf = w0mac[lev][2];
        const MultiFab& w0rcart_mf = w0r_cart[lev];
        const MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];
        const MultiFab& rho_Hext_mf = rho_Hext[lev];
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
                      BL_TO_FORTRAN_3D(rho_Hnuc_mf[mfi]),
                      BL_TO_FORTRAN_3D(rho_Hext_mf[mfi]),
                      BL_TO_FORTRAN_3D(uin_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macx_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0macz_mf[mfi]),
                      BL_TO_FORTRAN_3D(w0rcart_mf[mfi]),
                      dx,
                      BL_TO_FORTRAN_3D(normal_mf[mfi]),
                      &T_max_local, coord_Tmax_local.dataPtr(), vel_Tmax_local.dataPtr(),
                      &enuc_max_local, coord_enucmax_local.dataPtr(), vel_enucmax_local.dataPtr(),
                      &kin_ener_level, &int_ener_level, &nuc_ener_level,
                      &U_max_level, &Mach_max_level,
                      &ncenter_level, &T_center_level, vel_center_level.dataPtr(),
                      BL_TO_FORTRAN_3D(mask[mfi]), &use_mask);
        }

        // sum quantities over all processors
        ParallelDescriptor::ReduceRealSum(T_center_level);
        ParallelDescriptor::ReduceRealSum(vel_center_level.dataPtr(),AMREX_SPACEDIM);
        ParallelDescriptor::ReduceRealSum(kin_ener_level);
        ParallelDescriptor::ReduceRealSum(int_ener_level);
        ParallelDescriptor::ReduceRealSum(nuc_ener_level);
        ParallelDescriptor::ReduceIntSum(ncenter_level);

        // find the largest U and Mach number over all processors
        ParallelDescriptor::ReduceRealMax(U_max_level);
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

        // for enuc_max, we also want to know where the hot spot is, so
        // we do the same gather procedure as with the temperature
        Vector<Real> enuc_max_data(nprocs);

        if (nprocs == 1) {
            enuc_max_data[0] = enuc_max_local;
        } else {
            ParallelDescriptor::Gather(&enuc_max_local, 1, &enuc_max_data[0], 1, ioproc);
        }

        // determine index of max global enuc
        index_max = 0;
        Real enuc_max_level = 0.0;
        for (int ip=0; ip<nprocs; ++ip) {
            if (enuc_max_data[ip] > enuc_max_level) {
                enuc_max_level = enuc_max_data[ip];
                index_max = ip;
            }
        }

        Vector<Real> enuc_max_coords(2*AMREX_SPACEDIM);
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            enuc_max_coords[i] = coord_enucmax_local[i];
            enuc_max_coords[i+3] = vel_enucmax_local[i];
        }

        Vector<Real> enuc_max_coords_level(2*AMREX_SPACEDIM*nprocs);
        if (nprocs == 1) {
            for (int i=0; i<2*AMREX_SPACEDIM; ++i) {
                enuc_max_coords_level[i] = enuc_max_coords[i];
            }
        } else {
            ParallelDescriptor::Gather(&enuc_max_coords[0], 2*AMREX_SPACEDIM,
                                       &enuc_max_coords_level[0], 2*AMREX_SPACEDIM, ioproc);
        }

        // initialize global variables
        Vector<Real> coord_enucmax_level(AMREX_SPACEDIM);
        Vector<Real> vel_enucmax_level(AMREX_SPACEDIM);

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            coord_enucmax_level[i] = enuc_max_coords_level[2*AMREX_SPACEDIM*index_max+i];
            vel_enucmax_level[i] = enuc_max_coords_level[2*AMREX_SPACEDIM*index_max+i+3];
        }


        //
        // reduce the current level's data with the global data
        //
        if (ParallelDescriptor::IOProcessor()) {

            kin_ener = kin_ener + kin_ener_level;
            int_ener = int_ener + int_ener_level;
            nuc_ener = nuc_ener + nuc_ener_level;

            U_max = max(U_max, U_max_level);
            Mach_max = max(Mach_max, Mach_max_level);

            // compute center of domain
            Vector<Real> center(3, 0.0);
            const Real* probLo = geom[0].ProbLo();
            const Real* probHi = geom[0].ProbHi();

            for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                center[idim] = 0.5*(*(probLo+idim) + *(probHi+idim));
            }

            // if T_max_level is the new max, then copy the location as well
            if ( T_max_level > T_max ) {
                T_max = T_max_level;

                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    coord_Tmax[i] = coord_Tmax_level[i];
                    vel_Tmax[i] = vel_Tmax_level[i];
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

            // if enuc_max_level is the new max, then copy the location as well
            if ( enuc_max_level > enuc_max ) {
                enuc_max = enuc_max_level;

                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    coord_enucmax[i] = coord_enucmax_level[i];
                    vel_enucmax[i] = vel_enucmax_level[i];
                }

                // compute the radius of the bubble from the center
                Rloc_enucmax = sqrt( (coord_enucmax[0] - center[0])*(coord_enucmax[0] - center[0]) +
                                     (coord_enucmax[1] - center[1])*(coord_enucmax[1] - center[1]) +
                                     (coord_enucmax[2] - center[2])*(coord_enucmax[2] - center[2]) );

                // use the coordinates of the hot spot and the velocity components
                // to compute the radial velocity at the hotspot
                vr_enucmax = ((coord_enucmax[0] - center[0])/Rloc_enucmax)*vel_enucmax[0] +
                             ((coord_enucmax[1] - center[1])/Rloc_enucmax)*vel_enucmax[1] +
                             ((coord_enucmax[2] - center[2])/Rloc_enucmax)*vel_enucmax[2];

            }

            T_center = T_center + T_center_level;
            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                vel_center[i] = vel_center[i] + vel_center_level[i];
            }
            ncenter = ncenter + ncenter_level;
        }
    }

    // compute the graviational potential energy too
    Real grav_ener=0.0;
    diag_grav_energy(&grav_ener, rho0_in.dataPtr(), r_cc_loc.dataPtr(), r_edge_loc.dataPtr());

    // normalize
    if (ParallelDescriptor::IOProcessor()) {

        // the volume we normalize with is that of a single coarse-level
        // zone.  This is because the weight used in the loop over cells
        // was with reference to the coarse level
        const Real* dx = geom[0].CellSize();
        for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
            kin_ener *= *(dx+idim);
            int_ener *= *(dx+idim);
            nuc_ener *= *(dx+idim);
        }

        // for a full star ncenter should be 8 -- there are only 8 zones
        // that have a vertex at the center of the star.  For an octant,
        // ncenter should be 1
        if ( !((ncenter == 8 && !octant) ||
               (ncenter == 1 && octant)) ) {
            Abort("ERROR: ncenter invalid in Diag()");
        } else {
            T_center = T_center/ncenter;
            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                vel_center[i] = vel_center[i]/ncenter;
            }
        }
    }

    // write out diagnosis data if at initialization
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename1 = "diag_temp.out";
        std::ofstream diagfile1;
        const std::string& diagfilename2 = "diag_enuc.out";
        std::ofstream diagfile2;
        const std::string& diagfilename3 = "diag_vel.out";
        std::ofstream diagfile3;

        if (step == 0) {
            // create file after initialization
            diagfile1.open(diagfilename1, std::ofstream::out |
                           std::ofstream::trunc | std::ofstream::binary);

            // diag_temp.out
            // write variable names
            diagfile1 << std::setw(20) << std::left << "time";
            diagfile1 << std::setw(20) << std::left << "max{T}";
            diagfile1 << std::setw(20) << std::left << "x(max{T})";
            diagfile1 << std::setw(20) << std::left << "y(max{T})";
            diagfile1 << std::setw(20) << std::left << "z(max{T})";
            diagfile1 << std::setw(20) << std::left << "vx(max{T})";
            diagfile1 << std::setw(20) << std::left << "vy(max{T})";
            diagfile1 << std::setw(20) << std::left << "vz(max{T})";
            diagfile1 << std::setw(20) << std::left << "R(max{T})";
            diagfile1 << std::setw(20) << std::left << "vr(max{T})";
            diagfile1 << std::setw(20) << std::left << "T_center" << std::endl;

            // write data
            diagfile1.precision(10);
            diagfile1 << std::scientific;
            diagfile1 << std::setw(20) << std::left << t_in;
            diagfile1 << std::setw(20) << std::left << T_max;
            diagfile1 << std::setw(20) << std::left << coord_Tmax[0];
            diagfile1 << std::setw(20) << std::left << coord_Tmax[1];
            diagfile1 << std::setw(20) << std::left << coord_Tmax[2];
            diagfile1 << std::setw(20) << std::left << vel_Tmax[0];
            diagfile1 << std::setw(20) << std::left << vel_Tmax[1];
            diagfile1 << std::setw(20) << std::left << vel_Tmax[2];
            diagfile1 << std::setw(20) << std::left << Rloc_Tmax;
            diagfile1 << std::setw(20) << std::left << vr_Tmax;
            diagfile1 << std::setw(20) << std::left << T_center << std::endl;

            // close files
            diagfile1.close();

            // diag_enuc.out
            diagfile2.open(diagfilename2, std::ofstream::out |
                           std::ofstream::trunc | std::ofstream::binary);
            // write variable names
            diagfile2 << std::setw(20) << std::left << "time";
            diagfile2 << std::setw(20) << std::left << "max{enuc}";
            diagfile2 << std::setw(20) << std::left << "x(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "y(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "z(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "vx(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "vy(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "vz(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "R(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "vr(max{enuc})";
            diagfile2 << std::setw(20) << std::left << "tot nuc ener(erg/s)" << std::endl;

            // write data
            diagfile2.precision(10);
            diagfile2 << std::scientific;
            diagfile2 << std::setw(20) << std::left << t_in;
            diagfile2 << std::setw(20) << std::left << enuc_max;
            diagfile2 << std::setw(20) << std::left << coord_enucmax[0];
            diagfile2 << std::setw(20) << std::left << coord_enucmax[1];
            diagfile2 << std::setw(20) << std::left << coord_enucmax[2];
            diagfile2 << std::setw(20) << std::left << vel_enucmax[0];
            diagfile2 << std::setw(20) << std::left << vel_enucmax[1];
            diagfile2 << std::setw(20) << std::left << vel_enucmax[2];
            diagfile2 << std::setw(20) << std::left << Rloc_enucmax;
            diagfile2 << std::setw(20) << std::left << vr_enucmax;
            diagfile2 << std::setw(20) << std::left << nuc_ener << std::endl;

            // close file
            diagfile2.close();

            // diag_vel.out
            diagfile3.open(diagfilename3, std::ofstream::out |
                           std::ofstream::trunc | std::ofstream::binary);
            // write variable names
            diagfile3 << std::setw(20) << std::left << "time";
            diagfile3 << std::setw(20) << std::left << "max{U}";
            diagfile3 << std::setw(20) << std::left << "max{Mach}";
            diagfile3 << std::setw(20) << std::left << "tot kin energy";
            diagfile3 << std::setw(20) << std::left << "tot grav energy";
            diagfile3 << std::setw(20) << std::left << "tot int energy";
            diagfile3 << std::setw(20) << std::left << "velx_center";
            diagfile3 << std::setw(20) << std::left << "vely_center";
            diagfile3 << std::setw(20) << std::left << "velz_center";
            diagfile3 << std::setw(20) << std::left << "dt" << std::endl;

            // write data
            diagfile3.precision(10);
            diagfile3 << std::scientific;
            diagfile3 << std::setw(20) << std::left << t_in;
            diagfile3 << std::setw(20) << std::left << U_max;
            diagfile3 << std::setw(20) << std::left << Mach_max;
            diagfile3 << std::setw(20) << std::left << kin_ener;
            diagfile3 << std::setw(20) << std::left << grav_ener;
            diagfile3 << std::setw(20) << std::left << int_ener;
            diagfile3 << std::setw(20) << std::left << vel_center[0];
            diagfile3 << std::setw(20) << std::left << vel_center[1];
            diagfile3 << std::setw(20) << std::left << vel_center[2];
            diagfile3 << std::setw(20) << std::left << dt << std::endl;

            // close file
            diagfile3.close();

        } else {

            // store variable values in data array to be written later

            // temp
            diagfile1_data[index*11  ] = t_in;
            diagfile1_data[index*11+1] = T_max;
            diagfile1_data[index*11+2] = coord_Tmax[0];
            diagfile1_data[index*11+3] = coord_Tmax[1];
            diagfile1_data[index*11+4] = coord_Tmax[2];
            diagfile1_data[index*11+5] = vel_Tmax[0];
            diagfile1_data[index*11+6] = vel_Tmax[1];
            diagfile1_data[index*11+7] = vel_Tmax[2];
            diagfile1_data[index*11+8] = Rloc_Tmax;
            diagfile1_data[index*11+9] = vr_Tmax;
            diagfile1_data[index*11+10] = T_center;

            // enuc
            diagfile2_data[index*11  ] = t_in;
            diagfile2_data[index*11+1] = enuc_max;
            diagfile2_data[index*11+2] = coord_enucmax[0];
            diagfile2_data[index*11+3] = coord_enucmax[1];
            diagfile2_data[index*11+4] = coord_enucmax[2];
            diagfile2_data[index*11+5] = vel_enucmax[0];
            diagfile2_data[index*11+6] = vel_enucmax[1];
            diagfile2_data[index*11+7] = vel_enucmax[2];
            diagfile2_data[index*11+8] = Rloc_enucmax;
            diagfile2_data[index*11+9] = vr_enucmax;
            diagfile2_data[index*11+10] = nuc_ener;

            // vel
            diagfile3_data[index*10  ] = t_in;
            diagfile3_data[index*10+1] = U_max;
            diagfile3_data[index*10+2] = Mach_max;
            diagfile3_data[index*10+3] = kin_ener;
            diagfile3_data[index*10+4] = grav_ener;
            diagfile3_data[index*10+5] = int_ener;
            diagfile3_data[index*10+6] = vel_center[0];
            diagfile3_data[index*10+7] = vel_center[1];
            diagfile3_data[index*10+8] = vel_center[2];
            diagfile3_data[index*10+9] = dt;

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

    // write out diagnosis data
    if (ParallelDescriptor::IOProcessor()) {
        const std::string& diagfilename1 = "diag_temp.out";
        std::ofstream diagfile1(diagfilename1, std::ofstream::out |
                                std::ofstream::app | std::ofstream::binary);
        // time
        // T_max
        // coord_Tmax (3)
        // vel_Tmax (3)
        // Rloc_Tmax
        // vr_Tmax
        // T_center
        // Mach_max

        diagfile1.precision(10);
        diagfile1 << std::scientific;
        for (int ii=0; ii<index; ++ii) {
            for (int icomp=0; icomp<11; ++icomp) {
                diagfile1 << std::setw(20) << std::left << diagfile1_data[ii*11+icomp];
            }
            diagfile1 << std::endl;
        }

        // close file
        diagfile1.close();


        const std::string& diagfilename2 = "diag_enuc.out";
        std::ofstream diagfile2(diagfilename2, std::ofstream::out |
                                std::ofstream::app | std::ofstream::binary);
        // time
        // enuc_max
        // coord_enucmax (3)
        // vel_enucmax (3)
        // Rloc_enucmax
        // vr_enucmax
        // nuc_ener

        diagfile2.precision(10);
        diagfile2 << std::scientific;
        for (int ii=0; ii<index; ++ii) {
            for (int icomp=0; icomp<11; ++icomp) {
                diagfile2 << std::setw(20) << std::left << diagfile2_data[ii*11+icomp];
            }
            diagfile2 << std::endl;
        }

        // close file
        diagfile2.close();


        const std::string& diagfilename3 = "diag_vel.out";
        std::ofstream diagfile3(diagfilename3, std::ofstream::out |
                                std::ofstream::app | std::ofstream::binary);
        // time
        // U_max
        // Mach_max
        // kin_ener
        // grav_ener
        // int_ener
        // vel_center (3)
        // dt

        diagfile3.precision(10);
        diagfile3 << std::scientific;
        for (int ii=0; ii<index; ++ii) {
            for (int icomp=0; icomp<10; ++icomp) {
                diagfile3 << std::setw(20) << std::left << diagfile3_data[ii*10+icomp];
            }
            diagfile3 << std::endl;
        }

        // close file
        diagfile3.close();

        // reset buffer array
        index = 0;
    }

}
