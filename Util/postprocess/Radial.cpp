#include <Radial.H>
#include <Maestro.H>
//#include <Maestro_F.H>
//#include <iterator>     // std::istream_iterator
//#include <unistd.h>     // getcwd

using namespace amrex;


// write radial plotfile to disk
void
WriteRadialFile (const std::string& plotfilename,
		 const BaseState<Real>& rho0_in,
		 const BaseState<Real>& rhoh0_in,
		 const BaseState<Real>& p0_in,
		 const Vector<MultiFab>& u_in,
		 const Vector<MultiFab>& w0_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WritePlotFile()", WritePlotFile);

    std::string radialfilename = "radial_" + plotfilename;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    // write out the cell-centered diagnostics
    if (ParallelDescriptor::IOProcessor()) {

        for (int lev=0; lev<=base_geom.max_radial_level; ++lev) {

            std::ofstream RadialFile;
            RadialFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            std::string levStr = std::to_string(lev);
            radialfilename.append(levStr);
            RadialFile.open(radialfilename.c_str(), std::ofstream::out   |
                            std::ofstream::trunc |
                            std::ofstream::binary);
            if(!RadialFile.good()) {
                amrex::FileOpenFailed(radialfilename);
            }

            RadialFile.precision(17);

            RadialFile << "r_cc  rho0  p0  convect_vel  |N| \n";

            for (int i=0; i<base_geom.nr(lev); ++i) {
                RadialFile << base_geom.r_cc_loc(lev,i) << " "
                           << rho0_in.array()(lev,i) << " "
                           << p0_in.array()(lev,i) /*<< " "
                           << convectvel.array()(lev,i) << " "
                           << Nfreq.array()(lev,i)*/ << "\n";
            }
        }
    }

}

void
MakeRadialNFreq (const BaseState<Real>& p0_s,
		 const BaseState<Real>& rho0_s, 
		 BaseState<Real>& freq0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRadialNFreq()",MakeRadialNFreq);

    BaseState<Real> entropy_s(1, base_geom.nr_fine);
    BaseState<Real> grav_s(1, base_geom.nr_fine);
    // need to compute gravity analogous to Maestro::MakeGravCell(grav_s,rho0_s);
    
    const auto& r_cc_loc = base_geom.r_cc_loc;
    const auto rho0 = rho0_s.const_array();
    const auto p0 = p0_s.const_array();
    const auto grav = grav_s.const_array();
    auto entropy = entropy_s.array();
    Real gamma;
    
    // dimensionless entropy = 1/(gam - 1)*( log(p) - gamma*log(rho) )
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
	eos_t eos_state;
	
	eos_state.rho = rho0(0,r);
	eos_state.p     = p0(0,r);
	for (auto comp = 0; comp < NumSpec; ++comp) {
	    eos_state.xn[comp] = 1.0;
	}
	eos(eos_input_rp, eos_state);
	gamma = eos_state.gam1;
    
	// printf("HACK, r = %d, p0 = %f, rho0 = %f\n", r, p0(0,r), rho0(0,r));
	entropy(0,r) = 1.0/(gamma-1.0) * (log(p0(0,r)) - gamma*log(rho0(0,r)));
    }
    
    auto freq = freq0.array();
    for (auto r = 0; r < base_geom.nr_fine; ++r) {
	if (r == 0) {
	    freq(0,r) = -(gamma-1.0)/gamma*grav(0,r)*(entropy(0,r+1) - entropy(0,r))/(r_cc_loc(0,r+1)-r_cc_loc(0,r)); 
	} else if (r < base_geom.base_cutoff_density_coord(0)-1) {
	    freq(0,r) = -(gamma-1.0)/gamma*grav(0,r)*(entropy(0,r+1) - entropy(0,r-1))/(r_cc_loc(0,r+1)-r_cc_loc(0,r-1));
	} else {
	    freq(0,r) = 0.0;
	}
	freq(0,r) = std::sqrt(fabs(freq(0,r)));
    }
}

/*
void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

    for (int lev=0; lev<=finest_level; ++lev) {

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> vel_arr = vel[lev].array(mfi);
            const Array4<Real> radvel_arr = rad_vel[lev].array(mfi);
            const Array4<Real> circvel_arr = circ_vel[lev].array(mfi);
            const Array4<const Real> w0rcart_arr = w0rcart[lev].array(mfi);
            const Array4<const Real> normal_arr = normal[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                circvel_arr(i,j,k) = 0.0;
                radvel_arr(i,j,k) = 0.0;

                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    radvel_arr(i,j,k) += vel_arr(i,j,k,n) * normal_arr(i,j,k,n);
                }
                
                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    Real circ_comp = vel_arr(i,j,k,n) - radvel_arr(i,j,k) * normal_arr(i,j,k,n);
                    circvel_arr(i,j,k) += circ_comp * circ_comp;
                }

                circvel_arr(i,j,k) = sqrt(circvel_arr(i,j,k));

                // add base state vel to get full radial velocity
                radvel_arr(i,j,k) += w0rcart_arr(i,j,k);
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(rad_vel, 0, 1);
    FillPatch(t_old, rad_vel, rad_vel, rad_vel, 0, 0, 1, 0, bcs_f);
    AverageDown(circ_vel, 0, 1);
    FillPatch(t_old, circ_vel, circ_vel, circ_vel, 0, 0, 1, 0, bcs_f);
}

void
Maestro::MakeVorticity (const Vector<MultiFab>& vel,
                        Vector<MultiFab>& vorticity)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeVorticity()",MakeVorticity);

    for (int lev=0; lev<=finest_level; ++lev) {

        // get references to the MultiFabs at level lev
        const MultiFab& vel_mf = vel[lev];

        const Real* dx = geom[lev].CellSize();
        const Box& domainBox = geom[lev].Domain();

        const Real hx = dx[0];
        const Real hy = dx[1];
#if (AMREX_SPACEDIM == 3)
        const Real hz = dx[2];
#endif
        const int ilo = domainBox.loVect()[0];
        const int ihi = domainBox.hiVect()[0];
        const int jlo = domainBox.loVect()[1];
        const int jhi = domainBox.hiVect()[1];
#if (AMREX_SPACEDIM == 3)
        const int klo = domainBox.loVect()[2];
        const int khi = domainBox.hiVect()[2];
#endif

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            Array4<const Real> const u = vel[lev].array(mfi);
            Array4<Real> const vort = vorticity[lev].array(mfi);
            GpuArray<int,AMREX_SPACEDIM*2> physbc;
            for (int n = 0; n < AMREX_SPACEDIM*2; ++n) {
                physbc[n] = phys_bc[n];
            } 

#if (AMREX_SPACEDIM == 2)

            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;

                if (i == ilo && 
                    (physbc[0] == Inflow || 
                     physbc[0] == SlipWall || 
                     physbc[0] == NoSlipWall)) 
                {
                    vx = (u(i+1,j,k,1) + 3.0*u(i,j,k,1) - 
                          4.0*u(i-1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;

                } else if (i == ihi+1 &&
                        (physbc[AMREX_SPACEDIM] == Inflow || 
                        physbc[AMREX_SPACEDIM] == SlipWall || 
                        physbc[AMREX_SPACEDIM] == NoSlipWall))
                {
                    vx = -(u(i-1,j,k,1) + 3.0*u(i,j,k,1) - 
                         4.0*u(i+1,j,k,1)) / hx;
                    uy = 0.5 * (u(i,j+1,k,0) - u(i,j-1,k,0)) / hy;
                }

                if (j == jlo &&
                    (physbc[1] == Inflow || 
                     physbc[1] == SlipWall || 
                     physbc[1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,0)) / hx;
                    uy = (u(i,j+1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j-1,k,0)) / hy;

                } else if (j == jhi+1 && 
                           (physbc[AMREX_SPACEDIM+1] == Inflow || 
                            physbc[AMREX_SPACEDIM+1] == SlipWall || 
                            physbc[AMREX_SPACEDIM+1] == NoSlipWall))
                {
                    vx = 0.5 * (u(i+1,j,k,1) - u(i-1,j,k,1)) / hx;
                    uy = -(u(i,j-1,k,0) + 3.0*u(i,j,k,0) - 
                         4.0*u(i,j+1,k,0)) / hy;
                }

                vort(i,j,k) = vx - uy;
            });

#else 
            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k,
            {
                Real uy = 0.5*(u(i,j+1,k,0)-u(i,j-1,k,0))/hy;
                Real uz = 0.5*(u(i,j,k+1,0)-u(i,j,k-1,0))/hz;
                Real vx = 0.5*(u(i+1,j,k,1)-u(i-1,j,k,1))/hx;
                Real vz = 0.5*(u(i,j,k+1,1)-u(i,j,k-1,1))/hz;
                Real wx = 0.5*(u(i+1,j,k,2)-u(i-1,j,k,2))/hx;
                Real wy = 0.5*(u(i,j+1,k,2)-u(i,j-1,k,2))/hy;

                bool fix_lo_x = (physbc[0] == Inflow || 
                                 physbc[0] == NoSlipWall);
                bool fix_hi_x = (physbc[AMREX_SPACEDIM] == Inflow || 
                                 physbc[AMREX_SPACEDIM] == NoSlipWall);

                bool fix_lo_y = (physbc[1] == Inflow || 
                                 physbc[1] == NoSlipWall);
                bool fix_hi_y = (physbc[AMREX_SPACEDIM+1] == Inflow ||
                                 physbc[AMREX_SPACEDIM+1] == NoSlipWall);

                bool fix_lo_z = (physbc[2] == Inflow || 
                                 physbc[2] == NoSlipWall);
                bool fix_hi_z = (physbc[AMREX_SPACEDIM+2] == Inflow ||
                                 physbc[AMREX_SPACEDIM+2] == NoSlipWall);

                // First do all the faces
                if (fix_lo_x && i == ilo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                } else if (fix_hi_x && i == ihi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                }

                if (fix_lo_y && j == jlo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                } else if (fix_hi_y && j == jhi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_z && k == klo) {
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_z && k == khi+1) {
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                // Next do all the edges
                if (fix_lo_x && fix_lo_y && i == ilo && j == jlo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_hi_x && fix_lo_y && i == ihi+1 && j == jlo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_hi_y && i == ilo && j == jhi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                }

                if (fix_lo_x && fix_lo_z && i == ilo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_z && i == ihi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_z && i == ilo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_z && i == ihi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_lo_z && j == jlo && k == klo) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_lo_z && j == jhi+1 && k == klo) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_y && fix_hi_z && j == jlo && k == khi+1) {
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_y && fix_hi_z && j == jhi+1 && k == khi+1) {
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }
                
                // Finally do all the corners
                if (fix_lo_x && fix_lo_y && fix_lo_z && 
                    i == ilo && j == jlo && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_lo_z &&
                    i == ihi+1 && j == jlo && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_lo_z &&
                    i == ilo && j == jhi+1 && k == klo) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_lo_z &&
                    i == ihi+1 && j == jhi+1 && k == klo) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = (u(i,j,k+1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k-1,0))/(3.0*hz);
                    vz = (u(i,j,k+1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k-1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_lo_y && fix_hi_z &&
                    i == ilo && j == jlo && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_lo_y && fix_hi_z &&
                    i == ihi+1 && j == jlo && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = (u(i,j+1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j-1,k,0))/(3.0*hy);
                    wy = (u(i,j+1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j-1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_lo_x && fix_hi_y && fix_hi_z &&
                    i == ilo && j == jhi+1 && k == khi+1) {
                    vx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    wx = (u(i+1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i-1,j,k,1))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                if (fix_hi_x && fix_hi_y && fix_hi_z &&
                    i == ihi+1 && j == jhi+1 && k == khi+1) {
                    vx = -(u(i-1,j,k,1)+3.0*u(i,j,k,1)-4.0*u(i+1,j,k,1))/(3.0*hx);
                    wx = -(u(i-1,j,k,2)+3.0*u(i,j,k,2)-4.0*u(i+1,j,k,2))/(3.0*hx);
                    uy = -(u(i,j-1,k,0)+3.0*u(i,j,k,0)-4.0*u(i,j+1,k,0))/(3.0*hy);
                    wy = -(u(i,j-1,k,2)+3.0*u(i,j,k,2)-4.0*u(i,j+1,k,2))/(3.0*hy);
                    uz = -(u(i,j,k-1,0)+3.0*u(i,j,k,0)-4.0*u(i,j,k+1,0))/(3.0*hz);
                    vz = -(u(i,j,k-1,1)+3.0*u(i,j,k,1)-4.0*u(i,j,k+1,1))/(3.0*hz);
                }

                vort(i,j,k) = sqrt((wy-vz)*(wy-vz)+
                    (uz-wx)*(uz-wx)+(vx-uy)*(vx-uy));
            });
#endif
        }
    }

    // average down and fill ghost cells
    AverageDown(vorticity,0,1);
    FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeEntropy (const Vector<MultiFab>& state,
                      Vector<MultiFab>& entropy)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeEntropy()",MakeEntropy);

    for (int lev=0; lev<=finest_level; ++lev) {

        // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> state_arr = state[lev].array(mfi);
            const Array4<Real> entropy_arr = entropy[lev].array(mfi);

            AMREX_PARALLEL_FOR_3D(tileBox, i, j ,k, {
                eos_t eos_state;

                eos_state.rho = state_arr(i,j,k,Rho);
                eos_state.T = state_arr(i,j,k,Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = state_arr(i,j,k,FirstSpec+comp) / state_arr(i,j,k,Rho);
                }

                eos(eos_input_rt, eos_state);

                entropy_arr(i,j,k) = eos_state.s;
            });
        }
    }

    // average down and fill ghost cells
    AverageDown(entropy,0,1);
    FillPatch(t_old,entropy,entropy,entropy,0,0,1,0,bcs_f);
}

void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeDivw0()", MakeDivw0);

    for (int lev=0; lev<=finest_level; ++lev) {

        if (!spherical) {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0_arr = w0_cart[lev].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
#if (AMREX_SPACEDIM == 2)
                    divw0_arr(i,j,k) = (w0_arr(i,j+1,k,1) - w0_arr(i,j,k,1)) / dx[1];
#else
                    divw0_arr(i,j,k) = (w0_arr(i,j,k+1,2) - w0_arr(i,j,k,2)) / dx[2];
#endif
                });
            }

        } else {

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(divw0[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                // Get the index space of the valid region
                const Box& tileBox = mfi.tilebox();
                const auto dx = geom[lev].CellSizeArray();
                
                const Array4<const Real> w0macx = w0mac[lev][0].array(mfi);
                const Array4<const Real> w0macy = w0mac[lev][1].array(mfi);
                const Array4<const Real> w0macz = w0mac[lev][2].array(mfi);
                const Array4<Real> divw0_arr = divw0[lev].array(mfi);

                AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                    divw0_arr(i,j,k) = (w0macx(i+1,j,k) - w0macx(i,j,k)) / dx[0] + 
                        (w0macy(i,j+1,k) - w0macy(i,j,k)) / dx[1] + 
                        (w0macz(i,j,k+1) - w0macz(i,j,k)) / dx[2];
                });
            }
        }
    }

    // average down and fill ghost cells
    AverageDown(divw0, 0, 1);
    FillPatch(t_old, divw0, divw0, divw0, 0, 0, 1, 0, bcs_f);
}
*/
