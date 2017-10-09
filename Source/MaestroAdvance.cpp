
#include <Maestro.H>

using namespace amrex;

// advance solution to final time
void
Maestro::Evolve ()
{
    Real cur_time = t_new;
    int last_plot_file_step = 0;

    for (istep = 1; istep <= max_step && cur_time < stop_time; ++istep)
    {

        if (regrid_int > 0)  // We may need to regrid
        {
            if ( (istep-1) % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                regrid(0, cur_time);
            }
        }
    
        // wallclock time
        const Real strt_total = ParallelDescriptor::second();

        // compute time step
        ComputeDt();

        Print() << "\nTimestep " << istep << " starts with TIME = " << cur_time 
                       << " DT = " << dt << std::endl << std::endl;

        AdvanceTimeStep(cur_time);

        cur_time += dt;

        Print() << "\nTimestep " << istep << " ends with TIME = " << cur_time 
                       << " DT = " << dt << std::endl;

        // wallclock time
        Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (Verbose()) {
            Print() << "Time to advance time step: " << end_total << '\n';
        }

        if (plot_int > 0 && istep % plot_int == 0)
        {
            Print() << "\nWriting plotfile " << istep << std::endl;
            last_plot_file_step = istep;
            WritePlotFile(istep);
        }

    }

    // write a final plotfile if we haven't already
    if (plot_int > 0 && istep > last_plot_file_step)
    {
        Print() << "\nWriting plotfile " << istep-1 << std::endl;
        WritePlotFile(istep-1);
    }
}

// advance a single level for a single time step, updates flux registers
void
Maestro::AdvanceTimeStep (Real time)
{
    constexpr int num_grow = 3;

    t_old = t_new;
    t_new += dt;

    for (int lev=0; lev<=finest_level; ++lev) 
    {

        if (Verbose())
        {
            Print() << "[Level " << lev << " step " << istep+1 << "] ";
            Print() << "BEGIN ADVANCE with " << CountCells(lev) << " cells" << std::endl;
        }

        std::swap(sold[lev], snew[lev]);

        MultiFab& S_new = *snew[lev];

        const Real old_time = t_old;
        const Real new_time = t_new;
        const Real ctr_time = 0.5*(old_time+new_time);

        const Real* dx = geom[lev].CellSize();
        const Real* prob_lo = geom[lev].ProbLo();

        MultiFab fluxes[BL_SPACEDIM];
        if (do_reflux)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                BoxArray ba = grids[lev];
                ba.surroundingNodes(i);
                fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
            }
        }

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], S_new.nComp(), num_grow);
        FillPatch(lev, time, Sborder, 0, Sborder.nComp(), bcs_s);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

            for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                const FArrayBox& statein = Sborder[mfi];
                FArrayBox& stateout      =   S_new[mfi];

                // Allocate fabs for fluxes and Godunov velocities.
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                    const Box& bxtmp = surroundingNodes(bx,i);
                    flux[i].resize(bxtmp,S_new.nComp());
                    uface[i].resize(grow(bxtmp,1),1);
                }

                // compute velocities on faces (prescribed function of space and time)
                get_face_velocity(lev, ctr_time,
                                  AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                                               BL_TO_FORTRAN(uface[1]),
                                               BL_TO_FORTRAN(uface[2])),
                                  dx, prob_lo);

                // compute new state (stateout) and fluxes.
                advect(time, bx.loVect(), bx.hiVect(),
                       BL_TO_FORTRAN_3D(statein), 
                       BL_TO_FORTRAN_3D(stateout),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
                                    BL_TO_FORTRAN_3D(uface[1]),
                                    BL_TO_FORTRAN_3D(uface[2])),
                       AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]), 
                                    BL_TO_FORTRAN_3D(flux[1]), 
                                    BL_TO_FORTRAN_3D(flux[2])), 
                       dx, dt);

                if (do_reflux) {
                    for (int i = 0; i < BL_SPACEDIM ; i++) {
                        fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));	  
                    }
                }
            }
        }

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes have already been scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (do_reflux) { 
            if (lev < finest_level)
            {
                for (int i = 0; i < BL_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)   
                    flux_reg_s[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
                }	
            }
            if (lev > 0)
            {
                for (int i = 0; i < BL_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev) 
                    flux_reg_s[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
                }
            }
        }

        if (Verbose())
        {
            Print() << "[Level " << lev << " step " << istep+1 << "] ";
            Print() << "END ADVANCE" << std::endl;
        }

    } // end loop over levels

    if (Verbose())
    {
        Print() << "Synchronizing all levels" << std::endl;
    }

    // synchronize by refluxing and averagig down, starting from the finest_level-1/finest_level pair
    for (int lev=finest_level-1; lev>=0; --lev)
    {
        if (do_reflux) {
            // update lev based on coarse-fine flux mismatch
            flux_reg_s[lev+1]->Reflux(*snew[lev], 1.0, 0, 0, snew[lev]->nComp(), geom[lev]);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }

}
