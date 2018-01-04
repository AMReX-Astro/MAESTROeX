
#include <Maestro.H>
#include <AMReX_VisMF.H>
using namespace amrex;


// initialize AMR data
// perform initial projection
// perform divu iters
// perform initial (pressure) iterations
void
Maestro::Init ()
{
    Print() << "Calling Init()" << endl;

    // fill in multifab and base state data
    InitData();

    if (spherical == 1) {
        // FIXME
        // MakeNormal();
    }

    // make gravity
    make_grav_cell(grav_cell_old.dataPtr(),
                   rho0_old.dataPtr(),
                   r_cc_loc.dataPtr(),
                   r_edge_loc.dataPtr());

    // compute gamma1bar
    MakeGamma1bar(sold,gamma1bar_old,p0_old);

    // compute beta0
    make_beta0(beta0_old.dataPtr(),
               rho0_old.dataPtr(),
               p0_old.dataPtr(),
               gamma1bar_old.dataPtr(),
               grav_cell_old.dataPtr());

    // initial projection
    if (do_initial_projection) {
        Print() << "Doing initial projection" << endl;
        InitProj();
    }

    // compute initial time step
    FirstDt();

    // divu iters
    for (int i=1; i<=init_divu_iter; ++i) {
        Print() << "Doing initial divu iteration #" << i << endl;
        DivuIter(i);
    }

    // copy S_cc_old into S_cc_new for the pressure iterations
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Copy(S_cc_new[lev],S_cc_old[lev],0,0,1,0);
    }

    // initial (pressure) iters
    for (int i=1; i<= init_iter; ++i) {
        Print() << "Doing initial pressure iteration #" << i << endl;
        InitIter();
    }
}

// fill in multifab and base state data
void
Maestro::InitData ()
{
    Print() << "Calling InitData()" << endl;

    // read in model file and fill in s0_init and p0_init for all levels
    init_base_state(s0_init.dataPtr(),p0_init.dataPtr(),rho0_old.dataPtr(),
                    rhoh0_old.dataPtr(),p0_old.dataPtr(),tempbar.dataPtr(),max_level+1);

    // calls AmrCore::InitFromScratch(), which calls a MakeNewGrids() function 
    // that repeatedly calls Maestro::MakeNewLevelFromScratch() to build and initialize
    InitFromScratch(t_old);

    // set finest_radial_level in fortran
    // compute numdisjointchunks, r_start_coord, r_end_coord
    init_multilevel(finest_level);

    // synchronize levels
    AverageDown(sold,0,Nscal);
    AverageDown(uold,0,AMREX_SPACEDIM);

    // free memory in s0_init and p0_init by swapping it
    // with an empty vector that will go out of scope
    Vector<Real> s0_swap, p0_swap;
    std::swap(s0_swap,s0_init);
    std::swap(p0_swap,p0_init);

    if (fix_base_state) {
        // compute cutoff coordinates
        compute_cutoff_coords(rho0_old.dataPtr());
        make_grav_cell(grav_cell_old.dataPtr(),
                       rho0_old.dataPtr(),
                       r_cc_loc.dataPtr(),
                       r_edge_loc.dataPtr());
    }
    else {
        if (do_smallscale) {
            // first compute cutoff coordinates using initial density profile
            compute_cutoff_coords(rho0_old.dataPtr());
            // set rho0_old = rhoh0_old = 0.
            std::fill(rho0_old.begin(),  rho0_old.end(),  0.);
            std::fill(rhoh0_old.begin(), rhoh0_old.end(), 0.);
        }
        else {
            // set rho0 to be the average
            Average(sold,rho0_old,Rho);
            compute_cutoff_coords(rho0_old.dataPtr());

            // compute gravity
            make_grav_cell(grav_cell_old.dataPtr(),
                           rho0_old.dataPtr(),
                           r_cc_loc.dataPtr(),
                           r_edge_loc.dataPtr());

            // compute p0 with HSE
            enforce_HSE(rho0_old.dataPtr(),
                        p0_old.dataPtr(),
                        grav_cell_old.dataPtr(),
                        r_edge_loc.dataPtr());

            // call eos with r,p as input to recompute T,h
            TfromRhoP(sold,p0_old,1);

            // set rhoh0 to be the average
            Average(sold,rhoh0_old,RhoH);
        }
    }

    if (plot_int > 0) {
        Print() << "\nWriting plotfile 0 after initialization" << endl;
        WritePlotFile(0,t_old,rho0_old,p0_old,uold,sold);
    }
}

// During initialization of a simulation, Maestro::InitData() calls 
// AmrCore::InitFromScratch(), which calls 
// a MakeNewGrids() function that repeatedly calls this function to build 
// and initialize finer levels.  This function creates a new fine
// level that did not exist before by interpolating from the coarser level
// overrides the pure virtual function in AmrCore
void Maestro::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				       const DistributionMapping& dm)
{
    sold    [lev].define(ba, dm,          Nscal, 0);
    snew    [lev].define(ba, dm,          Nscal, 0);
    uold    [lev].define(ba, dm, AMREX_SPACEDIM, 0);
    unew    [lev].define(ba, dm, AMREX_SPACEDIM, 0);
    S_cc_old[lev].define(ba, dm,              1, 0);
    S_cc_new[lev].define(ba, dm,              1, 0);
    gpi     [lev].define(ba, dm, AMREX_SPACEDIM, 0);
    dSdt    [lev].define(ba, dm,              1, 0);
    if (spherical == 1) {
        normal[lev].define(ba, dm, 1, 1);
    }

    if (lev > 0 && do_reflux) {
        flux_reg_s[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, Nscal));
        flux_reg_u[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, AMREX_SPACEDIM));
    }

    const Real* dx = geom[lev].CellSize();

    MultiFab& scal = sold[lev];
    MultiFab& vel = uold[lev];

    // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for (MFIter mfi(scal); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        initdata(lev, t_old, ARLIM_3D(lo), ARLIM_3D(hi),
                 BL_TO_FORTRAN_FAB(scal[mfi]), 
                 BL_TO_FORTRAN_FAB(vel[mfi]), 
                 s0_init.dataPtr(), p0_init.dataPtr(),
                 ZFILL(dx));
    }
}


void Maestro::InitProj ()
{

    Vector<MultiFab> rho_omegadot(finest_level+1);
    Vector<MultiFab>      thermal(finest_level+1);
    Vector<MultiFab>     rho_Hnuc(finest_level+1);
    Vector<MultiFab>     rho_Hext(finest_level+1);
    Vector<MultiFab>      rhohalf(finest_level+1);
    Vector<MultiFab>         rhcc(finest_level+1);

    Vector<Real> Sbar( (max_radial_level+1)*nr_fine );

    for (int lev=0; lev<=finest_level; ++lev) {
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        thermal     [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_Hext    [lev].define(grids[lev], dmap[lev],       1, 0);
        rhohalf     [lev].define(grids[lev], dmap[lev],       1, 1);
        rhcc        [lev].define(grids[lev], dmap[lev],       1, 1);

        // we don't have a legit timestep yet, so we set rho_omegadot,
        // rho_Hnuc, and rho_Hext to 0 
        rho_omegadot[lev].setVal(0.);
        rho_Hnuc[lev].setVal(0.);
        rho_Hext[lev].setVal(0.);

        // initial projection does not use density weighting
        rhohalf[lev].setVal(1.);
    }

    // compute thermal diffusion
    if (use_thermal_diffusion) {
        Abort("InitProj: use_thermal_diffusion not implemented");
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            thermal[lev].setVal(0.);
        }
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_old,sold,rho_omegadot,rho_Hnuc,rho_Hext,thermal);

    if (evolve_base_state) {
        // average S into Sbar
        Average(S_cc_old,Sbar,0);
    }

    // make the nodal rhs for projection
    MakeCCRHSforNodalProj(rhcc,S_cc_old,Sbar,beta0_old);

    // define epsilon for initial projection
    //

    // perform a nodal projection
    NodalProj(initial_projection_comp,rhcc);

}


void Maestro::DivuIter (int istep_divu_iter)
{

    Vector<MultiFab> stemp       (finest_level+1);
    Vector<MultiFab> rho_Hext    (finest_level+1);
    Vector<MultiFab> rho_omegadot(finest_level+1);
    Vector<MultiFab> rho_Hnuc    (finest_level+1);
    Vector<MultiFab> thermal     (finest_level+1);
    Vector<MultiFab> rhohalf     (finest_level+1);
    Vector<MultiFab> rhcc        (finest_level+1);

    Vector<Real> Sbar                ( (max_radial_level+1)*nr_fine );
    Vector<Real> w0_force            ( (max_radial_level+1)*nr_fine );
    Vector<Real> p0_minus_peosbar    ( (max_radial_level+1)*nr_fine );
    Vector<Real> delta_chi_w0        ( (max_radial_level+1)*nr_fine );

    std::fill(etarho_ec.begin(),            etarho_ec.end(),            0.);
    std::fill(Sbar.begin(),                 Sbar.end(),                 0.);
    std::fill(w0_force.begin(),             w0_force.end(),             0.);
    std::fill(psi.begin(),                  psi.end(),                  0.);
    std::fill(etarho_cc.begin(),            etarho_cc.end(),            0.);
    std::fill(p0_minus_peosbar.begin(),     p0_minus_peosbar.end(),     0.);


    for (int lev=0; lev<=finest_level; ++lev) {
        stemp       [lev].define(grids[lev], dmap[lev],   Nscal, 0);
        rho_Hext    [lev].define(grids[lev], dmap[lev],       1, 0);
        rho_omegadot[lev].define(grids[lev], dmap[lev], NumSpec, 0);
        rho_Hnuc    [lev].define(grids[lev], dmap[lev],       1, 0);
        thermal     [lev].define(grids[lev], dmap[lev],       1, 0);
        rhohalf     [lev].define(grids[lev], dmap[lev],       1, 1);
        rhcc        [lev].define(grids[lev], dmap[lev],       1, 0);

        // divu_iters do not use density weighting
        rhohalf[lev].setVal(1.);
    }

    React(sold,stemp,rho_Hext,rho_omegadot,rho_Hnuc,p0_old,0.5*dt);

    if (use_thermal_diffusion) {
        Abort("DivuIter: use_thermal_diffusion not implemented");
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            thermal[lev].setVal(0.);
        }        
    }

    // compute S at cell-centers
    Make_S_cc(S_cc_old,sold,rho_omegadot,rho_Hnuc,rho_Hext,thermal);

    if (evolve_base_state) {
        Average(S_cc_old,Sbar,0);

        make_w0(w0.dataPtr(), w0.dataPtr(), w0_force.dataPtr() ,Sbar.dataPtr(),
                rho0_old.dataPtr(), rho0_new.dataPtr(), p0_old.dataPtr(), 
                p0_new.dataPtr(), gamma1bar_old.dataPtr(), gamma1bar_new.dataPtr(),
                p0_minus_peosbar.dataPtr(), psi.dataPtr(), etarho_ec.dataPtr(),
                etarho_cc.dataPtr(), delta_chi_w0.dataPtr(), r_cc_loc.dataPtr(),
                r_edge_loc.dataPtr(), dt, dt, 0);
    }

    MakeCCRHSforNodalProj(rhcc,S_cc_old,Sbar,beta0_old);

    // define epsilon for divu_iters
    //

    // perform a nodal projection
    //

    // compute new dt



}


void Maestro::InitIter ()
{
    dtold = dt;

    // advance the solution by dt
    AdvanceTimeStep(true);
}
