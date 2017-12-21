#ifndef Maestro_H_
#define Maestro_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FMultiGrid.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <PhysBCFunctMaestro.H>
#include <Maestro_F.H>

class Maestro
    : public amrex::AmrCore
{
public:

    /*
      public member functions
    */

    // constructor and destructor
    Maestro ();
    virtual ~Maestro ();

    // in MaestroSetup.cpp
    // read in C++/F90 parameters
    // define global C++/F90 variables and initialize network
    // set up boundary conditions
    // initialize base state geometry parameters
    // set istep, t_new, t_old
    // allocate MultiFabs and base state arrays
    void Setup ();

    // in MaestroInit.cpp
    // initialize multifab and base state data
    // perform initial projection
    // perform divu iters
    // perform initial (pressure) iterations
    void Init ();

    // in MaestroEvolve.cpp
    // advance solution to final time
    void Evolve ();

private:

    /*
      private member functions
    */

    ////////////
    // MaestroAdvance.cpp functions

    // advance solution at all levels for a single time step
    void AdvanceTimeStep (bool is_initIter);

    // end MaestroAdvance.cpp functions
    ////////////

    ////////////
    // MaestroAdvection.cpp functions
    void AdvancePremac (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);
    ////////////

    ////////////
    // MaestroAverage.cpp functions

    // compute the radial average of a quantitiy
    void Average (const amrex::Vector<amrex::MultiFab>& mf,
                  amrex::Vector<amrex::Real>& phibar,
                  int comp);

    // end MaestroAverage.cpp functions
    ////////////

    ////////////
    // Time step computation in MaestroDt.cpp

    // compute the time step
    void EstDt ();

    // compute initial time step
    void FirstDt ();

    // end MaestroDt.cpp functions
    ////////////////////////

    ////////////////////////
    // MaestroFillData.cpp functions

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases 
    // (fill fine grid ghost by interpolating from coarse)
    // scomp of the source component
    // dcomp is the destination component AND the bc component
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int scomp, int dcomp, int ncomp,
                    amrex::Vector<amrex::BCRec> bcs);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    // scomp of the source component
    // dcomp is the destination component AND the bc component
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<amrex::MultiFab>& mf_old,
                          amrex::Vector<amrex::MultiFab>& mf_new,
                          int scomp, int dcomp, int ncomp,
                          amrex::Vector<amrex::BCRec> bcs);

    // utility to copy in data from mf_old and/or mf_new into mf
    // if time=t_old we copy mf_old into mf
    // if time=t_new we copy mf_new into mf
    // otherwise copy copy in both mf_old and mf_new into mf and the fillpatch
    // routines know to interpolate in time.  However in MAESTRO since we don't
    // subcycle I'm not sure if we need this capability?
    void GetData (int lev, amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& mf,
                  amrex::Vector<amrex::Real>& mftime,
                  amrex::Vector<amrex::MultiFab>& mf_old,
                  amrex::Vector<amrex::MultiFab>& mf_new);

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown (amrex::Vector<amrex::MultiFab>& mf,
                      int comp,
                      int ncomp);

    // more flexible version of AverageDown() that lets you average 
    // down across multiple levels
    void AverageDownTo (int crse_lev,
                        amrex::Vector<amrex::MultiFab>& mf,
                        int comp,
                        int ncomp);

    // fill in ONE ghost cell for all components of a face-centered (MAC) velocity
    // field behind physical boundaries.  Does not modify the velocities on the boundary
    void FillUmacGhost (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);

    // end MaestroFillData.cpp functions
    ////////////

    ////////////////////////
    // MaestroFill3dData.cpp functions
    
    void Put1dArrayOnCart (const amrex::Vector<amrex::Real>& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           amrex::Vector<amrex::BCRec>& bcs,
                           int is_input_edge_centered,
                           int is_output_a_vector);

    // end MaestroFill3dData.cpp functions
    ////////////

    ////////////////////////
    // MaestroForce.cpp functions
    void MakeVelForce (amrex::Vector<amrex::MultiFab>& vel_force,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
                       int is_final_update,
                       int do_add_utilde_force);
    ////////////////////////

    ////////////
    // MaestroGamma.cpp functions

    void MakeGamma1bar (const amrex::Vector<amrex::MultiFab>& scal,
                        amrex::Vector<amrex::Real>& gamma1bar,
                        const amrex::Vector<amrex::Real>& p0);

    // end MaestroGamma.cpp functions
    ////////////

    ////////////
    // MaestroInit.cpp functions

    // fill in multifab and base state data
    void InitData ();

    // During initialization of a simulation, Maestro::InitData() calls 
    // AmrCore::InitFromScratch(), which calls 
    // a MakeNewGrids() function that repeatedly calls this function to create
    // finer levels.  This function creates a new fine
    // level that did not exist before by interpolating from the coarser level
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev,
                                          amrex::Real time,
                                          const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // initial projection
    void InitProj ();

    // divu iter
    void DivuIter (int istep_divu_iter);

    // init iter to initialize gradpi
    void InitIter ();

    // end MaestroInit.cpp functions
    ////////////

    ////////////
    // MaestroMacProj.cpp functions

    void MacProj ();

    // end MaestroMacProj.cpp functions
    ////////////

    ////////////
    // MaestroMakeS.cpp functions

    void Make_S_cc (amrex::Vector<amrex::MultiFab>& S_cc,
                    const amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::MultiFab>& rho_omegadot,
                    const amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                    const amrex::Vector<amrex::MultiFab>& rho_Hext,
                    const amrex::Vector<amrex::MultiFab>& thermal);

    void Make_NodalRHS (const amrex::Vector<amrex::MultiFab>& S_cc,
                        amrex::Vector<amrex::MultiFab>& nodalrhs,
                        const amrex::Vector<amrex::Real>& Sbar,
                        const amrex::Vector<amrex::Real>& beta0);

    // end MaestroMakeS.cpp functions
    ////////////

    ////////////
    // MaestroNodalProj.cpp functions

    // perform a nodal projection
    void NodalProj (int proj_type,
                    amrex::Vector<amrex::MultiFab>& nodalrhs,
                    amrex::Vector<amrex::MultiFab>& rhohalf);

    // modify unew depending on proj_type
    // initial_projection_comp: leave unew alone
    // divu_iters_comp:         leave unew alone
    // pressure_iters_comp:     unew = (unew-uold)/dt
    // regular_timestep_comp:   unew = unew + dt*gpi/rhohalf
    void CreateUvecForProj (int proj_type,
                            amrex::Vector<amrex::MultiFab>& rhohalf);

    // end MaestroNodalProj.cpp functions
    ////////////

    ////////////
    // MaestroPlot.cpp functions

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF (const amrex::Vector<amrex::MultiFab>& p0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& rho0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& u_in,
                                                            amrex::Vector<amrex::MultiFab>& s_in,
                                                      const amrex::Vector<amrex::Real>& p0_in);

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile (const int step,
                        const amrex::Real dt_in,
                        const amrex::Vector<amrex::Real>& rho0_in,
                        const amrex::Vector<amrex::Real>& p0_in,
                        const amrex::Vector<amrex::MultiFab>& u_in,
                              amrex::Vector<amrex::MultiFab>& s_in);

    // end MaestroPlot.cpp functions
    ////////////

    ////////////
    // MaestroReact.cpp functions

    // compute heating term, rho_Hext, then
    // react the state over dt_react and update rho_omegadot, rho_Hnuc
    void React (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const amrex::Vector<amrex::Real>& p0,
                const amrex::Real dt_react);

    void Burner (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                const amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const amrex::Vector<amrex::Real>& p0,
                const amrex::Real dt_react);

    // compute heating term, rho_Hext
    void ComputeHeating (amrex::Vector<amrex::MultiFab>& rho_Hext);

    // end MaestroReact.cpp functions
    ////////////

    ////////////
    // regridding functions in MaestroRegrid.cpp

    // check to see if we need to regrid, then regrid
    void Regrid ();

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ng) override;

    // within a call to AmrCore::regrid, this function fills in data at a level
    // that existed before, using pre-existing fine and interpolated coarse data
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev,
                              amrex::Real time,
                              const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // within a call to AmrCore::regrid, this function fills in data at a level
    // that did NOT exist before, using interpolated coarse data
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev,
                                         amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // within a call to AmrCore::regrid, this function deletes all data
    // at a level of refinement that is no longer needed
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // end regridding functions
    ////////////

    ////////////
    // MaestroRhoHT.cpp functions
    void TfromRhoH (amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::Real>& p0);

    void TfromRhoP (amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::Real>& p0,
                    int updateRhoH = 0);
    ////////////

    ////////////
    // MaestroSetup.cpp functions

    // read in C++ parameters from inputs file
    void ReadParameters ();

    // define variable mappings (Rho, RhoH, ..., NSCAL, etc.)
    void VariableSetup ();

    // set up BCRec definitions for BC types
    void BCSetup ();

    // end MaestroSetup.cpp functions
    ////////////

    ////////////
    // MaestroSponge.cpp functions
    void MakeSponge (amrex::Vector<amrex::MultiFab>& sponge);
    ////////////

    /*
      private data members
    */
    
    ////////////////
    // static data members defined in _cpp_parameters
    // this file is in Source/param_includes/
#include <maestro_params.H>
    ////////////////

    // these are defined in VariableSetup()
    static int Rho, RhoH, FirstSpec, NumSpec, Temp, Pi, Nscal;

    // these are defined in Maestro.cpp
    // these represent the four types of nodal projections
    static int initial_projection_comp;
    static int divu_iters_comp;
    static int pressure_iters_comp;
    static int regular_timestep_comp;

    // these are defined in Maestro.cpp
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;

    // which step?
    int istep;

    // keep track of old time, new time, and time step at each level
    amrex::Real t_new;
    amrex::Real t_old;
    amrex::Real dt;
    amrex::Real dtold;

    // number of ghost cells needed for hyperbolic step
    int ng_adv;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap" to change pointers to old/new
    // to avoid a copy.  These MultiFabs go through the regridding process
    amrex::Vector<amrex::MultiFab> sold;
    amrex::Vector<amrex::MultiFab> snew;
    amrex::Vector<amrex::MultiFab> uold;
    amrex::Vector<amrex::MultiFab> unew;
    amrex::Vector<amrex::MultiFab> S_cc_old;
    amrex::Vector<amrex::MultiFab> S_cc_new;

    // these MultiFabs are persistant and need to be regridded as well
    amrex::Vector<amrex::MultiFab> gpi;
    amrex::Vector<amrex::MultiFab> dSdt;

    // spherical only - 
    // we make this persistent in that we only have to rebuild and
    // fill this after regridding
    amrex::Vector<amrex::MultiFab> normal;

    // stores domain boundary conditions
    // these muse be vectors (rather than arrays) so we can ParmParse them
    amrex::Vector<int> lo_bc;
    amrex::Vector<int> hi_bc;    

    // boundary condition objects needed for FillPatch routines
    // this is essentially an array (over components)
    // of 2*DIM integer arrays storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs_s;  // for sold/snew
    amrex::Vector<amrex::BCRec> bcs_u;  // for uold/unew
    amrex::Vector<amrex::BCRec> bcs_h;  // for rhoh
    amrex::Vector<amrex::BCRec> bcs_t;  // for temperature
    amrex::Vector<amrex::BCRec> bcs_f;  // "first-order extrap" used for S_cc, gpi, dSdt

    // boundary conditions objects needed for linear solves
    int lo_inflow[AMREX_SPACEDIM];
    int hi_inflow[AMREX_SPACEDIM];
    int mg_bcs_p[2*AMREX_SPACEDIM];
    int mg_bcs_h[2*AMREX_SPACEDIM];

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    amrex::Vector<amrex::Real> s0_init;
    amrex::Vector<amrex::Real> p0_init;
    amrex::Vector<amrex::Real> rho0_old;
    amrex::Vector<amrex::Real> rho0_new;
    amrex::Vector<amrex::Real> rhoh0_old;
    amrex::Vector<amrex::Real> rhoh0_new;
    amrex::Vector<amrex::Real> p0_old;
    amrex::Vector<amrex::Real> p0_new;
    amrex::Vector<amrex::Real> tempbar;
    amrex::Vector<amrex::Real> tempbar_init;
    amrex::Vector<amrex::Real> beta0_old;
    amrex::Vector<amrex::Real> beta0_new;
    amrex::Vector<amrex::Real> gamma1bar_old;
    amrex::Vector<amrex::Real> gamma1bar_new;
    amrex::Vector<amrex::Real> etarho_cc;
    amrex::Vector<amrex::Real> psi;
    amrex::Vector<amrex::Real> grav_cell_old;
    amrex::Vector<amrex::Real> grav_cell_new;
    amrex::Vector<amrex::Real> r_cc_loc;

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    amrex::Vector<amrex::Real> w0;
    amrex::Vector<amrex::Real> etarho_ec;
    amrex::Vector<amrex::Real> r_edge_loc;

    // information for the base state
    // we also store these (except for the large arrays) 
    // in base_state.f90 in a fortran module for convenience
    int max_radial_level;
    int finest_radial_level;
    int nr_fine;
    amrex::Real dr_fine;

    // these provide information about the multilevel base state configuration
    amrex::Vector<int> numdisjointchunks;
    amrex::Vector<int> r_start_coord;
    amrex::Vector<int> r_end_coord;

    // tagging criteria
    amrex::Vector<amrex::Real> temperr;

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "max_level+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[max_level] are never actually 
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_u;

};

#endif
