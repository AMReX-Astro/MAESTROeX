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

#include <MaestroPhysBCFunct.H>
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
    virtual ~Maestro();

    // read in C++/F90 parameters
    // define global C++/F90 variables and initialize network
    // set up boundary conditions
    // initialize base state geometry parameters
    // set istep, t_new, t_old
    // allocate MultiFabs and base state arrays
    void Setup ();

    // initialize multifab and base state data
    // perform initial projection
    // perform divu iters
    // perform initial (pressure) iterations
    void Init ();


    // advance solution to final time
    void Evolve ();

private:

    /*
      private member functions
    */

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
    void InitProj();

    // divu iter
    void DivuIter();

    // init iter to initialize gradpi
    void InitIter();

    // end MaestroInit.cpp functions
    ////////////

    ////////////
    // MaestroAdvance.cpp functions

    // advance solution at all levels for a single time step
    void AdvanceTimeStep (bool is_initIter);

    // end MaestroAdvance.cpp functions
    ////////////

    ////////////
    // Time step computation in MaestroDt.cpp

    // a wrapper for EstDtLevel
    void EstDt ();

    // compute dt at a level from CFL considerations
    amrex::Real EstDtLevel (int lev) const;

    // compute initial time step
    void FirstDt ();

    // end MaestroDt.cpp functions
    ////////////

    ////////////
    // MaestroFillData.cpp functions

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases 
    // (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int icomp, int ncomp, amrex::Vector<amrex::BCRec> bcs);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<amrex::MultiFab>& mf_old,
                          amrex::Vector<amrex::MultiFab>& mf_new,
                          int icomp, int ncomp,
                          amrex::Vector<amrex::BCRec> bcs);

    // utility to copy in data from mf_old and/or mf_new into mf
    // if time=t_old we copy mf_old into mf
    // if time=t_new we copy mf_new into mf
    // otherwise copy copy in both mf_old and mf_new into mf and the fillpatch
    // routines know to interpolate in time.  However in MAESTRO since we don't
    // subcycle I don't think we need this capability
    void GetData (int lev, amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& mf,
                  amrex::Vector<amrex::Real>& mftime,
                  amrex::Vector<amrex::MultiFab>& mf_old,
                  amrex::Vector<amrex::MultiFab>& mf_new);

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown (amrex::Vector<amrex::MultiFab>& mf);

    // more flexible version of AverageDown() that lets you average 
    // down across multiple levels
    void AverageDownTo (int crse_lev,
                        amrex::Vector<amrex::MultiFab>& mf);

    // end MaestroFillData.cpp functions
    ////////////

    ////////////
    // MaestroMakeS.cpp functions

    void Make_S_cc (amrex::Vector<amrex::MultiFab>& S_cc,
                    const amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::MultiFab>& rho_omegadot,
                    const amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                    const amrex::Vector<amrex::MultiFab>& rho_Hext,
                    const amrex::Vector<amrex::MultiFab>& thermal);

    void Make_S_nodal (const amrex::Vector<amrex::MultiFab>& S_cc,
                       amrex::Vector<amrex::MultiFab>& S_nodal,
                       amrex::Real Sbar,
                       const amrex::Vector<amrex::Real>& div_coeff);

    // end MaestroMakeS.cpp functions
    ////////////

    ////////////
    // MaestroNodalProj.cpp functions

    void NodalProj (amrex::Vector<amrex::MultiFab>& phi,
                    amrex::Vector<amrex::MultiFab>& vel,
                    amrex::Vector<amrex::MultiFab>& rhcc,
                    amrex::Vector<amrex::MultiFab>& rhnd,
                    amrex::Vector<amrex::MultiFab>& beta0,
                    int* bc,
                    amrex::Real rel_tol,
                    amrex::Real abs_tol);


    // end MaestroNodalProj.cpp functions
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

    // end MaestroAverage.cpp functions
    ////////////

    ////////////
    // MaestroPlot.cpp functions

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile (int step) const;

    // end MaestroPlot.cpp functions
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
    static int Nscal;
    static int Rho, RhoH, FirstSpec, NumSpec, Temp, Pi;

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
    amrex::Vector<amrex::Real> div_coeff_old;
    amrex::Vector<amrex::Real> div_coeff_new;
    amrex::Vector<amrex::Real> gamma1bar;
    amrex::Vector<amrex::Real> gamma1bar_init;
    amrex::Vector<amrex::Real> etarho_cc;
    amrex::Vector<amrex::Real> psi;
    amrex::Vector<amrex::Real> grav_cell;
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
