#ifndef Maestro_H_
#define Maestro_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Tuple.H>

#include <MaestroPhysBCFunct.H>
#include <Maestro_F.H>

class Maestro
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel vectors and data structures
    Maestro ();
    virtual ~Maestro();

    // read in some parameters from inputs file
    void ReadParameters ();

    // define variable mappings (Rho, RhoH, ..., NSCAL, etc.)
    void VariableSetup();

    // set up BCRec definitions for BC types
    void BCSetup();

    // initializes multilevel data
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


    // advance solution to final time
    void Evolve ();

    ////////////
    // Regridding functions

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ngrow) override;

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

    //
    ////////////

private:

    /*
      private member functions
    */

    // initial projection
    void InitProj();

    // divu iter
    void DivuIter();

    // init iter to initialize gradpi
    void InitIter();

    // advance solution at all levels for a single time step
    void AdvanceTimeStep (bool is_initIter);

    ////////////
    // Time step computation

    // a wrapper for ComputeDtLevel
    void ComputeDt ();

    // compute dt at a level from CFL considerations
    amrex::Real ComputeDtLevel (int lev) const;

    //
    ////////////

    ////////////
    // Functions that fill data in MultiFabs

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases 
    // (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_old,
                    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_new,
                    int icomp, int ncomp, amrex::Vector<amrex::BCRec> bcs);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_old,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_new,
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
                  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_old,
                  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf_new);

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf);

    // more flexible version of AverageDown() that lets you average 
    // down across multiple levels
    void AverageDownTo (int crse_lev,
                        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf);

    //
    ////////////

    void Make_S_cc (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& S_cc);

    ////////////
    // PlotFile routines

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile (int step) const;

    //
    ////////////

    // compute the number of cells at a level
    long CountCells (int lev);

    /*
      private data members
    */

    // which step?
    int istep;

    // keep track of old time, new time, and time step at each level
    amrex::Real t_new;
    amrex::Real t_old;
    amrex::Real dt;

    static int NSCAL;
    static int Rho, RhoH, FirstSpec, NumSpec, Temp, Pi;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap" to change pointers to old/new
    // to avoid a copy
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sold;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > snew;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > uold;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > unew;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_cc_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_cc_new;

    // these MultiFabs need to be regridded as well
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > gpi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > dSdt;

    // these are not persistent - after regridding these can be rebuilt
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > S_nodal;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > normal;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sponge;

    // these are not persistent, in the original MAESTRO they are 
    // kept around for restarting into matching plotfiles only
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > pi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_omegadot;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_Hnuc;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_Hext;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > thermal;

    // this is essentially an array (over components)
    // of 2*DIM integer arrays storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs_s;
    amrex::Vector<amrex::BCRec> bcs_u;

    // this contains the base state
    // stores as one long vector
    amrex::Vector<amrex::Real> rho0_old;
    amrex::Vector<amrex::Real> rho0_new;
    amrex::Vector<amrex::Real> rhoh0_old;
    amrex::Vector<amrex::Real> rhoh0_new;
    amrex::Vector<amrex::Real> p0_old;
    amrex::Vector<amrex::Real> p0_new;

    // this contains the base state for velocity
    amrex::Vector<amrex::Real> w0;

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually 
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_u;
    
    ////////////////
    // runtime parameters - see maestro_params.H

#include <maestro_params.H>

    //
    ////////////////

};

#endif
