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

    // advance solution to final time
    void Evolve ();

    // initializes multilevel data
    void InitData ();

    // define variable mappings (Rho, RhoH, ..., NUM_STATE, etc.)
    void VariableSetup();

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ngrow) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev,
                                          amrex::Real time,
                                          const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // Make a new level using provided BoxArray and DistributionMapping and 
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev,
                                         amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and 
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev,
                              amrex::Real time,
                              const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters ();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute the number of cells at a level
    long CountCells (int lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev,
                    amrex::Real time,
                    amrex::MultiFab& mf,
                    int icomp,
                    int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev,
                          amrex::Real time,
                          amrex::MultiFab& mf,
                          int icomp,
                          int ncomp);

    // utility to copy in data from sold and/or snew into another multifab
    void GetData (int lev,
                  amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // advance solution at all levels for a single time step
    void AdvanceTimeStep (amrex::Real time);

    // a wrapper for ComputeDtLevel
    void ComputeDt ();

    // compute dt at a level from CFL considerations
    amrex::Real ComputeDtLevel (int lev) const;

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile (int step) const;

    ////////////////
    // private data members

    // which step?
    int istep;

    // keep track of old time, new time, and time step at each level
    amrex::Real t_new;
    amrex::Real t_old;
    amrex::Real dt;

    static int NUM_STATE;
    static int Rho, RhoH;
    static int FirstSpec, NumSpec;
    static int Temp, Pi;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap" to change pointers to old/new
    // to avoid a copy
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > sold;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > snew;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > uold;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > unew;

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > hgrhs;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > gamma1;

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > pi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > gpi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > dSdt;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > Source_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > Source_new;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_omegadot2;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_Hnuc2;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rho_Hext;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > thermal2;

    // this is essentially an array (over components)
    // of 2*DIM integer arrays storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs;

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually 
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_u;
    
    ////////////////
    // runtime parameters

#include <maestro_params.H>

};

#endif
