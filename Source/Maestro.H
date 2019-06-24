#ifndef Maestro_H_
#define Maestro_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <PhysBCFunctMaestro.H>
#include <Maestro_F.H>

#ifdef AMREX_USE_CUDA
#include <AMReX_CudaAllocators.H>
#endif

// Define vector types for CUDA
#ifdef AMREX_USE_CUDA
typedef amrex::Vector< amrex::Real, amrex::CudaManagedAllocator<amrex::Real> > RealVector;
#else
typedef amrex::Vector< amrex::Real > RealVector;
#endif

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

    // for keeping track of the amount of CPU time used -- this will persist
    // after restarts
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

    static amrex::Real getCPUTime();

private:

    /*
      private member functions
    */

    ////////////
    // MaestroAdvance.cpp functions

    // advance solution at all levels for a single time step
    void AdvanceTimeStep (bool is_initIter);

    // advance solution for a single time step with irregular base state spacing
    void AdvanceTimeStepIrreg (bool is_initIter);

    // advance solution for a single time step with regular base state spacing
    // and new time-stepping scheme
    void AdvanceTimeStepAverage (bool is_initIter);

    // end MaestroAdvance.cpp functions
    ////////////

    ////////////
    // MaestroAdvection.cpp functions
    void AdvancePremac (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const RealVector& w0_force,
			const amrex::Vector<amrex::MultiFab>& w0_force_cart);

    void MakeUtrans (const amrex::Vector<amrex::MultiFab>& utilde,
                     const amrex::Vector<amrex::MultiFab>& ufull,
                     amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);

    void VelPred (const amrex::Vector<amrex::MultiFab>& utilde,
                  const amrex::Vector<amrex::MultiFab>& ufull,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                  amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const amrex::Vector<amrex::MultiFab>& force);

    void MakeEdgeScal (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                       const amrex::Vector<amrex::MultiFab>& force,
                       int is_vel, const amrex::Vector<amrex::BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative);

    void MakeRhoXFlux (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		       amrex::Vector<amrex::MultiFab>& etarhoflux,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
		       const RealVector& r0_old,
		       const RealVector& r0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_old,
		       const RealVector& r0_new,
		       const RealVector& r0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_new,
		       const RealVector& r0_predicted_edge,
		       int start_comp, int num_comp);

    void MakeRhoHFlux (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
		       const RealVector& r0_old,
		       const RealVector& r0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_old,
		       const RealVector& r0_new,
		       const RealVector& r0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_new,
		       const RealVector& rh0_old,
		       const RealVector& rh0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& rh0mac_old,
		       const RealVector& rh0_new,
		       const RealVector& rh0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& rh0mac_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& h0mac_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& h0mac_new);

    void UpdateScal (const amrex::Vector<amrex::MultiFab>& stateold,
		     amrex::Vector<amrex::MultiFab>& statenew,
		     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		     const amrex::Vector<amrex::MultiFab>& force,
		     int start_scomp, int num_comp,
		     const amrex::Real* p0 = NULL,
		     const amrex::Vector<amrex::MultiFab>& p0_cart = {});

    void UpdateVel (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		    const amrex::Vector<amrex::MultiFab>& force,
		    const amrex::Vector<amrex::MultiFab>& sponge,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);
    ////////////

    ////////////
    // MaestroAverage.cpp functions

    // compute the radial average of a quantitiy
    void Average (const amrex::Vector<amrex::MultiFab>& mf,
                  RealVector& phibar,
                  int comp);

    // end MaestroAverage.cpp functions
    ////////////

    ////////////
    // MaestroCheckpoint.cpp functions
    void WriteCheckPoint (int step);
    int ReadCheckPoint ();
    void GotoNextLine (std::istream& is);
    ////////////

    ////////////
    // MaestroConvert.cpp functions
    void PutInPertForm (amrex::Vector<amrex::MultiFab>& scal,
                        const RealVector& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);

    void PutInPertForm (int level, amrex::Vector<amrex::MultiFab>& scal,
                        const RealVector& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);

    void ConvertRhoXToX (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);

    void ConvertRhoHToH (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);
    ////////////

    ////////////////////////
    // MaestroDebug.cpp functions

    // print out the contents of a Vector of MultiFabs
    void PrintMF   (amrex::Vector<amrex::MultiFab>& CC);


    // print out the contents of a Vector of edge-based MultiFabs
    void PrintEdge (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir);

    // utility to write out a multilevel multifab to a plotfile
    void WriteMF   (amrex::Vector<amrex::MultiFab>& mf,
                    std::string name);
    ////////////////////////

    ////////////////////////
    // MaestroDensityAdvance.cpp functions
    void DensityAdvance (int which_step,
                         amrex::Vector<amrex::MultiFab>& scalold,
                         amrex::Vector<amrex::MultiFab>& scalnew,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
                         amrex::Vector<amrex::MultiFab>& scal_force,
			 amrex::Vector<amrex::MultiFab>& etarhoflux,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			 const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
			 const RealVector& rho0_predicted_edge);
    ////////////////////////

    ////////////
    // MaestroDiag.cpp functions
    // put together an array of multifabs for writing
    void WriteDiagFile (int& index);

    // write plotfile to disk
    void DiagFile (const int step,
		   const amrex::Real dt_in,
		   const RealVector& rho0_in,
		   const RealVector& p0_in,
		   const amrex::Vector<amrex::MultiFab>& u_in,
		   const amrex::Vector<amrex::MultiFab>& s_in,
		   int& index);
    // end MaestroDiag.cpp functions
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
    // MaestroEnthalpyAdvance.cpp functions
    void EnthalpyAdvance (int which_step,
                          amrex::Vector<amrex::MultiFab>& scalold,
                          amrex::Vector<amrex::MultiFab>& scalnew,
                          amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                          amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
                          amrex::Vector<amrex::MultiFab>& scal_force,
                          amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			  const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                          const amrex::Vector<amrex::MultiFab>& thermal);
    ////////////////////////

    ////////////////////////
    // MaestroFillData.cpp functions

    // call FillPatch for all levels
    void FillPatch (amrex::Real time,
                    amrex::Vector<amrex::MultiFab>& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in, int variable_type=0);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases
    // (fill fine grid ghost by interpolating from coarse)
    // scomp of the source component
    // dcomp is the destination component AND the bc component
    void FillPatch (int lev, amrex::Real time,
                    amrex::MultiFab& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in, int variable_type=0);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    // scomp of the source component
    // dcomp is the destination component AND the bc component
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<amrex::MultiFab>& mf_old,
                          amrex::Vector<amrex::MultiFab>& mf_new,
                          int srccomp, int destcomp, int ncomp,
                          const amrex::Vector<amrex::BCRec>& bcs, int variable_type=0);

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

    // set covered faces to be the average of overlying fine faces
    void AverageDownFaces (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge);

    // fill in ONE ghost cell for all components of a face-centered (MAC) velocity
    // field behind physical boundaries.  Does not modify the velocities on the boundary
    void FillUmacGhost (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int level=-1);

    // fill in all ghost cells for an edge-based MAC velocity field
    void FillPatchUedge(amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge);

    // end MaestroFillData.cpp functions
    ////////////

    ////////////////////////
    // MaestroFill3dData.cpp functions

    void Put1dArrayOnCart (const RealVector& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0, int variable_type=0);

    void Put1dArrayOnCart (int level, const RealVector& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0);

    void Addw0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                const amrex::Real& mult);

    void MakeW0mac (amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac);

    void MakeS0mac (const RealVector& s0,
		    amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& s0mac);

    void MakeNormal ();

    void PutDataOnFaces(const amrex::Vector<amrex::MultiFab>& s_cc,
			amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM >>& face,
			int harmonic_avg);

    void MakeCCtoRadii ();
    // end MaestroFill3dData.cpp functions
    ////////////

    ////////////////////////
    // MaestroForce.cpp functions
    void MakeVelForce (amrex::Vector<amrex::MultiFab>& vel_force,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
                       const amrex::Vector<amrex::MultiFab>& rho,
                       const RealVector& rho0,
                       const RealVector& grav,
                       const RealVector& w0_force,
		       const amrex::Vector<amrex::MultiFab>& w0_force_cart,
                       int do_add_utilde_force);

    void ModifyScalForce(amrex::Vector<amrex::MultiFab>& scal_force,
			 const amrex::Vector<amrex::MultiFab>& state,
                         const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                         const RealVector& s0,
                         const RealVector& s0_edge,
			 const amrex::Vector<amrex::MultiFab>& s0_cart,
                         int comp,
                         const amrex::Vector<amrex::BCRec>& bcs,
                         int fullform);

    void MakeRhoHForce (amrex::Vector<amrex::MultiFab>& scal_force,
                        int is_prediction,
                        const amrex::Vector<amrex::MultiFab>& thermal,
                        const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int add_thermal,
			const int &which_step);
    ////////////////////////

    ////////////
    // MaestroGamma.cpp functions

    void MakeGamma1bar (const amrex::Vector<amrex::MultiFab>& scal,
                        RealVector& gamma1bar,
                        const RealVector& p0);

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

    void MacProj (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  amrex::Vector<amrex::MultiFab>& macphi,
                  const amrex::Vector<amrex::MultiFab>& macrhs,
                  const RealVector& beta0,
                  const int& is_predictor);

    // multiply (or divide) face-data by beta0
    void MultFacesByBeta0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge,
                           const RealVector& beta0,
                           const RealVector& beta0_edge,
                           const int& mult_or_div);

    // compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    void ComputeMACSolverRHS (amrex::Vector<amrex::MultiFab>& solverrhs,
                              const amrex::Vector<amrex::MultiFab>& macrhs,
                              const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);

    // Average bcoefs at faces using inverse of rho
    void AvgFaceBcoeffsInv(amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM > >& facebcoef,
			   const amrex::Vector<amrex::MultiFab>& rhocc);

    // set boundaries for LABecLaplacian to solve -div(B grad) phi = RHS
    void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

    // end MaestroMacProj.cpp functions
    ////////////

    ////////////
    // MaestroMakeEta.cpp functions
    void MakeEtarho (RealVector& etarho_edge,
		     RealVector& etarho_cell,
		     const amrex::Vector<amrex::MultiFab>& etarho_flux);

    void MakeEtarhoSphr (const amrex::Vector<amrex::MultiFab>& scal_old,
			 const amrex::Vector<amrex::MultiFab>& scal_new,
			 const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			 const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
			 RealVector& etarho_edge,
			 RealVector& etarho_cell);
    ////////////

    ////////////
    // MaestroMakeS.cpp functions

    // compute S at cell-centers
    void Make_S_cc (amrex::Vector<amrex::MultiFab>& S_cc,
                    amrex::Vector<amrex::MultiFab>& delta_gamma1_term,
                    amrex::Vector<amrex::MultiFab>& delta_gamma1,
                    const amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::MultiFab>& u,
                    const amrex::Vector<amrex::MultiFab>& rho_omegadot,
                    const amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                    const amrex::Vector<amrex::MultiFab>& rho_Hext,
                    const amrex::Vector<amrex::MultiFab>& thermal,
                    const RealVector& p0,
                    const RealVector& gamma1bar,
                    RealVector& delta_gamma1_termbar,
                    const RealVector& psi_in);

    // compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
    void MakeRHCCforNodalProj (amrex::Vector<amrex::MultiFab>& rhcc,
                               const amrex::Vector<amrex::MultiFab>& S_cc,
                               const RealVector& Sbar,
                               const RealVector& beta0,
                               const amrex::Vector<amrex::MultiFab>& delta_gamma1_term);

    // compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
    void MakeRHCCforMacProj (amrex::Vector<amrex::MultiFab>& rhcc,
			                 const RealVector& rho0,
                             const amrex::Vector<amrex::MultiFab>& S_cc,
                             const RealVector& Sbar,
                             const RealVector& beta0,
                             const amrex::Vector<amrex::MultiFab>& delta_gamma1_term,
                             const RealVector& gamma1bar,
			                 const RealVector& p0,
                             const amrex::Vector<amrex::MultiFab>& delta_p_term,
			                 amrex::Vector<amrex::MultiFab>& delta_chi,
                             int is_predictor);

    void CorrectRHCCforNodalProj(amrex::Vector<amrex::MultiFab>& rhcc,
				 const RealVector& rho0,
				 const RealVector& beta0,
				 const RealVector& gamma1bar,
				 const RealVector& p0,
				 const amrex::Vector<amrex::MultiFab>& delta_p_term);
    // end MaestroMakeS.cpp functions
    ////////////

    ////////////
    // MaestroNodalProj.cpp functions

    // perform a nodal projection
    void NodalProj (int proj_type,
                    amrex::Vector<amrex::MultiFab>& rhcc,
                    int istep_divu_iter=0);

    // fill in Vproj
    // initial_projection_comp: Utilde^0                        -- uold
    // divu_iters_comp:         Utilde^0                        -- uold
    // pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
    // regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
    // sig contains rhohalf if proj_type == regular_timestep_comp
    void CreateUvecForProj (int proj_type,
                            amrex::Vector<amrex::MultiFab>& Vproj,
                            const amrex::Vector<amrex::MultiFab>& sig);

    void SetBoundaryVelocity(amrex::Vector<amrex::MultiFab>& vel);

    // given a nodal phi, compute grad(phi) at cell centers
    void ComputeGradPhi(amrex::Vector<amrex::MultiFab>& phi,
                        amrex::Vector<amrex::MultiFab>& gphi);

    // average nodal pi to cell-centers and put in the Pi component of snew
    void MakePiCC(const amrex::Vector<amrex::MultiFab>& beta0_cart);

    // end MaestroNodalProj.cpp functions
    ////////////

    ////////////
    // MaestroPlot.cpp functions

    // get plotfile name
    void PlotFileName (const int lev, std::string* plotfilename);

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF (const int nPlot,
                                                      const amrex::Real t_in,
                                                      const amrex::Real dt_in,
                                                      const amrex::Vector<amrex::MultiFab>& rho0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& rhoh0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& p0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& gamma1bar_cart,
                                                      const amrex::Vector<amrex::MultiFab>& u_in,
                                                      amrex::Vector<amrex::MultiFab>& s_in,
                                                      const RealVector& p0_in,
                                                      const RealVector& gamma1bar_in,
                                                      const amrex::Vector<amrex::MultiFab>& S_cc_in);

    amrex::Vector<const amrex::MultiFab*> SmallPlotFileMF(const int nPlot,
                           const int nSmallPlot,
                           amrex::Vector<const amrex::MultiFab*> mf,
                           const amrex::Vector<std::string> varnames,
                           const amrex::Vector<std::string> small_plot_varnames);

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames (int * nPlot) const;

    // set small plotfile variables names
    amrex::Vector<std::string> SmallPlotFileVarNames (int * nPlot, amrex::Vector<std::string> varnames) const;

    // write a small plotfile to disk
    void WriteSmallPlotFile (const int step,
                        const amrex::Real t_in,
                        const amrex::Real dt_in,
                        const RealVector& rho0_in,
                        const RealVector& rhoh0_in,
                        const RealVector& p0_in,
                        const RealVector& gamma1bar_in,
                        const amrex::Vector<amrex::MultiFab>& u_in,
                        amrex::Vector<amrex::MultiFab>& s_in,
                        const amrex::Vector<amrex::MultiFab>& S_cc_in);

    // write plotfile to disk
    void WritePlotFile (const int step,
                        const amrex::Real t_in,
                        const amrex::Real dt_in,
                        const RealVector& rho0_in,
                        const RealVector& rhoh0_in,
                        const RealVector& p0_in,
                        const RealVector& gamma1bar_in,
                        const amrex::Vector<amrex::MultiFab>& u_in,
                        amrex::Vector<amrex::MultiFab>& s_in,
                        const amrex::Vector<amrex::MultiFab>& S_cc_in,
                        const bool is_small = false);

    void WriteJobInfo (const std::string& dir) const;

    void MakeMagvel (const amrex::Vector<amrex::MultiFab>& vel,
                     amrex::Vector<amrex::MultiFab>& magvel);

    void MakeVelrc (const amrex::Vector<amrex::MultiFab>& vel,
                    const amrex::Vector<amrex::MultiFab>& w0rcart,
                    amrex::Vector<amrex::MultiFab>& rad_vel,
                    amrex::Vector<amrex::MultiFab>& circ_vel);

    void MakeAdExcess (const amrex::Vector<amrex::MultiFab>& state,
                       amrex::Vector<amrex::MultiFab>& ad_excess);

    void MakeVorticity (const amrex::Vector<amrex::MultiFab>& vel,
                        amrex::Vector<amrex::MultiFab>& vorticity);

    void MakeDeltaGamma (const amrex::Vector<amrex::MultiFab>& state,
                         const RealVector& p0,
                         const amrex::Vector<amrex::MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const amrex::Vector<amrex::MultiFab>& gamma1bar_cart,
                         amrex::Vector<amrex::MultiFab>& deltagamma);

    void MakeEntropy (const amrex::Vector<amrex::MultiFab>& state,
                      amrex::Vector<amrex::MultiFab>& entropy);

    void MakeDivw0 (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                    amrex::Vector<amrex::MultiFab>& divw0);

    void MakePiDivu (const amrex::Vector<amrex::MultiFab>& vel,
                   const amrex::Vector<amrex::MultiFab>& pi_cc,
                   amrex::Vector<amrex::MultiFab>& pidivu);

    // mass fractions of the species
    void MakeAbar(const amrex::Vector<amrex::MultiFab>& state,
                amrex::Vector<amrex::MultiFab>& abar);

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
                const RealVector& p0,
                const amrex::Real dt_react);

    void Burner (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                const amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const RealVector& p0,
                const amrex::Real dt_react);

    // compute heating term, rho_Hext
    void MakeHeating (amrex::Vector<amrex::MultiFab>& rho_Hext,
                      const amrex::Vector<amrex::MultiFab>& scal);

    // end MaestroReact.cpp functions
    ////////////

    ////////////
    // regridding functions in MaestroRegrid.cpp

    // check to see if we need to regrid, then regrid
    void Regrid ();

    // set tagging array to include buffer zones for multilevel
    void TagArray ();

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
                    const RealVector& p0);

    void TfromRhoP (amrex::Vector<amrex::MultiFab>& scal,
                    const RealVector& p0,
                    int updateRhoH = 0);

    void PfromRhoH (const amrex::Vector<amrex::MultiFab>& state,
                    const amrex::Vector<amrex::MultiFab>& s_old,
		    amrex::Vector<amrex::MultiFab>& peos);

    void MachfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                       const amrex::Vector<amrex::MultiFab>& vel,
                       const RealVector& p0,
                       amrex::Vector<amrex::MultiFab>& mach);

    void MachfromRhoHSphr (const amrex::Vector<amrex::MultiFab>& scal,
                           const amrex::Vector<amrex::MultiFab>& vel,
                           const RealVector& p0,
                           const amrex::Vector<amrex::MultiFab>& w0cart,
                           amrex::Vector<amrex::MultiFab>& mach);

    void CsfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                     const RealVector& p0,
                     const amrex::Vector<amrex::MultiFab>& p0cart,
                     amrex::Vector<amrex::MultiFab>& cs);
    ////////////

    ////////////
    // MaestroSetup.cpp functions

    void ExternInit ();

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

    ////////////////////////
    // MaestroThermal.cpp functions
    void MakeExplicitThermal (amrex::Vector<amrex::MultiFab>& thermal,
			      const amrex::Vector<amrex::MultiFab>& scal,
			      const amrex::Vector<amrex::MultiFab>& Tcoeff,
			      const amrex::Vector<amrex::MultiFab>& hcoeff,
			      const amrex::Vector<amrex::MultiFab>& Xkcoeff,
			      const amrex::Vector<amrex::MultiFab>& pcoeff,
			      const RealVector& p0,
			      int temp_formulation);

    void ApplyThermal (amrex::MLABecLaplacian& mlabec,
		       amrex::Vector<amrex::MultiFab>& thermalout,
		 const amrex::Vector<amrex::MultiFab>& coeff,
		       amrex::Vector<amrex::MultiFab>& phi,
		 const amrex::Vector<amrex::BCRec>& bcs,
		 int bccomp);

    void MakeThermalCoeffs(const amrex::Vector<amrex::MultiFab>& scal,
			   amrex::Vector<amrex::MultiFab>& Tcoeff,
			   amrex::Vector<amrex::MultiFab>& hcoeff,
			   amrex::Vector<amrex::MultiFab>& Xkcoeff,
			   amrex::Vector<amrex::MultiFab>& pcoeff);

    void ThermalConduct (const amrex::Vector<amrex::MultiFab>& s1,
			       amrex::Vector<amrex::MultiFab>& s2,
			 const amrex::Vector<amrex::MultiFab>& hcoeff1,
			 const amrex::Vector<amrex::MultiFab>& Xkcoeff1,
			 const amrex::Vector<amrex::MultiFab>& pcoeff1,
			 const amrex::Vector<amrex::MultiFab>& hcoeff2,
			 const amrex::Vector<amrex::MultiFab>& Xkcoeff2,
			 const amrex::Vector<amrex::MultiFab>& pcoeff2);
    ////////////////////////

    ////////////////////////
    // MaestroVelocityAdvance.cpp functions
    void VelocityAdvance (const amrex::Vector<amrex::MultiFab>& rhohalf,
                          amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			  const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                          const RealVector& w0_force,
			  const amrex::Vector<amrex::MultiFab>& w0_force_cart,
                          const RealVector& rho0_nph,
                          const RealVector& grav_cell_nph,
			  const amrex::Vector<amrex::MultiFab>& sponge);
    ////////////////////////


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
    static int ng_s; // number of ghost cells for sold/new and uold/new

    // these are defined in Maestro.cpp
    // these represent the four types of nodal projections
    static int initial_projection_comp;
    static int divu_iters_comp;
    static int pressure_iters_comp;
    static int regular_timestep_comp;

    // these are defined in Maestro.cpp
    // density and enthalpy edge state prediction
    static int predict_rhoprime_and_X, predict_rhoX, predict_rho_and_X;
    static int predict_rhoh, predict_rhohprime, predict_h;
    static int predict_T_then_rhohprime, predict_T_then_h;
    static int predict_hprime, predict_Tprime_then_h;

    // these are defined in Maestro.cpp
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;

    // which step?
    int istep;
    int start_step;

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
    amrex::Vector<amrex::MultiFab> pi;   // nodal

    // this only needs to persist leading into the initial pressure iters
    // since we project (beta0^nph S^1 - beta0 S^0) / dt
    // during a regular time step we overwrite this
    amrex::Vector<amrex::MultiFab> rhcc_for_nodalproj;

    // spherical only -
    // we make this persistent in that we only have to rebuild and
    // fill this after regridding
    amrex::Vector<amrex::MultiFab> normal;
    amrex::Vector<amrex::MultiFab> cell_cc_to_r;

    // stores domain boundary conditions
    // these muse be vectors (rather than arrays) so we can ParmParse them
    amrex::Vector<int> phys_bc;

    // boundary condition objects needed for FillPatch routines
    // this is essentially an array (over components)
    // of 2*DIM integer arrays storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs_s;  // for sold/snew
    amrex::Vector<amrex::BCRec> bcs_u;  // for uold/unew
    amrex::Vector<amrex::BCRec> bcs_f;  // "first-order extrap" used for S_cc, gpi, dSdt

    // vectors store the multilevel 1D states as one very long array
    // these are cell-centered
    RealVector s0_init;
    RealVector p0_init;
    RealVector rho0_old;
    RealVector rho0_new;
    RealVector rhoh0_old;
    RealVector rhoh0_new;
    RealVector rhoX0_old;
    RealVector rhoX0_new;
    RealVector p0_old;
    RealVector p0_new;
    RealVector tempbar;
    RealVector tempbar_init;
    RealVector beta0_old;
    RealVector beta0_new;
    RealVector gamma1bar_old;
    RealVector gamma1bar_new;
    RealVector etarho_cc;
    RealVector psi;
    RealVector grav_cell_old;
    RealVector grav_cell_new;
    RealVector r_cc_loc;
    // p0^{n-1} needed to compute d(p0)/dt for nonuniform grid spacing
    RealVector p0_nm1;
    RealVector beta0_nm1;

    // vectors store the multilevel 1D states as one very long array
    // these are edge-centered
    RealVector w0;
    RealVector etarho_ec;
    RealVector r_edge_loc;

    // information for the base state
    // we also store these (except for the large arrays)
    // in base_state.f90 in a fortran module for convenience
    int max_radial_level;
    int finest_radial_level;
    int nr_fine;
    amrex::Real dr_fine;
    int nr_irreg;

    // these provide information about the multilevel base state configuration
    amrex::Vector<int> numdisjointchunks;
    amrex::Vector<int> r_start_coord;
    amrex::Vector<int> r_end_coord;

    // array of tagged boxes (planar)
    amrex::Vector<int> tag_array;

    // diag file array buffers
    amrex::Vector<amrex::Real> diagfile1_data;
    amrex::Vector<amrex::Real> diagfile2_data;
    amrex::Vector<amrex::Real> diagfile3_data;

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "max_level+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[max_level] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;

    // flag for writing plotfiles
    enum plotfile_flag {plotInitData = -9999999, plotInitProj = -9999998, plotDivuIter = -9999997};

};

#endif
