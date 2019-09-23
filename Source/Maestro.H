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
// #include <AMReX_CudaAllocators.H>
#include <AMReX_CudaContainers.H>
#endif

/// Define Real vector types for CUDA-compatability. If `AMREX_USE_CUDA`, then
/// this will be stored in CUDA managed memory.
#ifdef AMREX_USE_CUDA
typedef amrex::Gpu::ManagedVector<amrex::Real> RealVector;
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

    /// constructor
    Maestro ();
    /// destructor
    virtual ~Maestro ();

    // in `MaestroSetup.cpp`
    /// Setup the simulation.
    /// - read in C++/F90 parameters
    /// - define global C++/F90 variables and initialize network
    /// - set up boundary conditions
    /// - initialize base state geometry parameters
    /// - set `istep`, `t_new`, `t_old`
    /// - allocate MultiFabs and base state arrays
    void Setup ();

    // in `MaestroInit.cpp`
    /// Initialize the simulation.
    /// - initialize multifab and base state data
    /// - perform initial projection
    /// - perform divu iters
    /// - perform initial (pressure) iterations
    void Init ();

    // in `MaestroEvolve.cpp`
    /// advance solution to final time
    void Evolve ();

    // for keeping track of the amount of CPU time used. This will persist
    // after restarts
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

    /// Get the amount of CPU time used. This will persist
    /// after restarts.
    static amrex::Real getCPUTime();

    /// Dump build info
    static void WriteBuildInfo ();

private:

    /*
      private member functions
    */

    ////////////
    // MaestroAdvance.cpp functions

    /// Advance solution at all levels for a single time step. This uses the old
    /// temporal integration scheme.
    ///
    /// @param is_initIter is it the initial iteration?
    void AdvanceTimeStep (bool is_initIter);

    /// Advance solution for a single time step with irregular base state spacing
    ///
    /// @param is_initIter is it the initial iteration?
    void AdvanceTimeStepIrreg (bool is_initIter);

    /// Advance solution for a single time step with regular base state spacing
    /// and new time-stepping scheme
    ///
    /// @param is_initIter is it the initial iteration?
    void AdvanceTimeStepAverage (bool is_initIter);

    // end MaestroAdvance.cpp functions
    ////////////

    ////////////
    // MaestroAdvection.cpp functions

    /// Compute unprojected mac velocities
    ///
    /// @param umac             MAC velocity
    /// @param w0mac            MAC base-state velocity
    /// @param w0_force         base-state-velocity force
    /// @param w0_force_cart    base-state-velocity force on cartesian grid
    void AdvancePremac (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const RealVector& w0_force,
			const amrex::Vector<amrex::MultiFab>& w0_force_cart);

    /// Create `utrans`, the transverse velocity
    ///
    /// @param utilde           perturbed velocity
    /// @param ufull            full velocity
    /// @param utrans           transverse velocity
    /// @param w0mac            MAC base-state velocity
    void MakeUtrans (const amrex::Vector<amrex::MultiFab>& utilde,
                     const amrex::Vector<amrex::MultiFab>& ufull,
                     amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);

    /// Predict the normal velocities to the interfaces.  We don't care about the
    /// transverse velocities here.  The prediction is done piecewise linear (for now)
    ///
    /// @param utilde           perturbed velocity
    /// @param ufull            full velocity
    /// @param utrans           transverse velocity
    /// @param umac             MAC velocity
    /// @param w0mac            MAC base-state velocity
    /// @param force            velocity force
    void VelPred (const amrex::Vector<amrex::MultiFab>& utilde,
                  const amrex::Vector<amrex::MultiFab>& ufull,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                  amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const amrex::Vector<amrex::MultiFab>& force);

    /// Calculate scalars on cell edges given cell-centered state
    /// constructs the edge state of a scalar, using a
    /// second-order Taylor expansion in space (through `dx/2`) and time
    /// (though `dt/2`) (if `ppm_type = 0`) or using PPM (for `ppm_type = 1,2`).
    ///
    /// We use only MAC-projected edge velocities in this prediction.
    ///
    /// We are computing all edge states for each variable.  This is what is
    /// done for the final updates of the state variables and velocity.  For
    /// velocity, we should set `is_vel = true`
    ///
    /// @param state            cell-centered scalars
    /// @param sedge            edge state of scalars
    /// @param umac             MAC velocity
    /// @param force            velocity force
    /// @param is_vel           set to true if state is a velocity
    /// @param bcs              boundary conditions
    /// @param nbccomp          number of components of `bcs`
    /// @param start_scomp      index of component of `state` to begin with
    /// @param start_bccomp     index of component of `bcs` to begin with
    /// @param num_comp         number of components to perform calculation for
    /// @param is_conservative  are these conserved quantities?
    void MakeEdgeScal (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                       const amrex::Vector<amrex::MultiFab>& force,
                       int is_vel, const amrex::Vector<amrex::BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative);

    /// Calculate `rhoX` flux
    ///
    /// Takes the predicted edges states of the scalars
    /// and the MAC velocities and computes the flux through the
    /// interfaces.
    ///
    /// The construction of the fluxes depends on
    /// what form the incoming edge states take.  This depends on
    /// `species_pred_type`:
    ///
    /// - `predict_rhoprime_and_X`:
    ///   We have rho' and X, and need a edge-centered base state to
    ///   make the final fluxes
    ///
    /// - `predict_rhoX`:
    ///   We use the (rho X) edge state directly to compute the fluxes.
    ///   No base state input needed.
    ///
    /// - `predict_rho_and_X`:
    ///   The fluxes are computed from the product of the rho and X
    ///   edge states, again, no base state input needed.
    ///
    /// @param state            cell-centered scalars
    /// @param sflux            scalar flux
    /// @param etarhoflux       `eta_rho` flux
    /// @param sedge            edge state of scalars
    /// @param umac, w0mac      MAC velocity of full and base-state velocity
    /// @param r0_old           old base-state density
    /// @param r0_edge_old      old base-state density on cell-edges
    /// @param r0mac_old        old MAC-projected base-state density
    /// @param r0_new           new base-state density
    /// @param r0_edge_new      new base-state density on cell-edges
    /// @param r0mac_new        new MAC-projected base-state density
    /// @param r0_predicted_edge  new base-state density on cell edges
    /// @param start_comp       index of component of `state` to begin with
    /// @param num_comp         number of components to perform calculation for
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

    /// Calculate `rhoh` flux
    ///
    /// Takes the predicted edges states of the scalars
    /// and the MAC velocities and computes the flux through the
    /// interfaces.
    ///
    /// @param state            cell-centered scalars
    /// @param sflux            scalar flux
    /// @param sedge            edge state of scalars
    /// @param umac, w0mac      MAC velocity of full and base-state velocity
    /// @param r0_old           old base-state density
    /// @param r0_edge_old      old base-state density on cell-edges
    /// @param r0mac_old        old MAC-projected base-state density
    /// @param r0_new           new base-state density
    /// @param r0_edge_new      new base-state density on cell-edges
    /// @param r0mac_new        new MAC-projected base-state density
    /// @param rh0_old          old base-state conserved enthalpy
    /// @param rh0_edge_old     old base-state conserved enthalpy on cell-edges
    /// @param rh0mac_old       old MAC-projected base-state conserved enthalpy
    /// @param rh0_new          new base-state conserved enthalpy
    /// @param rh0_edge_new     new base-state conserved enthalpy on cell-edges
    /// @param rh0mac_new       new MAC-projected base-state conserved enthalpy
    /// @param h0mac_old, h0mac_new base-state primitive enthalpy
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

    /// Given scalar fluxes, update scalars
    ///
    /// @param stateold         cell-centered scalars
    /// @param statenew         new scalar flux
    /// @param sflux            scalar fluxes
    /// @param force            velocity force
    /// @param start_scomp      index of component of `state` to begin with
    /// @param num_comp         number of components to perform calculation for
    /// @param p0_cart          base state pressure on cartesian grid
    void UpdateScal (const amrex::Vector<amrex::MultiFab>& stateold,
		     amrex::Vector<amrex::MultiFab>& statenew,
		     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		     const amrex::Vector<amrex::MultiFab>& force,
		     int start_scomp, int num_comp,
		     const amrex::Vector<amrex::MultiFab>& p0_cart);

    /// Update velocity
    ///
    /// @param umac             MAC velocity
    /// @param uedge            edge-based velocity
    /// @param force            velocity force
    /// @param sponge
    /// @param w0mac            base state MAC velocity
    void UpdateVel (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		    const amrex::Vector<amrex::MultiFab>& force,
		    const amrex::Vector<amrex::MultiFab>& sponge,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);
    ////////////

    ////////////
    // MaestroAverage.cpp functions

    /// Compute the radial average of a quantitiy
    ///
    /// @param mf       MultiFab containing quantity to be averaged
    /// @param phibar   Averaged quantity
    /// @param comp     Index of component of `mf` to average
    void Average (const amrex::Vector<amrex::MultiFab>& mf,
                  RealVector& phibar,
                  int comp);

    // end MaestroAverage.cpp functions
    ////////////

    ////////////
    // MaestroCheckpoint.cpp functions

    /// Write a checkpoint at timestep `step`
    void WriteCheckPoint (int step);
    int ReadCheckPoint ();
    void GotoNextLine (std::istream& is);
    ////////////

    ////////////
    // MaestroConvert.cpp functions

    /// If `flag`, subtract the base state, returning the perturbed quantity.
    /// Otherwise, add the base state, returning the full quantity.
    /// This version iterates over all levels.
    ///
    /// @param scal     full/perturbed scalar quantity to subtract/add base state to
    /// @param s0       base state scalar
    /// @param comp     component of `scal` to perform calculation on
    /// @param bccomp   component of `bcs_in` to use to enforce boundary conditions
    /// @param bcs_in   boundary conditions
    /// @param flag     determines whether base state is subtracted (true) or added (false)
    void PutInPertForm (amrex::Vector<amrex::MultiFab>& scal,
                        const RealVector& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);

    /// If `flag`, subtract the base state, returning the perturbed quantity.
    /// Otherwise, add the base state, returning the full quantity.
    /// This version operates only on a single level.
    ///
    /// @param level    level to perform calculation on
    /// @param scal     full/perturbed scalar quantity to subtract/add base state to
    /// @param s0       base state scalar
    /// @param comp     component of `scal` to perform calculation on
    /// @param bccomp   component of `bcs_in` to use to enforce boundary conditions
    /// @param bcs_in   boundary conditions
    /// @param flag     determines whether base state is subtracted (true) or added (false)
    void PutInPertForm (int level, amrex::Vector<amrex::MultiFab>& scal,
                        const RealVector& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);

    /// If `flag`, returns species mass fraction `X` given the conserved variable `rhoX`.
    /// Otherwise, performs inverse operation
    void ConvertRhoXToX (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);

    /// If `flag`, return enthalpy `h` given the conserved variable `rhoh`.
    /// Otherwise, performs inverse operation
    void ConvertRhoHToH (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);
    ////////////

    ////////////////////////
    // MaestroDebug.cpp functions

    /// Print out the contents of a Vector of MultiFabs
    void PrintMF   (const amrex::Vector<amrex::MultiFab>& CC);


    /// Print out the contents of a Vector of edge-based MultiFabs
    void PrintEdge (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir);

    /// Utility to write out a multilevel multifab to a plotfile
    void WriteMF   (const amrex::Vector<amrex::MultiFab>& mf,
                    std::string name);
    ////////////////////////

    ////////////////////////
    // MaestroDensityAdvance.cpp functions

    /// Advance the density
    ///
    /// @param which_step       Is this the predictor (1) or corrector (2) step?
    /// @param scalold          old cell-centered scalars
    /// @param scalnew          new scalars
    /// @param sedge            edge-based scalars
    /// @param sflux            scalar fluxes
    /// @param scal_force       scalar force
    /// @param etarhoflux       `eta_rho` flux
    /// @param umac             MAC velocity
    /// @param w0mac            MAC base state velocity
    /// @param rho0_predicted_edge base state density predicted to cell edges
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

    /// Put together an array of multifabs for writing
    void WriteDiagFile (int& index);

    /// Write plotfile to disk
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
    // Time step computation in `MaestroDt.cpp`

    /// Compute the time step
    void EstDt ();

    /// Compute initial time step
    void FirstDt ();

    // end MaestroDt.cpp functions
    ////////////////////////

    ////////////////////////
    // MaestroEnthalpyAdvance.cpp functions

    /// Advance the enthalpy
    ///
    /// @param which_step       Is this the predictor (1) or corrector (2) step?
    /// @param scalold          old cell-centered scalars
    /// @param scalnew          new scalars
    /// @param sedge            edge-based scalars
    /// @param sflux            scalar fluxes
    /// @param scal_force       scalar force
    /// @param umac             MAC velocity
    /// @param w0mac            MAC base state velocity
    /// @param thermal          thermal term
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

    /// Call `FillPatch` for all levels
    void FillPatch (amrex::Real time,
                    amrex::Vector<amrex::MultiFab>& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in, int variable_type=0);

    /// Compute a new multifab by coping in phi from valid region and filling ghost cells
    /// - works for single level and 2-level cases
    /// (fill fine grid ghost by interpolating from coarse)
    /// - `srccomp` is the source component
    /// - `destcomp` is the destination component AND the bc component
    void FillPatch (int lev, amrex::Real time,
                    amrex::MultiFab& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in, int variable_type=0);

    /// Fill an entire multifab by interpolating from the coarser level
    /// - this comes into play when a new level of refinement appears
    /// - `srccomp` is the source component
    /// - `destcomp` is the destination component AND the bc component
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<amrex::MultiFab>& mf_old,
                          amrex::Vector<amrex::MultiFab>& mf_new,
                          int srccomp, int destcomp, int ncomp,
                          const amrex::Vector<amrex::BCRec>& bcs, int variable_type=0);

    /// Utility to copy in data from `mf_old` and/or `mf_new` into `mf`
    /// - if `time=t_old` we copy `mf_old` into` mf`
    /// - if `time=t_new` we copy `mf_new` into `mf`
    /// - otherwise copy in both `mf_old` and `mf_new` into `mf` and the `fillpatch`
    /// routines know to interpolate in time.
    // However in MAESTRO since we don't
    // subcycle I'm not sure if we need this capability?
    void GetData (int lev, amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& mf,
                  amrex::Vector<amrex::Real>& mftime,
                  amrex::Vector<amrex::MultiFab>& mf_old,
                  amrex::Vector<amrex::MultiFab>& mf_new);

    /// Set covered coarse cells to be the average of overlying fine cells
    ///
    /// @param mf       MultiFab to average
    /// @param comp     Index of first component to average
    /// @param ncomp    Number of components to average
    void AverageDown (amrex::Vector<amrex::MultiFab>& mf,
                      int comp,
                      int ncomp);

    /// Set covered faces to be the average of overlying fine faces
    void AverageDownFaces (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge);

    /// Fill in ONE ghost cell for all components of a face-centered (MAC) velocity
    /// field behind physical boundaries.  Does not modify the velocities on the boundary
    void FillUmacGhost (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int level=-1);

    /// Fill in all ghost cells for an edge-based MAC velocity field
    void FillPatchUedge(amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge);

    // end MaestroFillData.cpp functions
    ////////////

    ////////////////////////
    // MaestroFill3dData.cpp functions

    /// Maps 1d arrays onto multi-D cartesian MultiFabs
    ///
    /// @param s0           1d base state
    /// @param s0_cart      base state mapped to multi-d cartesian MultiFab
    /// @param is_input_edge_centered   is the input edge-centered?
    /// @param is_output_a_vector       is the output a vector?
    /// @param bcs          boundary conditions
    /// @param sbccomp      start boundary conditions component
    void Put1dArrayOnCart (const RealVector& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0, int variable_type=0);

    /// Maps 1d arrays onto multi-D cartesian MultiFabs
    ///
    /// @param level        AMR level to perform calculation on
    /// @param s0           1d base state
    /// @param s0_cart      base state mapped to multi-d cartesian MultiFab
    /// @param is_input_edge_centered   is the input edge-centered?
    /// @param is_output_a_vector       is the output a vector?
    /// @param bcs          boundary conditions
    /// @param sbccomp      start boundary conditions component
    void Put1dArrayOnCart (int level, const RealVector& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0);

    /// Add (`mult` times) the MAC-projected base state velocity to the edge-based
    /// velocity `uedge`
    ///
    /// @param uedge    edge based velocity
    /// @param w0mac    MAC base state velocity
    /// @param mult     multiplication factor
    void Addw0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                const amrex::Real& mult);

    /// MAC-project the base state velocity
    void MakeW0mac (amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac);

    /// MAC-project the base state scalar `s0`
    void MakeS0mac (const RealVector& s0,
		    amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& s0mac);

    /// Create the unit normal across the grids
    void MakeNormal ();

    /// Put the cell-centered data `s_cc` on faces by averaging adjacent cells
    void PutDataOnFaces(const amrex::Vector<amrex::MultiFab>& s_cc,
			amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM >>& face,
			int harmonic_avg);

    void MakeCCtoRadii ();
    // end MaestroFill3dData.cpp functions
    ////////////

    ////////////////////////
    // MaestroForce.cpp functions

    /// Calculate the velocity force term
    void MakeVelForce (amrex::Vector<amrex::MultiFab>& vel_force,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
                       const amrex::Vector<amrex::MultiFab>& rho,
                       const RealVector& rho0,
                       const RealVector& grav,
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

    /// Calculate the conserved enthalpy force term
    void MakeRhoHForce (amrex::Vector<amrex::MultiFab>& scal_force,
                        int is_prediction,
                        const amrex::Vector<amrex::MultiFab>& thermal,
                        const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int add_thermal,
			const int &which_step);
    ////////////////////////

    ////////////
    // MaestroGamma.cpp functions

    /// Calculate the horizontal average of \f$\Gamma_1\f$
    void MakeGamma1bar (const amrex::Vector<amrex::MultiFab>& scal,
                        RealVector& gamma1bar,
                        const RealVector& p0);

    // end MaestroGamma.cpp functions
    ////////////

    ////////////
    // MaestroInit.cpp functions

    /// fill in multifab and base state data
    void InitData ();

    /// During initialization of a simulation, `Maestro::InitData()` calls
    /// `AmrCore::InitFromScratch()`, which calls
    /// a `MakeNewGrids()` function that repeatedly calls this function to create
    /// finer levels.  This function creates a new fine
    /// level that did not exist before by interpolating from the coarser level
    /// overrides the pure virtual function in `AmrCore`
    virtual void MakeNewLevelFromScratch (int lev,
                                          amrex::Real time,
                                          const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    /// Performs the initial projection
    void InitProj ();

    /// Performs the divu iteration
    void DivuIter (int istep_divu_iter);

    /// Performs the initial iteration to initialize `gradpi`
    void InitIter ();

    // end MaestroInit.cpp functions
    ////////////

    ////////////
    // MaestroMacProj.cpp functions

    /// Do the MAC projection
    ///
    /// @param umac enters with face-centered, time-centered `Utilde^*` and should leave with `Utilde`
    /// @param macphi is the solution to the elliptic solve and
    ///   enters as either zero, or the solution to the predictor MAC projection
    /// @param macrhs enters as `beta0*(S-Sbar)`
    /// @param beta0 is a 1d cell-centered array
    void MacProj (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  amrex::Vector<amrex::MultiFab>& macphi,
                  const amrex::Vector<amrex::MultiFab>& macrhs,
                  const RealVector& beta0,
                  const int& is_predictor);

    /// Multiply (or divide) face-centered data by `beta0`
    ///
    /// @param edge         face-centered data
    /// @param beta0        cell-centered \f$\beta_0\f$
    /// @param beta0_edge   face-centered \f$\beta_0\f$
    /// @param mult_or_div  do we multiply or divide?
    void MultFacesByBeta0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge,
                           const RealVector& beta0,
                           const RealVector& beta0_edge,
                           const int& mult_or_div);

    /// Compute the RHS for the solve, `RHS = macrhs - div(beta0*umac)`
    ///
    /// @param solverrhs    RHS for the solve
    /// @param macrhs       `macrhs` term
    /// @param umac         MAC velocity
    void ComputeMACSolverRHS (amrex::Vector<amrex::MultiFab>& solverrhs,
                              const amrex::Vector<amrex::MultiFab>& macrhs,
                              const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);

    /// Average `bcoefs` at faces using inverse of `rho`
    ///
    /// @param facebcoef    face-centered `bcoefs`
    /// @param rhocc        cell-centered density
    void AvgFaceBcoeffsInv(amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM > >& facebcoef,
			   const amrex::Vector<amrex::MultiFab>& rhocc);

    /// Set boundaries for `LABecLaplacian` to solve `-div(B grad) phi = RHS`
    void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

    // end MaestroMacProj.cpp functions
    ////////////

    ////////////
    // MaestroMakeEta.cpp functions

    /// Compute `eta_rho` at edge- and cell-centers
    ///
    /// @param etarho_edge  face-centered \f$\eta_\rho\f$
    /// @param etarho_cell  cell-centered \f$\eta_\rho\f$
    /// @param etarho_flux  \f$\eta_\rho\f$ flux
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

    /// Compute S at cell-centers
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

    /// Compute `rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi`
    void MakeRHCCforNodalProj (amrex::Vector<amrex::MultiFab>& rhcc,
                               const amrex::Vector<amrex::MultiFab>& S_cc,
                               const RealVector& Sbar,
                               const RealVector& beta0,
                               const amrex::Vector<amrex::MultiFab>& delta_gamma1_term);

    /// Compute `rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi`
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

    /// Perform a nodal projection
    void NodalProj (int proj_type,
                    amrex::Vector<amrex::MultiFab>& rhcc,
                    int istep_divu_iter=0);

    /// Fill in `Vproj`
    ///
    /// `initial_projection_comp: Utilde^0                        -- uold`
    /// `divu_iters_comp:         Utilde^0                        -- uold`
    /// `pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt`
    /// `regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf`
    /// @param sig contains `rhohalf` if `proj_type == regular_timestep_comp`
    void CreateUvecForProj (int proj_type,
                            amrex::Vector<amrex::MultiFab>& Vproj,
                            const amrex::Vector<amrex::MultiFab>& sig);

    void SetBoundaryVelocity(amrex::Vector<amrex::MultiFab>& vel);

    /// Given a nodal `phi`, compute \f$\nabla(\phi)\f$ at cell centers
    void ComputeGradPhi(amrex::Vector<amrex::MultiFab>& phi,
                        amrex::Vector<amrex::MultiFab>& gphi);

    /// Average nodal `pi` to cell-centers and put in the `Pi` component of `snew`
    void MakePiCC(const amrex::Vector<amrex::MultiFab>& beta0_cart);

    // end MaestroNodalProj.cpp functions
    ////////////

    ////////////
    // MaestroPlot.cpp functions

    /// Get plotfile name
    void PlotFileName (const int lev, std::string* plotfilename);

    /// Put together an array of multifabs for writing
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

    /// Set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames (int * nPlot) const;

    /// Set small plotfile variables names
    amrex::Vector<std::string> SmallPlotFileVarNames (int * nPlot, amrex::Vector<std::string> varnames) const;

    /// Write a small plotfile to disk
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

    /// Write plotfile to disk
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

    /// Calculate the magnitude of the velocity
    void MakeMagvel (const amrex::Vector<amrex::MultiFab>& vel,
                     amrex::Vector<amrex::MultiFab>& magvel);

    /// Calculate the radial and circular components of the velocity
    void MakeVelrc (const amrex::Vector<amrex::MultiFab>& vel,
                    const amrex::Vector<amrex::MultiFab>& w0rcart,
                    amrex::Vector<amrex::MultiFab>& rad_vel,
                    amrex::Vector<amrex::MultiFab>& circ_vel);

    /// Calculate the adiabatic excess
    void MakeAdExcess (const amrex::Vector<amrex::MultiFab>& state,
                       amrex::Vector<amrex::MultiFab>& ad_excess);

    /// Calculate the gravitational acceleration
    void MakeGrav (const RealVector& rho0,
                   amrex::Vector<amrex::MultiFab>& grav);

    /// Calculate the vorticity
    void MakeVorticity (const amrex::Vector<amrex::MultiFab>& vel,
                        amrex::Vector<amrex::MultiFab>& vorticity);

    /// Calculate `deltagamma`
    void MakeDeltaGamma (const amrex::Vector<amrex::MultiFab>& state,
                         const RealVector& p0,
                         const amrex::Vector<amrex::MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const amrex::Vector<amrex::MultiFab>& gamma1bar_cart,
                         amrex::Vector<amrex::MultiFab>& deltagamma);

    /// Calculate the entropy
    void MakeEntropy (const amrex::Vector<amrex::MultiFab>& state,
                      amrex::Vector<amrex::MultiFab>& entropy);

    /// Calculate the divergenve of the base state velocity
    void MakeDivw0 (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                    amrex::Vector<amrex::MultiFab>& divw0);

    /// Calculate `pi` times the divergence of the velocity
    void MakePiDivu (const amrex::Vector<amrex::MultiFab>& vel,
                   const amrex::Vector<amrex::MultiFab>& pi_cc,
                   amrex::Vector<amrex::MultiFab>& pidivu);

    /// Mass fractions of the species
    void MakeAbar(const amrex::Vector<amrex::MultiFab>& state,
                amrex::Vector<amrex::MultiFab>& abar);

    // end MaestroPlot.cpp functions
    ////////////

    ////////////
    // MaestroReact.cpp functions

    /// Compute heating term, `rho_Hext`, then
    /// react the state over `dt_react` and update `rho_omegadot`, `rho_Hnuc`
    void React (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const RealVector& p0,
                const amrex::Real dt_react,
		const amrex::Real time_in);

    void Burner (const amrex::Vector<amrex::MultiFab>& s_in,
		 amrex::Vector<amrex::MultiFab>& s_out,
		 const amrex::Vector<amrex::MultiFab>& rho_Hext,
		 amrex::Vector<amrex::MultiFab>& rho_omegadot,
		 amrex::Vector<amrex::MultiFab>& rho_Hnuc,
		 const RealVector& p0,
		 const amrex::Real dt_react,
		 const amrex::Real time_in);

    /// Compute heating term, `rho_Hext`
    void MakeHeating (amrex::Vector<amrex::MultiFab>& rho_Hext,
                      const amrex::Vector<amrex::MultiFab>& scal);

    // end MaestroReact.cpp functions
    ////////////

    ////////////
    // regridding functions in MaestroRegrid.cpp

    /// Check to see if we need to regrid, then regrid
    void Regrid ();

    /// Set tagging array to include buffer zones for multilevel
    void TagArray ();

    /// Tag all cells for refinement
    ///
    /// Overrides the pure virtual function in `AmrCore`
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ng) override;

    /// Within a call to `AmrCore::regrid`, this function fills in data at a level
    /// that existed before, using pre-existing fine and interpolated coarse data
    ///
    /// overrides the pure virtual function in `AmrCore`
    virtual void RemakeLevel (int lev,
                              amrex::Real time,
                              const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    /// Within a call to `AmrCore::regrid`, this function fills in data at a level
    /// that did NOT exist before, using interpolated coarse data
    ///
    /// overrides the pure virtual function in `AmrCore`
    virtual void MakeNewLevelFromCoarse (int lev,
                                         amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    /// Within a call to `AmrCore::regrid`, this function deletes all data
    /// at a level of refinement that is no longer needed
    ///
    /// overrides the pure virtual function in `AmrCore`
    virtual void ClearLevel (int lev) override;

    // end regridding functions
    ////////////

    ////////////
    // MaestroRhoHT.cpp functions

    /// Calculate the temperature given the density and the enthalpy
    ///
    /// @param scal     scalars
    /// @param p0       base state pressure
    void TfromRhoH (amrex::Vector<amrex::MultiFab>& scal,
                    const RealVector& p0);

    /// Calculate the temperature given the density and the pressure
    ///
    /// @param scal     scalars
    /// @param p0       base state pressure
    void TfromRhoP (amrex::Vector<amrex::MultiFab>& scal,
                    const RealVector& p0,
                    int updateRhoH = 0);

    /// Calculate the pressure given the density and the enthalpy
    ///
    /// @param state    scalars
    /// @param s_old    scalars at old time step
    /// @param peos     pressure calculated from the equation of state
    void PfromRhoH (const amrex::Vector<amrex::MultiFab>& state,
                    const amrex::Vector<amrex::MultiFab>& s_old,
		    amrex::Vector<amrex::MultiFab>& peos);

    /// Calculate the Mach number given the density and the enthalpy
    ///
    /// @param scal     scalars
    /// @param vel      velocity
    /// @param p0       base state pressure
    /// @param mach     Mach number
    void MachfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                           const amrex::Vector<amrex::MultiFab>& vel,
                           const RealVector& p0,
                           const amrex::Vector<amrex::MultiFab>& w0cart,
                           amrex::Vector<amrex::MultiFab>& mach);

    /// Calculate the sound speed given the density and the enthalpy
    ///
    /// @param scal     scalars
    /// @param p0       base state pressure
    /// @param p0cart   base state pressure on cartesian grid
    /// @param cs       sound speed
    void CsfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                     const RealVector& p0,
                     const amrex::Vector<amrex::MultiFab>& p0cart,
                     amrex::Vector<amrex::MultiFab>& cs);
    ////////////

    ////////////
    // MaestroSetup.cpp functions

    void ExternInit ();

    /// Read in C++ parameters from inputs file
    void ReadParameters ();

    /// Define variable mappings (`Rho`, `RhoH`, ..., `Nscal`, etc.)
    void VariableSetup ();

    /// Set up `BCRec` definitions for BC types
    void BCSetup ();

    // end MaestroSetup.cpp functions
    ////////////

    ////////////
    // MaestroSponge.cpp functions

    void MakeSponge (amrex::Vector<amrex::MultiFab>& sponge);
    ////////////

    ////////////////////////
    // MaestroThermal.cpp functions

    /// Compute the quantity: `thermal = del dot kappa grad T`
    ///
    /// - if `temp_diffusion_formulation = 1`, then we compute this directly.
    /// - if `temp_diffusion_formulation = 2`, then we compute the algebraically
    ///     equivalent form with `grad h` - `grad X_k` - `grad p_0` formulation
    ///
    /// @param thermal  thermal term
    /// @param scal     scalars
    /// @param Tcoeff   temperature coefficient
    /// @param hcoeff   enthalpy coefficient
    /// @param Xkcoeff  species coefficients
    /// @param pcoeff   pressure coefficient
    /// @param p0       base state pressure
    void MakeExplicitThermal (amrex::Vector<amrex::MultiFab>& thermal,
			      const amrex::Vector<amrex::MultiFab>& scal,
			      const amrex::Vector<amrex::MultiFab>& Tcoeff,
			      const amrex::Vector<amrex::MultiFab>& hcoeff,
			      const amrex::Vector<amrex::MultiFab>& Xkcoeff,
			      const amrex::Vector<amrex::MultiFab>& pcoeff,
			      const RealVector& p0,
			      int temp_formulation);

    /// Use `apply()` to construct the form of the conduction term.
    /// `apply()` forms the generic quantity:
    ///
    ///   `(alpha * A - beta * div B grad) phi = RHS`
    void ApplyThermal (amrex::MLABecLaplacian& mlabec,
		       amrex::Vector<amrex::MultiFab>& thermalout,
		 const amrex::Vector<amrex::MultiFab>& coeff,
		       amrex::Vector<amrex::MultiFab>& phi,
		 const amrex::Vector<amrex::BCRec>& bcs,
		 int bccomp);

    /// create the coefficients for `grad{T}`, `grad{h}`, `grad{X_k}`, and `grad{p_0}`
    /// for the thermal diffusion term in the enthalpy equation.
    ///
    /// note: we explicitly fill the ghostcells by looping over them directly
    ///
    /// @param scal     scalars
    /// @param Tcoeff   temperature coefficient
    /// @param hcoeff   enthalpy coefficient
    /// @param Xkcoeff  species coefficients
    /// @param pcoeff   pressure coefficient
    void MakeThermalCoeffs(const amrex::Vector<amrex::MultiFab>& scal,
			   amrex::Vector<amrex::MultiFab>& Tcoeff,
			   amrex::Vector<amrex::MultiFab>& hcoeff,
			   amrex::Vector<amrex::MultiFab>& Xkcoeff,
			   amrex::Vector<amrex::MultiFab>& pcoeff);

    /// ThermalConduct implements thermal diffusion in the enthalpy equation.
    /// This is an implicit solve, using the multigrid solver.  This updates
    /// the enthalpy only.
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

    /// Advance the velocity
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
    // static data members defined in `_cpp_parameters`
    // this file is in Source/param_includes/
#include <maestro_params.H>
    ////////////////

    // these are defined in `VariableSetup()`
    static int Rho, RhoH, FirstSpec, NumSpec, Temp, Pi, Nscal;
    static int ng_s; //< number of ghost cells for `sold/new` and `uold/new`

    // these are defined in `Maestro.cpp`
    // these represent the four types of nodal projections
    static int initial_projection_comp;
    static int divu_iters_comp;
    static int pressure_iters_comp;
    static int regular_timestep_comp;

    // these are defined in `Maestro.cpp`
    // density and enthalpy edge state prediction
    static int predict_rhoprime_and_X, predict_rhoX, predict_rho_and_X;
    static int predict_rhoh, predict_rhohprime, predict_h;
    static int predict_T_then_rhohprime, predict_T_then_h;
    static int predict_hprime, predict_Tprime_then_h;

    // these are defined in `Maestro.cpp`
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;

    /// which step?
    int istep;
    int start_step;

    // keep track of old time, new time, and time step at each level
    amrex::Real t_new;
    amrex::Real t_old;
    amrex::Real dt;
    amrex::Real dtold;

    /// number of ghost cells needed for hyperbolic step
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
    amrex::Vector<amrex::MultiFab> pi;   //< nodal

    /// this doesn't have to be persistent, but we make it so that we avoid
    /// continually creating and filling temporaries
    /// saves on some flops and data movement (GPU)
    amrex::Vector<amrex::MultiFab> w0_cart;

    /// this only needs to persist leading into the initial pressure iters
    /// since we project `(beta0^nph S^1 - beta0 S^0) / dt`
    /// during a regular time step we overwrite this
    amrex::Vector<amrex::MultiFab> rhcc_for_nodalproj;

    /// spherical only -
    /// we make this persistent in that we only have to rebuild and
    /// fill this after regridding
    amrex::Vector<amrex::MultiFab> normal;
    amrex::Vector<amrex::MultiFab> cell_cc_to_r;

    /// Stores domain boundary conditions.
    /// These muse be vectors (rather than arrays) so we can ParmParse them
    amrex::Vector<int> phys_bc;

    /// Boundary condition objects needed for FillPatch routines.
    /// This is essentially an array (over components)
    /// of `2*DIM` integer arrays storing the physical boundary
    /// condition types at the `lo/hi` walls in each direction
    amrex::Vector<amrex::BCRec> bcs_s;  //< for `sold/snew`
    amrex::Vector<amrex::BCRec> bcs_u;  //< for `uold/unew`
    amrex::Vector<amrex::BCRec> bcs_f;  //< "first-order extrap" used for `S_cc`, `gpi`, `dSdt`

    // vectors store the multilevel 1D states as one very long array.
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
    /// `p0^{n-1}` needed to compute `d(p0)/dt` for nonuniform grid spacing
    RealVector p0_nm1;
    RealVector beta0_nm1;

    // Vectors store the multilevel 1D states as one very long array.
    // These are edge-centered
    RealVector w0;
    RealVector etarho_ec;
    RealVector r_edge_loc;

    // Information for the base state.
    // We also store these (except for the large arrays)
    // in `base_state.f90` in a fortran module for convenience
    int max_radial_level;
    int finest_radial_level;
    int nr_fine;
    amrex::Real dr_fine;
    int nr_irreg;

    // these provide information about the multilevel base state configuration
    amrex::Vector<int> numdisjointchunks;
    amrex::Vector<int> r_start_coord;
    amrex::Vector<int> r_end_coord;

    /// array of tagged boxes (planar)
    amrex::Vector<int> tag_array;

    // diag file array buffers
    amrex::Vector<amrex::Real> diagfile1_data;
    amrex::Vector<amrex::Real> diagfile2_data;
    amrex::Vector<amrex::Real> diagfile3_data;

    /// Stores fluxes at coarse-fine interface for synchronization.
    /// This will be sized `max_level+1`
    ///
    /// NOTE: the flux register associated with `flux_reg[lev]` is associated
    /// with the `lev/lev-1` interface (and has grid spacing associated with `lev-1`)
    /// therefore `flux_reg[0]` and `flux_reg[max_level]` are never actually
    /// used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;

    /// flag for writing plotfiles
    enum plotfile_flag {plotInitData = -9999999, plotInitProj = -9999998, plotDivuIter = -9999997};

};

#endif
