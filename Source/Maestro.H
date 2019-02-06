#ifndef Maestro_H_
#define Maestro_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FMultiGrid.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <PhysBCFunctMaestro.H>
#include <Maestro_F.H>

    ///
    /// @class Maestro
    ///
    /// @brief
    ///
class Maestro
    : public amrex::AmrCore
{
public:

    /*
      public member functions
    */

    ///
    /// constructor
    ///
    Maestro ();

    ///
    /// destructor
    ///
    virtual ~Maestro ();

    ///
    /// in ``MaestroSetup.cpp``
    ///     - read in C++/F90 parameters
    ///     - define global C++/F90 variables and initialize network
    ///     - set up boundary conditions
    ///     - initialize base state geometry parameters
    ///     - set ``istep``, ``t_new``, ``t_old``
    ///     - allocate MultiFabs and base state arrays
    ///
    void Setup ();

    ///
    /// in ``MaestroInit.cpp``
    ///     - initialize multifab and base state data
    ///     - perform initial projection
    ///     - perform divu iters
    ///     - perform initial (pressure) iterations
    ///
    void Init ();

    ///
    /// in ``MaestroEvolve.cpp``
    ///
    /// advance solution to final time
    ///
    void Evolve ();


    ///
    /// for keeping track of the amount of CPU time used - this will persist
    /// after restarts
    ///
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

    static amrex::Real getCPUTime();

private:

    /*
      private member functions
    */

    ////////////
    // MaestroAdvance.cpp functions

    ///
    /// advance solution at all levels for a single time step
    ///
    /// @param is_initIter
    ///
    void AdvanceTimeStep (bool is_initIter);

    ///
    /// advance solution for a single time step with irregular base state spacing
    ///
    /// @param is_initIter
    ///
    void AdvanceTimeStepIrreg (bool is_initIter);

    ///
    /// advance solution for a single time step with regular base state spacing
    /// and new time-stepping scheme
    ///
    /// @param is_initIter
    ///
    void AdvanceTimeStepAverage (bool is_initIter);

    // end MaestroAdvance.cpp functions
    ////////////

    ////////////
    // MaestroAdvection.cpp functions

    ///
    /// compute unprojected mac velocities
    ///
    /// @param umac
    /// @param w0mac            MAC base state velocity
    /// @param w0_force
    /// @param w0_force_cart
    ///
    void AdvancePremac (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const amrex::Vector<amrex::Real>& w0_force,
			const amrex::Vector<amrex::MultiFab>& w0_force_cart,
			const amrex::Vector<amrex::Real>& beta0,
			const int is_predictor);


    ///
    /// Compute the transverse velocity components.
    ///
    /// @param utilde       Perturbational velocity
    /// @param ufull        Full velocity
    /// @param utrans       Transverse velocity
    /// @param w0mac        Base state MAC velocity
    ///
    void MakeUtrans (const amrex::Vector<amrex::MultiFab>& utilde,
                     const amrex::Vector<amrex::MultiFab>& ufull,
                     amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);


    ///
    /// Predict the normal velocities to the interfaces.
    ///
    /// @param utilde   Perturbational velocity
    /// @param ufull    Full velocity
    /// @param utrans   Transverse velocity
    /// @param umac     MAC velocity
    /// @param w0mac    Base state MAC velocity
    /// @param force
    ///
    void VelPred (const amrex::Vector<amrex::MultiFab>& utilde,
                  const amrex::Vector<amrex::MultiFab>& ufull,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                  amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const amrex::Vector<amrex::MultiFab>& force);


    ///
    /// Predict scalars to edges
    ///
    /// @param state            Cell-centered state
    /// @param sedge            Edge-based scalars
    /// @param umac             MAC velocity
    /// @param force
    /// @param is_vel
    /// @param bcs              Boundary conditions
    /// @param nbccomp          Number of components of `bcs` to use
    /// @param start_scomp      Component of `state` to start with
    /// @param start_bccomp     Component of `bcs` to start with
    /// @param num_comp         Number of `state` components to predict
    /// @param is_conservative
    ///
    void MakeEdgeScal (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                       const amrex::Vector<amrex::MultiFab>& force,
                       int is_vel, const amrex::Vector<amrex::BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative);


    ///
    /// Make the conservative flux for the species evolution equation
    ///
    /// @param state        State
    /// @param sflux        Flux
    /// @param etarhoflux
    /// @param sedge        Edge-based state
    /// @param umac         MAC velocity
    /// @param w0mac        MAC base state velocity
    /// @param r0_old
    /// @param r0_edge_old
    /// @param r0mac_old
    /// @param r0_new
    /// @param r0_edge_new
    /// @param r0mac_new
    /// @param r0_predicted_edge
    /// @param start_comp   Component of state to begin with
    /// @param num_comp     Number of components to use
    ///
    void MakeRhoXFlux (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		       amrex::Vector<amrex::MultiFab>& etarhoflux,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
		       const amrex::Vector<amrex::Real>& r0_old,
		       const amrex::Vector<amrex::Real>& r0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_old,
		       const amrex::Vector<amrex::Real>& r0_new,
		       const amrex::Vector<amrex::Real>& r0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_new,
		       const amrex::Vector<amrex::Real>& r0_predicted_edge,
		       int start_comp, int num_comp);


    ///
    /// Make the conservative flux for the enthalpy equation
    ///
    /// @param state            State
    /// @param sflux            Flux
    /// @param sedge            Edge-based state
    /// @param umac             MAC velocity
    /// @param w0mac            MAC base state velocity
    /// @param r0_old
    /// @param r0_edge_old
    /// @param r0mac_old
    /// @param r0_new
    /// @param r0_edge_new
    /// @param r0mac_new
    /// @param rh0_old
    /// @param rh0_edge_old
    /// @param rh0mac_old
    /// @param rh0_new
    /// @param rh0_edge_new
    /// @param rh0mac_new
    /// @param h0mac_old
    /// @param h0mac_new
    ///
    void MakeRhoHFlux (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
		       const amrex::Vector<amrex::Real>& r0_old,
		       const amrex::Vector<amrex::Real>& r0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_old,
		       const amrex::Vector<amrex::Real>& r0_new,
		       const amrex::Vector<amrex::Real>& r0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& r0mac_new,
		       const amrex::Vector<amrex::Real>& rh0_old,
		       const amrex::Vector<amrex::Real>& rh0_edge_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& rh0mac_old,
		       const amrex::Vector<amrex::Real>& rh0_new,
		       const amrex::Vector<amrex::Real>& rh0_edge_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& rh0mac_new,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& h0mac_old,
		       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& h0mac_new);

    ///
    /// Advance scalars to new timestep
    ///
    /// @param stateold     State at old timestep
    /// @param statenew     State at new timestep
    /// @param sflux        Flux
    /// @param force
    /// @param start_scomp  Index of first scalar component
    /// @param num_comp     Number of scalar components
    /// @param p0           Base state pressure
    /// @param p0_cart      Base state pressure on cartesian grid
    ///
    void UpdateScal (const amrex::Vector<amrex::MultiFab>& stateold,
		     amrex::Vector<amrex::MultiFab>& statenew,
		     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		     const amrex::Vector<amrex::MultiFab>& force,
		     int start_scomp, int num_comp,
		     const amrex::Real* p0 = NULL,
		     const amrex::Vector<amrex::MultiFab>& p0_cart = {});


    ///
    /// Advance velocity to new timestep
    ///
    /// @param umac     MAC velocity
    /// @param uedge    Edge-based velocity
    /// @param force
    /// @param sponge
    /// @param w0mac    MAC base state velocity
    ///
    void UpdateVel (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		    const amrex::Vector<amrex::MultiFab>& force,
		    const amrex::Vector<amrex::MultiFab>& sponge,
		    const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);
    ////////////

    ////////////
    // MaestroAverage.cpp functions

    ///
    /// Given a multifab of data (``mf``), average down to a base state quantity, ``phibar``.
    /// If we are in plane-parallel, the averaging is at constant height.
    /// If we are spherical, then the averaging is done at constant radius.
    ///
    /// @param mf       Multifab of data to be averaged
    /// @param phibar   Averaged, base state quantity
    /// @param comp     Component of ``mf`` to be averaged
    ///
    void Average (const amrex::Vector<amrex::MultiFab>& mf,
                  amrex::Vector<amrex::Real>& phibar,
                  int comp);

    // end MaestroAverage.cpp functions
    ////////////

    ////////////
    // MaestroCheckpoint.cpp functions

    ///
    /// @param step
    ///
    void WriteCheckPoint (int step);
    int ReadCheckPoint ();

    ///
    /// utility to skip to next line in Header
    ///
    /// @param is   std::istream object
    ///
    void GotoNextLine (std::istream& is);
    ////////////

    ////////////
    // MaestroConvert.cpp functions

    ///
    /// Given a scalar, subtract/add the base state (``flag == true/false``) to
    /// return the perturbed/full state.
    ///
    /// @param scal     MultiFab containing scalar
    /// @param s0       Base state scalar
    /// @param comp     Component of ``scal`` to use
    /// @param bccomp   Component of ``bcs_in`` to use
    /// @param bcs_in   Boundary conditions
    /// @param flag     Should we subtract (true) or add (false) the base state?
    ///
    void PutInPertForm (amrex::Vector<amrex::MultiFab>& scal,
                        const amrex::Vector<amrex::Real>& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);

    void PutInPertForm (int level, amrex::Vector<amrex::MultiFab>& scal,
                        const amrex::Vector<amrex::Real>& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);


    ///
    /// Given the conserved species mass fraction, \f$ \rho X \f$, return the
    /// primitive species fraction \f$ X \f$. If flag is set to ``false``, then this
    /// does the inverse, multiplying \f$ X \f$ by the density to return \f$ \rho X \f$
    ///
    /// @param scal     MultiFab of scalars
    /// @param flag     Should we divide (true) or multiply (false) by rho?
    ///
    void ConvertRhoXToX (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);


    ///
    /// Given the conserved variable \f$ \rho h \f$, divides by the density to return
    /// the primitive enthalpy \f$ h \f$. If the flag is set to ``false``, it does
    /// the inverse, instead multiplying by the density.
    ///
    /// @param scal     MultiFab of scalars
    /// @param flag     Should we divide (true) or multiply (false) by rho?
    ///
    void ConvertRhoHToH (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);
    ////////////

    ////////////////////////
    // MaestroDebug.cpp functions

    ///
    /// print out the contents of a Vector of MultiFabs. Assumes the data is cell-centered.
    ///
    /// @param CC   Vector of MultiFabs to print
    ///
    void PrintMF   (amrex::Vector<amrex::MultiFab>& CC);

    ///
    /// print out the contents of a Vector of edge-based MultiFabs
    ///
    /// @param EDGE     Vector of edge-based MultiFabs
    /// @param dir      direction along which data has been predicted to edges
    ///
    void PrintEdge (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir);

    ///
    /// utility to write out a multilevel multifab to a plotfile
    ///
    /// @param mf       MultiFab to print
    /// @param name     directory in which to print
    ///
    void WriteMF   (amrex::Vector<amrex::MultiFab>& mf,
                    std::string name);
    ////////////////////////

    ////////////////////////
    // MaestroDensityAdvance.cpp functions

    ///
    /// Advance the density by ``dt``
    ///
    /// @param which_step           is this is predictor (1) or corrector (2) step?
    /// @param scalold              Scalars at old timestep
    /// @param scalnew              Scalars at new timestep
    /// @param sedge                Edge based scalars
    /// @param sflux                Scalar flux
    /// @param scal_force
    /// @param etarhoflux
    /// @param umac                 MAC velocity
    /// @param w0mac                MAC base velocity
    /// @param rho0_predicted_edge  base state density predicted to edges
    ///
    void DensityAdvance (int which_step,
                         amrex::Vector<amrex::MultiFab>& scalold,
                         amrex::Vector<amrex::MultiFab>& scalnew,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
                         amrex::Vector<amrex::MultiFab>& scal_force,
			 amrex::Vector<amrex::MultiFab>& etarhoflux,
                         amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			 const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
			 const amrex::Vector<amrex::Real>& rho0_predicted_edge);
    ////////////////////////

    ////////////
    // MaestroDiag.cpp functions
    ///
    /// put together an array of multifabs for writing
    ///
    /// @param index
    ///
    void WriteDiagFile (int& index);

    ///
    /// write plotfile to disk
    ///
    /// @param step      Index of current timestep
    /// @param dt_in     Timestep
    /// @param rho0_in   Base state density
    /// @param p0        Base state pressure
    /// @param u_in      Velocity
    /// @param s_in      Scalars
    /// @param index
    ///
    void DiagFile (const int step,
		   const amrex::Real dt_in,
		   const amrex::Vector<amrex::Real>& rho0_in,
		   const amrex::Vector<amrex::Real>& p0_in,
		   const amrex::Vector<amrex::MultiFab>& u_in,
		   const amrex::Vector<amrex::MultiFab>& s_in,
		   int& index);
    // end MaestroDiag.cpp functions
    ////////////

    ////////////
    // Time step computation in MaestroDt.cpp

    ///
    /// compute the time step
    ///
    void EstDt ();

    ///
    /// compute initial time step
    ///
    void FirstDt ();

    // end MaestroDt.cpp functions
    ////////////////////////

    ////////////////////////
    // MaestroEnthalpyAdvance.cpp functions

    ///
    /// Advance enthalpy by ``dt`` to new time
    ///
    /// @param which_step   Is the this predictor step (1) or the corrector step (2)?
    /// @param scalold      Old scalars
    /// @param scalnew      New scalars
    /// @param sedge        Edge-based scalars
    /// @param sflux        Scalar flux
    /// @param scal_force   Scalar force terms
    /// @param umac         MAC velocity
    /// @param w0mac        MAC base velocity
    /// @param thermal      Thermal terms
    ///
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

    ///
    /// call FillPatch for all levels
    ///
    /// @param time         Current time
    /// @param mf
    /// @param mf_old
    /// @param mf_new
    /// @param srccomp      the source component
    /// @param destcomp     the destination component AND the bc component
    /// @param ncomp        Number of components
    /// @param startbccomp
    /// @param bcs_in       boundary conditions
    ///
    void FillPatch (amrex::Real time,
                    amrex::Vector<amrex::MultiFab>& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in);

    ///
    /// Compute a new multifab by coping in phi from valid region and filling ghost cells.
    /// Works for single level and 2-level cases,
    /// (fill fine grid ghost by interpolating from coarse)
    ///
    /// @param lev
    /// @param time         Current time
    /// @param mf
    /// @param mf_old
    /// @param mf_new
    /// @param srccomp      the source component
    /// @param destcomp     the destination component AND the bc component
    /// @param ncomp        Number of components
    /// @param startbccomp
    /// @param bcs_in       boundary conditions
    ///
    void FillPatch (int lev, amrex::Real time,
                    amrex::MultiFab& mf,
                    amrex::Vector<amrex::MultiFab>& mf_old,
                    amrex::Vector<amrex::MultiFab>& mf_new,
                    int srccomp, int destcomp, int ncomp, int startbccomp,
                    const amrex::Vector<amrex::BCRec>& bcs_in);

    ///
    /// fill an entire multifab by interpolating from the coarser level.
    /// This comes into play when a new level of refinement appears
    ///
    /// @param lev
    /// @param time         Current time
    /// @param mf
    /// @param mf_old
    /// @param mf_new
    /// @param srccomp      the source component
    /// @param destcomp     the destination component AND the bc component
    /// @param ncomp        Number of components
    /// @param bcs          boundary conditions
    ///
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf,
                          amrex::Vector<amrex::MultiFab>& mf_old,
                          amrex::Vector<amrex::MultiFab>& mf_new,
                          int srccomp, int destcomp, int ncomp,
                          const amrex::Vector<amrex::BCRec>& bcs);

    ///
    /// Utility to copy in data from ``mf_old`` and/or ``mf_new`` into ``mf``
    ///     - if ``time``=``t_old``, we copy ``mf_old`` into ``mf``
    ///     - if ``time``=``t_new``, we copy ``mf_new`` into ``mf``
    ///     - otherwise copy both ``mf_old`` and ``mf_new`` into ``mf`` and the fillpatch
    ///       routines know to interpolate in time.  However in MAESTRO since we don't
    ///       subcycle I'm not sure if we need this capability?
    ///
    /// @param lev
    /// @param time    Current time
    /// @param mf
    /// @param mftime
    /// @param mf_old
    /// @param mf_new
    ///
    void GetData (int lev, amrex::Real time,
                  amrex::Vector<amrex::MultiFab*>& mf,
                  amrex::Vector<amrex::Real>& mftime,
                  amrex::Vector<amrex::MultiFab>& mf_old,
                  amrex::Vector<amrex::MultiFab>& mf_new);

    ///
    /// set covered coarse cells to be the average of overlying fine cells
    ///
    /// @param mf       MultiFab to average
    /// @param comp     Index of first component to average down
    /// @param ncomp    Number of components to average down
    ///
    void AverageDown (amrex::Vector<amrex::MultiFab>& mf,
                      int comp,
                      int ncomp);

    ///
    /// set covered faces to be the average of overlying fine faces
    ///
    /// @param edge     MultiFab of edge-based data to average down
    ///
    void AverageDownFaces (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge);

    ///
    /// fill in ONE ghost cell for all components of a face-centered (MAC) velocity
    /// field behind physical boundaries.  Does not modify the velocities on the boundary
    ///
    /// @param umac     MAC velocity
    /// @param level
    ///
    void FillUmacGhost (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int level=-1);

    ///
    /// fill in all ghost cells for an edge-based MAC velocity field
    ///
    /// @param uedge        edge-based MAC velocity field
    ///
    void FillPatchUedge(amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge);

    // end MaestroFillData.cpp functions
    ////////////

    ////////////////////////
    // MaestroFill3dData.cpp functions

    void Put1dArrayOnCart (const amrex::Vector<amrex::Real>& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0);

    void Put1dArrayOnCart (int level, const amrex::Vector<amrex::Real>& s0,
                           amrex::Vector<amrex::MultiFab>& s0_cart,
                           int is_input_edge_centered,
                           int is_output_a_vector,
                           const amrex::Vector<amrex::BCRec>& bcs = amrex::Vector<amrex::BCRec>(),
                           int sbccomp = 0);


    ///
    /// @param uedge    edge-based MAC velocity field
    /// @param w0mac    MAC base state velocity
    /// @param mult
    ///
    void Addw0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                const amrex::Real& mult);


    ///
    /// @param w0mac    MAC base state velocity
    ///
    void MakeW0mac (amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac);


    ///
    /// @param s0       Base state scalars
    /// @param s0mac    MAC base state scalars
    ///
    void MakeS0mac (const amrex::Vector<amrex::Real>& s0,
		    amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& s0mac);

    /// Computes ``normal`` on all levels.
    ///
    void MakeNormal ();


    ///
    /// Predict cell-centered data to faces
    ///
    /// @param s_cc         Cell-centered data
    /// @param face         Edge-based data
    /// @param harmonic_avg If (1), use harmonic average on faces, otherwise use mean
    ///
    void PutDataOnFaces(const amrex::Vector<amrex::MultiFab>& s_cc,
			amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM >>& face,
			int harmonic_avg);

    /// map cell centers to base state indices
    ///
    void MakeCCtoRadii ();
    // end MaestroFill3dData.cpp functions
    ////////////

    ////////////////////////
    // MaestroForce.cpp functions

    ///
    /// Make force terms for momentum equation
    ///
    /// @param vel_force            Velocity force terms
    /// @param uedge                edge-based MAC velocity field
    /// @param rho                  Density
    /// @param rho0                 Base state density
    /// @param grav                 Gravitational acceleration
    /// @param w0_force
    /// @param w0_force_cart
    /// @param do_add_utilde_force
    ///
    void MakeVelForce (amrex::Vector<amrex::MultiFab>& vel_force,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
                       const amrex::Vector<amrex::MultiFab>& rho,
                       const amrex::Vector<amrex::Real>& rho0,
                       const amrex::Vector<amrex::Real>& grav,
                       const amrex::Vector<amrex::Real>& w0_force,
		       const amrex::Vector<amrex::MultiFab>& w0_force_cart,
		       const amrex::Vector<amrex::Real>& beta0,
		       const int is_predictor,
                       int do_add_utilde_force);


    ///
    /// Add velocity divergence terms to scalar force
    ///
    /// @param scal_force       Scalar force terms
    /// @param state            Scalar
    /// @param umac             MAC velocity
    /// @param s0               Base state
    /// @param s0_edge          Base state predicted to edges
    /// @param s0_cart          Base state on cartesian grid
    /// @param comp             Component of state/scalar force to use
    /// @param bcs              Boundary conditions
    /// @param fullform         Use full scalar (1) or do perturbed and base components separately
    ///
    void ModifyScalForce(amrex::Vector<amrex::MultiFab>& scal_force,
			 const amrex::Vector<amrex::MultiFab>& state,
                         const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                         const amrex::Vector<amrex::Real>& s0,
                         const amrex::Vector<amrex::Real>& s0_edge,
			 const amrex::Vector<amrex::MultiFab>& s0_cart,
                         int comp,
                         const amrex::Vector<amrex::BCRec>& bcs,
                         int fullform);


    ///
    /// Compute the source terms for the non-reactive part of the enthalpy equation \f$ w dp_0/dr \f$
    ///
    /// @param scal_force       Force terms
    /// @param is_prediction    Predictor step (1) or final update (0)?
    /// @param thermal          Thermal term
    /// @param umac             MAC velocity
    /// @param add_thermal      Do we add the thermal term (1)?
    /// @param which_step
    ///
    void MakeRhoHForce (amrex::Vector<amrex::MultiFab>& scal_force,
                        int is_prediction,
                        const amrex::Vector<amrex::MultiFab>& thermal,
                        const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int add_thermal,
			const int &which_step);
    ////////////////////////

    ////////////
    // MaestroGamma.cpp functions


    ///
    /// Compute gamma1 for the full state then average to create gamma1bar
    ///
    /// @param scal         Scalar state
    /// @param gamma1bar    Gamma1bar
    /// @param p0           Base state pressure
    ///
    void MakeGamma1bar (const amrex::Vector<amrex::MultiFab>& scal,
                        amrex::Vector<amrex::Real>& gamma1bar,
                        const amrex::Vector<amrex::Real>& p0);

    // end MaestroGamma.cpp functions
    ////////////

    ////////////
    // MaestroInit.cpp functions

    /// fill in multifab and base state data
    ///
    void InitData ();

    ///
    /// During initialization of a simulation, Maestro::InitData() calls
    /// AmrCore::InitFromScratch(), which calls
    /// a MakeNewGrids() function that repeatedly calls this function to create
    /// finer levels.  This function creates a new fine
    /// level that did not exist before by interpolating from the coarser level
    /// overrides the pure virtual function in AmrCore
    ///
    /// @param lev     AMR level index
    /// @param time    Current time
    /// @param ba      ``BoxArray``
    /// @param dm      ``DistributionMapping``
    ///
    virtual void MakeNewLevelFromScratch (int lev,
                                          amrex::Real time,
                                          const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    /// initial projection
    ///
    void InitProj ();

    ///
    /// divu iter
    ///
    /// @param istep_divu_iter
    ///
    void DivuIter (int istep_divu_iter);

    /// initial iteration to initialize gradpi
    ///
    void InitIter ();

    // end MaestroInit.cpp functions
    ////////////

    ////////////
    // MaestroMacProj.cpp functions


    ///
    /// ``umac`` enters with face-centered, time-centered Utilde^* and should leave with Utilde
    /// ``macphi`` is the solution to the elliptic solve and
//   enters as either zero, or the solution to the predictor MAC projection
    /// ``macrhs`` enters as beta0*(S-Sbar)
    ///
    /// @param umac         MAC velocity
    /// @param macphi       MAC phi
    /// @param macrhs       MAC rhs
    /// @param beta0        a 1d cell-centered array
    /// @param is_predictor Predictor step (1) or corrector step (0)?
    ///
    void MacProj (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  amrex::Vector<amrex::MultiFab>& macphi,
                  const amrex::Vector<amrex::MultiFab>& macrhs,
                  const amrex::Vector<amrex::Real>& beta0,
                  const int& is_predictor);

    ///
    /// Multiply (or divide) face-data by ``beta0``
    ///
    /// @param edge             Edge-based (face) data
    /// @param beta0            cell-centered beta0
    /// @param beta0_edge       edge-centered beta0
    /// @param mult_or_div      Multiply (1) or divide (0)?
    ///
    void MultFacesByBeta0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge,
                           const amrex::Vector<amrex::Real>& beta0,
                           const amrex::Vector<amrex::Real>& beta0_edge,
                           const int& mult_or_div);

    ///
    /// Compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
    ///
    /// @param solverrhs    RHS
    /// @param macrhs       MAC rhs
    /// @param umac         MAC velocity
    ///
    void ComputeMACSolverRHS (amrex::Vector<amrex::MultiFab>& solverrhs,
                              const amrex::Vector<amrex::MultiFab>& macrhs,
                              const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);

    ///
    /// Average bcoefs at faces using inverse of rho
    ///
    /// @param facebcoef    ``bcoefs`` averaged at faces
    /// @param rhocc        cell-centered rho
    ///
    void AvgFaceBcoeffsInv(amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM > >& facebcoef,
			   const amrex::Vector<amrex::MultiFab>& rhocc);

    ///
    /// Set boundaries for ``LABecLaplacian`` to solve \f$-\nabla\cdot(B \nabla) \phi = \mathrm{RHS}\f$
    ///
    /// @param mlabec   pointer to ``amrex::MLABecLaplacian``
    ///
    void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

    // end MaestroMacProj.cpp functions
    ////////////

    ////////////
    // MaestroMakeEta.cpp functions

    ///
    /// Compute ``eta_rho`` at edge- and cell-centers
    ///
    /// @param etarho_edge      ``eta_rho`` at edges
    /// @param etarho_cell      ``eta_rho`` at cell-centers
    /// @param etarho_flux      ``eta_rho`` flux
    ///
    void MakeEtarho (amrex::Vector<amrex::Real>& etarho_edge,
		     amrex::Vector<amrex::Real>& etarho_cell,
		     const amrex::Vector<amrex::MultiFab>& etarho_flux);


    ///
    /// Compute ``eta_rho`` at edge- and cell-centers in spherical coordinates
    ///
    /// @param scal_old     Old scalars
    /// @param scal_new     New scalars
    /// @param umac         MAC velocity
    /// @param w0mac        MAC base state velocity
    /// @param etarho_edge  ``eta_rho`` at edges
    /// @param etarho_cell  ``eta_rho`` at cell-centers
    ///
    void MakeEtarhoSphr (const amrex::Vector<amrex::MultiFab>& scal_old,
			 const amrex::Vector<amrex::MultiFab>& scal_new,
			 const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			 const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
			 amrex::Vector<amrex::Real>& etarho_edge,
			 amrex::Vector<amrex::Real>& etarho_cell);
    ////////////

    ////////////
    // MaestroMakeS.cpp functions

    ///
    /// compute ``S`` at cell-centers
    ///
    /// @param S_cc                 cell-centered S
    /// @param delta_gamma1_term    delta gamma1 term
    /// @param delta_gamma1         delta gamma1
    /// @param scal                 Scalars
    /// @param u                    Velocity
    /// @param rho_omegadot         product of density and species creation rates
    /// @param rho_Hnuc             product of density and nuclear energy release
    /// @param rho_Hext             product of density and external energy source
    /// @param thermal              Thermal source term
    /// @param p0                   Base state pressure
    /// @param gamma1bar            gamma1bar
    /// @param delta_gamma1_termbar delta gamma1bar term
    /// @param psi                  psi
    ///
    void Make_S_cc (amrex::Vector<amrex::MultiFab>& S_cc,
                    amrex::Vector<amrex::MultiFab>& delta_gamma1_term,
                    amrex::Vector<amrex::MultiFab>& delta_gamma1,
                    const amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::MultiFab>& u,
                    const amrex::Vector<amrex::MultiFab>& rho_omegadot,
                    const amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                    const amrex::Vector<amrex::MultiFab>& rho_Hext,
                    const amrex::Vector<amrex::MultiFab>& thermal,
                    const amrex::Vector<amrex::Real>& p0,
                    const amrex::Vector<amrex::Real>& gamma1bar,
                    amrex::Vector<amrex::Real>& delta_gamma1_termbar,
                    const amrex::Vector<amrex::Real>& psi);

    ///
    /// Compute \f[
    /// \textrm{rhcc} = \beta_0(S_{cc}-\overline{S}) + \beta_0\delta\chi
    /// \f]
    ///
    /// @param rhcc                 rhs cell-centered source term
    /// @param S_cc                 cell-centered source term
    /// @param Sbar                 average source term
    /// @param beta0                beta0
    /// @param delta_gamma1_term    delta gamma1 term
    ///
    void MakeRHCCforNodalProj (amrex::Vector<amrex::MultiFab>& rhcc,
                               const amrex::Vector<amrex::MultiFab>& S_cc,
                               const amrex::Vector<amrex::Real>& Sbar,
                               const amrex::Vector<amrex::Real>& beta0,
                               const amrex::Vector<amrex::MultiFab>& delta_gamma1_term);

    ///
    /// Compute \f[
    /// \textrm{rhcc} = \beta_0(S_{cc}-\overline{S}) + \beta_0\delta\chi
    /// \f]
    ///
    /// @param rhcc                 rhs cell-centered source term
    /// @param rho0                 Base state density
    /// @param S_cc                 cell-centered source term
    /// @param Sbar                 average source term
    /// @param beta0                beta0
    /// @param delta_gamma1_term    delta gamma1 term
    /// @param gamma1bar            gamma1bar
    /// @param p0                   Base state pressure
    /// @param delta_p_term         delta p term
    /// @param delta_chi            delta chi
    /// @param is_predictor         Is this the predictor step (1) or corrector step (0)?
    ///
    void MakeRHCCforMacProj (amrex::Vector<amrex::MultiFab>& rhcc,
			                       const amrex::Vector<amrex::Real>& rho0,
                             const amrex::Vector<amrex::MultiFab>& S_cc,
                             const amrex::Vector<amrex::Real>& Sbar,
                             const amrex::Vector<amrex::Real>& beta0,
                             const amrex::Vector<amrex::MultiFab>& delta_gamma1_term,
                             const amrex::Vector<amrex::Real>& gamma1bar,
			                       const amrex::Vector<amrex::Real>& p0,
                             const amrex::Vector<amrex::MultiFab>& delta_p_term,
			                       amrex::Vector<amrex::MultiFab>& delta_chi,
                             int is_predictor);


    ///
    /// @param rhcc         cell-centered rhs
    /// @param rho0         base state density
    /// @param beta0
    /// @param gamma1bar
    /// @param p0           Base state pressure
    /// @param delta_p_term
    ///
    void CorrectRHCCforNodalProj(amrex::Vector<amrex::MultiFab>& rhcc,
				 const amrex::Vector<amrex::Real>& rho0,
				 const amrex::Vector<amrex::Real>& beta0,
				 const amrex::Vector<amrex::Real>& gamma1bar,
				 const amrex::Vector<amrex::Real>& p0,
				 const amrex::Vector<amrex::MultiFab>& delta_p_term);
    // end MaestroMakeS.cpp functions
    ////////////

    ////////////
    // MaestroNodalProj.cpp functions

    ///
    /// Perform a nodal projection
    ///
    /// Given a cell-centered velocity field ``Vproj`` \f$ V^\star \f$ (assembled in CreateUvecForProj()),
    /// \f$ V^\star \f$ can decomposed into \f[
    /// V^\star = \tilde{U} + \sigma \nabla \Phi,
    /// \f]
    /// where ``Utilde`` \f$ \tilde{U} \f$ satisfies the constraint \f[
    /// \nabla \cdot \left(\beta_0\tilde{U}\right) = \beta_0(S-\overline{S}).
    /// \f]
    /// Depending on ``proj_type`` we use different values of \f$ V^\star \f$,\f$ \sigma \f$, and \f$ \beta_0 \f$
    /// to solve for \f$ \Phi \f$ we use \f[
    /// \nabla \cdot \sigma \nabla \Phi = \nabla \cdot(\beta_0 V^\star) - \beta_0(S-\overline{S})
    /// \f]
    /// where \f$ \sigma = \beta_0 \f$ or \f$ \beta_0/\rho \f$ depending on ``proj_type``
    /// then solve for\f$ \tilde{U} \f$, \f$ \pi \f$, and \f$ \nabla \pi \f$ based on ``proj_type``.
    /// rhcc should enter as \f$ \beta_0(S-\overline{S}) \f$ so we need to multiply by -1 in this routine
    /// the projection (done below)
    ///
    /// @param proj_type        Projection type
    /// @param rhcc             rhcc
    /// @param istep_divu_iter
    ///
    void NodalProj (int proj_type,
                    amrex::Vector<amrex::MultiFab>& rhcc,
                    int istep_divu_iter=0);
    ///
    /// Fill in ``Vproj``
    ///
    /// - ``initial_projection_comp``: \f$ \tilde{U}^0 \f$, ``uold``
    /// - ``divu_iters_comp``:         \f$ \tilde{U}^0 \f$, ``uold``
    /// - ``pressure_iters_comp``:     \f$ (\tilde{U}^{n+1,*} - \tilde{U}^n)/dt\f$, ``(unew-uold)/dt``
    /// - ``regular_timestep_comp``:  \f$ (\tilde{U}^{n+1,*} + dt \nabla \pi/\rho^{n+1/2})\f$,  ``u_new + dt*gpi/rhohalf``
    ///
    /// @param proj_type        Projection type
    /// @param Vproj
    /// @param sig              contains rhohalf if ``proj_type == regular_timestep_comp``
    ///
    void CreateUvecForProj (int proj_type,
                            amrex::Vector<amrex::MultiFab>& Vproj,
                            const amrex::Vector<amrex::MultiFab>& sig);


    ///
    /// Sets the velocity at the boundaries
    ///
    /// 1) At non-inflow faces, the normal component of velocity will be completely zero'd
    ///
    /// 2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces
    ///                                will be zero'd outside of Neumann boundaries
    ///                                (slipWall, noSlipWall, Symmetry)
    ///                                BUT will retain non-zero values at periodic corners
    ///
    /// @param vel  Velocity
    ///
    void SetBoundaryVelocity(amrex::Vector<amrex::MultiFab>& vel);

    ///
    /// Given a nodal ``phi``, compute \f$ \nabla \Phi \f$ at cell centers
    ///
    /// @param phi      phi
    /// @param gphi     gradphi
    ///
    void ComputeGradPhi(amrex::Vector<amrex::MultiFab>& phi,
                        amrex::Vector<amrex::MultiFab>& gphi);

    ///
    /// Average nodal \f$ \pi \f$ to cell-centers and put in the Pi component of ``snew``
    ///
    /// @param beta0_cart       beta0 on cartesian grid
    ///
    void MakePiCC(const amrex::Vector<amrex::MultiFab>& beta0_cart);

    // end MaestroNodalProj.cpp functions
    ////////////

    ////////////
    // MaestroPlot.cpp functions

    ///
    /// get plotfile name
    ///
    /// @param lev
    ///
    std::string PlotFileName (int lev) const;


    ///
    /// Put together an array of multifabs for writing
    ///
    /// @param nPlot
    /// @param t_in
    /// @param dt_in
    /// @param rho0_cart
    /// @param rhoh0_cart
    /// @param p0    Base state pressure_cart
    /// @param gamma1bar_cart
    /// @param u_in
    /// @param s_in
    /// @param p0    Base state pressure_in
    /// @param gamma1bar_in
    /// @param S_cc_in
    ///

    amrex::Vector<const amrex::MultiFab*> PlotFileMF (const int nPlot,
                                                      const amrex::Real t_in,
                                                      const amrex::Real dt_in,
                                                      const amrex::Vector<amrex::MultiFab>& rho0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& rhoh0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& p0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& gamma1bar_cart,
                                                      const amrex::Vector<amrex::MultiFab>& u_in,
                                                      amrex::Vector<amrex::MultiFab>& s_in,
                                                      const amrex::Vector<amrex::Real>& p0_in,
                                                      const amrex::Vector<amrex::Real>& gamma1bar_in,
                                                      const amrex::Vector<amrex::MultiFab>& S_cc_in);

    ///
    /// Set plotfile variables' names
    ///
    /// @param nPlot
    ///
    amrex::Vector<std::string> PlotFileVarNames (int * nPlot) const;

    ///
    /// Write plotfile to disk
    ///
    /// @param step
    /// @param dt_in
    /// @param rho0_in
    /// @param rhoh0_in
    /// @param p0    Base state pressure_in
    /// @param u_in
    /// @param s_in
    ///
    void WritePlotFile (const int step,
                        const amrex::Real t_in,
                        const amrex::Real dt_in,
                        const amrex::Vector<amrex::Real>& rho0_in,
                        const amrex::Vector<amrex::Real>& rhoh0_in,
                        const amrex::Vector<amrex::Real>& p0_in,
                        const amrex::Vector<amrex::Real>& gamma1bar_in,
                        const amrex::Vector<amrex::MultiFab>& u_in,
                        amrex::Vector<amrex::MultiFab>& s_in,
                        const amrex::Vector<amrex::MultiFab>& S_cc_in);


    ///
    /// @param dir
    ///
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
                         const amrex::Vector<amrex::Real>& p0,
                         const amrex::Vector<amrex::MultiFab>& p0_cart,
                         const amrex::Vector<amrex::Real>& gamma1bar,
                         const amrex::Vector<amrex::MultiFab>& gamma1bar_cart,
                         amrex::Vector<amrex::MultiFab>& deltagamma);

    void MakeDivw0 (const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                    amrex::Vector<amrex::MultiFab>& divw0);

    void MakePiDivu (const amrex::Vector<amrex::MultiFab>& vel,
                   const amrex::Vector<amrex::MultiFab>& pi_cc,
                   amrex::Vector<amrex::MultiFab>& pidivu);

    // end MaestroPlot.cpp functions
    ////////////

    ////////////
    // MaestroReact.cpp functions

    ///
    /// compute heating term, ``rho_Hext``, then
    /// react the state over ``dt_react`` and update ``rho_omegadot``, ``rho_Hnuc``
    ///
    /// @param s_in
    /// @param s_out
    /// @param rho_Hext
    /// @param rho_omegadot
    /// @param rho_Hnuc
    /// @param p0    Base state pressure
    /// @param dt_react
    ///
    void React (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const amrex::Vector<amrex::Real>& p0,
                const amrex::Real dt_react);


    ///
    /// @param s_in
    /// @param s_out
    /// @param rho_Hext
    /// @param rho_omegadot
    /// @param rho_Hnuc
    /// @param p0    Base state pressure
    /// @param dt_react
    ///
    void Burner (const amrex::Vector<amrex::MultiFab>& s_in,
                amrex::Vector<amrex::MultiFab>& s_out,
                const amrex::Vector<amrex::MultiFab>& rho_Hext,
                amrex::Vector<amrex::MultiFab>& rho_omegadot,
                amrex::Vector<amrex::MultiFab>& rho_Hnuc,
                const amrex::Vector<amrex::Real>& p0,
                const amrex::Real dt_react);

    ///
    /// Compute heating term, ``rho_Hext``
    ///
    /// @param rho_Hext
    /// @param scal
    ///
    void MakeHeating (amrex::Vector<amrex::MultiFab>& rho_Hext,
                      const amrex::Vector<amrex::MultiFab>& scal);

    // end MaestroReact.cpp functions
    ////////////

    ////////////
    // regridding functions in MaestroRegrid.cpp

    /// check to see if we need to regrid, then regrid
    ///
    void Regrid ();

    /// set tagging array to include buffer zones for multilevel
    ///
    void TagArray ();

    ///
    /// tag all cells for refinement,
    /// overrides the pure virtual function in AmrCore
    ///
    /// @param lev
    /// @param tags
    /// @param time     Current time
    /// @param ng       Number of ghost cells
    ///
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ng) override;

    ///
    /// Within a call to ``AmrCore::regrid``, this function fills in data at a level
    /// that existed before, using pre-existing fine and interpolated coarse data
    /// overrides the pure virtual function in ``AmrCore``
    ///
    /// @param lev
    /// @param time    Current time
    /// @param ba
    /// @param dm
    ///
    virtual void RemakeLevel (int lev,
                              amrex::Real time,
                              const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    ///
    /// Within a call to ``AmrCore::regrid``, this function fills in data at a level
    /// that did NOT exist before, using interpolated coarse data
    /// overrides the pure virtual function in ``AmrCore``
    ///
    /// @param lev
    /// @param time    Current time
    /// @param ba
    /// @param dm
    ///
    virtual void MakeNewLevelFromCoarse (int lev,
                                         amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    ///
    /// within a call to ``AmrCore::regrid``, this function deletes all data
    /// at a level of refinement that is no longer needed
    /// overrides the pure virtual function in ``AmrCore``
    ///
    /// @param lev
    ///
    virtual void ClearLevel (int lev) override;

    // end regridding functions
    ////////////

    ////////////
    // MaestroRhoHT.cpp functions

    ///
    /// Calculate the temperature T given the density and enthalpy
    ///
    /// @param scal
    /// @param p0    Base state pressure
    ///
    void TfromRhoH (amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::Real>& p0);


    ///
    /// Calculate the temperature T given the density and pressure
    ///
    /// @param scal
    /// @param p0    Base state pressure
    /// @param updateRhoH
    ///
    void TfromRhoP (amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::Real>& p0,
                    int updateRhoH = 0);


    ///
    /// Calculate the pressure given the density and enthalpy
    ///
    /// @param state
    /// @param s_old
    /// @param peos
    ///
    void PfromRhoH (const amrex::Vector<amrex::MultiFab>& state,
                    const amrex::Vector<amrex::MultiFab>& s_old,
		    amrex::Vector<amrex::MultiFab>& peos);


    ///
    /// Calculate the Mach number given the density and enthalpy
    ///
    /// @param scal
    /// @param vel
    /// @param p0    Base state pressure
    /// @param mach
    ///
    void MachfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                       const amrex::Vector<amrex::MultiFab>& vel,
                       const amrex::Vector<amrex::Real>& p0,
                       amrex::Vector<amrex::MultiFab>& mach);


    ///
    /// @param scal
    /// @param vel
    /// @param p0    Base state pressure
    /// @param w0cart
    /// @param mach
    ///
    void MachfromRhoHSphr (const amrex::Vector<amrex::MultiFab>& scal,
                           const amrex::Vector<amrex::MultiFab>& vel,
                           const amrex::Vector<amrex::Real>& p0,
                           const amrex::Vector<amrex::MultiFab>& w0cart,
                           amrex::Vector<amrex::MultiFab>& mach);

    void CsfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                     const amrex::Vector<amrex::Real>& p0,
                     const amrex::Vector<amrex::MultiFab>& p0cart,
                     amrex::Vector<amrex::MultiFab>& cs);
    ////////////

    ////////////
    // MaestroSetup.cpp functions

    void ExternInit ();

    /// read in C++ parameters from inputs file
    ///
    void ReadParameters ();

    /// define variable mappings (Rho, RhoH, ..., NSCAL, etc.)
    ///
    void VariableSetup ();

    /// set up ``BCRec`` definitions for BC types
    ///
    void BCSetup ();

    // end MaestroSetup.cpp functions
    ////////////

    ////////////
    // MaestroSponge.cpp functions

    ///
    /// @param sponge
    ///
    void MakeSponge (amrex::Vector<amrex::MultiFab>& sponge);
    ////////////

    ////////////////////////
    // MaestroThermal.cpp functions

    ///
    /// Compute the quantity: ``thermal = del dot kappa grad T``
    ///
    /// - if ``temp_diffusion_formulation = 1``, then we compute this directly.
    /// - if ``temp_diffusion_formulation = 2``, then we compute the algebraically
    ///     equivalent form with ``grad h - grad X_k - grad p_0`` formulation
    ///
    /// @param thermal
    /// @param scal
    /// @param Tcoeff
    /// @param hcoeff
    /// @param Xkcoeff
    /// @param pcoeff
    /// @param p0    Base state pressure
    /// @param temp_formulation
    ///
    void MakeExplicitThermal (amrex::Vector<amrex::MultiFab>& thermal,
			      const amrex::Vector<amrex::MultiFab>& scal,
			      const amrex::Vector<amrex::MultiFab>& Tcoeff,
			      const amrex::Vector<amrex::MultiFab>& hcoeff,
			      const amrex::Vector<amrex::MultiFab>& Xkcoeff,
			      const amrex::Vector<amrex::MultiFab>& pcoeff,
			      const amrex::Vector<amrex::Real>& p0,
			      int temp_formulation);


    ///
    /// Use apply() to construct the form of the conduction term.
    /// apply() forms the generic quantity:
    /// \f[
    ///   (\alpha A - \beta \nabla \cdot B \nabla) \phi = \textrm{RHS}.
    /// \f]
    /// @param mlabec
    /// @param thermalout
    /// @param coeff
    /// @param phi
    /// @param bcs
    /// @param bccomp
    ///
    void ApplyThermal (amrex::MLABecLaplacian& mlabec,
		       amrex::Vector<amrex::MultiFab>& thermalout,
		 const amrex::Vector<amrex::MultiFab>& coeff,
		       amrex::Vector<amrex::MultiFab>& phi,
		 const amrex::Vector<amrex::BCRec>& bcs,
		 int bccomp);


    ///
    /// create the coefficients for grad{T}, grad{h}, grad{X_k}, and grad{p_0}
    /// for the thermal diffusion term in the enthalpy equation.
    ///
    /// @param scal
    /// @param Tcoeff
    /// @param hcoeff
    /// @param Xkcoeff
    /// @param pcoeff
    ///
    void MakeThermalCoeffs(const amrex::Vector<amrex::MultiFab>& scal,
			   amrex::Vector<amrex::MultiFab>& Tcoeff,
			   amrex::Vector<amrex::MultiFab>& hcoeff,
			   amrex::Vector<amrex::MultiFab>& Xkcoeff,
			   amrex::Vector<amrex::MultiFab>& pcoeff);


    ///
    /// ThermalConduct implements thermal diffusion in the enthalpy equation.
    /// This is an implicit solve, using the multigrid solver.  This updates
    /// the enthalpy only.
    ///
    /// @param s1
    /// @param s2
    /// @param hcoeff1
    /// @param Xkcoeff1
    /// @param pcoeff1
    /// @param hcoeff2
    /// @param Xkcoeff2
    /// @param pcoeff2
    ///
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

    ///
    /// @param rhohalf
    /// @param umac
    /// @param w0mac            MAC base state velocity
    /// @param w0_force
    /// @param w0_force_cart
    /// @param rho0_nph
    /// @param grav_cell_nph
    /// @param sponge
    ///
    void VelocityAdvance (const amrex::Vector<amrex::MultiFab>& rhohalf,
                          amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			  const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                          const amrex::Vector<amrex::Real>& w0_force,
			  const amrex::Vector<amrex::MultiFab>& w0_force_cart,
                          const amrex::Vector<amrex::Real>& rho0_nph,
                          const amrex::Vector<amrex::Real>& grav_cell_nph,
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

    ///
    /// @name Variables
    /// these are defined in VariableSetup()
    /// @{
    static int Rho, RhoH, FirstSpec, NumSpec, Temp, Pi, Nscal;
    static int ng_s; //< number of ghost cells for sold/new and uold/new
    /// @}
    ///

    ///
    /// @name Nodal projections
    /// these are defined in ``Maestro.cpp``,
    /// these represent the four types of nodal projections
    /// @{
    static int initial_projection_comp;
    static int divu_iters_comp;
    static int pressure_iters_comp;
    static int regular_timestep_comp;
    /// @}
    ///


    ///
    /// @name Edge state prediction
    /// these are defined in ``Maestro.cpp``,
    /// density and enthalpy edge state prediction
    /// @{
    static int predict_rhoprime_and_X, predict_rhoX, predict_rho_and_X;
    static int predict_rhoh, predict_rhohprime, predict_h;
    static int predict_T_then_rhohprime, predict_T_then_h;
    static int predict_hprime, predict_Tprime_then_h;
    /// @}
    ///


    ///
    /// @name Nodal flags
    /// these are defined in Maestro.cpp
    /// @{
    static amrex::IntVect nodal_flag;
    static amrex::IntVect nodal_flag_x;
    static amrex::IntVect nodal_flag_y;
    static amrex::IntVect nodal_flag_z;
    /// @}
    ///


    ///
    /// @name which step?
    /// @{
    int istep;
    int start_step;
    /// @}
    ///


    ///
    /// @name time
    /// keep track of old time, new time, and time step at each level
    /// @{
    amrex::Real t_new;
    amrex::Real t_old;
    amrex::Real dt;
    amrex::Real dtold;
    /// @}
    ///


    ///
    /// number of ghost cells needed for hyperbolic step
    ///
    int ng_adv;


    ///
    /// @name  array of multifabs to store the solution at each level of refinement
    /// after advancing a level we use "swap" to change pointers to old/new
    /// to avoid a copy.  These MultiFabs go through the regridding process
    /// @{
    amrex::Vector<amrex::MultiFab> sold;
    amrex::Vector<amrex::MultiFab> snew;
    amrex::Vector<amrex::MultiFab> uold;
    amrex::Vector<amrex::MultiFab> unew;
    amrex::Vector<amrex::MultiFab> S_cc_old;
    amrex::Vector<amrex::MultiFab> S_cc_new;
    /// @}
    ///


    ///
    /// these `MultiFab`s are persistant and need to be regridded as well
    ///
    amrex::Vector<amrex::MultiFab> gpi;
    amrex::Vector<amrex::MultiFab> dSdt;
    amrex::Vector<amrex::MultiFab> pi;   // nodal


    ///
    /// this only needs to persist leading into the initial pressure iters
    /// since we project ``(beta0^nph S^1 - beta0 S^0) / dt``
    /// during a regular time step we overwrite this
    ///
    amrex::Vector<amrex::MultiFab> rhcc_for_nodalproj;


    ///
    /// spherical only -
    /// we make this persistent in that we only have to rebuild and
    /// fill this after regridding
    ///
    amrex::Vector<amrex::MultiFab> normal;
    amrex::Vector<amrex::MultiFab> cell_cc_to_r;


    ///
    /// stores domain boundary conditions.
    /// These muse be vectors (rather than arrays) so we can ParmParse them
    ///
    amrex::Vector<int> phys_bc;


    ///
    /// @name boundary condition objects needed for FillPatch routines
    ///
    /// this is essentially an array (over components)
    /// of ``2*DIM`` integer arrays storing the physical boundary
    /// condition types at the lo/hi walls in each direction
    /// @{
    amrex::Vector<amrex::BCRec> bcs_s;  // for sold/snew
    amrex::Vector<amrex::BCRec> bcs_u;  // for uold/unew
    amrex::Vector<amrex::BCRec> bcs_f;  // "first-order extrap" used for S_cc, gpi, dSdt
    /// @}
    ///


    ///
    /// @name Multilevel 1d states
    /// vectors store the multilevel 1D states as one very long array.
    /// These are cell-centered
    /// @{
    amrex::Vector<amrex::Real> s0_init;
    amrex::Vector<amrex::Real> p0_init;
    amrex::Vector<amrex::Real> rho0_old;
    amrex::Vector<amrex::Real> rho0_new;
    amrex::Vector<amrex::Real> rhoh0_old;
    amrex::Vector<amrex::Real> rhoh0_new;
    amrex::Vector<amrex::Real> rhoX0_old;
    amrex::Vector<amrex::Real> rhoX0_new;
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
    /// @}
    ///

    ///
    /// \f$ p_0^{n-1} \f$ needed to compute \f$ d(p_0)/dt \f$ for nonuniform grid spacing
    ///
    amrex::Vector<amrex::Real> p0_nm1;
    amrex::Vector<amrex::Real> beta0_nm1;


    ///
    /// @name 1D states
    /// vectors store the multilevel 1D states as one very long array.
    /// These are edge-centered
    /// @{
    amrex::Vector<amrex::Real> w0;
    amrex::Vector<amrex::Real> etarho_ec;
    amrex::Vector<amrex::Real> r_edge_loc;
    /// @}
    ///


    ///
    /// @name information for the base state
    /// we also store these (except for the large arrays)
    /// in ``base_state.f90`` in a fortran module for convenience
    /// @{
    int max_radial_level;
    int finest_radial_level;
    int nr_fine;
    amrex::Real dr_fine;
    int nr_irreg;
    /// @}
    ///


    ///
    /// these provide information about the multilevel base state configuration
    ///
    amrex::Vector<int> numdisjointchunks;
    amrex::Vector<int> r_start_coord;
    amrex::Vector<int> r_end_coord;


    ///
    /// @name tagging criteria
    /// @{
    amrex::Vector<amrex::Real> temperr;
    amrex::Vector<amrex::Real> denserr;
    /// @}

    ///
    /// currently has 2 possible criterias
    ///
    amrex::Vector< amrex::Vector<amrex::Real> > tag_err;

    ///
    /// array of tagged boxes (planar)
    ///
    amrex::Vector<int> tag_array;

    ///
    /// diag file array buffers
    ///
    amrex::Vector<amrex::Real> diagfile_data;


    ///
    /// @name interface fluxes
    /// stores fluxes at coarse-fine interface for synchronization
    /// this will be sized "max_level+1"
    ///
    /// @note the flux register associated with ``flux_reg[lev]`` is associated
    /// with the ``lev/lev-1`` interface (and has grid spacing associated with ``lev-1``)
    /// therefore ``flux_reg[0]`` and ``flux_reg[max_level]`` are never actually
    /// used in the reflux operation
    /// @{
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_u;
    /// @}

};

#endif
