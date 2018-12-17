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
/// @param w0mac
/// @param w0_force
/// @param w0_force_cart
///
    void AdvancePremac (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
			const amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac,
                        const amrex::Vector<amrex::Real>& w0_force,
			const amrex::Vector<amrex::MultiFab>& w0_force_cart);


///
/// @param utilde
/// @param ufull
/// @param utrans
/// @param w0mac
///
    void MakeUtrans (const amrex::Vector<amrex::MultiFab>& utilde,
                     const amrex::Vector<amrex::MultiFab>& ufull,
                     amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac);


///
/// @param utilde
/// @param ufull
/// @param utrans
/// @param umac
/// @param w0mac
/// @param force
///
    void VelPred (const amrex::Vector<amrex::MultiFab>& utilde,
                  const amrex::Vector<amrex::MultiFab>& ufull,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& utrans,
                  amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                  const amrex::Vector<amrex::MultiFab>& force);


///
/// @param state
/// @param sedge
/// @param umac
/// @param force
/// @param is_vel
/// @param bcs
/// @param nbccomp
/// @param start_scomp
/// @param start_bccomp
/// @param num_comp
/// @param is_conservative
///
    void MakeEdgeScal (const amrex::Vector<amrex::MultiFab>& state,
		       amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sedge,
                       const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                       const amrex::Vector<amrex::MultiFab>& force,
                       int is_vel, const amrex::Vector<amrex::BCRec>& bcs, int nbccomp,
                       int start_scomp, int start_bccomp, int num_comp, int is_conservative);


///
/// @param state
/// @param sflux
/// @param etarhoflux
/// @param sedge
/// @param umac
/// @param w0mac
/// @param r0_old
/// @param r0_edge_old
/// @param r0mac_old
/// @param r0_new
/// @param r0_edge_new
/// @param r0mac_new
/// @param r0_predicted_edge
/// @param start_comp
/// @param num_comp
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
/// @param state
/// @param sflux
/// @param sedge
/// @param umac
/// @param w0mac
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

    void UpdateScal (const amrex::Vector<amrex::MultiFab>& stateold,
		     amrex::Vector<amrex::MultiFab>& statenew,
		     const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& sflux,
		     const amrex::Vector<amrex::MultiFab>& force,
		     int start_scomp, int num_comp,
		     const amrex::Real* p0 = NULL,
		     const amrex::Vector<amrex::MultiFab>& p0_cart = {});


///
/// @param umac
/// @param uedge
/// @param force
/// @param sponge
/// @param w0mac
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
/// @param is
///
    void GotoNextLine (std::istream& is);
    ////////////

    ////////////
    // MaestroConvert.cpp functions

///
/// @param scal
/// @param s0
/// @param comp
/// @param bccomp
/// @param bcs_in
/// @param flag
///
    void PutInPertForm (amrex::Vector<amrex::MultiFab>& scal,
                        const amrex::Vector<amrex::Real>& s0,
                        int comp, int bccomp,
                        const amrex::Vector<amrex::BCRec>& bcs_in,
                        bool flag);


///
/// @param scal
/// @param flag
///
    void ConvertRhoXToX (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);


///
/// @param scal
/// @param flag
///
    void ConvertRhoHToH (amrex::Vector<amrex::MultiFab>& scal,
                         bool flag);
    ////////////

    ////////////////////////
    // MaestroDebug.cpp functions

///
/// print out the contents of a Vector of MultiFabs
///
/// @param CC
///
    void PrintMF   (amrex::Vector<amrex::MultiFab>& CC);

///
/// print out the contents of a Vector of edge-based MultiFabs
///
/// @param EDGE
/// @param dir
///
    void PrintEdge (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& EDGE,
                    int dir);

///
/// utility to write out a multilevel multifab to a plotfile
///
/// @param mf
/// @param name
///
    void WriteMF   (amrex::Vector<amrex::MultiFab>& mf,
                    std::string name);
    ////////////////////////

    ////////////////////////
    // MaestroDensityAdvance.cpp functions

///
/// @param which_step
/// @param scalold
/// @param scalnew
/// @param sedge
/// @param sflux
/// @param scal_force
/// @param etarhoflux
/// @param umac
/// @param w0mac
/// @param rho0_predicted_edge
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
/// @param step
/// @param dt_in
/// @param rho0_in
/// @param p0_in
/// @param u_in
/// @param s_in
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
/// @param which_step
/// @param scalold
/// @param scalnew
/// @param sedge
/// @param sflux
/// @param scal_force
/// @param umac
/// @param w0mac
/// @param thermal
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
/// @param time
/// @param mf
/// @param mf_old
/// @param mf_new
/// @param srccomp
/// @param destcomp
/// @param ncomp
/// @param startbccomp
/// @param bcs_in
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
/// @param time
/// @param mf
/// @param mf_old
/// @param mf_new
/// @param srccomp      the source component
/// @param destcomp     the destination component AND the bc component
/// @param ncomp
/// @param startbccomp
/// @param bcs_in
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
/// @param time
/// @param mf
/// @param mf_old
/// @param mf_new
/// @param srccomp      the source component
/// @param destcomp     the destination component AND the bc component
/// @param ncomp
/// @param bcs
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
/// @param time
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
/// @param mf
/// @param comp
/// @param ncomp
///
    void AverageDown (amrex::Vector<amrex::MultiFab>& mf,
                      int comp,
                      int ncomp);

///
/// set covered faces to be the average of overlying fine faces
///
/// @param edge
///
    void AverageDownFaces (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge);

///
/// fill in ONE ghost cell for all components of a face-centered (MAC) velocity
/// field behind physical boundaries.  Does not modify the velocities on the boundary
///
/// @param umac
/// @param level
///
    void FillUmacGhost (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                        int level=-1);

///
/// fill in all ghost cells for an edge-based MAC velocity field
///
/// @param uedge
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


///
/// @param uedge
/// @param w0mac
/// @param mult
///
    void Addw0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& uedge,
		const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& w0mac,
                const amrex::Real& mult);


///
/// @param w0mac
///
    void MakeW0mac (amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& w0mac);


///
/// @param s0
/// @param s0mac
///
    void MakeS0mac (const amrex::Vector<amrex::Real>& s0,
		    amrex::Vector<std::array< amrex::MultiFab,AMREX_SPACEDIM > >& s0mac);

    void MakeNormal ();


///
/// @param s_cc
/// @param face
/// @param harmonic_avg
///
    void PutDataOnFaces(const amrex::Vector<amrex::MultiFab>& s_cc,
			amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM >>& face,
			int harmonic_avg);

    void MakeCCtoRadii ();
    // end MaestroFill3dData.cpp functions
    ////////////

    ////////////////////////
    // MaestroForce.cpp functions

///
/// @param vel_force
/// @param uedge
/// @param rho
/// @param rho0
/// @param grav
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
                       int do_add_utilde_force);


///
/// @param scal_force
/// @param state
/// @param umac
/// @param s0
/// @param s0_edge
/// @param s0_cart
/// @param comp
/// @param bcs
/// @param fullform
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
/// @param scal_force
/// @param is_prediction
/// @param thermal
/// @param umac
/// @param add_thermal
/// @param &which_step
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
/// @param scal
/// @param gamma1bar
/// @param p0
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
/// @param lev
/// @param time
/// @param ba
/// @param dm
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

    /// init iter to initialize gradpi
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
/// @param umac
/// @param macphi
/// @param macrhs
/// @param beta0        a 1d cell-centered array
/// @param is_predictor
///
    void MacProj (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac,
                  amrex::Vector<amrex::MultiFab>& macphi,
                  const amrex::Vector<amrex::MultiFab>& macrhs,
                  const amrex::Vector<amrex::Real>& beta0,
                  const int& is_predictor);

///
/// multiply (or divide) face-data by ``beta0``
///
/// @param edge
/// @param beta0
/// @param beta0_edge
/// @param mult_or_div
///
    void MultFacesByBeta0 (amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& edge,
                           const amrex::Vector<amrex::Real>& beta0,
                           const amrex::Vector<amrex::Real>& beta0_edge,
                           const int& mult_or_div);

///
/// compute the RHS for the solve, RHS = macrhs - div(beta0*umac)
///
/// @param solverrhs
/// @param macrhs
/// @param umac
///
    void ComputeMACSolverRHS (amrex::Vector<amrex::MultiFab>& solverrhs,
                              const amrex::Vector<amrex::MultiFab>& macrhs,
                              const amrex::Vector<std::array< amrex::MultiFab, AMREX_SPACEDIM > >& umac);

///
/// Average bcoefs at faces using inverse of rho
///
/// @param facebcoef
/// @param rhocc
///
    void AvgFaceBcoeffsInv(amrex::Vector<std::array<amrex::MultiFab, AMREX_SPACEDIM > >& facebcoef,
			   const amrex::Vector<amrex::MultiFab>& rhocc);

///
/// set boundaries for ``LABecLaplacian`` to solve \f$-\nabla\cdot(B \nabla) \phi = RHS\f$
///
/// @param mlabec
///
    void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

    // end MaestroMacProj.cpp functions
    ////////////

    ////////////
    // MaestroMakeEta.cpp functions

///
/// compute eta_rho at edge- and cell-centers
///
/// @param etarho_edge
/// @param etarho_cell
/// @param etarho_flux
///
    void MakeEtarho (amrex::Vector<amrex::Real>& etarho_edge,
		     amrex::Vector<amrex::Real>& etarho_cell,
		     const amrex::Vector<amrex::MultiFab>& etarho_flux);


///
/// @param scal_old
/// @param scal_new
/// @param umac
/// @param w0mac
/// @param etarho_edge
/// @param etarho_cell
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
/// @param S_cc
/// @param delta_gamma1_term
/// @param delta_gamma1
/// @param scal
/// @param u
/// @param rho_omegadot
/// @param rho_Hnuc
/// @param rho_Hext
/// @param thermal
/// @param p0
/// @param gamma1bar
/// @param delta_gamma1_termbar
/// @param psi
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
/// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
///
/// @param rhcc
/// @param S_cc
/// @param Sbar
/// @param beta0
/// @param delta_gamma1_term
///
    void MakeRHCCforNodalProj (amrex::Vector<amrex::MultiFab>& rhcc,
                               const amrex::Vector<amrex::MultiFab>& S_cc,
                               const amrex::Vector<amrex::Real>& Sbar,
                               const amrex::Vector<amrex::Real>& beta0,
                               const amrex::Vector<amrex::MultiFab>& delta_gamma1_term);

///
/// compute rhcc = beta0*(S_cc-Sbar) + beta0*delta_chi
///
/// @param rhcc
/// @param rho0
/// @param S_cc
/// @param Sbar
/// @param beta0
/// @param delta_gamma1_term
/// @param gamma1bar
/// @param p0
/// @param delta_p_term
/// @param delta_chi
/// @param is_predictor
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
/// @param rhcc
/// @param rho0
/// @param beta0
/// @param gamma1bar
/// @param p0
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
/// perform a nodal projection
///
/// Given a cell-centered velocity field Vproj (assembled in CreateUvecForProj),
/// Vproj can decomposed into Vproj = Utilde + sig grad phi,
/// where Utilde satisfies the constraint div(beta0*Utilde) = beta0*(S-Sbar).
/// Depending on proj_type we use different values of Vproj, sig, and beta0
/// to solve for phi we use div sigma grad phi = div(beta*Vproj) - beta0*(S-Sbar)
/// where sigma = beta0 or beta0/rho depending on proj_type
/// then solve for Utilde, pi, and grad(pi) based on proj_type.
/// rhcc should enter as beta0*(S-Sbar) so we need to multiply by -1 in this routine
/// the projection (done below)
///
/// @param proj_type
/// @param rhcc
/// @param istep_divu_iter
///
    void NodalProj (int proj_type,
                    amrex::Vector<amrex::MultiFab>& rhcc,
                    int istep_divu_iter=0);
///
/// fill in Vproj
///
/// - ``initial_projection_comp``: Utilde^0                        -- uold
/// - ``divu_iters_comp``:         Utilde^0                        -- uold
/// - ``pressure_iters_comp``:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
/// - ``regular_timestep_comp``:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
///
/// @param proj_type
/// @param Vproj
/// @param sig      contains rhohalf if proj_type == regular_timestep_comp
///
    void CreateUvecForProj (int proj_type,
                            amrex::Vector<amrex::MultiFab>& Vproj,
                            const amrex::Vector<amrex::MultiFab>& sig);


///
/// 1) At non-inflow faces, the normal component of velocity will be completely zero'd
///
/// 2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces
///                                will be zero'd outside of Neumann boundaries
///                                (slipWall, noSlipWall, Symmetry)
///                                BUT will retain non-zero values at periodic corners
///
/// @param vel
///
    void SetBoundaryVelocity(amrex::Vector<amrex::MultiFab>& vel);

///
/// given a nodal ``phi``, compute grad(phi) at cell centers
///
/// @param phi
/// @param gphi
///
    void ComputeGradPhi(amrex::Vector<amrex::MultiFab>& phi,
                        amrex::Vector<amrex::MultiFab>& gphi);

///
/// average nodal pi to cell-centers and put in the Pi component of ``snew``
///
/// @param beta0_cart
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
/// put together an array of multifabs for writing
///
/// @param rho0_cart
/// @param rhoh0_cart
/// @param p0_cart
/// @param u_in
/// @param s_in
/// @param p0_in
///
    amrex::Vector<const amrex::MultiFab*> PlotFileMF (const amrex::Vector<amrex::MultiFab>& rho0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& rhoh0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& p0_cart,
                                                      const amrex::Vector<amrex::MultiFab>& u_in,
                                                      amrex::Vector<amrex::MultiFab>& s_in,
                                                      const amrex::Vector<amrex::Real>& p0_in);

    /// set plotfile variables names
    ///
    amrex::Vector<std::string> PlotFileVarNames () const;

///
/// write plotfile to disk
///
/// @param step
/// @param dt_in
/// @param rho0_in
/// @param rhoh0_in
/// @param p0_in
/// @param u_in
/// @param s_in
///
    void WritePlotFile (const int step,
                        const amrex::Real dt_in,
                        const amrex::Vector<amrex::Real>& rho0_in,
                        const amrex::Vector<amrex::Real>& rhoh0_in,
                        const amrex::Vector<amrex::Real>& p0_in,
                        const amrex::Vector<amrex::MultiFab>& u_in,
                        amrex::Vector<amrex::MultiFab>& s_in);


///
/// @param dir
///
    void WriteJobInfo (const std::string& dir) const;

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
/// @param p0
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
/// @param p0
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
/// compute heating term, ``rho_Hext``
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

///
/// tag all cells for refinement,
/// overrides the pure virtual function in AmrCore
///
/// @param lev
/// @param tags
/// @param time
/// @param ng
///
    virtual void ErrorEst (int lev,
                           amrex::TagBoxArray& tags,
                           amrex::Real time,
                           int ng) override;

///
/// Within a call to ``AmrCore::regrid``, this function fills in data at a level
/// that existed before, using pre-existing fine and interpolated coarse data
/// overrides the pure virtual function in AmrCore
///
/// @param lev
/// @param time
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
/// overrides the pure virtual function in AmrCore
///
/// @param lev
/// @param time
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
/// overrides the pure virtual function in AmrCore
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
/// @param p0
///
    void TfromRhoH (amrex::Vector<amrex::MultiFab>& scal,
                    const amrex::Vector<amrex::Real>& p0);


///
/// Calculate the temperature T given the density and pressure
///
/// @param scal
/// @param p0
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
/// @param p0
/// @param mach
///
    void MachfromRhoH (const amrex::Vector<amrex::MultiFab>& scal,
                       const amrex::Vector<amrex::MultiFab>& vel,
                       const amrex::Vector<amrex::Real>& p0,
                       amrex::Vector<amrex::MultiFab>& mach);


///
/// @param scal
/// @param vel
/// @param p0
/// @param w0cart
/// @param mach
///
    void MachfromRhoHSphr (const amrex::Vector<amrex::MultiFab>& scal,
                           const amrex::Vector<amrex::MultiFab>& vel,
                           const amrex::Vector<amrex::Real>& p0,
                           const amrex::Vector<amrex::MultiFab>& w0cart,
                           amrex::Vector<amrex::MultiFab>& mach);
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

    /// set up BCRec definitions for BC types
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
/// Compute the quantity: thermal = del dot kappa grad T
///
/// - if temp_diffusion_formulation = 1, then we compute this directly.
/// - if temp_diffusion_formulation = 2, then we compute the algebraically
///     equivalent form with grad h - grad X_k - grad p_0 formulation
///
/// @param thermal
/// @param scal
/// @param Tcoeff
/// @param hcoeff
/// @param Xkcoeff
/// @param pcoeff
/// @param p0
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
///   (\alpha A - \beta \nabla \cdot B \nabla) \phi = RHS
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
/// @param w0mac
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
    static int ng_s; // number of ghost cells for sold/new and uold/new
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
/// these MultiFabs are persistant and need to be regridded as well
///
    amrex::Vector<amrex::MultiFab> gpi;
    amrex::Vector<amrex::MultiFab> dSdt;
    amrex::Vector<amrex::MultiFab> pi;   // nodal


///
/// this only needs to persist leading into the initial pressure iters
/// since we project (beta0^nph S^1 - beta0 S^0) / dt
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
/// of 2*DIM integer arrays storing the physical boundary
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
/// p0^{n-1} needed to compute d(p0)/dt for nonuniform grid spacing
///
    amrex::Vector<amrex::Real> p0_nm1;


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
/// diag file array buffers
///
    amrex::Vector<amrex::Real> diagfile_data;


///
/// @name interface fluxes
/// stores fluxes at coarse-fine interface for synchronization
/// this will be sized "max_level+1"
///
/// @note the flux register associated with flux_reg[lev] is associated
/// with the lev/lev-1 interface (and has grid spacing associated with lev-1)
/// therefore flux_reg[0] and flux_reg[max_level] are never actually
/// used in the reflux operation
/// @{
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_s;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg_u;
/// @}

};

#endif
