
#ifndef _Projection_H_
#define _Projection_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <RegType.H>
#include <SyncRegister.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFabUtil.H>

class Projection 
{
protected:
  static const int maxlev = 10;
  static const amrex::Real BogusValue;
  static const int  LEVEL_PROJ    = 1001;
  static const int  INITIAL_VEL   = 1002;
  static const int  INITIAL_PRESS = 1003;
  static const int  INITIAL_SYNC  = 1004;
  static const int  SYNC_PROJ     = 1005;
  static const int  FILTER_P      = 1006;

public:

    Projection (amrex::Amr*   _parent,
                amrex::BCRec* _phys_bc, 
                int    _do_sync_proj,
                int    /*_finest_level*/, 
                int    _radius_grow );

    ~Projection ();

    void install_anelastic_coefficient (int                        level,
                                        amrex::Real                **anel_coeff);

    //
    // The initial velocity projection in post_init.
    // This function ensures that the velocities are nondivergent.
    //
    void initialVelocityProject (
                                 amrex::Real divu_time,
                                 int  have_divu);

    //
    // Multiply/divide a amrex::MultiFab by anel_coeff.
    //
    void AnelCoeffMult (int       level,
                        amrex::MultiFab& mf,
                        int       comp);
    void AnelCoeffDiv (int       level,
                       amrex::MultiFab& mf,
                       int       comp);

    void put_divu_in_cc_rhs (amrex::MultiFab&       rhs,
                             int             level,
                             amrex::Real            time);

    //
    // Helper function for NavierStokes::getGradP().
    //
    static void getGradP (amrex::FArrayBox& p_fab,
                          amrex::FArrayBox& gp,
                          const amrex::Box& gpbox_to_fill,
                          const amrex::Real* dx);
protected:
    //
    // Read parameters from the input file.
    //
    static void Initialize ();

    static void Finalize ();

    //
    // Convert U from an Accleration like quantity to a velocity
    // Unew = Uold + alpha*Unew
    //
    void UnConvertUnew (amrex::MultiFab&       Uold,
                        amrex::Real            alpha,
                        amrex::MultiFab&       Unew, 
                        const amrex::BoxArray& grids );
    //
    // Convert U from an Accleration like quantity to a velocity
    // Unew = Uold + alpha*Unew
    //
    void UnConvertUnew (amrex::FArrayBox& Uold,
                        amrex::Real       alpha,
                        amrex::FArrayBox& Unew,
                        const amrex::Box& grd);
    //
    // Convert U to an Accleration like quantity
    // Unew = (Unew - Uold)/alpha
    //
    void ConvertUnew (amrex::MultiFab&       Unew,
                      amrex::MultiFab&       Uold,
                      amrex::Real            alpha,
                      const amrex::BoxArray& grids );
    //
    // Convert U to an Accleration like quantity.
    // Unew = (Unew - Uold)/alpha
    //
    void ConvertUnew (amrex::FArrayBox& Unew,
                      amrex::FArrayBox& Uold,
                      amrex::Real       alpha,
                      const amrex::Box& grd);
    //
    // Update a quantity U using the formula
    // Unew = Unew + alpha*Uold
    //
    void UpdateArg1 (amrex::MultiFab&       Unew,
                     amrex::Real            alpha,
                     amrex::MultiFab&       Uold,
                     int             nvar,
                     const amrex::BoxArray& grids,
                     int ngrow);
    //
    // Update a quantity U using the formula
    // Unew = Unew + alpha*Uold
    //
    void UpdateArg1 (amrex::FArrayBox& Unew,
                     amrex::Real       alpha,
                     amrex::FArrayBox& Uold,
                     int        nvar,
                     const amrex::Box& grd,
                     int        ngrow);
    //
    // Add phi to P.
    //
    void AddPhi (amrex::MultiFab&       p,
                 amrex::MultiFab&       phi);
    //
    // Convert phi into p^n+1/2.
    //
    void incrPress (int  level,
                    amrex::Real dt);
    //
    // This function scales variables at the start of a projection.
    //
    void scaleVar (int             which_call,
                   amrex::MultiFab*       sig,
                   int             sig_nghosts,
                   amrex::MultiFab*       vel,
                   int             level);
    //
    // This function rescales variables at the end of a projection.
    //
    void rescaleVar (int             which_call,
                     amrex::MultiFab*       sig,
                     int             sig_nghosts,
                     amrex::MultiFab*       vel,
                     int             level);

    void set_outflow_bcs (int        which_call,
                          const amrex::Vector<amrex::MultiFab*>& phi,
                          const amrex::Vector<amrex::MultiFab*>& Vel_in,
                          const amrex::Vector<amrex::MultiFab*>& Divu_in,
                          const amrex::Vector<amrex::MultiFab*>& Sig_in,
                          int        c_lev,
                          int        f_lev,
                          int        have_divu);

    void set_outflow_bcs_at_level (int          which_call,
                                   int          lev,
                                   int          c_lev,
                                   amrex::Box*         state_strip,
                                   amrex::Orientation* outFacesAtThisLevel,
                                   int          numOutFlowFaces,
                                   const amrex::Vector<amrex::MultiFab*>& phi,
                                   amrex::MultiFab*    Vel_in,
                                   amrex::MultiFab*    Divu_in,
                                   amrex::MultiFab*    Sig_in,
                                   int          have_divu,
                                   amrex::Real         gravity);

protected:

    void putDown(const amrex::Vector<amrex::MultiFab*>& phi, amrex::FArrayBox* phi_fine_strip,
		 int c_lev, int f_lev, const amrex::Orientation* outFaces,
		 int numOutFlowFaces, int ncStripWidth);
    //
    // Pointers to amrlevel and amr.
    //
    amrex::Amr*             parent;
    amrex::Vector<amrex::AmrLevel*> LevelData;
    //
    // The array of Radii, and number of radii ghost cells.
    //
    int                    radius_grow;
    amrex::Vector< amrex::Vector< amrex::Vector<amrex::Real> > * > radius;
    //
    // Anelastic coefficient array.
    //
    amrex::Vector<amrex::Real**> anel_coeff;

    //
    // Boundary objects.
    //
    amrex::BCRec*                         phys_bc;

    int do_sync_proj;

    static int  verbose;
    static int  P_code;
    static amrex::Real proj_tol;
    static amrex::Real sync_tol;
    static amrex::Real proj_abs_tol;
    static int  proj_2;
    static int  filter_u;
    static int  rho_wgt_vel_proj;
    static int  do_outflow_bcs;
    static int  make_sync_solvable;
    static amrex::Real divu_minus_s_factor;
    static int  add_vort_proj;

    void doMLMGNodalProjection(const amrex::Vector<amrex::MultiFab*>& vel, 
                               const amrex::Vector<amrex::MultiFab*>& phi,
                               const amrex::Vector<amrex::MultiFab*>& sig,
                               const amrex::Vector<amrex::MultiFab*>& rhs_cc,
                               const amrex::Vector<amrex::MultiFab*>& rhnd, 
                               amrex::Real rel_tol, amrex::Real abs_tol,
                               amrex::MultiFab* sync_resid_crse=nullptr,
                               amrex::MultiFab* sync_resid_fine=nullptr,
                               bool doing_initial_velproj=false);
    
    // set velocity in ghost cells to zero except for inflow
    void set_boundary_velocity(int c_lev, int nlevel, const amrex::Vector<amrex::MultiFab*>& vel, 
                               bool doing_initial_velproj, bool inflowCorner);

};
#endif
