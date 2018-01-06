

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <NavierStokesBase.H>
#include <NS_BC.H>
#include <AMReX_BLProfiler.H>
#include <Projection.H>
#include <PROJECTION_F.H>
#include <NAVIERSTOKES_F.H>
#include <ProjOutFlowBC.H>

#include <AMReX_MGT_Solver.H>
#include <mg_cpp_f.h>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_BOX_LIMITS(box,boxlo,boxhi)   \
const int* boxlo = (box).loVect();           \
const int* boxhi = (box).hiVect();

const Real Projection::BogusValue = 1.e200;

int  Projection::P_code              = 0;
int  Projection::proj_2              = 1;
int  Projection::verbose             = 0;
Real Projection::proj_tol            = 1.0e-12;
Real Projection::proj_abs_tol        = 1.e-16;
int  Projection::do_outflow_bcs      = 1;

namespace
{
    bool initialized = false;
}


void
Projection::Initialize ()
{
    if (initialized) return;

    ParmParse pp("proj");

    pp.query("v",                   verbose);
    pp.query("Pcode",               P_code);
    pp.query("proj_2",              proj_2);
    pp.query("proj_tol",            proj_tol);
    pp.query("proj_abs_tol",        proj_abs_tol);
    pp.query("do_outflow_bcs",      do_outflow_bcs);

    if (!proj_2) 
	amrex::Error("With new gravity and outflow stuff, must use proj_2");

    amrex::ExecOnFinalize(Projection::Finalize);

    initialized = true;
}

void
Projection::Finalize ()
{
    initialized = false;
}

Projection::Projection (Amr*   _parent,
                        BCRec* _phys_bc, 
                        int    /*_finest_level*/, 
                        int    _radius_grow )
   :
    parent(_parent),
    LevelData(_parent->finestLevel()+1),
    radius_grow(_radius_grow), 
    radius(_parent->finestLevel()+1),
    phys_bc(_phys_bc), 
{

    BL_ASSERT ( parent->finestLevel()+1 <= maxlev );

    Initialize();

    if (verbose) amrex::Print() << "Creating projector\n";

}

Projection::~Projection ()
{
  if (verbose) amrex::Print() << "Deleting projector\n";
}


//
// The initial velocity projection in post_init.
// this function ensures that the velocities are nondivergent
//


// there are four types of nodal projections we are interested in
// all of them are multilevel composite
// initial projection
// divu iter
// init iter
// regular timestep
void
Projection::initialVelocityProject (Real cur_divu_time)
{
    int c_lev = 0;
    int f_lev = finest_level;

    Vector<MultiFab*> vel(maxlev, nullptr);
    Vector<MultiFab*> phi(maxlev, nullptr);
    Vector<std::unique_ptr<MultiFab> > sig(maxlev);

    // fixme
    // set pressure to zero

    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev] = &(LevelData[lev]->get_new_data(State_Type));
        phi[lev] = &(LevelData[lev]->get_old_data(Press_Type));

        const int       nghost = 1;
        const BoxArray& grids  = LevelData[lev]->boxArray();
        const DistributionMapping& dmap = LevelData[lev]->DistributionMap();
        sig[lev].reset(new MultiFab(grids,dmap,1,nghost));


        // fixme
        // for initial projection and divu iters, we want "1" weighting
        // for pressure iterations and regular time step we want "rho" weighting
        if (rho_wgt_vel_proj) 
        {
            LevelData[lev]->get_new_data(State_Type).setBndry(BogusValue,Density,1);

	    AmrLevel& amr_level = parent->getLevel(lev);
	    
	    MultiFab& S_new = amr_level.get_new_data(State_Type);
	    
            Real curr_time = amr_level.get_state_data(State_Type).curTime();
	    
            for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
	      amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,
					      Density,Density,1);
            }

            MultiFab::Copy(*sig[lev],
                           LevelData[lev]->get_new_data(State_Type),
                           Density,
                           0,
                           1,
                           nghost);
        }
        else 
        {
            sig[lev]->setVal(1,nghost);
        }
    }

    Vector<std::unique_ptr<MultiFab> > rhcc(maxlev);
    const int nghost = 1; 

    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
        vel[lev]->setBndry(BogusValue,Xvel,BL_SPACEDIM);
        //
        // Set the physical boundary values.
        //
        AmrLevel& amr_level = parent->getLevel(lev);

        MultiFab& S_new = amr_level.get_new_data(State_Type);

        Real curr_time = amr_level.get_state_data(State_Type).curTime();

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            amr_level.setPhysBoundaryValues(S_new[mfi],State_Type,curr_time,Xvel,Xvel,BL_SPACEDIM);
        }

        // divu
        int Divu_Type, Divu;
        if (!LevelData[lev]->isStateVariable("divu", Divu_Type, Divu)) 
            amrex::Error("Projection::initialVelocityProject(): Divu not found");
        //
        // Make sure ghost cells are properly filled.
        //
        MultiFab& divu_new = amr_level.get_new_data(Divu_Type);
        divu_new.FillBoundary();
        
        Real curr_divu_time_int lev = amr_level.get_state_data(Divu_Type).curTime();
        
        for (MFIter mfi(divu_new); mfi.isValid(); ++mfi)
        {
            amr_level.setPhysBoundaryValues(divu_new[mfi],Divu_Type,curr_divu_time_lev,0,0,1);
        }
	
        const BoxArray& grids     = amr_level.boxArray();
        const DistributionMapping& dmap = amr_level.DistributionMap();
        rhcc[lev].reset(new MultiFab(grids,dmap,1,nghost));
        put_divu_in_cc_rhs(*rhcc[lev],lev,cur_divu_time);
    }

    if (OutFlowBC::HasOutFlowBC(phys_bc) && do_outflow_bcs)
       set_outflow_bcs(INITIAL_VEL,phi,vel,
                       amrex::GetVecOfPtrs(rhcc),
                       amrex::GetVecOfPtrs(sig),
                       c_lev,f_lev); 

    //
    // Project
    //
    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
        rhcc[lev]->mult(-1.0,0,1,nghost);
    }

    doMLMGNodalProjection(vel,
                          phi,
                          amrex::GetVecOfPtrs(sig),
                          amrex::GetVecOfPtrs(rhcc),
                          {},
                          proj_tol, proj_abs_tol)

    // fixme
    // update pressure
    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
        LevelData[lev]->get_old_data(Press_Type).setVal(0);
        LevelData[lev]->get_new_data(Press_Type).setVal(0);
    }
}



//
// Put S in the rhs of the projector--cell based version.
//

void
Projection::put_divu_in_cc_rhs (MultiFab&       rhs,
                                int             level,
                                Real            time)
{
    rhs.setVal(0);

    NavierStokesBase* ns = dynamic_cast<NavierStokesBase*>(&parent->getLevel(level));

    BL_ASSERT(!(ns == 0));

    std::unique_ptr<MultiFab> divu (ns->getDivCond(1,time));

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        rhs[mfi].copy((*divu)[mfi]);
    }
}

void 
Projection::putDown (const Vector<MultiFab*>& phi,
                     FArrayBox*         phi_fine_strip,
                     int                c_lev,
                     int                f_lev,
                     const Orientation* outFaces,
                     int                numOutFlowFaces,
                     int                ncStripWidth)
{
    BL_PROFILE("Projection::putDown()");
    //
    // Put down to coarser levels.
    //
    const int nCompPhi = 1; // phi_fine_strip.nComp();
    const int nGrow    = 0; // phi_fine_strip.nGrow();
    IntVect ratio      = IntVect::TheUnitVector();

    for (int lev = f_lev-1; lev >= c_lev; lev--)
    {
        ratio *= parent->refRatio(lev);
        const Box& domainC = parent->Geom(lev).Domain();

        for (int iface = 0; iface < numOutFlowFaces; iface++) 
        {
            Box phiC_strip = 
                amrex::surroundingNodes(amrex::bdryNode(domainC, outFaces[iface], ncStripWidth));
            phiC_strip.grow(nGrow);
            BoxArray ba(phiC_strip);
            DistributionMapping dm{ba};
            MultiFab phi_crse_strip(ba, dm, nCompPhi, 0);
            phi_crse_strip.setVal(0);

            for (MFIter mfi(phi_crse_strip); mfi.isValid(); ++mfi)
            {
                Box ovlp = amrex::coarsen(phi_fine_strip[iface].box(),ratio) & mfi.validbox();

                if (ovlp.ok())
                {
                    FArrayBox& cfab = phi_crse_strip[mfi];
                    FORT_PUTDOWN (BL_TO_FORTRAN(cfab),
                                  BL_TO_FORTRAN(phi_fine_strip[iface]),
                                  ovlp.loVect(), ovlp.hiVect(), ratio.getVect());
                }
            }

            phi[lev]->copy(phi_crse_strip);
        }
    }
}

//
// Given a nodal pressure P compute the pressure gradient at the
// contained cell centers.

void
Projection::getGradP (FArrayBox& p_fab,
                      FArrayBox& gp,
                      const Box& gpbox_to_fill,
                      const Real* dx)
{
    BL_PROFILE("Projection::getGradP()");
    //
    // Test to see if p_fab contains gpbox_to_fill
    //
    BL_ASSERT(amrex::enclosedCells(p_fab.box()).contains(gpbox_to_fill));

    const int*  plo    = p_fab.loVect();
    const int*  phi    = p_fab.hiVect();
    const int*  glo    = gp.box().loVect();
    const int*  ghi    = gp.box().hiVect();
    const int*   lo    = gpbox_to_fill.loVect();
    const int*   hi    = gpbox_to_fill.hiVect();
    const Real* p_dat  = p_fab.dataPtr();
    const Real* gp_dat = gp.dataPtr();

#if (BL_SPACEDIM == 2)
    int is_full = 0;
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx,
               &is_full);
#elif (BL_SPACEDIM == 3)
    FORT_GRADP(p_dat,ARLIM(plo),ARLIM(phi),gp_dat,ARLIM(glo),ARLIM(ghi),lo,hi,dx);
#endif
}

void
Projection::set_outflow_bcs (int        which_call,
                             const Vector<MultiFab*>& phi,
                             const Vector<MultiFab*>& Vel_in,
                             const Vector<MultiFab*>& Divu_in,
                             const Vector<MultiFab*>& Sig_in,
                             int        c_lev,
                             int        f_lev)
{
    BL_ASSERT((which_call == INITIAL_VEL  ) || 
              (which_call == INITIAL_PRESS) || 
              (which_call == LEVEL_PROJ   ) );

    if (which_call != LEVEL_PROJ)
      BL_ASSERT(c_lev == 0);

    if (verbose)
      amrex::Print() << "...setting outflow bcs for the nodal projection ... " << '\n';

    bool        hasOutFlow;
    Orientation outFaces[2*BL_SPACEDIM];
    Orientation outFacesAtThisLevel[maxlev][2*BL_SPACEDIM];

    int fine_level[2*BL_SPACEDIM];

    int numOutFlowFacesAtAllLevels;
    int numOutFlowFaces[maxlev];
    OutFlowBC::GetOutFlowFaces(hasOutFlow,outFaces,phys_bc,numOutFlowFacesAtAllLevels);

    //
    // Get 2-wide cc box, state_strip, along interior of top. 
    // Get 1-wide nc box, phi_strip  , along top.
    //
    const int ccStripWidth = 2;

//    const int nCompPhi    = 1;
//    const int srcCompVel  = Xvel;
//    const int srcCompDivu = 0;
//    const int   nCompVel  = BL_SPACEDIM;
//    const int   nCompDivu = 1;

    //
    // Determine the finest level such that the entire outflow face is covered
    // by boxes at this level (skip if doesnt touch, and bomb if only partially
    // covered).
    //
    Box state_strip[maxlev][2*BL_SPACEDIM];

    int icount[maxlev];
    for (int i=0; i < maxlev; i++) icount[i] = 0;

    //
    // This loop is only to define the number of outflow faces at each level.
    //
    Box temp_state_strip;
    for (int iface = 0; iface < numOutFlowFacesAtAllLevels; iface++) 
    {
      const int outDir    = outFaces[iface].coordDir();

      fine_level[iface] = -1;
      for (int lev = f_lev; lev >= c_lev; lev--)
      {
        Box domain = parent->Geom(lev).Domain();

        if (outFaces[iface].faceDir() == Orientation::high)
        {
            temp_state_strip = amrex::adjCellHi(domain,outDir,ccStripWidth);
            temp_state_strip.shift(outDir,-ccStripWidth);
        }
        else
        {
            temp_state_strip = amrex::adjCellLo(domain,outDir,ccStripWidth);
            temp_state_strip.shift(outDir,ccStripWidth);
        }
        // Grow the box by one tangentially in order to get velocity bc's.
        for (int dir = 0; dir < BL_SPACEDIM; dir++) 
          if (dir != outDir) temp_state_strip.grow(dir,1);

        const BoxArray& Lgrids               = parent->getLevel(lev).boxArray();
        const Box&      valid_state_strip    = temp_state_strip & domain;
        const BoxArray  uncovered_outflow_ba = amrex::complementIn(valid_state_strip,Lgrids);

        BL_ASSERT( !(uncovered_outflow_ba.size() &&
                     amrex::intersect(Lgrids,valid_state_strip).size()) );

        if ( !(uncovered_outflow_ba.size()) && fine_level[iface] == -1) {
            int ii = icount[lev];
            outFacesAtThisLevel[lev][ii] = outFaces[iface];
            state_strip[lev][ii] = temp_state_strip;
            fine_level[iface] = lev;
            icount[lev]++;
        }
      }
    }

    for (int lev = f_lev; lev >= c_lev; lev--) {
      numOutFlowFaces[lev] = icount[lev];
    }

    NavierStokesBase* ns0 = dynamic_cast<NavierStokesBase*>(LevelData[c_lev]);
    BL_ASSERT(!(ns0 == 0));
   
    int Divu_Type, Divu;
    Real gravity = 0;

    if (which_call == INITIAL_VEL)
    {
      gravity = 0;
      if (!LevelData[c_lev]->isStateVariable("divu", Divu_Type, Divu))
        amrex::Error("Projection::set_outflow_bcs: No divu.");
    }

    if (which_call == INITIAL_PRESS || which_call == LEVEL_PROJ)
    {
      gravity = ns0->getGravity();
      if (!LevelData[c_lev]->isStateVariable("divu", Divu_Type, Divu) &&
          (gravity == 0) )
        amrex::Error("Projection::set_outflow_bcs: No divu or gravity.");
    }

    for (int lev = c_lev; lev <= f_lev; lev++) 
    {
      if (numOutFlowFaces[lev] > 0) 
        set_outflow_bcs_at_level (which_call,lev,c_lev,
                                  state_strip[lev],
                                  outFacesAtThisLevel[lev],
                                  numOutFlowFaces[lev],
                                  phi,
                                  Vel_in[lev],
                                  Divu_in[lev],
                                  Sig_in[lev],
                                  gravity);
                                  
    }

}

void
Projection::set_outflow_bcs_at_level (int          which_call,
                                      int          lev,
                                      int          c_lev,
                                      Box*         state_strip,
                                      Orientation* outFacesAtThisLevel,
                                      int          numOutFlowFaces,
                                      const Vector<MultiFab*>&  phi, 
                                      MultiFab*    Vel_in,
                                      MultiFab*    Divu_in,
                                      MultiFab*    Sig_in,
                                      Real         gravity)
{
    BL_ASSERT(dynamic_cast<NavierStokesBase*>(LevelData[lev]) != nullptr);

    Box domain = parent->Geom(lev).Domain();

    const int ncStripWidth = 1;

    FArrayBox  rho[2*BL_SPACEDIM];
    FArrayBox dsdt[2*BL_SPACEDIM];
    FArrayBox dudt[1][2*BL_SPACEDIM];
    FArrayBox phi_fine_strip[2*BL_SPACEDIM];

    const int ngrow = 1;

    for (int iface = 0; iface < numOutFlowFaces; iface++)
    {
        dsdt[iface].resize(state_strip[iface],1);
        dudt[0][iface].resize(state_strip[iface],BL_SPACEDIM);

        rho[iface].resize(state_strip[iface],1);

        (*Sig_in).copyTo(rho[iface],0,0,1,ngrow);

        Box phi_strip = 
            amrex::surroundingNodes(amrex::bdryNode(domain,
                                                      outFacesAtThisLevel[iface],
                                                      ncStripWidth));
        phi_fine_strip[iface].resize(phi_strip,1);
        phi_fine_strip[iface].setVal(0.);
    }

    ProjOutFlowBC projBC;
    if (which_call == INITIAL_PRESS) 
    {

        const int*      lo_bc = phys_bc->lo();
        const int*      hi_bc = phys_bc->hi();
        projBC.computeRhoG(rho,phi_fine_strip,
                           parent->Geom(lev),
                           outFacesAtThisLevel,numOutFlowFaces,gravity,
                           lo_bc,hi_bc);
    }
    else
    {
        Vel_in->FillBoundary();

	for (int iface = 0; iface < numOutFlowFaces; iface++) 
	    (*Vel_in).copyTo(dudt[0][iface],0,0,BL_SPACEDIM,1);

        // since we have divu
        for (int iface = 0; iface < numOutFlowFaces; iface++) 
            (*Divu_in).copyTo(dsdt[iface],0,0,1,1);

        const int*      lo_bc = phys_bc->lo();
        const int*      hi_bc = phys_bc->hi();
        projBC.computeBC(dudt, dsdt, rho, phi_fine_strip,
                         parent->Geom(lev),
                         outFacesAtThisLevel,
                         numOutFlowFaces, lo_bc, hi_bc, gravity);
    }

    for (int i = 0; i < 2*BL_SPACEDIM; i++)
    {
        rho[i].clear();
        dsdt[i].clear();
        dudt[0][i].clear();
    }

    for ( int iface = 0; iface < numOutFlowFaces; iface++)
    {
        BoxArray phi_fine_strip_ba(phi_fine_strip[iface].box());
        DistributionMapping dm {phi_fine_strip_ba};
        MultiFab phi_fine_strip_mf(phi_fine_strip_ba,dm,1,0);

        for (MFIter mfi(phi_fine_strip_mf); mfi.isValid(); ++mfi) {
            phi_fine_strip_mf[mfi].copy(phi_fine_strip[iface]);
        }

        phi[lev]->copy(phi_fine_strip_mf);
    }

    if (lev > c_lev) 
    {
        putDown(phi, phi_fine_strip, c_lev, lev, outFacesAtThisLevel,
                numOutFlowFaces, ncStripWidth);
    }
}


//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void Projection::doMLMGNodalProjection (const Vector<MultiFab*>& vel, 
                                        const Vector<MultiFab*>& phi,
                                        const Vector<MultiFab*>& sig,
                                        const Vector<MultiFab*>& rhcc,
                                        const Vector<MultiFab*>& rhnd,
                                        Real rel_tol, Real abs_tol)
{
    BL_PROFILE("Projection:::doMLMGNodalProjection()");

    int c_lev = 0;
    int f_lev = finest_level;
    int nlevel = finest_level+1;

    BL_ASSERT(vel[c_lev]->nGrow() == 1);
    BL_ASSERT(vel[f_lev]->nGrow() == 1);
    BL_ASSERT(phi[c_lev]->nGrow() == 1);
    BL_ASSERT(phi[f_lev]->nGrow() == 1);
    BL_ASSERT(sig[c_lev]->nGrow() == 1);
    BL_ASSERT(sig[f_lev]->nGrow() == 1);
    
    if (rhcc[c_lev]) {
        BL_ASSERT(rhcc[c_lev]->nGrow() == 1);
        BL_ASSERT(rhcc[f_lev]->nGrow() == 1);
    }

    // Set velocity in ghost cells to zero except for inflow
    set_boundary_velocity(vel);

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geometry::isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            if (phys_bc->lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc->hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc->hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    Vector<Geometry> mg_geom(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_geom[lev] = parent->Geom(lev+c_lev);
    }  

    Vector<BoxArray> mg_grids(nlevel);
    for (int lev = 0; lev < nlevel; lev++) {
        mg_grids[lev] = parent->boxArray(lev+c_lev);
    }

    Vector<DistributionMapping> mg_dmap(nlevel);
    for (int lev=0; lev < nlevel; lev++ ) {
        mg_dmap[lev] = LevelData[lev+c_lev]->get_new_data(State_Type).DistributionMap();
    }

    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMetricTerm(false);

    MLNodeLaplacian mlndlap(mg_geom, mg_grids, mg_dmap, info);
    mlndlap.setGaussSeidel(true);
    mlndlap.setHarmonicAverage(false);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);
  
    for (int ilev = 0; ilev < nlevel; ++ilev) {
        mlndlap.setSigma(ilev, *sig[c_lev+ilev]);
    }

    Vector<MultiFab> rhs(nlevel);
    for (int ilev = 0; ilev < nlevel; ++ilev)
    {
        const auto& ba = amrex::convert(mg_grids[ilev], IntVect::TheNodeVector());
        rhs[ilev].define(ba, mg_dmap[ilev], 1, 0);
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), vel, rhnd, rhcc);

    MLMG mlmg(mlndlap);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(P_code);

    Real mlmg_err = mlmg.solve(phi, amrex::GetVecOfConstPtrs(rhs), rel_tol, abs_tol);

    mlndlap.updateVelocity(vel, amrex::GetVecOfConstPtrs(phi));
}


// Set velocity in ghost cells to zero except for inflow
void Projection::set_boundary_velocity(const Vector<MultiFab*>& vel)
{
    const int* lo_bc = phys_bc->lo();
    const int* hi_bc = phys_bc->hi();

    int c_lev = 0;
    int nlevel = finest_level+1;

    // 1) At non-inflow faces, the normal component of velocity will be completely zero'd 
    // 2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces 
    //                                will be zero'd outside of Neumann boundaries 
    //                                (slipWall, noSlipWall, Symmetry) 
    //                                BUT will retain non-zero values at periodic corners

    for (int lev=c_lev; lev < c_lev+nlevel; lev++) {
        const BoxArray& grids = parent->boxArray(lev);
        const Box& domainBox = parent->Geom(lev).Domain();

        const Geometry& geom = parent->Geom(lev);

        for (int idir=0; idir<BL_SPACEDIM; idir++) {
            if (lo_bc[idir] != Inflow && hi_bc[idir] != Inflow) {
                vel[lev]->setBndry(0.0, Xvel+idir, 1);
            }
            else {
                for (MFIter mfi(*vel[lev]); mfi.isValid(); ++mfi) {
                    int i = mfi.index();
                    FArrayBox& v_fab = (*vel[lev])[mfi];

                    const Box& reg = grids[i];
                    const Box& bxg1 = amrex::grow(reg, 1);

                    BoxList bxlist(reg);

                    if (lo_bc[idir] == Inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
                        Box bx; // bx is the region we *protect* from zero'ing

                        bx = amrex::adjCellLo(reg, idir);

                        for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                            if (odir != idir) {
                                if (geom.isPeriodic(odir)) bx.grow(odir,1);
                                if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                                if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
                            }
                        }
                        bxlist.push_back(bx);
                    }

                    if (hi_bc[idir] == Inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
                        Box bx; // bx is the region we *protect* from zero'ing

                        bx = amrex::adjCellHi(reg, idir);

                        for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                            if (odir != idir) {
                                if (geom.isPeriodic(odir)) bx.grow(odir,1);
                                if (reg.bigEnd  (odir) != domainBox.bigEnd  (odir) ) bx.growHi(odir,1);
                                if (reg.smallEnd(odir) != domainBox.smallEnd(odir) ) bx.growLo(odir,1);
                            }
                        }

                        bxlist.push_back(bx);
                    }

                    BoxList bxlist2 = amrex::complementIn(bxg1, bxlist); 
 
                    for (BoxList::iterator it=bxlist2.begin(); it != bxlist2.end(); ++it) {
                        Box ovlp = *it & v_fab.box();
                        if (ovlp.ok()) {
                            v_fab.setVal(0.0, ovlp, Xvel+idir, 1);
                        }
                    }

                } // end loop over grids
            } // end if/else logic for inflow
        } // end loop over direction
    } // end loop over levels
}
