
#include <Maestro.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

using namespace amrex;


// Perform a nodal projection.
// Given a cell-centered velocity field Vproj (assembled in CreateUvecForProj),
// Vproj can decomposed into Vproj = Utilde + sig grad phi,
// where Utilde satisfies the constraint div(beta0*Utilde) = beta0*(S-Sbar)
// Depending on proj_type we use different values of Vproj, sig, and beta0
// to solve for phi we use div sigma grad phi = div(beta*Vproj) - beta0*(S-Sbar)
// where sigma = beta0 or beta0/rho depending on proj_type
// then solve for Utilde, pi, and grad(pi) based on proj_type.
// rhcc should enter as beta0*(S-Sbar) so we need to multiply by -1 before
// the projection (done below)
void
Maestro::NodalProj (int proj_type,
                    Vector<MultiFab>& rhcc,
                    int istep_divu_iter)
{
    AMREX_ASSERT(rhcc[0].nGrow() == 1);

    if ( !(proj_type == initial_projection_comp ||
           proj_type == divu_iters_comp ||
           proj_type == pressure_iters_comp ||
           proj_type == regular_timestep_comp) ) {
        amrex::Abort("Maestro::NodalProj - invalid proj_type");
    }

    // build a multifab sig with 1 ghost cell
    Vector<MultiFab> sig(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        sig[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    // fill in sig
    // initial_projection_comp: 1
    // divu_iters_comp:         1
    // pressure_iters_comp:     rho^1/2   -- (rhoold+rhonew)/2
    // regular_timestep_comp:   rho^n+1/2 -- (rhoold+rhonew)/2
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        for (int lev=0; lev<=finest_level; ++lev) {
            sig[lev].setVal(1.,1);
        }
    }
    else if (proj_type == pressure_iters_comp || proj_type == regular_timestep_comp) {
        for (int lev=0; lev<=finest_level; ++lev) {
            FillPatch(lev, 0.5*(t_old+t_new), sig[lev], sold, snew, 0, 0, 1, bcs_s);
        }
    }

    // build a multifab Vproj (the quantity that will be projected) with 1 ghost cell
    Vector<MultiFab> Vproj(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        Vproj[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
    }

    // fill in Vproj
    // initial_projection_comp: Utilde^0                        -- uold
    // divu_iters_comp:         Utilde^0                        -- uold
    // pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
    // regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
    // note sig is only used for regular_timestep_comp, and currently holds rhohalf
    CreateUvecForProj(proj_type,Vproj,sig);

    // build a multifab to store a Cartesian version of beta0 with 1 ghost cell
    Vector<MultiFab> beta0_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        beta0_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    // convert beta0 to multi-D MultiFab
    // initial_projection_comp: beta0_old
    // divu_iters_comp:         beta0_old
    // pressure_iters_comp:     (beta0_old+beta0_new)/2
    // regular_timestep_comp:   (beta0_old+beta0_new)/2
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        Put1dArrayOnCart(beta0_old,beta0_cart,bcs_f,0,0);
    }
    else {
        Vector<Real> beta0_nph( (max_radial_level+1)*nr_fine );
        for(int i=0; i<beta0_nph.size(); ++i) {
            beta0_nph[i] = 0.5*(beta0_old[i]+beta0_new[i]);
        }
        Put1dArrayOnCart(beta0_nph,beta0_cart,bcs_f,0,0);
    }

    // convert Vproj to beta0*Vproj
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            MultiFab::Multiply(Vproj[lev],beta0_cart[lev],0,dir,1,1);
        }
    }

    // invert sig then multiply by beta0 so sig now holds:
    // initial_projection_comp: beta0
    // divu_iters_comp:         beta0
    // pressure_iters_comp:     beta0/rho
    // regular_timestep_comp:   beta0/rho
    for (int lev=0; lev<=finest_level; ++lev) {
        sig[lev].invert(1.0,1);
        MultiFab::Multiply(sig[lev],beta0_cart[lev],0,0,1,1);
    }


/*
  set_outflow_bcs() modifies boundary conditions on phi at outflow.
  Typically these are homogeneous Dirichlet but in some cases you need
  to modify this.  First, if div(u) is "large enough" at outflow you need
  to modify the boundary conditions.  This often happens in combustion
  flame systems (PeleLM),  Second if you have hydrostatic effects
  due to gravity you need to modify the boundary conditions.  This happens
  in IAMR_type runs with strong gravity.  MAESTRO
  requires neither of these since we do not run problems with strong
  div(u) at outflow, and our pi (pressure) does not need to capture
  stratification since we have put this pressure in perturbational form
  by adding (rho-rho0)*g to the RHS of the velocity equation.
*/
/*
    if (OutFlowBC::HasOutFlowBC(phys_bc))
       set_outflow_bcs(INITIAL_VEL,phi,vel,
                       amrex::GetVecOfPtrs(rhcc),
                       amrex::GetVecOfPtrs(sig),
                       0,finest_level); 
*/

/*
  set_boundary_velocity() is a simplified version of the function in IAMR/Projection.cpp
  basically this sets velocity in ghost cells to zero except at inflow
  right now, the nodal solver expects ghost cells to be filled this way
  in slightly more detail
  1) At non-inflow faces, the normal component of velocity will be completely zero'd 
  2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces 
     will be zero'd outside of Neumann boundaries (slipWall, noSlipWall, Symmetry) 
     BUT will retain non-zero values at periodic corners
*/
    set_boundary_velocity(Vproj);

    // 
    for (int lev=0; lev<=finest_level; ++lev) {
        rhcc[lev].mult(-1.0,0,1,1);
    }

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else {
            if (lo_bc[idim] == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (lo_bc[idim] == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (hi_bc[idim] == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (hi_bc[idim] == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
        }
    }

    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMetricTerm(false);

    MLNodeLaplacian mlndlap(geom, grids, dmap, info);
    mlndlap.setGaussSeidel(true);
    mlndlap.setHarmonicAverage(false);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);
  
    // set sig in the MLNodeLaplacian object
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    // define a nodal rhs and set it to zero (everything we want is in rhcc)
    Vector<MultiFab> rhnd(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rhnd[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 0);
        rhnd[lev].setVal(0.);
    }

    // rhstotal is nodal and will contain the complete right-hand-side (rhnd + rhcc)
    Vector<MultiFab> rhstotal(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        rhstotal[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 0);
    }

    // phi is nodal and will contain the solution (need 1 ghost cell)
    Vector<MultiFab> phi(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        phi[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 1);
    }

    // Assemble the nodal RHS as the sum of the cell-centered RHS averaged to nodes 
    // plus div (beta0*Vproj) on nodes
    mlndlap.compRHS(amrex::GetVecOfPtrs(rhstotal),
                    amrex::GetVecOfPtrs(Vproj),
                    amrex::GetVecOfConstPtrs(rhnd),
                    amrex::GetVecOfPtrs(rhcc));

    MLMG mlmg(mlndlap);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(10);

    Real abs_tol = -1.; // disable absolute tolerance
    Real rel_tol;

    // logic for choosing multigrid tolerance
    // parameters are defined in Maestro.cpp
    // we change the tolerance depending on what part of the algorithm we
    // are in, and other factors including planar vs. spherical, the number
    // of AMR levels, etc.
    if (proj_type == initial_projection_comp) {
        rel_tol = (spherical==1) ? eps_init_proj_sph : eps_init_proj_cart;
    }
    else if (proj_type == divu_iters_comp) {
        if (spherical == 1) {
            if (istep_divu_iter == init_divu_iter) {
                rel_tol = eps_divu_sph;
            }
            else if (istep_divu_iter == init_divu_iter-1) {
                rel_tol = eps_divu_sph*divu_iter_factor;
            }
            else if (istep_divu_iter <= init_divu_iter-2) {
                rel_tol = eps_divu_sph*pow(divu_iter_factor,2);
            }
        }
        else {
            if (istep_divu_iter == init_divu_iter) {
                rel_tol = std::min(eps_divu_cart*pow(divu_level_factor,finest_radial_level), 
                                   eps_divu_cart*pow(divu_level_factor,2));
            }
            else if (istep_divu_iter == init_divu_iter-1) {
                rel_tol = std::min(eps_divu_cart*divu_iter_factor*pow(divu_level_factor,finest_radial_level),
                                   eps_divu_cart*divu_iter_factor*pow(divu_level_factor,2));
            }
            else if (istep_divu_iter <= init_divu_iter-2) {
                rel_tol = std::min(eps_divu_cart*pow(divu_iter_factor,2)*pow(divu_level_factor,finest_radial_level),
                                   eps_divu_cart*pow(divu_iter_factor,2)*pow(divu_level_factor,2));
            }
        }
    }
    else if (proj_type == pressure_iters_comp || proj_type == regular_timestep_comp) {
        rel_tol = std::min( eps_hg_max, eps_hg*pow(hg_level_factor,finest_level) );
    }

    // solve for phi
    Print() << "Calling nodal solver" << endl;
    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi),
                               amrex::GetVecOfConstPtrs(rhstotal),
                               rel_tol, abs_tol);
    Print() << "Done calling nodal solver" << endl;

    // compute a cell-centered grad(phi) from nodal phi
    // fixme need to write a routine for this


    // convert beta0*Vproj back to Vproj
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            MultiFab::Divide(Vproj[lev],beta0_cart[lev],0,dir,1,1);
        }
    }

    // divide sig by beta0 so it now holds
    // initial_projection_comp: 1
    // divu_iters_comp:         1
    // pressure_iters_comp:     1/rho
    // regular_timestep_comp:   1/rho
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Divide(sig[lev],beta0_cart[lev],0,0,1,1);
    }
  
    // reset sig in the MLNodeLaplacian object so the call to updateVelocity
    // works properly
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    // update velocity
    // initial_projection_comp: Utilde^0   = V - sig*grad(phi)
    // divu_iters_comp:         Utilde^0   = V - sig*grad(phi)
    // pressure_iters_comp:     Utilde^n+1 = Utilde^n + dt(V-sig*grad(phi))
    // regular_timestep_comp:   Utilde^n+1 = V - sig*grad(phi)
    if (proj_type == initial_projection_comp || 
        proj_type == divu_iters_comp) {
        // Vproj = Vproj - (1/sig)*grad(phi)
        mlndlap.updateVelocity(amrex::GetVecOfPtrs(Vproj), amrex::GetVecOfConstPtrs(phi));
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(uold[lev],Vproj[lev],0,0,AMREX_SPACEDIM,0);
        }
    }
    else if (proj_type == pressure_iters_comp) {
        // Vproj = Vproj - (1/sig)*grad(phi)
        mlndlap.updateVelocity(amrex::GetVecOfPtrs(Vproj), amrex::GetVecOfConstPtrs(phi));
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(unew[lev],Vproj[lev],0,0,AMREX_SPACEDIM,0);
            // multiply by dt
            unew[lev].mult(dt);
            MultiFab::Add(unew[lev],uold[lev],0,0,AMREX_SPACEDIM,0);
        }

    }
    else if (proj_type == regular_timestep_comp) {
        // Vproj = Vproj - (1/sig)*grad(phi)
        mlndlap.updateVelocity(amrex::GetVecOfPtrs(Vproj), amrex::GetVecOfConstPtrs(phi));
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(unew[lev],Vproj[lev],0,0,AMREX_SPACEDIM,0);
        }
    }

    // update pi and grad(pi)
    // initial_projection_comp: pi = 0         grad(pi) = 0
    // divu_iters_comp:         pi = 0         grad(pi) = 0
    // pressure_iters_comp:     pi = pi + phi  grad(pi) = grad(pi) + grad(phi)
    // regular_timestep_comp:   pi = phi/dt    grad(pi) = grad(phi)/dt
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {

    for (int lev=0; lev<=finest_level; ++lev) {
        sold[lev].setVal(0.); // fixme only want to set Pi component to zero
        gpi[lev].setVal(0.);
    }
    }
    else if (proj_type == pressure_iters_comp) {
        // fixme
    }
    else if (proj_type == regular_timestep_comp) {
        // fixme
    }




}

// fill in Vproj
// initial_projection_comp: Utilde^0                        -- uold
// divu_iters_comp:         Utilde^0                        -- uold
// pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
// regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
    // note sig is only used for regular_timestep_comp, and currently holds rhohalf
void
Maestro::CreateUvecForProj (int proj_type,
                            Vector<MultiFab>& Vproj,
                            const Vector<MultiFab>& sig) {

    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        // Vproj = uold
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(Vproj[lev],uold[lev],0,0,AMREX_SPACEDIM,0);
        }
    }
    else if (proj_type == pressure_iters_comp) {
        // Vproj = (unew-uold)/dt
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(Vproj[lev],unew[lev],0,0,AMREX_SPACEDIM,0);
            MultiFab::Saxpy(Vproj[lev],-1.0,uold[lev],0,0,AMREX_SPACEDIM,0);
            Vproj[lev].mult(1/dt);
        }

    }
    else if (proj_type == regular_timestep_comp) {
        // Vproj = unew + dt*gpi/rhohalf
        for (int lev=0; lev<=finest_level; ++lev) {
            // Vproj = unew
            MultiFab::Copy(Vproj[lev],unew[lev],0,0,AMREX_SPACEDIM,0);
            // convert gpi to gpi/rhohalf
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                MultiFab::Divide(gpi[lev],sig[lev],0,dir,1,0);
            }
            // Vproj = unew + dt*gpi/rhohalf
            MultiFab::Saxpy(Vproj[lev],dt,gpi[lev],0,0,AMREX_SPACEDIM,0);
            // revert gpi/rhohalf back to gpi
            for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                MultiFab::Multiply(gpi[lev],sig[lev],0,dir,1,0);
            }
        }
    }
    else {
        amrex::Abort("MaestroNodalProj: invalid proj_type");
    }

    Real time;
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        time = t_old;
    }
    else  {
        time = t_new;
    }

    // fill ghost cells
    for (int lev=0; lev<=finest_level; ++lev) {
        FillPatch(lev, time, Vproj[lev], Vproj, Vproj, 0, 0, AMREX_SPACEDIM, bcs_u);
    }
}


void Maestro::set_boundary_velocity(Vector<MultiFab>& vel)
{
    // 1) At non-inflow faces, the normal component of velocity will be completely zero'd 
    // 2) If a face is an inflow face, then the normal velocity at corners just outside inflow faces 
    //                                will be zero'd outside of Neumann boundaries 
    //                                (slipWall, noSlipWall, Symmetry) 
    //                                BUT will retain non-zero values at periodic corners

    for (int lev=0; lev <= finest_level; lev++) {
        const BoxArray& grids_lev = grids[lev];
        const Box& domainBox = geom[lev].Domain();

        for (int idir=0; idir<BL_SPACEDIM; idir++) {
            if (lo_bc[idir] != Inflow && hi_bc[idir] != Inflow) {
                vel[lev].setBndry(0.0, idir, 1);
            }
            else {
                for (MFIter mfi(vel[lev]); mfi.isValid(); ++mfi) {
                    int i = mfi.index();
                    FArrayBox& v_fab = (vel[lev])[mfi];

                    const Box& reg = grids_lev[i];
                    const Box& bxg1 = amrex::grow(reg, 1);

                    BoxList bxlist(reg);

                    if (lo_bc[idir] == Inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
                        Box bx; // bx is the region we *protect* from zero'ing

                        bx = amrex::adjCellLo(reg, idir);

                        for (int odir = 0; odir < BL_SPACEDIM; odir++) {
                            if (odir != idir) {
                                if (geom[lev].isPeriodic(odir)) bx.grow(odir,1);
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
                                if (geom[lev].isPeriodic(odir)) bx.grow(odir,1);
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
                            v_fab.setVal(0.0, ovlp, idir, 1);
                        }
                    }

                } // end loop over grids
            } // end if/else logic for inflow
        } // end loop over direction
    } // end loop over levels
}

/*
void Maestro::set_outflow_bcs (int        which_call,
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

void Maestro::set_outflow_bcs_at_level (int          which_call,
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
*/
