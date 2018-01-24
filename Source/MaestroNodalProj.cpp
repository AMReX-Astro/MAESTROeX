
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
// rhcc should enter as beta0*(S-Sbar) so we need to multiply by -1 in this routine
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
  set_outflow_bcs() was used in IAMR/Projection.cpp
  This routine modifies boundary conditions on phi at outflow.
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
    set_outflow_bcs()
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
            if (phys_bc[idim] == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc[idim] == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

            if (phys_bc[AMREX_SPACEDIM+idim] == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (phys_bc[AMREX_SPACEDIM+idim] == Inflow) {
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
        phi[lev].setVal(0.);
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
            pi[lev].setVal(0.);
            gpi[lev].setVal(0.);
        }
    }
    else if (proj_type == pressure_iters_comp) {
        // fixme
    }
    else if (proj_type == regular_timestep_comp) {
        // fixme
    }

    // average pi from nodes to cell-centers and store in the Pi component of s
    if (proj_type == initial_projection_comp || proj_type == divu_iters_comp) {
        for (int lev=0; lev<=finest_level; ++lev) {
            sold[lev].setVal(0.,Pi,1,0);
        }
    } else if (proj_type == pressure_iters_comp || proj_type == regular_timestep_comp) {
        // fixme need a new routine
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
            if (phys_bc[idir] != Inflow && phys_bc[AMREX_SPACEDIM+idir] != Inflow) {
                vel[lev].setBndry(0.0, idir, 1);
            }
            else {
                for (MFIter mfi(vel[lev]); mfi.isValid(); ++mfi) {
                    int i = mfi.index();
                    FArrayBox& v_fab = (vel[lev])[mfi];

                    const Box& reg = grids_lev[i];
                    const Box& bxg1 = amrex::grow(reg, 1);

                    BoxList bxlist(reg);

                    if (phys_bc[idir] == Inflow && reg.smallEnd(idir) == domainBox.smallEnd(idir)) {
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

                    if (phys_bc[AMREX_SPACEDIM+idir] == Inflow && reg.bigEnd(idir) == domainBox.bigEnd(idir)) {
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
