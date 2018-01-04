
#include <Maestro.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

using namespace amrex;


// Perform a nodal projection.
// Given a cell-centered velocity field Vproj (assembled in CreateUvecForProj),
// Vproj can decomposed into Vproj = Utilde + (1/sig) grad phi,
// where Utilde satisfies the constraint div(beta0*Utilde) = beta0*(S-Sbar)
// Depending on proj_type we use different values of Vproj, sig, and beta0
// to solve for phi we use div (beta0/sig) grad phi = div(beta*Vproj) - beta0*(S-Sbar)
// then solve for Utilde, pi, and grad(pi) based on proj_type.
void
Maestro::NodalProj (int proj_type,
                    Vector<MultiFab>& rhcc)
{
    // build a multifab "sig" div (1/sig) grad phi = RHS with 1 ghost cell
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
            sig[lev].setVal(1.);
        }
    }
    else {
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

    // convert Vproj to beta0*Vproj in valid region
    for (int lev=0; lev<=finest_level; ++lev) {
        for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
            MultiFab::Multiply(Vproj[lev],beta0_cart[lev],0,dir,1,1);
        }
    }

    // divide sig by beta0
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Divide(sig[lev],beta0_cart[lev],0,0,1,1);
    }

    // create a unew with a filled ghost cell
    Vector<MultiFab> unew_ghost(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        unew_ghost[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        FillPatch(lev, t_new, unew_ghost[lev], unew, unew, 0, 0, AMREX_SPACEDIM, bcs_u);
    }

    // Assemble the nodal RHS as the sum of the cell-centered RHS averaged to nodes 
    // plus div (beta0*Vproj) on nodes
    // fix up the boundary condition ghost cells for Vproj
    // solve for phi
    

    // make grad phi



    // convert sig/beta0 back to sig
    for (int lev=0; lev<=finest_level; ++lev) {
        MultiFab::Multiply(sig[lev],beta0_cart[lev],0,0,1,1);
    }

    // update velocity
    // initial_projection_comp: Utilde^0   = V - (1/sig)*grad(phi)
    // divu_iters_comp:         Utilde^0   = V - (1/sig)*grad(phi)
    // pressure_iters_comp:     Utilde^n+1 = Utilde^n + dt(V-(1/sig)*grad(phi))
    // regular_timestep_comp:   Utilde^n+1 = V - (1/sig)*grad(phi)




    // update pi and grad(pi)
    // initial_projection_comp: pi = 0         grad(pi) = 0
    // divu_iters_comp:         pi = 0         grad(pi) = 0
    // pressure_iters_comp:     pi = pi + phi  grad(pi) = grad(pi) + grad(phi)
    // regular_timestep_comp:   pi = phi/dt    grad(pi) = grad(phi)/dt




}

// fill in Vproj
// initial_projection_comp: Utilde^0                        -- uold
// divu_iters_comp:         Utilde^0                        -- uold
// pressure_iters_comp:     (Utilde^n+1,* - Utilde^n)/dt    -- (unew-uold)/dt
// regular_timestep_comp:   (Utilde^n+1,* + dt*gpi/rhohalf) -- unew + dt*gpi/rhohalf
// sig contains rhohalf if proj_type == regular_timestep_comp
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



void
Maestro::doNodalProj (int proj_type,
                      Vector<MultiFab>& Vproj,
                      Vector<MultiFab>& phi,
                      Vector<MultiFab>& sig,
                      Vector<MultiFab>& rhcc,
                      Real rel_tol, Real abs_tol)
{
    BL_PROFILE("Projection:::doMLMGNodalProjection()");

    BL_ASSERT(Vproj[0].nGrow() == 1);
    BL_ASSERT(Vproj[0].nComp() == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(Vproj[0].boxArray().ixType().cellCentered());

    BL_ASSERT(phi[0].nGrow() == 1);
    BL_ASSERT(phi[0].nComp() == 1);
    AMREX_ALWAYS_ASSERT(phi[0].boxArray().ixType().nodeCentered());

    BL_ASSERT(sig[0].nGrow() == 1);
    BL_ASSERT(sig[0].nComp() == 1);
    AMREX_ALWAYS_ASSERT(sig[0].boxArray().ixType().cellCentered());

    BL_ASSERT(rhcc[0].nGrow() == 1);
    BL_ASSERT(rhcc[0].nComp() == 1);
    AMREX_ALWAYS_ASSERT(rhcc[0].boxArray().ixType().cellCentered());

    set_boundary_velocity(Vproj);

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
  
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    Vector<MultiFab> rhs(finest_level+1);
    for (int ilev = 0; ilev <= finest_level; ++ilev)
    {
        const auto& ba = amrex::convert(grids[ilev], IntVect::TheNodeVector());
        rhs[ilev].define(ba, dmap[ilev], 1, 0);
    }

        // define a nodal rhs and set it to zero (everything we want is in rhcc)
    Vector<MultiFab> rhnd(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rhnd[lev].define(convert(grids[lev],nodal_flag), dmap[lev], 1, 0);
        rhnd[lev].setVal(0.);
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs),
                    amrex::GetVecOfPtrs(Vproj),
                    amrex::GetVecOfConstPtrs(rhnd),
                    amrex::GetVecOfPtrs(rhcc));

    MLMG mlmg(mlndlap);
    mlmg.setMaxFmgIter(0);
    mlmg.setVerbose(0);

    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi),
                               amrex::GetVecOfConstPtrs(rhs),
                               rel_tol, abs_tol);

    // FIXME - update velocity, not Vproj
    // mlndlap.updateVelocity(Vproj, amrex::GetVecOfConstPtrs(phi));

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
