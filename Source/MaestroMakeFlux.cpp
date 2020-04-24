
#include <Maestro.H>

using namespace amrex;

const int pred_rhoh             = 0;
const int pred_rhohprime        = 1;
const int pred_h                = 2;
const int pred_T_then_rhohprime = 3;
const int pred_T_then_h         = 4;
const int pred_hprime           = 5;
const int pred_Tprime_then_h    = 6;

const int pred_rhoprime_and_X   = 1;
const int pred_rhoX             = 2;
const int pred_rho_and_X        = 3;

void
Maestro::MakeRhoXFlux (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                       Vector<MultiFab>& etarhoflux,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                       const RealVector& r0_old,
                       const RealVector& r0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_old,
                       const RealVector& r0_new,
                       const RealVector& r0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_new,
                       const RealVector& r0_predicted_edge,
                       int start_comp, int num_comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoXFlux()", MakeRhoXFlux);

    const int rho_comp = Rho;
    const int spec_comp = FirstSpec;
    const int nspec = NumSpec;
    const int max_lev = base_geom.max_radial_level;
    const int species_pred_type_loc = species_pred_type;
    const bool use_exact_base_state_loc = use_exact_base_state;
    const bool evolve_base_state_loc = evolve_base_state;

    for (int lev=0; lev<=finest_level; ++lev) {
   
#if (AMREX_SPACEDIM == 3)
        MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;

        if (spherical == 1) {
            rho0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            rho0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
            rho0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
            MultiFab::LinComb(rho0mac_edgex,0.5,r0mac_old[lev][0],0,0.5,r0mac_new[lev][0],0,0,1,1);
            MultiFab::LinComb(rho0mac_edgey,0.5,r0mac_old[lev][1],0,0.5,r0mac_new[lev][1],0,0,1,1);
            MultiFab::LinComb(rho0mac_edgez,0.5,r0mac_old[lev][2],0,0.5,r0mac_new[lev][2],0,0,1,1);
        }
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<Real> etarhoflux_arr = etarhoflux[lev].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
            const Array4<Real> sfluxy = sflux[lev][1].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
            const Array4<Real> sfluxz = sflux[lev][2].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
#endif

            Real * AMREX_RESTRICT w0_p = w0.dataPtr();
            const Real * AMREX_RESTRICT rho0_old_p = r0_old.dataPtr();
            const Real * AMREX_RESTRICT rho0_new_p = r0_new.dataPtr();
            const Real * AMREX_RESTRICT rho0_edge_old_p = r0_edge_old.dataPtr();
            const Real * AMREX_RESTRICT rho0_edge_new_p = r0_edge_new.dataPtr();
            const Real * AMREX_RESTRICT rho0_predicted_edge_p = r0_predicted_edge.dataPtr();

#if (AMREX_SPACEDIM == 2)

            // x-direction
            AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                int comp = n+start_comp;

                // reset density flux
                if (n == 0) {
                    sfluxx(i,j,k,0) = 0.0;
                }

                Real rho0_edge = 0.5*(rho0_old_p[lev+j*(max_lev+1)]+rho0_new_p[lev+j*(max_lev+1)]);

                if (species_pred_type_loc == pred_rhoprime_and_X) {
                    // edge states are rho' and X.  To make the (rho X) flux,
                    // we need the edge state of rho0
                    sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                        (rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rhoX) {
                    // edge states are (rho X)
                    sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rho_and_X) {
                    // edge states are rho and X
                    sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                }

                // compute the density fluxes by summing the species fluxes
                sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
            });

            // y-direction
            AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                int comp = n+start_comp;

                // reset density flux
                if (n == 0) {
                    sfluxy(i,j,k,0) = 0.0;
                }

                Real rho0_edge = 0.5*(rho0_edge_old_p[lev+j*(max_lev+1)]+rho0_edge_new_p[lev+j*(max_lev+1)]);

                if (species_pred_type_loc == pred_rhoprime_and_X) {
                    //   ! edge states are rho' and X.  To make the (rho X) flux,
                    //   ! we need the edge state of rho0
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rhoX) {
                    // ! edge states are (rho X)
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*sedgey(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rho_and_X) {
                    // ! edge state are rho and X
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                }

                if (evolve_base_state_loc && !use_exact_base_state_loc) {
                    if (comp >= spec_comp && comp <= spec_comp+nspec-1) {
                        etarhoflux_arr(i,j,k) += sfluxy(i,j,k,comp);
                    }

                    if (comp==spec_comp+nspec-1) {
                        etarhoflux_arr(i,j,k) -= w0_p[lev+j*(max_lev+1)]*rho0_predicted_edge_p[lev+j*(max_lev+1)];
                    }
                } 

                // compute the density fluxes by summing the species fluxes
                sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
            });

#elif (AMREX_SPACEDIM == 3)

            if (spherical == 0) {
                // x-direction
                AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxx(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            (rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
                });

                // y-direction
                AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxy(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            (rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxy(i,j,k,comp) = vmac(i,j,k)*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
                });

                // z-direction
                AMREX_PARALLEL_FOR_4D(zbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxz(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_edge_old_p[lev+k*(max_lev+1)]+rho0_edge_new_p[lev+k*(max_lev+1)]);

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        //   ! edge states are rho' and X.  To make the (rho X) flux,
                        //   ! we need the edge state of rho0
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*(rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // ! edge states are (rho X)
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // ! edge state are rho and X
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp);
                    }

                    if (evolve_base_state_loc && !use_exact_base_state_loc) {
                        if (comp >= spec_comp && comp <= spec_comp+nspec-1) {
                            etarhoflux_arr(i,j,k) += sfluxz(i,j,k,comp);
                        }

                        if (comp == spec_comp+nspec-1) {
                            etarhoflux_arr(i,j,k) -= w0_p[lev+k*(max_lev+1)]*rho0_predicted_edge_p[lev+k*(max_lev+1)];
                        }
                    } 

                    // compute the density fluxes by summing the species fluxes
                    sfluxz(i,j,k,0) += sfluxz(i,j,k,comp);
                });
            } else {

                const Array4<const Real> rho0_edgex = rho0mac_edgex.array(mfi);
                const Array4<const Real> rho0_edgey = rho0mac_edgey.array(mfi);
                const Array4<const Real> rho0_edgez = rho0mac_edgez.array(mfi);

                // x-direction
                AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxx(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            (rho0_edgex(i,j,k)+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
                });

                // y-direction
                AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxy(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            (rho0_edgey(i,j,k)+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxy(i,j,k,comp) = vmac(i,j,k)*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
                });

                // z-direction
                AMREX_PARALLEL_FOR_4D(zbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxz(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*(rho0_edgez(i,j,k)+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // ! edge states are (rho X)
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // ! edge state are rho and X
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxz(i,j,k,0) += sfluxz(i,j,k,comp);
                });
            } // end spherical
#endif
        } // end MFIter loop

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {

            // Get the grid size
            const Real* dx = geom[lev].CellSize();
            // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
            const Real area[3] = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};
#else
            const Real area[2] = {dx[1], dx[0]};
#endif

            if (flux_reg_s[lev+1])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,start_comp,start_comp,num_comp, -1.0*dt*area[i]);
                    // also include density flux
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,Rho,Rho,1, -1.0*dt*area[i]);
                }
            }
            if (flux_reg_s[lev])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,start_comp,start_comp,num_comp, 1.0*dt*area[i]);
                    // also include density flux
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,Rho,Rho,1, 1.0*dt*area[i]);
                }
            }

            if (spherical == 0) {
                // need edge_restrict for etarhoflux
            }
        }
    } // end loop over levels

    // average down fluxes
    if (reflux_type == 1) {
        AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}

void
Maestro::MakeRhoXFlux (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                       Vector<MultiFab>& etarhoflux,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                       const RealVector& r0_old,
                       const BaseState<Real>& rho0_edge_old_state,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_old,
                       const RealVector& r0_new,
                       const BaseState<Real>& rho0_edge_new_state,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_new,
                       const BaseState<Real>& rho0_predicted_edge_state,
                       int start_comp, int num_comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoXFlux()", MakeRhoXFlux);

    const int rho_comp = Rho;
    const int spec_comp = FirstSpec;
    const int nspec = NumSpec;
    const int max_lev = base_geom.max_radial_level;
    const int species_pred_type_loc = species_pred_type;
    const bool use_exact_base_state_loc = use_exact_base_state;
    const bool evolve_base_state_loc = evolve_base_state;

    auto rho0_edge_old = rho0_edge_old_state.array();
    auto rho0_edge_new = rho0_edge_new_state.array();
    auto rho0_predicted_edge = rho0_predicted_edge_state.array();

    for (int lev=0; lev<=finest_level; ++lev) {
   
#if (AMREX_SPACEDIM == 3)
        MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;

        if (spherical == 1) {
            rho0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
            rho0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
            rho0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
            MultiFab::LinComb(rho0mac_edgex,0.5,r0mac_old[lev][0],0,0.5,r0mac_new[lev][0],0,0,1,1);
            MultiFab::LinComb(rho0mac_edgey,0.5,r0mac_old[lev][1],0,0.5,r0mac_new[lev][1],0,0,1,1);
            MultiFab::LinComb(rho0mac_edgez,0.5,r0mac_old[lev][2],0,0.5,r0mac_new[lev][2],0,0,1,1);
        }
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif

            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<Real> etarhoflux_arr = etarhoflux[lev].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
            const Array4<Real> sfluxy = sflux[lev][1].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
            const Array4<Real> sfluxz = sflux[lev][2].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
#endif

            Real * AMREX_RESTRICT w0_p = w0.dataPtr();
            const Real * AMREX_RESTRICT rho0_old_p = r0_old.dataPtr();
            const Real * AMREX_RESTRICT rho0_new_p = r0_new.dataPtr();
            
#if (AMREX_SPACEDIM == 2)

            // x-direction
            AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                int comp = n+start_comp;

                // reset density flux
                if (n == 0) {
                    sfluxx(i,j,k,0) = 0.0;
                }

                Real rho0_edge = 0.5*(rho0_old_p[lev+j*(max_lev+1)]+rho0_new_p[lev+j*(max_lev+1)]);

                if (species_pred_type_loc == pred_rhoprime_and_X) {
                    // edge states are rho' and X.  To make the (rho X) flux,
                    // we need the edge state of rho0
                    sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                        (rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rhoX) {
                    // edge states are (rho X)
                    sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rho_and_X) {
                    // edge states are rho and X
                    sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                        sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                }

                // compute the density fluxes by summing the species fluxes
                sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
            });

            // y-direction
            AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                int comp = n+start_comp;

                // reset density flux
                if (n == 0) {
                    sfluxy(i,j,k,0) = 0.0;
                }

                Real rho0_edge = 0.5*(rho0_edge_old(lev,j)+rho0_edge_new(lev,j));

                if (species_pred_type_loc == pred_rhoprime_and_X) {
                    //   ! edge states are rho' and X.  To make the (rho X) flux,
                    //   ! we need the edge state of rho0
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rhoX) {
                    // ! edge states are (rho X)
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*sedgey(i,j,k,comp);

                } else if (species_pred_type_loc == pred_rho_and_X) {
                    // ! edge state are rho and X
                    sfluxy(i,j,k,comp) = 
                        vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                }

                if (evolve_base_state_loc && !use_exact_base_state_loc) {
                    if (comp >= spec_comp && comp <= spec_comp+nspec-1) {
                        etarhoflux_arr(i,j,k) += sfluxy(i,j,k,comp);
                    }

                    if (comp==spec_comp+nspec-1) {
                        etarhoflux_arr(i,j,k) -= w0_p[lev+j*(max_lev+1)]*rho0_predicted_edge(lev,j);
                    }
                } 

                // compute the density fluxes by summing the species fluxes
                sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
            });

#elif (AMREX_SPACEDIM == 3)

            if (spherical == 0) {
                // x-direction
                AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxx(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            (rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
                });

                // y-direction
                AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxy(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            (rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxy(i,j,k,comp) = vmac(i,j,k)*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
                });

                // z-direction
                AMREX_PARALLEL_FOR_4D(zbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxz(i,j,k,0) = 0.0;
                    }

                    Real rho0_edge = 0.5*(rho0_edge_old(lev,k)+rho0_edge_new(lev,k));

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        //   ! edge states are rho' and X.  To make the (rho X) flux,
                        //   ! we need the edge state of rho0
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*(rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // ! edge states are (rho X)
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // ! edge state are rho and X
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp);
                    }

                    if (evolve_base_state_loc && !use_exact_base_state_loc) {
                        if (comp >= spec_comp && comp <= spec_comp+nspec-1) {
                            etarhoflux_arr(i,j,k) += sfluxz(i,j,k,comp);
                        }

                        if (comp == spec_comp+nspec-1) {
                            etarhoflux_arr(i,j,k) -= w0_p[lev+k*(max_lev+1)]*rho0_predicted_edge(lev,k);
                        }
                    } 

                    // compute the density fluxes by summing the species fluxes
                    sfluxz(i,j,k,0) += sfluxz(i,j,k,comp);
                });
            } else {

                const Array4<const Real> rho0_edgex = rho0mac_edgex.array(mfi);
                const Array4<const Real> rho0_edgey = rho0mac_edgey.array(mfi);
                const Array4<const Real> rho0_edgez = rho0mac_edgez.array(mfi);

                // x-direction
                AMREX_PARALLEL_FOR_4D(xbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxx(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            (rho0_edgex(i,j,k)+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxx(i,j,k,comp) = umacx(i,j,k)*sedgex(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxx(i,j,k,comp) = umacx(i,j,k)* 
                            sedgex(i,j,k,rho_comp)*sedgex(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxx(i,j,k,0) += sfluxx(i,j,k,comp);
                });

                // y-direction
                AMREX_PARALLEL_FOR_4D(ybx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxy(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            (rho0_edgey(i,j,k)+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // edge states are (rho X)
                        sfluxy(i,j,k,comp) = vmac(i,j,k)*sedgey(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // edge states are rho and X
                        sfluxy(i,j,k,comp) = vmac(i,j,k)* 
                            sedgey(i,j,k,rho_comp)*sedgey(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxy(i,j,k,0) += sfluxy(i,j,k,comp);
                });

                // z-direction
                AMREX_PARALLEL_FOR_4D(zbx, num_comp, i, j, k, n, {
                    int comp = n + start_comp;

                    // reset density flux
                    if (n == 0) {
                        sfluxz(i,j,k,0) = 0.0;
                    }

                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // edge states are rho' and X.  To make the (rho X) flux,
                        // we need the edge state of rho0
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*(rho0_edgez(i,j,k)+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rhoX) {
                        // ! edge states are (rho X)
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,comp);

                    } else if (species_pred_type_loc == pred_rho_and_X) {
                        // ! edge state are rho and X
                        sfluxz(i,j,k,comp) = 
                            wmac(i,j,k)*sedgez(i,j,k,rho_comp)*sedgez(i,j,k,comp);
                    }

                    // compute the density fluxes by summing the species fluxes
                    sfluxz(i,j,k,0) += sfluxz(i,j,k,comp);
                });
            } // end spherical
#endif
        } // end MFIter loop

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {

            // Get the grid size
            const Real* dx = geom[lev].CellSize();
            // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
            const Real area[3] = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};
#else
            const Real area[2] = {dx[1], dx[0]};
#endif

            if (flux_reg_s[lev+1])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,start_comp,start_comp,num_comp, -1.0*dt*area[i]);
                    // also include density flux
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,Rho,Rho,1, -1.0*dt*area[i]);
                }
            }
            if (flux_reg_s[lev])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,start_comp,start_comp,num_comp, 1.0*dt*area[i]);
                    // also include density flux
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,Rho,Rho,1, 1.0*dt*area[i]);
                }
            }

            if (spherical == 0) {
                // need edge_restrict for etarhoflux
            }
        }
    } // end loop over levels

    // average down fluxes
    if (reflux_type == 1) {
        AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}

void
Maestro::MakeRhoHFlux (const Vector<MultiFab>& state,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                       Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& w0mac,
                       const RealVector& r0_old,
                       const BaseState<Real>& rho0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_old,
                       const RealVector& r0_new,
                       const BaseState<Real>& rho0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& r0mac_new,
                       const BaseState<Real>& rhoh0_old_in,
                       const BaseState<Real>& rhoh0_edge_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& rh0mac_old,
                       const BaseState<Real>& rhoh0_new_in,
                       const BaseState<Real>& rhoh0_edge_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& rh0mac_new,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& h0mac_old,
                       const Vector<std::array< MultiFab, AMREX_SPACEDIM > >& h0mac_new)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeRhoHFlux()", MakeRhoHFlux);

    const bool have_h = enthalpy_pred_type == pred_h ||
                        enthalpy_pred_type == pred_T_then_h || 
                        enthalpy_pred_type == pred_Tprime_then_h;
#if (AMREX_SPACEDIM == 3)
    const bool have_hprime = enthalpy_pred_type == pred_hprime;
#endif
    const bool have_rhoh = enthalpy_pred_type == pred_rhoh;

    const int rho_comp = Rho;
    const int rhoh_comp = RhoH;
    const int max_lev = base_geom.max_radial_level;
    const int species_pred_type_loc = species_pred_type;
    const int enthalpy_pred_type_loc = enthalpy_pred_type;

    for (int lev=0; lev<=finest_level; ++lev) {

#if (AMREX_SPACEDIM == 3)
        MultiFab rho0mac_edgex, rho0mac_edgey, rho0mac_edgez;
        MultiFab h0mac_edgex, h0mac_edgey, h0mac_edgez;
        MultiFab rhoh0mac_edgex, rhoh0mac_edgey, rhoh0mac_edgez;
        
        rho0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
        rho0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
        rho0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);
        h0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
        h0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
        h0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);
        rhoh0mac_edgex.define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);
        rhoh0mac_edgey.define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);
        rhoh0mac_edgez.define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0);

        rho0mac_edgex.setVal(0.);
        rho0mac_edgey.setVal(0.);
        rho0mac_edgez.setVal(0.);

        h0mac_edgex.setVal(0.);
        h0mac_edgey.setVal(0.);
        h0mac_edgez.setVal(0.);

        rhoh0mac_edgex.setVal(0.);
        rhoh0mac_edgey.setVal(0.);
        rhoh0mac_edgez.setVal(0.);

        if (spherical == 1) {
            if (use_exact_base_state) {
                MultiFab::LinComb(rhoh0mac_edgex,0.5,rh0mac_old[lev][0],0,0.5,rh0mac_new[lev][0],0,0,1,0);
                MultiFab::LinComb(rhoh0mac_edgey,0.5,rh0mac_old[lev][1],0,0.5,rh0mac_new[lev][1],0,0,1,0);
                MultiFab::LinComb(rhoh0mac_edgez,0.5,rh0mac_old[lev][2],0,0.5,rh0mac_new[lev][2],0,0,1,0);
            } else {
                MultiFab::LinComb(rho0mac_edgex,0.5,r0mac_old[lev][0],0,0.5,r0mac_new[lev][0],0,0,1,0);
                MultiFab::LinComb(rho0mac_edgey,0.5,r0mac_old[lev][1],0,0.5,r0mac_new[lev][1],0,0,1,0);
                MultiFab::LinComb(rho0mac_edgez,0.5,r0mac_old[lev][2],0,0.5,r0mac_new[lev][2],0,0,1,0);
                MultiFab::LinComb(h0mac_edgex,0.5,h0mac_old[lev][0],0,0.5,h0mac_new[lev][0],0,0,1,0);
                MultiFab::LinComb(h0mac_edgey,0.5,h0mac_old[lev][1],0,0.5,h0mac_new[lev][1],0,0,1,0);
                MultiFab::LinComb(h0mac_edgez,0.5,h0mac_old[lev][2],0,0.5,h0mac_new[lev][2],0,0,1,0);
            }
        }
#endif

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Get the index space of the valid region
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif
            const Array4<Real> sedgex = sedge[lev][0].array(mfi);
            const Array4<Real> sfluxx = sflux[lev][0].array(mfi);
            const Array4<const Real> umacx = umac[lev][0].array(mfi);
            const Array4<Real> sedgey = sedge[lev][1].array(mfi);
            const Array4<Real> sfluxy = sflux[lev][1].array(mfi);
            const Array4<const Real> vmac = umac[lev][1].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const Array4<Real> sedgez = sedge[lev][2].array(mfi);
            const Array4<Real> sfluxz = sflux[lev][2].array(mfi);
            const Array4<const Real> wmac = umac[lev][2].array(mfi);
#endif

            const Real * AMREX_RESTRICT rho0_old_p = r0_old.dataPtr();
            const Real * AMREX_RESTRICT rho0_new_p = r0_new.dataPtr();

            const auto rho0_edge_old_arr = rho0_edge_old.const_array();
            const auto rho0_edge_new_arr = rho0_edge_new.const_array();

            const auto rhoh0_old_arr = rhoh0_old_in.const_array();
            const auto rhoh0_new_arr = rhoh0_new_in.const_array();
            const auto rhoh0_edge_old_arr = rhoh0_edge_old.const_array();
            const auto rhoh0_edge_new_arr = rhoh0_edge_new.const_array();

#if (AMREX_SPACEDIM == 2)
            AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                // create x-fluxes
                if (have_h) {
                    // enthalpy edge state is h
                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // density edge state is rho'
                        Real rho0_edge = 0.5*(rho0_old_p[lev+j*(max_lev+1)]+rho0_new_p[lev+j*(max_lev+1)]);

                        sfluxx(i,j,k,rhoh_comp) = 
                           umacx(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp);

                    } else if (species_pred_type_loc == pred_rho_and_X ||
                               species_pred_type_loc == pred_rhoX) {
                        // density edge state is rho
                        sfluxx(i,j,k,rhoh_comp) = 
                            umacx(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp);
                    }
                } else if (have_rhoh) {
                    sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*sedgex(i,j,k,rhoh_comp);

                } else if (enthalpy_pred_type_loc == pred_rhohprime || 
                    enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                    //   ! enthalpy edge state is (rho h)'
                    Real rhoh0_edge = 0.5*(rhoh0_old_arr(lev,j)+rhoh0_new_arr(lev,j));

                    sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*(rhoh0_edge+sedgex(i,j,k,rhoh_comp));
                }
            });

            AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                // create y-fluxes
                if (have_h) {
                    // enthalpy edge state is h
                    if (species_pred_type_loc == pred_rhoprime_and_X) {
                        // density edge state is rho'
                        Real rho0_edge = 0.5*(rho0_edge_old_arr(lev,j)+rho0_edge_new_arr(lev,j));

                        sfluxy(i,j,k,rhoh_comp) = 
                            vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp);

                    } else if (species_pred_type_loc == pred_rho_and_X || 
                               species_pred_type_loc == pred_rhoX) {
                        // density edge state is rho
                        sfluxy(i,j,k,rhoh_comp) = 
                            vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp);
                    }
                } else if (have_rhoh) {
                    sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp);

                } else if (enthalpy_pred_type_loc == pred_rhohprime || 
                    enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                    // enthalpy edge state is (rho h)'
                    Real rhoh0_edge = 0.5*(rhoh0_edge_old_arr(lev,j)+rhoh0_edge_new_arr(lev,j));
                    
                    sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*(sedgey(i,j,k,rhoh_comp)+rhoh0_edge);
                }
            });

#elif (AMREX_SPACEDIM == 3)

            if (!spherical) {
                AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                    // create x-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                            sfluxx(i,j,k,rhoh_comp) = 
                            umacx(i,j,k)*(rho0_edge+sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp);

                        } else if (species_pred_type_loc == pred_rho_and_X ||
                                species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxx(i,j,k,rhoh_comp) = 
                                umacx(i,j,k)*sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp);
                        }
                    } else if (have_rhoh) {
                        sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*sedgex(i,j,k,rhoh_comp);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime || 
                        enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                        //   ! enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5*(rhoh0_old_arr(lev,k)+rhoh0_new_arr(lev,k));

                        sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*(rhoh0_edge+sedgex(i,j,k,rhoh_comp));
                    }
                });

                AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                    // create y-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5*(rho0_old_p[lev+k*(max_lev+1)]+rho0_new_p[lev+k*(max_lev+1)]);

                            sfluxy(i,j,k,rhoh_comp) = 
                            vmac(i,j,k)*(rho0_edge+sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp);

                        } else if (species_pred_type_loc == pred_rho_and_X ||
                                species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxy(i,j,k,rhoh_comp) = 
                                vmac(i,j,k)*sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp);
                        }
                    } else if (have_rhoh) {
                        sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime ||
                        enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                        //   ! enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5*(rhoh0_old_arr(lev,k)+rhoh0_new_arr(lev,k));

                        sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*(rhoh0_edge+sedgey(i,j,k,rhoh_comp));
                    }
                });

                AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                    // create z-fluxes
                    if (have_h) {
                        // enthalpy edge state is h
                        if (species_pred_type_loc == pred_rhoprime_and_X) {
                            // density edge state is rho'
                            Real rho0_edge = 0.5*(rho0_edge_old_arr(lev,k)+rho0_edge_new_arr(lev,k));

                            sfluxz(i,j,k,rhoh_comp) = 
                                wmac(i,j,k)*(rho0_edge+sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp);

                        } else if (species_pred_type_loc == pred_rho_and_X || 
                                species_pred_type_loc == pred_rhoX) {
                            // density edge state is rho
                            sfluxz(i,j,k,rhoh_comp) = 
                                wmac(i,j,k)*sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp);
                        }
                    } else if (have_rhoh) {
                        sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*sedgez(i,j,k,rhoh_comp);

                    } else if (enthalpy_pred_type_loc == pred_rhohprime || 
                        enthalpy_pred_type_loc == pred_T_then_rhohprime) {
                        // enthalpy edge state is (rho h)'
                        Real rhoh0_edge = 0.5*(rhoh0_edge_old_arr(lev,k)+rhoh0_edge_new_arr(lev,k));
                        
                        sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*(sedgez(i,j,k,rhoh_comp)+rhoh0_edge);
                    }
                });
            } else {

                if (use_exact_base_state) {

                    const Array4<const Real> rhoh0_edgex = rhoh0mac_edgex.array(mfi);
                    const Array4<const Real> rhoh0_edgey = rhoh0mac_edgey.array(mfi);
                    const Array4<const Real> rhoh0_edgez = rhoh0mac_edgez.array(mfi);

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*sedgex(i,j,k,rhoh_comp);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0

                            sfluxx(i,j,k,rhoh_comp) = 
                                umacx(i,j,k)*(rhoh0_edgex(i,j,k)+sedgex(i,j,k,rhoh_comp));
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0
                            sfluxy(i,j,k,rhoh_comp) = 
                                vmac(i,j,k)*(rhoh0_edgey(i,j,k)+sedgey(i,j,k,rhoh_comp));
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {
                        if (have_h) {
                            // enthalpy edge state is h
                            // this is not supported on irregular-spaced base state
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            // this is not supported on irregular-spaced base state
                        } else if (have_rhoh) {
                            sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*sedgez(i,j,k,rhoh_comp);
                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + rhoh_0
                            sfluxz(i,j,k,rhoh_comp) = 
                                wmac(i,j,k)*(rhoh0_edgez(i,j,k)+sedgez(i,j,k,rhoh_comp));
                        }
                    });
                } else {
                    const Array4<const Real> rho0_edgex = rho0mac_edgex.array(mfi);
                    const Array4<const Real> rho0_edgey = rho0mac_edgey.array(mfi);
                    const Array4<const Real> rho0_edgez = rho0mac_edgez.array(mfi);

                    const Array4<const Real> h0_edgex = h0mac_edgex.array(mfi);
                    const Array4<const Real> h0_edgey = h0mac_edgey.array(mfi);
                    const Array4<const Real> h0_edgez = h0mac_edgez.array(mfi);

                    AMREX_PARALLEL_FOR_3D(xbx, i, j, k, {

                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k) * 
                                    (rho0_edgex(i,j,k) + sedgex(i,j,k,rho_comp))*sedgex(i,j,k,rhoh_comp);

                            } else if (species_pred_type_loc == pred_rho_and_X || 
                                species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k) * 
                                    sedgex(i,j,k,rho_comp)*sedgex(i,j,k,rhoh_comp);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k) * 
                                    (sedgex(i,j,k,rho_comp)+rho0_edgex(i,j,k)) * (sedgex(i,j,k,rhoh_comp)+h0_edgex(i,j,k));
                            }

                        } else if (have_rhoh) {

                            sfluxx(i,j,k,rhoh_comp) = umacx(i,j,k)*sedgex(i,j,k,rhoh_comp);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0

                            sfluxx(i,j,k,rhoh_comp) = 
                                umacx(i,j,k)*(rho0_edgex(i,j,k)*h0_edgex(i,j,k)+sedgex(i,j,k,rhoh_comp));
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(ybx, i, j, k, {

                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k) * 
                                    (rho0_edgey(i,j,k) + sedgey(i,j,k,rho_comp))*sedgey(i,j,k,rhoh_comp);

                            } else if (species_pred_type_loc == pred_rho_and_X || 
                                species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k) * 
                                    sedgey(i,j,k,rho_comp)*sedgey(i,j,k,rhoh_comp);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k) * 
                                    (sedgey(i,j,k,rho_comp)+rho0_edgey(i,j,k)) * (sedgey(i,j,k,rhoh_comp)+h0_edgey(i,j,k));
                            }
                        } else if (have_rhoh) {

                            sfluxy(i,j,k,rhoh_comp) = vmac(i,j,k)*sedgey(i,j,k,rhoh_comp);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0
                            sfluxy(i,j,k,rhoh_comp) = 
                                vmac(i,j,k)*(rho0_edgey(i,j,k)*h0_edgey(i,j,k)+sedgey(i,j,k,rhoh_comp));
                        }
                    });

                    AMREX_PARALLEL_FOR_3D(zbx, i, j, k, {

                        if (have_h) {
                            // enthalpy edge state is h
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'
                                sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k) * 
                                    (rho0_edgez(i,j,k) + sedgez(i,j,k,rho_comp))*sedgez(i,j,k,rhoh_comp);

                            } else if (species_pred_type_loc == pred_rho_and_X || 
                                species_pred_type_loc == pred_rhoX) {
                                // density edge state is rho
                                sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k) * 
                                    sedgez(i,j,k,rho_comp)*sedgez(i,j,k,rhoh_comp);
                            }
                        } else if (have_hprime) {
                            // enthalpy edge state is h'
                            if (species_pred_type_loc == pred_rhoprime_and_X) {
                                // density edge state is rho'

                                // (rho h)_edge = (h' + h_0) * (rho' + rho_0) where h0 is
                                // computed from (rho h)_0 / rho_0
                                // sfluxx = (umac(i,j,k)+w0macx(i,j,k)) * (rho h)_edge
                                sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k) * 
                                    (sedgez(i,j,k,rho_comp)+rho0_edgez(i,j,k)) * (sedgez(i,j,k,rhoh_comp)+h0_edgez(i,j,k));
                            }
                        } else if (have_rhoh) {

                            sfluxz(i,j,k,rhoh_comp) = wmac(i,j,k)*sedgez(i,j,k,rhoh_comp);

                        } else {
                            // enthalpy edge state is (rho h)'

                            // Average (rho h) onto edges by averaging rho and h
                            // separately onto edges.
                            //  (rho h)_edge = (rho h)' + (rho_0 * h_0)
                            // where h_0 is computed from (rho h)_0 / rho_0
                            sfluxz(i,j,k,rhoh_comp) = 
                                wmac(i,j,k)*(rho0_edgez(i,j,k)*h0_edgez(i,j,k)+sedgez(i,j,k,rhoh_comp));
                        }
                    });
                }
            }
#endif
        } // end MFIter loop

        // increment or decrement the flux registers by area and time-weighted fluxes
        // Note that the fluxes need to be scaled by dt and area
        // In this example we are solving s_t = -div(+F)
        // The fluxes contain, e.g., F_{i+1/2,j} = (s*u)_{i+1/2,j}
        // Keep this in mind when considering the different sign convention for updating
        // the flux registers from the coarse or fine grid perspective
        // NOTE: the flux register associated with flux_reg_s[lev] is associated
        // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
        if (reflux_type == 2) {

            // Get the grid size
            const Real* dx = geom[lev].CellSize();
            // NOTE: areas are different in DIM=2 and DIM=3
#if (AMREX_SPACEDIM == 3)
            const Real area[3] = {dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1]};
#else
            const Real area[2] = {dx[1], dx[0]};
#endif

            if (flux_reg_s[lev+1])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev+1/lev flux register (index lev+1)
                    flux_reg_s[lev+1]->CrseInit(sflux[lev][i],i,RhoH,RhoH,1, -1.0*dt*area[i]);
                }
            }
            if (flux_reg_s[lev])
            {
                for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                    // update the lev/lev-1 flux register (index lev)
                    flux_reg_s[lev]->FineAdd(sflux[lev][i],i,RhoH,RhoH,1, 1.0*dt*area[i]);
                }
            }
        }
    } // end loop over levels

    if (reflux_type == 1) {
        AverageDownFaces(sflux);
    }

    // Something analogous to edge_restriction is done in UpdateScal()
}
