
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::DensityAdvance(
    int which_step, Vector<MultiFab>& scalold, Vector<MultiFab>& scalnew,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sedge,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sflux,
    Vector<MultiFab>& scal_force, Vector<MultiFab>& etarhoflux,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const BaseState<Real>& rho0_predicted_edge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DensityAdvance()", DensityAdvance);

    BaseState<Real> rho0_edge_old(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine + 1);
    BaseState<Real> rho0_edge_new(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine + 1);

    if (!spherical) {
        // create edge-centered base state quantities.
        // Note: rho0_edge_{old,new}
        // contains edge-centered quantities created via spatial interpolation.
        // This is to be contrasted to rho0_predicted_edge which is the half-time
        // edge state created in advect_base.
        CelltoEdge(rho0_old, rho0_edge_old);
        CelltoEdge(rho0_new, rho0_edge_new);
    }

    //////////////////////////////////
    // Create source terms at time n
    //////////////////////////////////

    // source terms for X and for tracers are zero - do nothing
    for (int lev = 0; lev <= finest_level; ++lev) {
        scal_force[lev].setVal(0.);
    }

    Vector<MultiFab> rho0_old_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rho0_old_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        rho0_old_cart[lev].setVal(0.);
    }

    Put1dArrayOnCart(rho0_old, rho0_old_cart, false, false, bcs_s, Rho);

    /////////////////////////////////////////////////////////////////
    // Subtract w0 from MAC velocities (MAC velocities has w0 already).
    /////////////////////////////////////////////////////////////////

    Addw0(umac, w0mac, -1.);

    /////////////////////////////////////////////////////////////////
    // Compute source terms
    /////////////////////////////////////////////////////////////////

    // ** density source term **

    // Make source term for rho or rho'
    if (species_pred_type == predict_rhoprime_and_X) {
        // rho' source term
        // this is needed for pred_rhoprime_and_X
        ModifyScalForce(scal_force, scalold, umac, rho0_edge_old, rho0_old_cart,
                        Rho, bcs_s, 0);

    } else if (species_pred_type == predict_rho_and_X) {
        // rho source term
        ModifyScalForce(scal_force, scalold, umac, rho0_edge_old, rho0_old_cart,
                        Rho, bcs_s, 1);
    }

    // ** species source term **

    // for species_pred_types predict_rhoprime_and_X and
    // predict_rho_and_X, there is no force for X.

    // for predict_rhoX, we are predicting (rho X)
    // as a conservative equation, and there is no force.

    /////////////////////////////////////////////////////////////////
    // Add w0 back to MAC velocities (trans velocities already have w0).
    /////////////////////////////////////////////////////////////////

    Addw0(umac, w0mac, 1.);

    /////////////////////////////////////////////////////////////////
    // Create the edge states of (rho X)' or X and rho'
    /////////////////////////////////////////////////////////////////

    if ((species_pred_type == predict_rhoprime_and_X) ||
        (species_pred_type == predict_rho_and_X)) {
        // we are predicting X to the edges, so convert the scalar
        // data to those quantities

        // convert (rho X) --> X in scalold
        ConvertRhoXToX(scalold, true);
    }

    if (species_pred_type == predict_rhoprime_and_X) {
        // convert rho -> rho' in scalold
        //   . this is needed for predict_rhoprime_and_X
        PutInPertForm(scalold, rho0_old, Rho, 0, bcs_f, true);
    }

    // predict species at the edges -- note, either X or (rho X) will be
    // predicted here, depending on species_pred_type

    const auto is_vel = false;  // false
    if (species_pred_type == predict_rhoprime_and_X ||
        species_pred_type == predict_rho_and_X) {
        // we are predicting X to the edges, using the advective form of
        // the prediction
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     FirstSpec, FirstSpec, NumSpec, false);

    } else if (species_pred_type == predict_rhoX) {
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     FirstSpec, FirstSpec, NumSpec, true);
    }

    // predict rho or rho' at the edges (depending on species_pred_type)
    if (species_pred_type == predict_rhoprime_and_X ||
        species_pred_type == predict_rho_and_X) {
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     Rho, Rho, 1, false);

    } else if (species_pred_type == predict_rhoX) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Copy(sedge[lev][idim], sedge[lev][idim], FirstSpec,
                               Rho, 1, 0);
                for (int ispec = 1; ispec < NumSpec; ++ispec) {
                    MultiFab::Add(sedge[lev][idim], sedge[lev][idim],
                                  FirstSpec + ispec, Rho, 1, 0);
                }
            }
        }
    }

    if (species_pred_type == predict_rhoprime_and_X) {
        // convert rho' -> rho in scalold
        PutInPertForm(scalold, rho0_old, Rho, Rho, bcs_s, false);
    }

    if ((species_pred_type == predict_rhoprime_and_X) ||
        (species_pred_type == predict_rho_and_X)) {
        // convert X --> (rho X) in scalold
        ConvertRhoXToX(scalold, false);
    }

    /////////////////////////////////////////////////////////////////
    // Compute fluxes
    /////////////////////////////////////////////////////////////////

    if (which_step == 1) {
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_old(finest_level +
                                                                  1);
#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                AMREX_D_TERM(
                    rho0mac_old[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_old[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_old[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););
            }
            MakeS0mac(rho0_old, rho0mac_old);
        }
#endif

        // compute species fluxes
        MakeRhoXFlux(scalold, sflux, etarhoflux, sedge, umac, rho0_old,
                     rho0_edge_old, rho0mac_old, rho0_old, rho0_edge_old,
                     rho0mac_old, rho0_predicted_edge, FirstSpec, NumSpec);

    } else if (which_step == 2) {
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_old(finest_level +
                                                                  1);
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_new(finest_level +
                                                                  1);

#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                AMREX_D_TERM(
                    rho0mac_old[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_old[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_old[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););

                AMREX_D_TERM(
                    rho0mac_new[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_new[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_new[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););
            }
            MakeS0mac(rho0_old, rho0mac_old);
            MakeS0mac(rho0_new, rho0mac_new);
        }
#endif

        // compute species fluxes
        MakeRhoXFlux(scalold, sflux, etarhoflux, sedge, umac, rho0_old,
                     rho0_edge_old, rho0mac_old, rho0_new, rho0_edge_new,
                     rho0mac_new, rho0_predicted_edge, FirstSpec, NumSpec);
    }

    //**************************************************************************
    //     1) Set force for (rho X)_i at time n+1/2 = 0.
    //     2) Update (rho X)_i with conservative differencing.
    //     3) Define density as the sum of the (rho X)_i
    //     4) Update tracer with conservative differencing as well.
    //**************************************************************************

    for (int lev = 0; lev <= finest_level; ++lev) {
        scal_force[lev].setVal(0.);
    }

    Vector<MultiFab> p0_new_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_new_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_new, p0_new_cart, false, false, bcs_f, 0);

    // p0 only used in rhoh update so it's an optional parameter
    UpdateScal(scalold, scalnew, sflux, scal_force, FirstSpec, NumSpec,
               p0_new_cart);
}

// Density advance for SDC using intra(global var)
void Maestro::DensityAdvanceSDC(
    int which_step, Vector<MultiFab>& scalold, Vector<MultiFab>& scalnew,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sedge,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& sflux,
    Vector<MultiFab>& scal_force, Vector<MultiFab>& etarhoflux,
    Vector<std::array<MultiFab, AMREX_SPACEDIM> >& umac,
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
    const BaseState<Real>& rho0_predicted_edge) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DensityAdvanceSDC()", DensityAdvanceSDC);

    BaseState<Real> rho0_edge_old(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine + 1);
    BaseState<Real> rho0_edge_new(base_geom.max_radial_level + 1,
                                  base_geom.nr_fine + 1);

    if (!spherical) {
        // create edge-centered base state quantities.
        // Note: rho0_edge_{old,new}
        // contains edge-centered quantities created via spatial interpolation.
        // This is to be contrasted to rho0_predicted_edge which is the half-time
        // edge state created in advect_base.
        CelltoEdge(rho0_old, rho0_edge_old);
        CelltoEdge(rho0_new, rho0_edge_new);
    }

    //////////////////////////////////
    // Create source terms at time n
    //////////////////////////////////

    // source terms for X and for tracers include reaction forcing terms
    for (int lev = 0; lev <= finest_level; ++lev) {
        scal_force[lev].setVal(0.);
        MultiFab::Add(scal_force[lev], intra[lev], FirstSpec, FirstSpec,
                      NumSpec, 0);
    }

    if (finest_level == 0) {
        // fill periodic ghost cells
        for (int lev = 0; lev <= finest_level; ++lev) {
            scal_force[lev].FillBoundary(geom[lev].periodicity());
        }
    }
    // fill ghost cells behind physical boundaries
    // !!!!!! uncertain about this
    FillPatch(t_old, scal_force, scal_force, scal_force, FirstSpec, FirstSpec,
              NumSpec, FirstSpec, bcs_f);

    Vector<MultiFab> rho0_old_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        rho0_old_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(rho0_old, rho0_old_cart, false, false, bcs_s, Rho);

    /////////////////////////////////////////////////////////////////
    // Subtract w0 from MAC velocities (MAC velocities has w0 already).
    /////////////////////////////////////////////////////////////////

    Addw0(umac, w0mac, -1.);

    /////////////////////////////////////////////////////////////////
    // Compute source terms
    /////////////////////////////////////////////////////////////////

    // ** density source term **

    // Make source term for rho or rho'
    if (species_pred_type == predict_rhoprime_and_X) {
        // rho' source term
        // this is needed for pred_rhoprime_and_X
        ModifyScalForce(scal_force, scalold, umac, rho0_edge_old, rho0_old_cart,
                        Rho, bcs_s, 0);
    } else if (species_pred_type == predict_rho_and_X) {
        // rho source term
        ModifyScalForce(scal_force, scalold, umac, rho0_edge_old, rho0_old_cart,
                        Rho, bcs_s, 1);
    }

    // ** species source term **

    // for species_pred_types predict_rhoprime_and_X and
    // predict_rho_and_X, there is no force for X.

    // for predict_rhoX, we are predicting (rho X)
    // as a conservative equation, and there is no force.

    /////////////////////////////////////////////////////////////////
    // Add w0 to MAC velocities (trans velocities already have w0).
    /////////////////////////////////////////////////////////////////

    Addw0(umac, w0mac, 1.);

    /////////////////////////////////////////////////////////////////
    // Create the edge states of (rho X)' or X and rho'
    /////////////////////////////////////////////////////////////////

    if ((species_pred_type == predict_rhoprime_and_X) ||
        (species_pred_type == predict_rho_and_X)) {
        // we are predicting X to the edges, so convert the scalar
        // data to those quantities

        // convert (rho X) --> X in scalold
        ConvertRhoXToX(scalold, true);
    }

    if (species_pred_type == predict_rhoprime_and_X) {
        // convert rho -> rho' in scalold
        //   . this is needed for predict_rhoprime_and_X
        PutInPertForm(scalold, rho0_old, Rho, 0, bcs_f, true);
    }

    // predict species at the edges -- note, either X or (rho X) will be
    // predicted here, depending on species_pred_type

    const auto is_vel = false;  // false
    if (species_pred_type == predict_rhoprime_and_X ||
        species_pred_type == predict_rho_and_X) {
        // we are predicting X to the edges, using the advective form of
        // the prediction
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     FirstSpec, FirstSpec, NumSpec, false);

    } else if (species_pred_type == predict_rhoX) {
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     FirstSpec, FirstSpec, NumSpec, true);
    }

    // predict rho or rho' at the edges (depending on species_pred_type)
    if (species_pred_type == predict_rhoprime_and_X ||
        species_pred_type == predict_rho_and_X) {
        MakeEdgeScal(scalold, sedge, umac, scal_force, is_vel, bcs_s, Nscal,
                     Rho, Rho, 1, false);

    } else if (species_pred_type == predict_rhoX) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                MultiFab::Copy(sedge[lev][idim], sedge[lev][idim], FirstSpec,
                               Rho, 1, 0);
                for (int ispec = 1; ispec < NumSpec; ++ispec) {
                    MultiFab::Add(sedge[lev][idim], sedge[lev][idim],
                                  FirstSpec + ispec, Rho, 1, 0);
                }
            }
        }
    }

    if (species_pred_type == predict_rhoprime_and_X) {
        // convert rho' -> rho in scalold
        PutInPertForm(scalold, rho0_old, Rho, Rho, bcs_s, false);
    }

    if ((species_pred_type == predict_rhoprime_and_X) ||
        (species_pred_type == predict_rho_and_X)) {
        // convert X --> (rho X) in scalold
        ConvertRhoXToX(scalold, false);
    }

    /////////////////////////////////////////////////////////////////
    // Compute fluxes
    /////////////////////////////////////////////////////////////////

    if (which_step == 1) {
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_old(finest_level +
                                                                  1);

#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                AMREX_D_TERM(
                    rho0mac_old[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_old[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_old[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););
            }
            MakeS0mac(rho0_old, rho0mac_old);
        }
#endif

        // compute species fluxes
        MakeRhoXFlux(scalold, sflux, etarhoflux, sedge, umac, rho0_old,
                     rho0_edge_old, rho0mac_old, rho0_old, rho0_edge_old,
                     rho0mac_old, rho0_predicted_edge, FirstSpec, NumSpec);

    } else if (which_step == 2) {
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_old(finest_level +
                                                                  1);
        Vector<std::array<MultiFab, AMREX_SPACEDIM> > rho0mac_new(finest_level +
                                                                  1);

#if (AMREX_SPACEDIM == 3)
        if (spherical) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                AMREX_D_TERM(
                    rho0mac_old[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_old[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_old[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););

                AMREX_D_TERM(
                    rho0mac_new[lev][0].define(
                        convert(grids[lev], nodal_flag_x), dmap[lev], 1, 1);
                    , rho0mac_new[lev][1].define(
                          convert(grids[lev], nodal_flag_y), dmap[lev], 1, 1);
                    , rho0mac_new[lev][2].define(
                          convert(grids[lev], nodal_flag_z), dmap[lev], 1, 1););
            }
            MakeS0mac(rho0_old, rho0mac_old);
            MakeS0mac(rho0_new, rho0mac_new);
        }
#endif

        // compute species fluxes
        MakeRhoXFlux(scalold, sflux, etarhoflux, sedge, umac, rho0_old,
                     rho0_edge_old, rho0mac_old, rho0_new, rho0_edge_new,
                     rho0mac_new, rho0_predicted_edge, FirstSpec, NumSpec);
    }

    //**************************************************************************
    //     1) Set force for (rho X)_i at time n+1/2 = 0.
    //     2) Update (rho X)_i with conservative differencing.
    //     3) Define density as the sum of the (rho X)_i
    //     4) Update tracer with conservative differencing as well.
    //**************************************************************************

    // reaction forcing terms
    for (int lev = 0; lev <= finest_level; ++lev) {
        scal_force[lev].setVal(0.);
        MultiFab::Add(scal_force[lev], intra[lev], FirstSpec, FirstSpec,
                      NumSpec, 0);
    }

    Vector<MultiFab> p0_new_cart(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        p0_new_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_new, p0_new_cart, false, false, bcs_f, 0);

    // p0 only used in rhoh update so it's an optional parameter
    UpdateScal(scalold, scalnew, sflux, scal_force, FirstSpec, NumSpec,
               p0_new_cart);
}
