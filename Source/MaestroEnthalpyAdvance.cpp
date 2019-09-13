
#include <Maestro.H>

using namespace amrex;

void
Maestro::EnthalpyAdvance (int which_step,
                          Vector<MultiFab>& scalold,
                          Vector<MultiFab>& scalnew,
                          Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                          Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                          Vector<MultiFab>& scal_force,
                          Vector<std::array< MultiFab, AMREX_SPACEDIM > >& umac,
                          const Vector<std::array< MultiFab,AMREX_SPACEDIM > >& w0mac,
                          const Vector<MultiFab>& thermal)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::EnthalpyAdvance()",EnthalpyAdvance);

    // Create cell-centered base state quantity
    RealVector h0_old( (max_radial_level+1)*nr_fine );
    RealVector h0_new( (max_radial_level+1)*nr_fine );
    h0_old.shrink_to_fit();
    h0_new.shrink_to_fit();

    // Create edge-centered base state quantities.
    // Note: rho0_edge_{old,new} and rhoh0_edge_{old,new}
    // contain edge-centered quantities created via spatial interpolation.
    RealVector  rho0_edge_old( (max_radial_level+1)*(nr_fine+1) );
    RealVector  rho0_edge_new( (max_radial_level+1)*(nr_fine+1) );
    RealVector rhoh0_edge_old( (max_radial_level+1)*(nr_fine+1) );
    RealVector rhoh0_edge_new( (max_radial_level+1)*(nr_fine+1) );
    rho0_edge_old.shrink_to_fit();
    rho0_edge_new.shrink_to_fit();
    rhoh0_edge_old.shrink_to_fit();
    rhoh0_edge_new.shrink_to_fit();

    if (spherical == 0) {
        cell_to_edge( rho0_old.dataPtr(), rho0_edge_old.dataPtr());
        cell_to_edge( rho0_new.dataPtr(), rho0_edge_new.dataPtr());
        cell_to_edge(rhoh0_old.dataPtr(),rhoh0_edge_old.dataPtr());
        cell_to_edge(rhoh0_new.dataPtr(),rhoh0_edge_new.dataPtr());
    }

    if (enthalpy_pred_type == predict_h ||
        enthalpy_pred_type == predict_hprime) {
        // convert (rho h) -> h
        ConvertRhoHToH(scalold,true);
    }

    //////////////////////////////////
    // Create scalar source term at time n
    //////////////////////////////////

    for (int lev=0; lev<=finest_level; ++lev) {
        scal_force[lev].setVal(0.,RhoH,1,1);
    }

    Vector<MultiFab> rhoh0_old_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        rhoh0_old_cart[lev].define(grids[lev], dmap[lev], 1, 1);
        // needed to avoid NaNs in filling corner ghost cells with 2 physical boundaries
        rhoh0_old_cart[lev].setVal(0.);
    }

    // compute forcing terms
    if (enthalpy_pred_type == predict_rhohprime) {
        // make force for (rho h)'
        MakeRhoHForce(scal_force,1,thermal,umac,1,1);

        Put1dArrayOnCart(rhoh0_old,rhoh0_old_cart,0,0,bcs_s,RhoH);

        ModifyScalForce(scal_force,scalold,umac,rhoh0_old,rhoh0_edge_old,rhoh0_old_cart,RhoH,bcs_s,0);

    }
    else if (enthalpy_pred_type == predict_h ||
             enthalpy_pred_type == predict_rhoh) {
        // make force for (rho h)
        MakeRhoHForce(scal_force,1,thermal,umac,1,1);

        // make force for h by calling mkrhohforce then dividing by rho
        if (enthalpy_pred_type == predict_h) {
            for (int lev=0; lev<=finest_level; ++lev) {
                MultiFab::Divide(scal_force[lev],scalold[lev],RhoH,Rho,1,1);
            }
        }

    }
    else if (enthalpy_pred_type == predict_hprime) {
        // first compute h0_old
        // make force for hprime
        Abort("MaestroEnthalpyAdvance forcing");
    }
    else if (enthalpy_pred_type == predict_T_then_rhohprime ||
             enthalpy_pred_type == predict_T_then_h ||
             enthalpy_pred_type == predict_Tprime_then_h) {
        // make force for temperature
        Abort("MaestroEnthalpyAdvance forcing");
    }

    //////////////////////////////////
    // Add w0 to MAC velocities
    //////////////////////////////////

    Addw0(umac,w0mac,1.);

    //////////////////////////////////
    // Create the edge states of (rho h)' or h or T
    //////////////////////////////////

    if (enthalpy_pred_type == predict_rhohprime) {
        // convert (rho h) -> (rho h)'
        PutInPertForm(scalold, rhoh0_old, RhoH, 0, bcs_f, true);
    }

    if (enthalpy_pred_type == predict_hprime) {
        // convert h -> h'
        Abort("MaestroEnthalpyAdvance predict_hprime");
    }

    if (enthalpy_pred_type == predict_Tprime_then_h) {
        // convert T -> T'
        Abort("MaestroEnthalpyAdvance predict_Tprime_then_h");
    }

    // predict either T, h, or (rho h)' at the edges
    int pred_comp;
    if ( enthalpy_pred_type == predict_T_then_rhohprime ||
         enthalpy_pred_type == predict_T_then_h         ||
         enthalpy_pred_type == predict_Tprime_then_h)  {
        pred_comp = Temp;
    }
    else {
        pred_comp = RhoH;
    }

    if (enthalpy_pred_type == predict_rhoh) {
        // use the conservative form of the prediction
        MakeEdgeScal(scalold,sedge,umac,scal_force,0,bcs_s,Nscal,pred_comp,pred_comp,1,1);
    }
    else {
        // use the advective form of the prediction
        MakeEdgeScal(scalold,sedge,umac,scal_force,0,bcs_s,Nscal,pred_comp,pred_comp,1,0);
    }

    if (enthalpy_pred_type == predict_rhohprime) {
        // convert (rho h)' -> (rho h)
        PutInPertForm(scalold, rhoh0_old, RhoH, RhoH, bcs_s, false);
    }

    if (enthalpy_pred_type == predict_hprime) {
        // convert h' -> h
        Abort("MaestroEnthalpyAdavnce predict_hprime");
    }

    if (enthalpy_pred_type == predict_Tprime_then_h) {
        // convert T' -> T
        Abort("MaestroEnthalpyAdavnce predict_Tprime_then_h");
    }

    if (enthalpy_pred_type == predict_h ||
        enthalpy_pred_type == predict_hprime) {
        // convert (rho h) -> h
        ConvertRhoHToH(scalold,false);
    }

    // Compute enthalpy edge states if we were predicting temperature.  This
    // needs to be done after the state was returned to the full state.
    if ( (enthalpy_pred_type == predict_T_then_rhohprime) ||
         (enthalpy_pred_type == predict_T_then_h        ) ||
         (enthalpy_pred_type == predict_Tprime_then_h) ) {
        Abort("MaestroEnthalpyAdvance need makeHfromRhoT_edge");
    }

    //////////////////////////////////
    // Subtract w0 from MAC velocities
    //////////////////////////////////

    Addw0(umac,w0mac,-1.);

    //////////////////////////////////
    // Compute fluxes
    //////////////////////////////////

    // for which_step .eq. 1, we pass in only the old base state quantities
    // for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step == 1) {
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >  rho0mac_old(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > > rhoh0mac_old(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >    h0mac_old(finest_level+1);

        if (spherical == 1) {
            for (int i=0; i<h0_old.size(); ++i) {
                h0_old[i] = rhoh0_old[i] / rho0_old[i];
            }

            for (int lev=0; lev<=finest_level; ++lev) {
                AMREX_D_TERM(rho0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rho0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rho0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(rhoh0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rhoh0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rhoh0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(h0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             h0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             h0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
            }

            MakeS0mac(rho0_old,rho0mac_old);
            MakeS0mac(rhoh0_old,rhoh0mac_old);
            MakeS0mac(h0_old,h0mac_old);
        }

        // compute enthalpy fluxes
        MakeRhoHFlux(scalold, sflux, sedge, umac, w0mac,
                     rho0_old,rho0_edge_old,rho0mac_old,
                     rho0_old,rho0_edge_old,rho0mac_old,
                     rhoh0_old,rhoh0_edge_old,rhoh0mac_old,
                     rhoh0_old,rhoh0_edge_old,rhoh0mac_old,
                     h0mac_old,h0mac_old);
    }
    else if (which_step == 2) {
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >  rho0mac_old(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > > rhoh0mac_old(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >    h0mac_old(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >  rho0mac_new(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > > rhoh0mac_new(finest_level+1);
        Vector< std::array< MultiFab,AMREX_SPACEDIM > >    h0mac_new(finest_level+1);

        if (spherical == 1) {
            for (int i=0; i<h0_old.size(); ++i) {
                h0_old[i] = rhoh0_old[i] / rho0_old[i];
                h0_new[i] = rhoh0_new[i] / rho0_new[i];
            }

            for (int lev=0; lev<=finest_level; ++lev) {
                AMREX_D_TERM(rho0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rho0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rho0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(rhoh0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rhoh0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rhoh0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(h0mac_old[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             h0mac_old[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             h0mac_old[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(rho0mac_new[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rho0mac_new[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rho0mac_new[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(rhoh0mac_new[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             rhoh0mac_new[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             rhoh0mac_new[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );

                AMREX_D_TERM(h0mac_new[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
                             h0mac_new[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
                             h0mac_new[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
            }

            MakeS0mac(rho0_old,rho0mac_old);
            MakeS0mac(rhoh0_old,rhoh0mac_old);
            MakeS0mac(h0_old,h0mac_old);
            MakeS0mac(rho0_new,rho0mac_new);
            MakeS0mac(rhoh0_new,rhoh0mac_new);
            MakeS0mac(h0_new,h0mac_new);
        }

        // compute enthalpy fluxes
        MakeRhoHFlux(scalold, sflux, sedge, umac, w0mac,
                     rho0_old,rho0_edge_old,rho0mac_old,
                     rho0_new,rho0_edge_new,rho0mac_new,
                     rhoh0_old,rhoh0_edge_old,rhoh0mac_old,
                     rhoh0_new,rhoh0_edge_new,rhoh0mac_new,
                     h0mac_old,h0mac_new);
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        scal_force[lev].setVal(0.,RhoH,1,1);
    }

    //**************************************************************************
    //     1) Create (rho h)' force at time n+1/2.
    //          (NOTE: we don't worry about filling ghost cells of the scal_force
    //                 because we only need them in valid regions...)
    //     2) Update (rho h) with conservative differencing.
    //**************************************************************************

    MakeRhoHForce(scal_force,0,thermal,umac,0,which_step);

    Vector<MultiFab> p0_new_cart(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        p0_new_cart[lev].define(grids[lev], dmap[lev], 1, 1);
    }

    Put1dArrayOnCart(p0_new,p0_new_cart,0,0,bcs_f,0);

    UpdateScal(scalold, scalnew, sflux, scal_force, RhoH, 1, p0_new_cart);
}
