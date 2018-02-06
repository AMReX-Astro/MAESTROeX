
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
                          const Vector<MultiFab>& thermal)
{
    // Create edge-centered base state quantities.
    // Note: rho0_edge_{old,new} and rhoh0_edge_{old,new}
    // contain edge-centered quantities created via spatial interpolation.
    Vector<Real>  rho0_edge_old( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real>  rho0_edge_new( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real> rhoh0_edge_old( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real> rhoh0_edge_new( (max_radial_level+1)*(nr_fine+1) );
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
        ConvertRhoHToH(sold,true);
    }

    //////////////////////////////////
    // Create scalar source term at time n
    //////////////////////////////////

    for (int lev=0; lev<=finest_level; ++lev) {
        scal_force[lev].setVal(0.,RhoH,1,1);
    }

    // compute forcing terms    
    if (enthalpy_pred_type == predict_rhohprime) {
        // make force for (rho h)'
        MakeRhoHForce(scal_force,1,thermal,umac,1);

        ModifyScalForce(scal_force,umac,rhoh0_old,rhoh0_edge_old,RhoH,bcs_s,0);

    }
    else if (enthalpy_pred_type == predict_h ||
             enthalpy_pred_type == predict_rhoh) {
        // make force for h by calling mkrhohforce then dividing by rho
        Abort("MaestroEnthalpyAdvance forcing");
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

    Addw0(umac,1.);

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
        Abort("MaestroEnthalpyAdavnce predict_rhoh");
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
        ConvertRhoHToH(sold,false);
    }

    // Compute enthalpy edge states if we were predicting temperature.  This
    // needs to be done after the state was returned to the full state.
    if ( (enthalpy_pred_type == predict_T_then_rhohprime) ||
         (enthalpy_pred_type == predict_T_then_h        ) ||
         (enthalpy_pred_type == predict_Tprime_then_h) ) {
        Abort("MaestroEnthalpyAdavnce need makeHfromRhoT_edge");
    }
   
    //////////////////////////////////
    // Subtract w0 from MAC velocities
    //////////////////////////////////

    Addw0(umac,-1.);
   
    //////////////////////////////////
    // Compute fluxes
    //////////////////////////////////

    // for which_step .eq. 1, we pass in only the old base state quantities
    // for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step == 1) {

        if (spherical == 1) {
        }

        // compute enthalpy fluxes
	MakeRhoHFlux(scalold, sflux, sedge, umac,
		     rho0_old,rho0_edge_old,rho0_old,rho0_edge_old,
		     rhoh0_old,rhoh0_edge_old,rhoh0_old,rhoh0_edge_old);
    }
    else if (which_step == 2) {

        if (spherical == 1) {
        }

        // compute enthalpy fluxes
	MakeRhoHFlux(scalold, sflux, sedge, umac,
		     rho0_old,rho0_edge_old,rho0_new,rho0_edge_new,
		     rhoh0_old,rhoh0_edge_old,rhoh0_new,rhoh0_edge_new);
    }

    for (int lev=0; lev<=finest_level; ++lev) {
        scal_force[lev].setVal(0.,RhoH,1,1);
    }

    MakeRhoHForce(scal_force,0,thermal,umac,0);

    if (spherical == 1) {
    }

    /*
    call update_scal(mla,rhoh_comp,rhoh_comp,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level)
    */
}
