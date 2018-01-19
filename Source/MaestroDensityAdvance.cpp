
#include <Maestro.H>

using namespace amrex;

void
Maestro::DensityAdvance (bool is_predictor,
                         Vector<MultiFab>& scalold,
                         Vector<MultiFab>& scalnew,
                         Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sedge,
                         Vector<std::array< MultiFab, AMREX_SPACEDIM > >& sflux,
                         Vector<MultiFab>& scal_force)
{
    Vector<Real> rho0_edge_old( (max_radial_level+1)*(nr_fine+1) );
    Vector<Real> rho0_edge_new( (max_radial_level+1)*(nr_fine+1) );

    if (spherical == 0) {
        // create edge-centered base state quantities.
        // Note: rho0_edge_{old,new} 
        // contains edge-centered quantities created via spatial interpolation.
        // This is to be contrasted to rho0_predicted_edge which is the half-time
        // edge state created in advect_base. 
        cell_to_edge(rho0_old.dataPtr(),rho0_edge_old.dataPtr());
        cell_to_edge(rho0_new.dataPtr(),rho0_edge_new.dataPtr());
    }



    //////////////////////////////////
    // Create source terms at time n
    //////////////////////////////////

    // source terms for X and for tracers are zero - do nothing
    for (int lev=0; lev<=finest_level; ++lev) {
        scal_force[lev].setVal(0.);
    }

    if (spherical == 1) {

    }

    // ** density source term **

    // Make source term for rho or rho' 
    if (species_pred_type == predict_rhoprime_and_X) {
        // rho' souce term
        // this is needed for pred_rhoprime_and_X


    }

    else if (species_pred_type == predict_rho_and_X) {
        // rho source term

    }



    // ** species source term **

    // for species_pred_types predict_rhoprime_and_X and
    // predict_rho_and_X, there is no force for X.

    // for predict_rhoX, we are predicting (rho X)
    // as a conservative equation, and there is no force.

}
