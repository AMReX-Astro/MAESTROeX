
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext, then
// react the state over dt_in and update rho_omegadot, rho_Hnuc
void
Maestro::React (const Vector<MultiFab>& s_in,
                Vector<MultiFab>& s_out,
                Vector<MultiFab>& rho_Hext,
                Vector<MultiFab>& rho_omegadot,
                Vector<MultiFab>& rho_Hnuc,
                const Vector<Real>& p0,
                const Real dt_in)
{

    // external heating
    if (do_heating) {

        // computing heating term
        ComputeHeating(rho_Hext);

        // if we aren't burning, then we should just copy the old state to the
        // new and only update the rhoh component with the heating term
        if (!do_burning) {
            for (int lev=0; lev<=finest_level; ++lev) {
                // copy s_in to s_out
                MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);

                // add in the heating term, s_out += dt_in * rho_Hext
                MultiFab::Saxpy(s_out[lev],dt_in,rho_Hext[lev],0,RhoH,1,0);
            }
        }
    }
    else {
        // not heating, so we zero rho_Hext
        for (int lev=0; lev<=finest_level; ++lev) {
            rho_Hext[lev].setVal(0.);
        }        
    }

    // apply burning term
    if (do_burning) {

        // do the burning, update rho_omegadot and rho_Hnuc
        // we pass in rho_Hext so that we can add it to rhoh in case we applied heating
        Burner(s_in,s_out,rho_Hext,rho_omegadot,rho_Hnuc,p0,dt_in);

        // pass temperature through for seeding the temperature update eos call
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],Temp,Temp,1,0);
        }
    }
    else {
        // not burning, so we zero rho_omegadot and rho_Hnuc
        for (int lev=0; lev<=finest_level; ++lev) {
            rho_omegadot[lev].setVal(0.);
            rho_Hnuc[lev].setVal(0.);
        }        

    }
    

    // if we aren't doing any heating/burning, then just copy the old to the new
    if (!do_heating && !do_burning) {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,0);
        }
    }


    // average fine data onto coarser cells
    AverageDown(s_out);
    AverageDown(rho_Hext);
    AverageDown(rho_omegadot);
    AverageDown(rho_Hnuc);

/* FIXME
      ! now update temperature
      if (use_tfromp) then
         call makeTfromRhoP(snew,p0,mla,the_bc_level,dx)
      else
         call makeTfromRhoH(snew,p0,mla,the_bc_level,dx)
      endif
*/



}

void Maestro::Burner(const Vector<MultiFab>& s_in,
                     Vector<MultiFab>& s_out,
                     const Vector<MultiFab>& rho_Hext,
                     Vector<MultiFab>& rho_omegadot,
                     Vector<MultiFab>& rho_Hnuc,
                     const Vector<Real>& p0,
                     const Real dt_in) {

    // get references to the MultiFabs at level lev
    for (int lev=0; lev<=finest_level; ++lev) {

        const MultiFab& s_in_mf = s_in[lev];
        MultiFab& s_out_mf = s_out[lev];
        const MultiFab& rho_Hext_mf = rho_Hext[lev];
        MultiFab& rho_omegadot_mf = rho_omegadot[lev];
        MultiFab& rho_Hnuc_mf = rho_Hnuc[lev];

        // loop over boxes
        for ( MFIter mfi(s_in_mf); mfi.isValid(); ++mfi ) {

            // get references to the FABs, each containing data and the valid+ghost box
            const FArrayBox& s_in_fab = s_in_mf[mfi];
            FArrayBox& s_out_fab = s_out_mf[mfi];
            const FArrayBox& rho_Hext_fab = rho_Hext_mf[mfi];
            FArrayBox& rho_omegadot_fab = rho_omegadot_mf[mfi];
            FArrayBox& rho_Hnuc_fab = rho_Hnuc_mf[mfi];

            // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

            // Get the index space of the valid+ghost region for each FAB
            // Note each of these boxes may contain ghost cells, and thus are
            // larger than or equal to mfi.validbox().
            const Box & s_in_box = s_in_fab.box();
            const Box & s_out_box = s_out_fab.box();
            const Box & rho_Hext_box = rho_Hext_fab.box();
            const Box & rho_omegadot_box = rho_omegadot_fab.box();
            const Box & rho_Hnuc_box = rho_Hnuc_fab.box();

            // We can now pass the information to a Fortran routine,
            // e.g. s_in_fab.dataPtr() gives a double*, which is reshaped into
            // a multi-dimensional array with dimensions specified by
            // the information in "s_in_box". We will also pass "box",
            // which specifies our "work" region .
            burner_loop(lev,ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                        s_in_fab.dataPtr(),
                        ARLIM_3D(s_in_box.loVect()), ARLIM_3D(s_in_box.hiVect()),
                        s_in_fab.nComp(),
                        s_out_fab.dataPtr(),
                        ARLIM_3D(s_out_box.loVect()), ARLIM_3D(s_out_box.hiVect()),
                        s_out_fab.nComp(),
                        rho_Hext_fab.dataPtr(),
                        ARLIM_3D(rho_Hext_box.loVect()), ARLIM_3D(rho_Hext_box.hiVect()),
                        rho_Hext_fab.nComp(),
                        rho_omegadot_fab.dataPtr(),
                        ARLIM_3D(rho_omegadot_box.loVect()), ARLIM_3D(rho_omegadot_box.hiVect()),
                        rho_omegadot_fab.nComp(),
                        rho_Hnuc_fab.dataPtr(),
                        ARLIM_3D(rho_Hnuc_box.loVect()), ARLIM_3D(rho_Hnuc_box.hiVect()),
                        rho_Hnuc_fab.nComp(),
                        tempbar_init.dataPtr(), dt);



        }


    }




}

// compute heating term, rho_Hext
void
Maestro::ComputeHeating (Vector<MultiFab>& rho_Hext) {

    // FIXME
    for (int lev=0; lev<=finest_level; ++lev) {
        rho_Hext[lev].setVal(0.);
    }        

}
