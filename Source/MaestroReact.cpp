
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

        // we pass in rho_Hext so that we can add it to rhoh in case we applied heating
        // FIXME
        // BurnerLoop();

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
            MultiFab::Copy(s_out[lev],s_in[lev],0,0,Nscal,s_out[lev].nGrow());
        }
    }


    // average fine data onto coarser cells
    AverageDown(rho_Hext);
    AverageDown(rho_omegadot);
    AverageDown(rho_Hnuc);


/*

      ! restrict data and fill all ghost cells
      call ml_restrict_and_fill(nlevs,snew,mla%mba%rr,the_bc_level, &
                                icomp=rho_comp, &
                                bcomp=dm+rho_comp, &
                                nc=nscal, &
                                ng=snew(1)%ng)

      ! the loop over nlevs must count backwards to make sure the finer grids are done first
      do n=nlevs,2,-1
         ! set level n-1 data to be the average of the level n data covering it
         call ml_cc_restriction(rho_omegadot(n-1),rho_omegadot(n),mla%mba%rr(n-1,:))
         call ml_cc_restriction(rho_Hext(n-1)    ,rho_Hext(n)    ,mla%mba%rr(n-1,:))
         call ml_cc_restriction(rho_Hnuc(n-1)    ,rho_Hnuc(n)    ,mla%mba%rr(n-1,:))
      enddo

      ! now update temperature
      if (use_tfromp) then
         call makeTfromRhoP(snew,p0,mla,the_bc_level,dx)
      else
         call makeTfromRhoH(snew,p0,mla,the_bc_level,dx)
      endif
*/



}

// compute heating term, rho_Hext
void
Maestro::ComputeHeating (Vector<MultiFab>& rho_Hext) {

    // FIXME
    for (int lev=0; lev<=finest_level; ++lev) {
        rho_Hext[lev].setVal(0.);
    }        

}
