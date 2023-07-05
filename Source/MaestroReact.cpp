
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext, then
// react the state over dt_in and update rho_omegadot, rho_Hnuc
void Maestro::React(const Vector<MultiFab>& s_in, Vector<MultiFab>& s_out,
                    Vector<MultiFab>& rho_Hext, Vector<MultiFab>& rho_omegadot,
                    Vector<MultiFab>& rho_Hnuc, const BaseState<Real>& p0,
                    const Real dt_in, [[maybe_unused]] const Real time_in) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::React()", React);

    // external heating
    if (do_heating) {
        // computing heating term
        MakeHeating(rho_Hext, s_in);

        // if we aren't burning, then we should just copy the old state to the
        // new and only update the rhoh component with the heating term
        if (!do_burning) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                // copy s_in to s_out
                MultiFab::Copy(s_out[lev], s_in[lev], 0, 0, Nscal, 0);

                // add in the heating term, s_out += dt_in * rho_Hext
                MultiFab::Saxpy(s_out[lev], dt_in, rho_Hext[lev], 0, RhoH, 1,
                                0);
            }
        }
    } else {
        // not heating, so we zero rho_Hext
        for (int lev = 0; lev <= finest_level; ++lev) {
            rho_Hext[lev].setVal(0.);
        }
    }

    // apply burning term
    if (do_burning) {
        // do the burning, update rho_omegadot and rho_Hnuc
        // we pass in rho_Hext so that we can add it to rhoh in case we applied heating
        Burner(s_in, s_out, rho_Hext, rho_omegadot, rho_Hnuc, p0, dt_in,
               time_in);
        // pass temperature through for seeding the temperature update eos call
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(s_out[lev], s_in[lev], Temp, Temp, 1, 0);
        }
    } else {
        // not burning, so we zero rho_omegadot and rho_Hnuc
        for (int lev = 0; lev <= finest_level; ++lev) {
            rho_omegadot[lev].setVal(0.);
            rho_Hnuc[lev].setVal(0.);
        }
    }

    // if we aren't doing any heating/burning, then just copy the old to the new
    if (!do_heating && !do_burning) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(s_out[lev], s_in[lev], 0, 0, Nscal, 0);
        }
    }

    // average down and fill ghost cells
    AverageDown(s_out, 0, Nscal);
    FillPatch(t_old, s_out, s_out, s_out, 0, 0, Nscal, 0, bcs_s);

    // average down (no ghost cells)
    AverageDown(rho_Hext, 0, 1);
    AverageDown(rho_omegadot, 0, NumSpec);
    AverageDown(rho_Hnuc, 0, 1);

    // now update temperature
    if (use_tfromp) {
        TfromRhoP(s_out, p0);
    } else {
        TfromRhoH(s_out, p0);
    }
}
