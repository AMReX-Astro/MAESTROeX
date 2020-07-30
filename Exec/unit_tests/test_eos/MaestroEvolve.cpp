
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    Print() << "Running tests..." << std::endl;

    Vector<MultiFab> error(finest_level + 1);

    // number of EoS tests that do_tests performs
    const int n_tests = 5;

    // setup error MultiFab
    for (int lev = 0; lev <= finest_level; ++lev) {
        error[lev].define(grids[lev], dmap[lev], n_tests, ng_s);
        error[lev].setVal(0.);
    }

    // do the tests
    for (int lev = 0; lev <= finest_level; ++lev) {
        const Box& domainBox = geom[lev].Domain();

        const auto dom_lo = domainBox.loVect3d();
        const auto dom_hi = domainBox.hiVect3d();

        const auto ih1 = network_spec_index("hydrogen-1");
        const auto ihe4 = network_spec_index("helium-4");

        const auto dlogrho = (std::log10(dens_max) - std::log10(dens_min)) /
                             (dom_hi[0] - dom_lo[0]);
        const auto dlogT = (std::log10(temp_max) - std::log10(temp_min)) /
                           (dom_hi[1] - dom_lo[1]);
        const auto dmetal = metalicity_max / (dom_hi[2] - dom_lo[2]);

        const auto temp_min_l = temp_min;
        const auto dens_min_l = dens_min;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sold[lev], true); mfi.isValid(); ++mfi) {
            const Box& tileBox = mfi.tilebox();

            const Array4<Real> scal = sold[lev].array(mfi);
            const Array4<Real> error_arr = error[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // set the composition -- approximately solar
                const auto metalicity = Real(k) * dmetal;

                Real xn_zone[NumSpec];

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    xn_zone[comp] =
                        metalicity / (NumSpec - 2);  // all but H, He
                }

                xn_zone[ih1] = 0.75 - 0.5 * metalicity;
                xn_zone[ihe4] = 0.25 - 0.5 * metalicity;

                // set the temperature
                const auto temp_zone =
                    std::pow(10.0, std::log10(temp_min_l) + Real(j) * dlogT);

                // set the density
                const auto dens_zone =
                    std::pow(10.0, std::log10(dens_min_l) + Real(i) * dlogrho);

                eos_t eos_state;

                // call the EOS -- rho, T, X directly
                // input: eos_input_rt
                eos_state.T = temp_zone;
                eos_state.rho = dens_zone;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                eos(eos_input_rt, eos_state);

                // store the thermodynamic state
                scal(i, j, k, Rho) = dens_zone;
                scal(i, j, k, RhoH) = dens_zone * eos_state.h;
                scal(i, j, k, Temp) = temp_zone;

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    scal(i, j, k, FirstSpec + comp) = xn_zone[comp];
                }

                auto h = eos_state.h;
                auto p = eos_state.p;
                auto e = eos_state.e;
                auto s = eos_state.s;

                // call the EOS with rho, h
                // input: eos_input_rh

                // change initial T guess to make the root find do some work
                eos_state.T = 0.5 * temp_zone;
                eos_state.rho = dens_zone;
                eos_state.h = h;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                eos(eos_input_rh, eos_state);

                // store the thermodynamic state
                const auto tfromrh = eos_state.T;
                error_arr(i, j, k, 0) += (eos_state.T - temp_zone) / temp_zone;

                // call the EOS with T, p
                // input: eos_input_tp

                // change initial rho guess to make the root find do some work
                eos_state.T = temp_zone;
                eos_state.rho = dens_zone / 3.0;
                eos_state.p = p;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                eos(eos_input_tp, eos_state);

                // store the thermodynamic state
                const auto rfromtp = eos_state.rho;
                error_arr(i, j, k, 1) +=
                    (eos_state.rho - dens_zone) / dens_zone;

                // call the EOS with rho, p
                // input: eos_input_rp

                // change initial T guess to make the root find do some work
                eos_state.T = 0.5 * temp_zone;
                eos_state.rho = dens_zone;
                eos_state.p = p;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                eos(eos_input_rp, eos_state);

                // store the thermodynamic state
                const auto tfromrp = eos_state.T;
                error_arr(i, j, k, 2) += (eos_state.T - temp_zone) / temp_zone;

                // call the EOS with rho, e
                // input: eos_input_re

                // change initial T guess to make the root find do some work
                eos_state.T = 0.5 * temp_zone;
                eos_state.rho = dens_zone;
                eos_state.e = e;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                eos(eos_input_re, eos_state);

                // store the thermodynamic state
                const auto tfromre = eos_state.T;
                error_arr(i, j, k, 3) += (eos_state.T - temp_zone) / temp_zone;

                // call the EOS with p, s
                // input: eos_input_ps

                // change initial T guess to make the root find do some work
                eos_state.T = 0.5 * temp_zone;
                eos_state.rho = 0.5 * dens_zone;
                eos_state.p = p;
                eos_state.s = s;
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] = xn_zone[comp];
                }

                // some EOSes don't have physically valid treatments
                // of entropy throughout the entire rho-T plane
                if (eos_state.s > 0.0) {
                    eos(eos_input_ps, eos_state);

                    // store the thermodynamic state
                    const auto tfromps = eos_state.T;
                    error_arr(i, j, k, 4) +=
                        (eos_state.T - temp_zone) / temp_zone;
                }
            });

            // const int* lo  = tilebox.loVect();
            // const int* hi  = tilebox.hiVect();

            // do_tests(ARLIM_3D(lo), ARLIM_3D(hi),
            //          BL_TO_FORTRAN_3D(scal[mfi]),
            //          BL_TO_FORTRAN_FAB(error_mf[mfi]),
            //          domainBox.loVect(), domainBox.hiVect());
        }
    }

    // print the errors
    {
        int lev = finest_level;

        Print() << "\nError in T   from rho, h = " << error[lev].norm2(0)
                << std::endl;
        Print() << "Error in rho from T  , p = " << error[lev].norm2(1)
                << std::endl;
        Print() << "Error in T   from rho, p = " << error[lev].norm2(2)
                << std::endl;
        Print() << "Error in T   from rho, e = " << error[lev].norm2(3)
                << std::endl;
        Print() << "Error in T   from p  , s = " << error[lev].norm2(4)
                << std::endl;
    }
}
