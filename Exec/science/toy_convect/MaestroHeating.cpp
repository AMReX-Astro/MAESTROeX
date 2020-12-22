
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    static int ih1, ic12, in14, io16;
    static bool firstCall = true;

    if (firstCall) {
        ih1 = network_spec_index("hydrogen-1");
        ic12 = network_spec_index("carbon-12");
        in14 = network_spec_index("nitrogen-14");
        io16 = network_spec_index("oxygen-16");

        firstCall = false;
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(rho_Hext[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();

            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real rho = scal_arr(i, j, k, Rho);
                Real temp;
                // if (drive_initial_convection)
                //    temp = tempbar_init(j)
                // else
                temp = scal_arr(i, j, k, Temp);
                // endif

                Real T6 = temp / 1.e6;
                Real T613 = std::pow(T6, 1.0 / 3.0);

                // total CNO abundance
                Real X_CNO = (scal_arr(i, j, k, FirstSpec + ic12) +
                              scal_arr(i, j, k, FirstSpec + in14) +
                              scal_arr(i, j, k, FirstSpec + io16)) /
                             rho;

                // H abundance
                Real X_1 = scal_arr(i, j, k, FirstSpec + ih1) / rho;

                // CNO heating from Kippenhahn & Weigert, Eq. 18.65
                Real g14 =
                    1.0 + 2.7e-3 * T613 - 7.78e-3 * T613 * T613 - 1.49e-4 * T6;
                Real eps_CNO = 8.67e27 * g14 * X_CNO * X_1 * rho *
                               exp(-152.28e0 / T613) / (T613 * T613);

                rho_Hext_arr(i, j, k) = scal_arr(i, j, k, Rho) * eps_CNO;
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
