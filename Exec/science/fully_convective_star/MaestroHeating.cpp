
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    static bool model_init = true;
    if (model_init) {
        const int model_file_length = model_file.length();
        Vector<int> model_file_name(model_file_length);
        for (int i = 0; i < model_file_length; i++)
            model_file_name[i] = model_file[i];
        ca_read_model_file(model_file_name.dataPtr(), &model_file_length);

        model_init = false;
    }

    const auto& center_p = center;

    for (int lev = 0; lev <= finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        MultiFab& rho_Hext_mf = rho_Hext[lev];
        const MultiFab& scal_mf = scal[lev];

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
            const auto dx = geom[lev].CellSizeArray();
            const auto prob_lo = geom[lev].ProbLoArray();

            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real starting_rad =
                    (spherical) ? 0.0 : prob_lo[AMREX_SPACEDIM - 1];

                Real xloc[3];
                xloc[0] = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
                xloc[1] = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
                xloc[2] = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

                Real rloc = (spherical) ? std::sqrt(xloc[0] * xloc[0] +
                                                    xloc[1] * xloc[1] +
                                                    xloc[2] * xloc[2])
                                        : xloc[AMREX_SPACEDIM - 1];

                Real T9 = scal_arr(i, j, k, Temp) / 1.e9;
                Real rho = scal_arr(i, j, k, Rho);
                rho_Hext_arr(i, j, k) =
                    rho * 2.4e4 * rho / std::pow(T9, 2.0 / 3.0) *
                    std::exp(-3.38 / std::pow(T9, 1.0 / 3.0));
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
