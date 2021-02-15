
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

// compute cp and xi
void Maestro::MakeIntraCoeffs(const Vector<MultiFab>& scal1,
                              const Vector<MultiFab>& scal2,
                              Vector<MultiFab>& cp, Vector<MultiFab>& xi) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeIntraCoeffs()", MakeIntraCoeffs);

    for (int lev = 0; lev <= finest_level; ++lev) {
        Print() << "... Level " << lev << " create intra coeffs:" << std::endl;

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(scal1[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& gtbx = mfi.tilebox();

            const Array4<const Real> scalold = scal1[lev].array(mfi);
            const Array4<const Real> scalnew = scal2[lev].array(mfi);
            const Array4<Real> cp_arr = cp[lev].array(mfi);
            const Array4<Real> xi_arr = xi[lev].array(mfi);

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // old state first
                eos_t eos_state;

                eos_state.rho = scalold(i, j, k, Rho);
                eos_state.T = scalold(i, j, k, Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] =
                        scalold(i, j, k, FirstSpec + comp) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto comp = 0; comp < NumAux; ++comp) {
                    eos_state.aux[comp] =
                        scalold(i, j, k, FirstAux + comp) / eos_state.rho;
                }
#endif

                // dens, temp, and xmass are inputs
                eos(eos_input_rt, eos_state);

                eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

                cp_arr(i, j, k) = eos_state.cp;

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    xi_arr(i, j, k, comp) = eos_xderivs.dhdX[comp];
                }

                // new state now -- average results
                eos_state.rho = scalnew(i, j, k, Rho);
                eos_state.T = scalnew(i, j, k, Temp);
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    eos_state.xn[comp] =
                        scalnew(i, j, k, FirstSpec + comp) / eos_state.rho;
                }
#if NAUX_NET > 0
                for (auto comp = 0; comp < NumAux; ++comp) {
                    eos_state.aux[comp] =
                        scalnew(i, j, k, FirstAux + comp) / eos_state.rho;
                }
#endif

                // dens, temp, and xmass are inputs
                eos(eos_input_rt, eos_state);

                eos_xderivs = composition_derivatives(eos_state);

                cp_arr(i, j, k) = 0.5 * (eos_state.cp + cp_arr(i, j, k));

                for (auto comp = 0; comp < NumSpec; ++comp) {
                    xi_arr(i, j, k, comp) =
                        0.5 * (eos_xderivs.dhdX[comp] + xi_arr(i, j, k, comp));
                }
            });
        }
    }

    // average down
    AverageDown(cp, 0, 1);
    AverageDown(xi, 0, 1);
}
