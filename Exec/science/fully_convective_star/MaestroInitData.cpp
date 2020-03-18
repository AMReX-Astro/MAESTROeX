
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel, 
                       const Real* s0_p, 
                       const Real* p0_p)
{
    Abort("Planar InitLevelData not implemented");
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
                       MultiFab& scal, MultiFab& vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);


    const int max_lev = max_radial_level + 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const auto tileBox = mfi.tilebox();

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);

        // set velocity to zero 
        AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
            vel_arr(i,j,k,n) = 0.0;
        });

        AMREX_PARALLEL_FOR_4D(tileBox, Nscal, i, j, k, n, {
            scal_arr(i,j,k,n) = 0.0;
        });

    }

    // if we are spherical, we want to make sure that p0 is good, since that is
    // what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    // then initialize h from rho, X, and p0.
    Vector<MultiFab> p0_cart(finest_level+1);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal 
    Vector<MultiFab> temp_mf(finest_level+1);

    for (auto i = 0; i <= finest_level; ++i) {
        p0_cart[i].define(grids[i], dmap[i], 1, 0);
        temp_mf[i].define(grids[i], dmap[i], Nscal, ng_s);

        temp_mf[i].setVal(0.0);
    }

    RealVector temp_vec(max_lev*nr_fine);

    // initialize temperature 
    for (auto l = 0; l < max_lev; ++l) {
        for (auto i = 0; i < nr_fine; ++i)
        temp_vec[l+max_lev*i] = s0_init[l+max_lev*(i+nr_fine*Temp)];
    }
    Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, Temp);
    // MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, ng_s);
    
    // initialize p0_cart
    Put1dArrayOnCart(p0_init, p0_cart, 0, 0, bcs_f, 0);

    // initialize species 
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto i = 0; i < max_lev*nr_fine; ++i) {
            temp_vec[i] = s0_init[i*(FirstSpec+comp)];
        }
        Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec+comp);
        // MultiFab::Copy(scal, temp_mf[lev], 0, FirstSpec+comp, 1, scal.nGrow());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const auto tileBox = mfi.tilebox();
        const int max_lev = max_radial_level + 1;

        const Array4<Real> vel_arr = vel.array(mfi);
        const Array4<Real> scal_arr = scal.array(mfi);
        const Array4<Real> temp_arr = temp_mf[lev].array(mfi);

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            scal_arr(i,j,k,Temp) = temp_arr(i,j,k,Temp);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i,j,k,FirstSpec+comp) = temp_arr(i,j,k,FirstSpec+comp);
            }
        });

    const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

//     // initialize rho as sum of partial densities rho*X_i
//     AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
//         for (auto comp = 0; comp < NumSpec; ++comp) {
//             scal_arr(i,j,k,Rho) += scal_arr(i,j,k,FirstSpec+comp);
//         }

//         // initialize (rho h) and T using the EOS
//         eos_t eos_state;
//         eos_state.T     = scal_arr(i,j,k,Temp);
//         eos_state.p     = p0_arr(i,j,k);
//         eos_state.rho   = scal_arr(i,j,k,Rho);
//         for (auto comp = 0; comp < NumSpec; ++comp) {
//             eos_state.xn[comp] = scal_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
//         }

//         eos(eos_input_rp, eos_state);

//         scal_arr(i,j,k,RhoH) = eos_state.rho*eos_state.h;
//         scal_arr(i,j,k,Temp) = eos_state.T;
//     });

//     if (perturb_model) {

//         // set the radial bounds of the perturbation
//         const auto velpert_r_outer = velpert_radius;
//         const auto velpert_r_inner = 0.0;

//         // random numbers between -1 and 1
//         GpuArray<Real,27> alpha;
//         GpuArray<Real,27> beta;
//         GpuArray<Real,27> gamma;

//         // random numbers between 0 and 2*pi
//         GpuArray<Real,27> phix;
//         GpuArray<Real,27> phiy;
//         GpuArray<Real,27> phiz;

//         for (auto i = 0; i < 27; ++i) {
//             alpha[i] = 2.0 * amrex::Random() - 1.0;
//             beta[i] = 2.0 * amrex::Random() - 1.0;
//             gamma[i] = 2.0 * amrex::Random() - 1.0;

//             phix[i] = 2.0 * M_PI * amrex::Random();
//             phiy[i] = 2.0 * M_PI * amrex::Random();
//             phiz[i] = 2.0 * M_PI * amrex::Random();
//         }

//         // compute the norm of k
//         GpuArray<Real,27> normk;
//         for (auto k = 1; k <= 3; ++k) {
//             for (auto j = 1; j <= 3; ++j) {
//                 for (auto i = 1; i <= 3; ++i) {
//                     normk[i + 3 * (j + 3 * k)] = std::sqrt(Real(i*i) + Real(j*j) + Real(k*k));
//                 }
//             }
//         }

//         const auto prob_lo = geom[lev].ProbLoArray();
//         const auto prob_hi = geom[lev].ProbHiArray();
//         const auto dx = geom[lev].CellSizeArray();

//         const Real * AMREX_RESTRICT center_p = center.dataPtr();

//         const auto velpert_scale_loc = velpert_scale;
//         const auto velpert_amplitude_loc = velpert_amplitude;
//         const auto velpert_steep_loc = velpert_steep;

//         // now do the big loop over all the points in the domain
//         AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
//             const Real x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - center_p[0];
//             const Real y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - center_p[1];
//             const Real z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - center_p[2];

//             // set perturbational velocity to zero
//             Real vpert[3];

//             for (auto n = 0; n < 3; ++n) {
//                 vpert[n] = 0.0;
//             }

//             // compute distance to the center of the star
//             Real rloc = std::sqrt(x*x + y*y + z*z);

//             // loop over the 27 combinations of Fourier components 
//             for (auto kk = 1; kk <= 3; ++kk) {
//                 for (auto jj = 1; jj <= 3; ++jj) {
//                     for (auto ii = 1; ii <= 3; ++ii) {
//                         // array index
//                         const int n = ii-1 + 3 * (jj-1 + 3 * (kk-1));
//                         // compute cosines and sines
//                         Real cx = std::cos(2.0 * M_PI * Real(ii) * x / velpert_scale_loc + phix[n]);
//                         Real cy = std::cos(2.0 * M_PI * Real(jj) * y / velpert_scale_loc + phiy[n]);
//                         Real cz = std::cos(2.0 * M_PI * Real(kk) * z / velpert_scale_loc + phiz[n]);

//                         Real sx = std::sin(2.0 * M_PI * Real(ii) * x / velpert_scale_loc + phix[n]);
//                         Real sy = std::sin(2.0 * M_PI * Real(jj) * y / velpert_scale_loc + phiy[n]);
//                         Real sz = std::sin(2.0 * M_PI * Real(kk) * z / velpert_scale_loc + phiz[n]);

//                         // compute contribution from perturbation velocity from each mode
//                         vpert[0] += (-gamma[n] * Real(jj) * cx*cz*sy + beta[n] * Real(kk) * cx*cy*sz) / normk[n];

//                         vpert[1] += (gamma[n] * Real(ii) * cy*cz*sx - alpha[n] * Real(kk) * cx*cy*sz) / normk[n];

//                         vpert[2] += (-beta[n] * Real(ii) * cy*cz*sx + alpha[n] * Real(jj) * cx*cz*sy) / normk[n];
//                     }
//                 }
//             }

//             // apply the cutoff function to the perturbational velocity
//             // add perturbational velocity to background velocity
//             for (auto n = 0; n < 3; ++n) {
//                 vpert[n] *= velpert_amplitude_loc * 0.5 * \
//                     (1.0 + std::tanh((velpert_r_outer - velpert_steep_loc - rloc) / velpert_steep_loc)) * 0.5 * \
//                     (1.0 + std::tanh((rloc - velpert_r_inner - velpert_steep_loc) / velpert_steep_loc));

//                 vel_arr(i,j,k,n) += vpert[n];
//             }
//         });
//     }
    }
}
