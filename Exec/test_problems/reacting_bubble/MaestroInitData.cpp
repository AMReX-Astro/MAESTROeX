
#include <Maestro.H>
using namespace amrex;

// prototype for pertubation function to be called on the 
// device (if USE_CUDA=TRUE)
AMREX_GPU_DEVICE
void Perturb(const Real p0, 
             const Real* s0,
             Real* perturbations,  
             const Real x, const Real y, const Real z,
             const Real pert_rad_factor,
             const Real pert_temp_factor,
             const bool do_small_domain,
             bool spherical=false);

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, 
                       const Array4<Real> scal, 
                       const Array4<Real> vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();
    const int max_lev = base_geom.max_radial_level + 1;
    const int nrf = base_geom.nr_fine;

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

    const Real * AMREX_RESTRICT s0_p = s0_init.dataPtr();
    const auto& p0_p = p0_init;

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i,j,k,Rho) = s0_p[lev+max_lev*(r+nrf*Rho)];
        scal(i,j,k,RhoH) = s0_p[lev+max_lev*(r+nrf*RhoH)];
        scal(i,j,k,Temp) = s0_p[lev+max_lev*(r+nrf*Temp)];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = s0_p[lev+max_lev*(r+nrf*(FirstSpec+comp))];
        }
        // initialize pi to zero for now
        scal(i,j,k,Pi) = 0.0;
    });    
    
    // add an optional perturbation
    if (perturb_model) {

        const auto prob_lo = geom[lev].ProbLoArray();
        const auto dx = geom[lev].CellSizeArray();

        const auto pert_rad_factor_loc = pert_rad_factor;
        const auto pert_temp_factor_loc = pert_temp_factor;
        const auto do_small_domain_loc = do_small_domain;

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0];
            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1];
            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2];

            Real perturbations[Nscal];
            Real s0[Nscal];

            for (auto n = 0; n < Nscal; ++n) {
                s0[n] = s0_p[lev+max_lev*(r+nrf*n)];
            }

            Perturb(p0_p(lev,r), s0, perturbations, 
                    x, y, z, 
                    pert_rad_factor_loc, 
                    pert_temp_factor_loc, 
                    do_small_domain_loc);

            scal(i,j,k,Rho) = perturbations[Rho];
            scal(i,j,k,RhoH) = perturbations[RhoH];
            scal(i,j,k,Temp) = perturbations[Temp];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i,j,k,FirstSpec+comp) = perturbations[FirstSpec+comp];
            }
        });
    }
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
                       MultiFab& scal, MultiFab& vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelDataSphr()", InitLevelDataSphr);
    const int max_lev = base_geom.max_radial_level + 1;
    
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
    Vector<MultiFab> p0_cart(finest_level);

    // make a temporary MultiFab and RealVector to hold the cartesian data then copy it back to scal 
    Vector<MultiFab> temp_mf(finest_level);

    for (auto i = 0; i <= finest_level; ++i) {
        p0_cart[i].define(grids[lev], dmap[lev], 1, ng_s);
        temp_mf[i].define(grids[lev], dmap[lev], 1, ng_s);
    }

    RealVector temp_vec(max_lev*base_geom.nr_fine);

    // initialize temperature 
    for (auto i = 0; i < max_lev*base_geom.nr_fine; ++i) {
        temp_vec[i] = s0_init[i*Temp];
    }
    Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, Temp);
    MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());
    
    // initialize p0_cart
    Put1dArrayOnCart(p0_init, p0_cart, 0, 0, bcs_f, 0);

    // initialize species 
    for (auto comp = 0; comp < NumSpec; ++comp) {
        for (auto i = 0; i < max_lev*base_geom.nr_fine; ++i) {
            temp_vec[i] = s0_init[i*(FirstSpec+comp)];
        }
        Put1dArrayOnCart(temp_vec, temp_mf, 0, 0, bcs_s, FirstSpec+comp);
        MultiFab::Copy(scal, temp_mf[lev], 0, Temp, 1, scal.nGrow());
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(scal, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto tileBox = mfi.tilebox();
        
        const Array4<Real> scal_arr = scal.array(mfi);
        const Array4<const Real> p0_arr = p0_cart[lev].array(mfi);

        // initialize rho as sum of partial densities rho*X_i
        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal_arr(i,j,k,Rho) += scal_arr(i,j,k,FirstSpec+comp);
            }

            // initialize (rho h) and T using the EOS
            eos_t eos_state;
            eos_state.T     = scal_arr(i,j,k,Temp);
            eos_state.p     = p0_arr(i,j,k);
            eos_state.rho   = scal_arr(i,j,k,Rho);
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = scal_arr(i,j,k,FirstSpec+comp)/eos_state.rho;
            }

            eos(eos_input_rp, eos_state);

            scal_arr(i,j,k,RhoH) = eos_state.rho*eos_state.h;
            scal_arr(i,j,k,Temp) = eos_state.T;
        });

        if (perturb_model) {
            const auto prob_lo = geom[lev].ProbLoArray();
            const auto dx = geom[lev].CellSizeArray();

            const auto pert_rad_factor_loc = pert_rad_factor;
            const auto pert_temp_factor_loc = pert_temp_factor;
            const auto do_small_domain_loc = do_small_domain;

            // add an optional perturbation
            AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
                Real x = prob_lo[0] + (Real(i)+0.5) * dx[0];
                Real y = prob_lo[1] + (Real(j)+0.5) * dx[1];
                Real z = prob_lo[2] + (Real(k)+0.5) * dx[2];

                Real perturbations[Nscal];
                Real s0[Nscal];
                
                for (auto n = 0; n < Nscal; ++n) {
                    s0[n] = scal_arr(i,j,k,n);
                }

                Perturb(p0_arr(i,j,k), s0, perturbations, 
                    x, y, z, 
                    pert_rad_factor_loc, 
                    pert_temp_factor_loc, 
                    do_small_domain_loc, true);

                scal_arr(i,j,k,Rho) = perturbations[Rho];
                scal_arr(i,j,k,RhoH) = perturbations[RhoH];
                scal_arr(i,j,k,Temp) = perturbations[Temp];
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    scal_arr(i,j,k,FirstSpec+comp) = perturbations[FirstSpec+comp];
                }
            });
        }
    }
}

void Perturb(const Real p0_init, 
             const Real* s0_init,
             Real* perturbations,  
             const Real x, const Real y, const Real z,
             const Real pert_rad_factor,
             const Real pert_temp_factor,
             const bool do_small_domain,
             bool spherical)
{
    Real t0 = s0_init[Temp];

#if (AMREX_SPACEDIM == 2)

    Real x1 = 5.0e7;
    Real y1 = 6.5e7;
    Real r1 = std::sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)) / (2.5e6*pert_rad_factor);

    Real x2 = 1.2e8;
    Real y2 = 8.5e7;
    Real r2 = std::sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)) / (2.5e6*pert_rad_factor);

    Real x3 = 2.0e8;
    Real y3 = 7.5e7;
    Real r3 = std::sqrt((x-x3)*(x-x3)+(y-y3)*(y-y3)) / (2.5e6*pert_rad_factor);

    // this is a tiny bubble for inputs_2d_smalldomain
    Real x4 = 0.5;
    Real y4 = 85218750.25;
    Real r4 = std::sqrt((x-x4)*(x-x4)+(y-y4)*(y-y4)) / (2.5e-2*pert_rad_factor);

    Real temp = 0.0;

    if (do_small_domain) {
        temp = t0 * (1.0 + pert_temp_factor * 
                (0.150 * (1.0 + std::tanh(2.0-r1)) + 
                0.300 * (1.0 + std::tanh(2.0-r2)) + 
                0.225 * (1.0 + std::tanh(2.0-r3)) + 
                0.300 * (1.0 + std::tanh(2.0-r4))));
    } else {
        temp = t0 * (1.0 + pert_temp_factor * 
                (0.150 * (1.0 + std::tanh(2.0-r1)) + 
                0.300 * (1.0 + std::tanh(2.0-r2)) + 
                0.225 * (1.0 + std::tanh(2.0-r3))));
    }

#else 

    Real temp = 0.0;

    if (!spherical) {
        Real x0 = 1.8e7;
        Real y0 = 1.8e7;
        Real z0 = 8.5e7;

        Real r0 = std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) / 2.5e6;

        temp = t0 * (1.0 + 2.0 * (0.15 * (1.0 + std::tanh((2.0-r0)))));
    } else {
        // center of the star is a 2.5e8
        Real x0 = 2.5e8;
        Real y0 = 2.5e8;
        Real z0 = 2.5e8;

        Real r0 = std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0)) / 2.5e6;

        // note extra factor of 0.5 for spherical case
        temp = t0 * (1.0 + 2.0 * (0.15 * 0.5 * (1.0 + std::tanh((2.0-r0)))));
    }

#endif

    eos_t eos_state;
    
    // Use the EOS to make this temperature perturbation occur at constant
    // pressure
    eos_state.T     = temp;
    eos_state.p     = p0_init;
    eos_state.rho   = s0_init[Rho];
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = s0_init[FirstSpec+comp]/s0_init[Rho];
    }

    eos(eos_input_tp, eos_state);

    perturbations[Rho] = eos_state.rho;
    perturbations[RhoH] = eos_state.rho * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        perturbations[FirstSpec+comp] = eos_state.rho*eos_state.xn[comp];
    }

    perturbations[Temp] = temp;
}
