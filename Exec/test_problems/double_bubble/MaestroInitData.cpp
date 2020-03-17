
#include <Maestro.H>
using namespace amrex;

// prototype for pertubation function to be called on the 
// device (if USE_CUDA=TRUE)
AMREX_GPU_DEVICE
void Perturb(const Real p0_init, 
             const Real* s0_init,
             Real* perturbations,  
             const GpuArray<Real,AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real,AMREX_SPACEDIM> prob_hi,
             const Real x, const Real y, 
             const Real y_pert_center,
             const Real pert_width,
             const bool single);

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, 
                       const Array4<Real> scal, const Array4<Real> vel, 
                       const Real* s0_init, 
                       const Real* p0_init)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitData()", InitData);

    const auto tileBox = mfi.tilebox();
    const auto max_lev = max_radial_level + 1;
    const auto nrf = nr_fine;

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

    AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
        int r = AMREX_SPACEDIM == 2 ? j : k;

        // set the scalars using s0
        scal(i,j,k,Rho) = s0_init[lev+max_lev*(r+nr_fine*Rho)];
        scal(i,j,k,RhoH) = s0_init[lev+max_lev*(r+nr_fine*RhoH)];
        scal(i,j,k,Temp) = s0_init[lev+max_lev*(r+nr_fine*Temp)];
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i,j,k,FirstSpec+comp) = s0_init[lev+max_lev*(r+nr_fine*(FirstSpec+comp))];
        }
        // initialize pi to zero for now
        scal(i,j,k,Pi) = 0.0;
    });     
    
    // add an optional perturbation
    if (perturb_model) {

        const auto prob_lo = geom[lev].ProbLoArray();
        const auto prob_hi = geom[lev].ProbHiArray();
        const auto dx = geom[lev].CellSizeArray();

        const auto y_pert_center_loc = y_pert_center;
        const auto pert_width_loc = pert_width;
        const auto single_loc = single;

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            Real x = prob_lo[0] + (Real(i)+0.5) * dx[0];
            Real y = prob_lo[1] + (Real(j)+0.5) * dx[1];
            Real z = prob_lo[2] + (Real(k)+0.5) * dx[2];

            Real perturbations[Nscal];
            Real s0[Nscal];

            for (auto n = 0; n < Nscal; ++n) {
                s0[n] = s0_init[lev+max_lev*(r+nrf*n)];
            }

            Perturb(p0_init[lev+max_lev*r], s0, perturbations, 
                    prob_lo, prob_hi,
                    x, y,
                    y_pert_center_loc, 
                    pert_width_loc, 
                    single_loc);

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
                       const MFIter& mfi, MultiFab& scal, MultiFab& vel, 
                       const RealVector& s0_init, 
                       const RealVector& p0_init)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitDataSphr()", InitDataSphr);

#ifndef AMREX_USE_GPU
    Abort("Error: InitLevelDataSphr not implemented");
#endif
}

void Perturb(const Real p0_init, 
             const Real* s0_init,
             Real* perturbations,  
             const GpuArray<Real,AMREX_SPACEDIM> prob_lo,
             const GpuArray<Real,AMREX_SPACEDIM> prob_hi,
             const Real x, const Real y, 
             const Real y_pert_center,
             const Real pert_width,
             const bool single)
{
    // apply an optional perturbation to the initial temperature field
    // to see some bubbles

    eos_t eos_state;

#if (AMREX_SPACEDIM == 2)

    if (!single) {

        Real x1 = prob_lo[0] + (prob_hi[0]-prob_lo[0])/3.0;
        Real x2 = prob_lo[0] + 2.0*(prob_hi[0]-prob_lo[0])/3.0;

        Real y1 = y_pert_center;
        Real y2 = y_pert_center;

        Real r1 = std::sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)) / pert_width;
        Real r2 = std::sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)) / pert_width;

        if (r1 < 2.0) {
            eos_state.rho = s0_init[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0-r1))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[1] = 1.0;
        } else if (r2 < 2.0) {
            eos_state.rho = s0_init[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0-r2))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[2] = 1.0;
        } else {
            eos_state.rho = s0_init[Rho];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = s0_init[FirstSpec+comp]/s0_init[Rho];
            }
        } 
    } else {
        
        Real x1 = prob_lo[0] + 0.5*(prob_hi[0]-prob_lo[0]);

        Real y1 = y_pert_center;

        Real r1 = std::sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)) / pert_width;

        if (r1 < 2.0) {
            eos_state.rho = s0_init[Rho] * (1.0 - (pert_factor * (1.0 + std::tanh(2.0-r1))));
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = 0.0;
            }
            eos_state.xn[1] = 1.0;
        } else {
            eos_state.rho = s0_init[Rho];
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = s0_init[FirstSpec+comp]/s0_init[Rho];
            }
        } 
    }
    
    // Use the EOS to make this temperature perturbation occur at constant
    // pressure
    eos_state.T     = 10000.0; // guess
    eos_state.p     = p0_init;

    eos(eos_input_rp, eos_state);

    perturbations[Rho] = eos_state.rho;
    perturbations[RhoH] = eos_state.rho * eos_state.h;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        perturbations[FirstSpec+comp] = eos_state.rho*eos_state.xn[comp];
    }

    perturbations[Temp] = eos_state.T;

#else 
#ifndef AMREX_USE_GPU
    Abort("Error: Perturb not implemented for 3d");
#endif
#endif

}