
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void
Maestro::InitLevelData(const int lev, const Real time, 
                       const MFIter& mfi, const Array4<Real> scal, const Array4<Real> vel, 
                       const Real* s0_p, 
                       const Real* p0_p)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();
    const int max_lev = base_geom.max_radial_level + 1;
    const auto nrf = base_geom.nr_fine;

    // set velocity to zero 
    AMREX_PARALLEL_FOR_4D(tileBox, AMREX_SPACEDIM, i, j, k, n, {
        vel(i,j,k,n) = 0.0;
    });

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

    const auto prob_lo = geom[lev].ProbLoArray();
    const auto prob_hi = geom[lev].ProbHiArray();
    const auto dx = geom[lev].CellSizeArray();

    if (perturb_model) {
        const auto xcen = center[0];
        const auto ycen = AMREX_SPACEDIM == 2 ? xrb_pert_height : center[1];
        const auto zcen = AMREX_SPACEDIM == 2 ? 0.0 : xrb_pert_height;

        const auto rad_pert = -xrb_pert_size*xrb_pert_size / (4.0 * std::log(0.5));

        const bool perturb_temp_true = xrb_pert_type == 1;

        const auto xrb_pert_factor_loc = xrb_pert_factor;
        const auto rad_pert_loc = rad_pert;

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            int r = AMREX_SPACEDIM == 2 ? j : k;

            const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0] - xcen;
            const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1] - ycen;
            const auto z = prob_lo[2] + (Real(k) + 0.5) * dx[2] - zcen;

            const auto dist = std::sqrt(x*x + y*y + z*z);

            Real temp = 0.0;
            Real dens = 0.0;
            auto eos_input_flag = eos_input_tp;

            if (perturb_temp_true) {
                Real t0 = s0_init[lev+max_lev*(r + nrf*Temp)];
                temp = t0 * (1.0 + xrb_pert_factor_loc * std::exp(-dist*dist / rad_pert_loc));
                dens = s0_init[lev+max_lev*(r + nrf*Rho)];
                eos_input_flag = eos_input_tp;
            } else {
                Real d0 = s0_init[lev+max_lev*(r + nrf*Rho)];
                dens = d0 * (1.0 + xrb_pert_factor_loc * std::exp(-dist*dist / rad_pert_loc));
                temp = s0_init[lev+max_lev*(r + nrf*Temp)];
                eos_input_flag = eos_input_rp;
            }

            eos_t eos_state;

            eos_state.T = temp;
            eos_state.p = p0_init[lev+max_lev*r];
            eos_state.rho = dens;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = s0_init[lev+max_lev*(r+nrf*(FirstSpec+comp))] / s0_init[lev+max_lev*(r + nrf*Rho)];
            }

            eos(eos_input_flag, eos_state);

            scal(i,j,k,Rho) = eos_state.rho;
            scal(i,j,k,RhoH) = eos_state.rho * eos_state.h;
            scal(i,j,k,Temp) = eos_state.T;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i,j,k,FirstSpec+comp) = eos_state.rho * eos_state.xn[comp];
            }
        });
    }

    if (apply_vel_field) {
        const auto velpert_amplitude_loc = velpert_amplitude;
        const auto velpert_scale_loc = velpert_scale;
        const auto num_vortices_loc = num_vortices;
        const auto velpert_height = velpert_height_loc;

        const Real offset = (prob_hi[0] - prob_lo[0]) / (num_vortices + 1);

        // vortex x-coords
        RealVector vortices_xloc(num_vortices);
        for (auto i = 0; i < num_vortices; ++i) {
            vortices_xloc[i] = (Real(i) + 0.5) * offset;
        }

        const Real * vortices_xloc_p = vortices_xloc.dataPtr();

        AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {
            const auto x = prob_lo[0] + (Real(i) + 0.5) * dx[0];
            const auto y = prob_lo[1] + (Real(j) + 0.5) * dx[1];

            Real ydist = y - velpert_height;

            Real upert = 0.0;
            Real vpert = 0.0;

            for (auto vortex = 0; vortex < num_vortices_loc; ++vortex) {
                Real xdist = x - vortices_xloc_p[vortex];

                Real rad = std::sqrt(x*x + y*y);

                // e.g. Calder et al. ApJSS 143, 201-229 (2002)
                // we set things up so that every other vortex has the same
                // orientation
                upert -= ydist/velpert_scale_loc * velpert_amplitude_loc * std::exp(-rad*rad / (2.0 * velpert_scale_loc*velpert_scale_loc)) * pow(-1.0, vortex+1);

                vpert += xdist/velpert_scale_loc * velpert_amplitude_loc * std::exp(-rad*rad / (2.0 * velpert_scale_loc*velpert_scale_loc)) * pow(-1.0, vortex+1);
            }

            vel(i,j,k,0) += upert;
            vel(i,j,k,1) += vpert;
        });
    }
}

void
Maestro::InitLevelDataSphr(const int lev, const Real time, 
                           MultiFab& scal, MultiFab& vel)
{
    Abort("InitLevelDataSphr not implemented.");
}
