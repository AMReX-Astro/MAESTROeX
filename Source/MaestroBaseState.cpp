#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& rho0, RealVector& rhoh0, 
                       RealVector& p0, 
                       const int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    if (use_exact_base_state && !spherical) {
        Abort("Irregular base state not valid for planar");
    }

    if (!input_model.model_initialized) {
        // read model file 
        input_model.ReadFile(model_file);
    }

    const int npts_model = input_model.npts_model;
    const Real TINY = 1.e-10;
    const int max_lev = max_radial_level + 1;
    const int n = lev;

    Real base_cutoff_density_loc = 1.e99;
    Real model_dr = (input_model.model_r[npts_model-1] - input_model.model_r[0]) / Real(npts_model - 1);
    Real rmax = input_model.model_r[npts_model-1];

    // irregular base state dr 
    RealVector model_dr_irreg;
    if (use_exact_base_state) {
        model_dr_irreg.resize(npts_model);

        model_dr_irreg[0] = input_model.model_r[0];
        for (auto i = 1; i < npts_model; ++i) {
            model_dr_irreg[i] = input_model.model_r[i] - input_model.model_r[i-1];
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        if (!spherical) {
            log_file.Log("model file mapping, level: ", n);
        } else {
            log_file.Log("model file mapping (spherical base state)");
        }

        log_file.Log("dr of MAESTRO base state =                            ", dr.array()(n));
        log_file.Log("dr of input file data =                               ", model_dr);
        log_file.Log(" ");
        log_file.Log("maximum radius (cell-centered) of input model =       ", rmax);

        Real mod_dr = 0.0;
        if (use_exact_base_state) {
            mod_dr = r_cc_loc_b.array()(lev,0) < model_dr_irreg[0] ? remainder(model_dr_irreg[0], r_cc_loc_b.array()(lev,0)) : remainder(r_cc_loc_b.array()(lev,0), model_dr_irreg[0]);
        } else {
            mod_dr = dr.array()(n) < model_dr ? 
                remainder(model_dr, dr.array()(n)) : remainder(dr.array()(n), model_dr);
        }

        if (mod_dr > TINY) {
            log_file.Log(" ");
            log_file.Log("WARNING: resolution of base state array is not an integer");
            log_file.Log("         multiple of the initial model's resolution.     ");
            log_file.Log("         make sure this is a desired property as this    ");
            log_file.Log("         could lead to aliasing when performing the      ");
            log_file.Log("         interpolation.                                  ");
            log_file.Log(" ");
            log_file.Log("modulus = ", mod_dr);
        }
    }

    Real starting_rad = spherical ? 0.0 : geom[0].ProbLo(AMREX_SPACEDIM-1);

    Real rho_above_cutoff = s0_init[n+max_lev*nr_fine*Rho];
    Real rhoh_above_cutoff = s0_init[n+max_lev*nr_fine*RhoH];
    RealVector spec_above_cutoff(NumSpec);
    for (auto comp = 0; comp < NumSpec; ++comp) {
        spec_above_cutoff[comp] = s0_init[n+max_lev*nr_fine*(FirstSpec+comp)];
    } 
    Real temp_above_cutoff = s0_init[n+max_lev*nr_fine*Temp];
    Real p_above_cutoff = p0_init[n];

    for (auto r = 0; r < nr.array()(n); ++r) {

        Real rloc = starting_rad + (Real(r) + 0.5)*dr.array()(n);

        // here we account for r > rmax of the model.hse array, assuming
        // that the state stays constant beyond rmax
        rloc = min(rloc, rmax);

        // also, if we've falled below the cutoff density, just keep the
        // model constant
        if (rloc > base_cutoff_density_loc) {

            s0_init[n+max_lev*(r+nr_fine*Rho)] = rho_above_cutoff;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = rhoh_above_cutoff;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = spec_above_cutoff[comp];
            }
            p0_init[n+max_lev*r] = p_above_cutoff;
            s0_init[n+max_lev*(r+nr_fine*Temp)] = temp_above_cutoff;

        } else {

            Real d_ambient = input_model.Interpolate(rloc, input_model.idens_model);
            Real t_ambient = input_model.Interpolate(rloc, input_model.itemp_model);
            Real p_ambient = input_model.Interpolate(rloc, input_model.ipres_model);

            RealVector xn_ambient(NumSpec);

            Real sumX = 0.0;

            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = max(0.0, min(1.0, 
                    input_model.Interpolate(rloc, input_model.ispec_model+comp)));
                sumX += xn_ambient[comp];
            }

            for (auto comp = 0; comp < NumSpec; ++comp) {
                if (sumX != 0.0) xn_ambient[comp] /= sumX;
            }

            eos_t eos_state;

            // use the EOS to make the state consistent
            eos_state.T     = t_ambient;
            eos_state.rho   = d_ambient;
            eos_state.p     = p_ambient;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                eos_state.xn[comp] = xn_ambient[comp];
            }

            // (rho,T) --> p,h
            eos(eos_input_rt, eos_state);

            s0_init[n+max_lev*(r+nr_fine*Rho)] = d_ambient;
            s0_init[n+max_lev*(r+nr_fine*RhoH)] = d_ambient * eos_state.h;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))] = 
                    d_ambient * xn_ambient[comp];
            }
            p0_init[n+max_lev*r] = eos_state.p; // p_ambient !
            s0_init[n+max_lev*(r+nr_fine*Temp)] = t_ambient;

            // keep track of the height where we drop below the cutoff density
            if (s0_init[n+max_lev*(r+nr_fine*Rho)] <= base_cutoff_density && 
                base_cutoff_density_loc == 1.e99) {

                Print() << ' ' << std::endl;
                Print() << "setting r_cutoff to " << r << std::endl;
                Print() << "radius at r_cutoff " << rloc << std::endl;

                base_cutoff_density_loc = rloc;

                rho_above_cutoff = s0_init[n+max_lev*(r+nr_fine*Rho)];
                rhoh_above_cutoff = s0_init[n+max_lev*(r+nr_fine*RhoH)];
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    spec_above_cutoff[comp] = 
                        s0_init[n+max_lev*(r+nr_fine*(FirstSpec+comp))];
                }
                temp_above_cutoff = s0_init[n+max_lev*(r+nr_fine*Temp)];
                p_above_cutoff = p0_init[n+max_lev*r];
            }
        }
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    for (auto i = 0; i < nr_fine; ++i) {
        rho0[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Rho)];
        rhoh0[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*RhoH)];
        tempbar[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        tempbar_init[lev+max_lev*i] = s0_init[lev+max_lev*(i+nr_fine*Temp)];
        p0[lev+max_lev*i] = p0_init[lev+max_lev*i];
    }

    // check whether we are in HSE

    Real mencl = 0.0;

    if (use_exact_base_state) {
        Real dr_irreg = r_edge_loc_b.array()(n,1) - r_edge_loc_b.array()(n,0); // edge-to-edge

        if (spherical || do_2d_planar_octant) {
            mencl = 4.0/3.0 * M_PI * dr_irreg*dr_irreg*dr_irreg * s0_init[n+max_lev*nr_fine*Rho];
        }
    } else {
        if (spherical || do_2d_planar_octant) {
            mencl = 4.0/3.0 * M_PI * dr.array()(n)*dr.array()(n)*dr.array()(n) * s0_init[n+max_lev*nr_fine*Rho];
        }
    }

    Real max_hse_error = -1.e30;

    for (auto r = 1; r < nr.array()(n); ++r) {

        Real rloc = starting_rad + (Real(r) + 0.5) * dr.array()(n);
        rloc = min(rloc, rmax);

        if (rloc < base_cutoff_density_loc) {

            Real r_r = starting_rad;
            r_r += use_exact_base_state ? r_edge_loc_b.array()(n,r+1) : Real(r+1) * dr.array()(n);
            Real r_l = starting_rad;
            r_l += use_exact_base_state ? r_edge_loc_b.array()(n,r) : Real(r) * dr.array()(n);

            Real g = 0.0;

            if (spherical || do_2d_planar_octant) {
                g = -Gconst * mencl / (r_l*r_l);
                mencl += 4.0/3.0 * M_PI * dr.array()(n) * 
                    (r_l*r_l+r_l*r_r+r_r*r_r) * 
                    s0_init[n+max_lev*(r+nr_fine*Rho)];
            } else {
                if (!do_planar_invsq_grav) {
                    g = grav_const;
                } else {
                    g = -Gconst*planar_invsq_mass / (r_l*r_l);
                }
            }
            
            Real dpdr = 0.0;
            Real rhog = 0.0;
            if (use_exact_base_state) {
                Real dr_irreg = r_cc_loc_b.array()(n,r) - r_cc_loc_b.array()(n,r-1);
                dpdr = (p0_init[n+max_lev*r] - p0_init[n+max_lev*(r-1)])/dr_irreg;

                Real rfrac = (r_edge_loc_b.array()(n,r) - r_cc_loc_b.array()(n,r-1)) / dr_irreg; 
                rhog = ((1.0-rfrac) * s0_init[n+max_lev*(r+nr_fine*Rho)] +  
                        rfrac * s0_init[n+max_lev*(r-1+nr_fine*Rho)]) * g;
            } else {
                dpdr = (p0_init[n+max_lev*r] - p0_init[n+max_lev*(r-1)])/dr.array()(n);
                rhog = 0.5 * (s0_init[n+max_lev*(r+nr_fine*Rho)] +  
                              s0_init[n+max_lev*(r-1+nr_fine*Rho)]) * g;
            }

            if (print_init_hse_diag) {
                Print() << "r, dpdr, rhog, err: " << rloc << ", " 
                        << dpdr  << ", " << rhog  << ", " 
                        << fabs(dpdr - rhog)/fabs(rhog) << std::endl;
            }

            max_hse_error = max(max_hse_error, fabs(dpdr - rhog)/fabs(rhog));
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        log_file.Log(" ");
        log_file.Log("Maximum HSE Error = ", max_hse_error);
        log_file.Log("   (after putting initial model into base state arrays, and");
        log_file.Log("    for density < base_cutoff_density)");
        log_file.Log(" ");
    }

    // initialize any inlet BC parameters
    SetInletBCs();
}
