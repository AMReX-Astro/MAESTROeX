#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void 
Maestro::InitBaseState(RealVector& s0_init, RealVector& p0_init, 
                       RealVector& rho0, RealVector& rhoh0, 
                       RealVector& p0 ,RealVector& tempbar, 
                       RealVector& tempbar_init)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitBaseState()", InitBaseState); 

    if (!input_model.model_initialized) {
        // read model file 
        input_model.ReadFile(model_file);
    }

    Real base_cutoff_density_loc = 1.e99;
    Real model_dr = (input_model.model_r[npts_model-1] - input_model.model_r[0]) / Real(npts_model - 1);
    Real rmax = input_model.model_r[npts_model-1];
    const int max_lev = max_radial_level + 1;

    int n = lev;

    if ( parallel_IOProcessor() ) {
        write (*,887)
        if (!spherical) {
            call log('model file mapping, level: ', n)
        } else {
            call log('model file mapping (spherical base state)')
        }

        call log('dr of MAESTRO base state =                            ', dr[n])
        call log('dr of input file data =                               ', model_dr)
        call log(' ')
        call log('maximum radius (cell-centered) of input model =       ', rmax)

        if (dr[n] > model_dr) {
            mod_dr = mod(model_dr,dr[n])
        } else {
            mod_dr = mod(dr[n],model_dr)
        }

        if (mod_dr > TINY) {
            call log(' ')
            call log("WARNING: resolution of base state array is not an integer")
            call log("         multiple of the initial model's resolution.     ")
            call log("         make sure this is a desired property as this    ")
            call log("         could lead to aliasing when performing the      ")
            call log("         interpolation.                                  ")
            call log(" ")
            call log("modulus = ", mod_dr)
        }
    }

    if (!spherical) {
        starting_rad = prob_lo[AMREX_SPACEDIM-1];
    } else {
        starting_rad = 0.0
    }

    // do r=0,nr(n)-1
    for (auto r = 0; r < nr[n]; ++r) {

        Real rloc = starting_rad + (Real(r) + 0.5)*dr[n];

        // here we account for r > rmax of the model.hse array, assuming
        // that the state stays constant beyond rmax
        rloc = min(rloc, rmax);

        // also, if we've falled below the cutoff density, just keep the
        // model constant
        if (rloc > base_cutoff_density_loc) {

            s0_init(n,r,Rho) = rho_above_cutoff
            s0_init(n,r,RhoH) = rhoh_above_cutoff
            s0_init(n,r,FirstSpec:FirstSpec+NumSpec-1) = spec_above_cutoff(1:NumSpec)
            p0_init[n+max_lev*r] = p_above_cutoff
            s0_init(n,r,Temp) = temp_above_cutoff

        } else {

            d_ambient = input_model.Interpolate(rloc, idens_model)
            t_ambient = input_model.Interpolate(rloc, itemp_model)
            p_ambient = input_model.Interpolate(rloc, ipres_model)

            Real sumX = 0.0;
            // do comp = 1, NumSpec
            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] = max(0.0, min(1.0, 
                        input_model.Interpolate(rloc, ispec_model+comp)))
                sumX += xn_ambient[comp]''
            }

            for (auto comp = 0; comp < NumSpec; ++comp) {
                xn_ambient[comp] /= sumX;
            }

            eos_t eos_state;

            // use the EOS to make the state consistent
            eos_state.T     = t_ambient;
            eos_state.rho   = d_ambient;
            eos_state.p     = p_ambient;
            eos_state.xn(:) = xn_ambient(:)

            // (rho,T) --> p,h
            eos(eos_input_rt, eos_state)

            s0_init(n,r, Rho ) = d_ambient
            s0_init(n,r,RhoH ) = d_ambient * eos_state.h
            s0_init(n,r,FirstSpec:FirstSpec+NumSpec-1) = d_ambient * xn_ambient(1:NumSpec)
            p0_init[n+max_lev*r] = eos_state.p; // p_ambient !
            s0_init(n,r,Temp) = t_ambient

            // keep track of the height where we drop below the cutoff density
            if (s0_init(n,r,Rho) <= base_cutoff_density && 
                base_cutoff_density_loc == 1.e99) {

                if ( parallel_IOProcessor() ) {
                    Print() << ' ' << std::endl;
                    Print() << 'setting r_cutoff to ' << r << std::endl;
                    Print() << 'radius at r_cutoff ' << rloc << std::endl;
                }

                base_cutoff_density_loc = rloc;

                rho_above_cutoff = s0_init(n,r,Rho)
                rhoh_above_cutoff = s0_init(n,r,RhoH)
                spec_above_cutoff(1:NumSpec) = s0_init(n,r,FirstSpec:FirstSpec+NumSpec-1)
                temp_above_cutoff = s0_init(n,r,Temp)
                p_above_cutoff = p0_init[n+max_lev*r];

            }
        }
    }

    // copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
    rho0 = s0_init(:,:,Rho)
    rhoh0 = s0_init(:,:,RhoH)
    tempbar = s0_init(:,:,Temp)
    tempbar_init = s0_init(:,:,Temp)
    p0 = p0_init

    // check whether we are in HSE

    Real mencl = 0.0;

    if (spherical || do_2d_planar_octant) {
        mencl = four3rd*M_PI*dr[n]**3*s0_init(n,0,Rho)
    }

    Real max_hse_error = -1.e30;

    // do r=1,nr(n)-1
    for (auto r = 1; r < nr[n]; ++r) {

        Real rloc = starting_rad + (Real(r) + 0.5)*dr[n];
        rloc = min(rloc, rmax);

        if (rloc > base_cutoff_density_loc) {

            Real r_r = starting_rad + Real(r+1)*dr[n];
            Real r_l = starting_rad + Real(r)*dr[n];

            if (spherical || do_2d_planar_octant) {
                g = -Gconst*mencl/(r_l*r_l)
                mencl += 4.0/3.0*M_PI*dr[n]*(r_l*r_l+r_l*r_r+r_r*r_r)*s0_init(n,r,Rho)
            } else {
                if (!do_planar_invsq_grav) {
                    g = grav_const;
                } else {
                    g = -Gconst*planar_invsq_mass / (r_l*r_l);
                }
            }

            dpdr = (p0_init[n+max_lev*r] - p0_init[n+max_lev*(r-1)])/dr[n];
            rhog = 0.5*(s0_init(n,r,Rho) + s0_init(n,r-1,Rho))*g

            if (print_init_hse_diag) {
                if ( parallel_IOProcessor() ) {
                    print *, 'r, dpdr, rhog, err: ', rloc, dpdr, rhog, 
                        fabs(dpdr - rhog)/fabs(rhog)
                }
            }

            max_hse_error = max(max_hse_error, fabs(dpdr - rhog)/fabs(rhog));
        }
    }

    if ( parallel_IOProcessor() ) {
        call log(' ')
        call log('Maximum HSE Error = ', max_hse_error)
        call log('   (after putting initial model into base state arrays, and')
        call log('    for density < base_cutoff_density)')
        write (*,887)
        call log(' ')
    }

    // initialize any inlet BC parameters
    call set_inlet_bcs()

}