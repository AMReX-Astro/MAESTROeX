//  Adapted from original AMReX-Astro code
//  (github.com/AMReX-Astro/initial_models/tree/main/toy_atm/init_1d.H)
//  by Simon Guichandut (2023).

//  Create a 1-d hydrostatic, atmosphere with an isothermal region
//  (T_star) representing the NS, a hyperbolic tangent rise to a
//  peak temperature (T_base) representing the base of an accreted
//  layer at r=H_star, an isentropic profile representing the convective
//  region at ignition, then a radiative profile from r=H_rad. The 
//  composition transitions from pure ash, to pure "fuel 2" (helium), 
//  to a mixture of fuels given with mass fractions specified in the 
//  inputs file.

//  This can serve as an initial model for a He-triggered XRB with an
//  overlying hydrogen rich shell
//
//  The temperature profile is:
//
//         T
//         ^
//         |
//  T_base +        /\
//         |       /  \
//         |      /  . \
//  T_star +-----+      \  nabla=nabla_ad
//         |     .   .   \
//         |              \
//  T_lo   +     .   .     +_
//         |     .   .     .  `-
//         |     .   .     .     `-  nabla=nabla_rad
//         |     .   .     .        `-
//         |     .   .     .           `-
//         |     .   .     .              `-
//         +-----+---+-----+---------------------> r
//         |      \  /
//         |      delta
//         |< H_star>|
//         |---<H_rad>----| 
//
//  The composition profile is:
//
//         X
//         ^
//         |  star    layer1      layer2 
//    1.0  +--------*----------                 
//         |  ash   .          \
//         |        .           \    +--------- fuel1_frac
//         |        .            \  /
//         |        .             \/
//         |        .             /\
//         |        .            /  \
//         |        .           /    +--------- fuel2_frac
//         |        .          /     +--------- fuel3_frac
//    0.0  +--------+---------+-----------------> r
//         |<H_star>|    
//         |---<H_rad>----| 
//         |------<H_sig>-----|     
//         |
//
//  The transition from ash to fuel2 is a hyperbolic tangent, parametrized
//  by delta (the half-width of the tanh).  The transition from layer1 to
//  layer2 is a sigmoid, parametrized by alpha_sig and beta_sig.  H_star and
//  H_sig controls the location of these transitions.
//
//  The fuel and ash compositions are specified by the fuel?_name,
//  fuel?_frac and ash?_name, ash?_frac parameters (name of the species
//  and mass fraction).  Where ? = 1,2,3.
//
//  The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
//
//  This will be iterated over in tandem with the EOS call,
//  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1))
//

#include <AMReX_Array.H>

#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <model_util.H>

using namespace amrex;

AMREX_INLINE void init_1d()
{

    // define convenient indices for the scalars

    const int nvar = 3 + NumSpec;
    const int idens = 0;
    const int itemp = 1;
    const int ipres = 2;
    const int ispec = 3;

    const Real TOL = 1.e-10_rt;

    const int MAX_ITER = 250;

    const Real smallx = 1.e-10_rt;

    // get the species indices

    bool species_defined{true};

    int ifuel1 = network_spec_index(fuel1_name);
    if (ifuel1 < 0)
        species_defined = false;

    int ifuel2{-1};
    if (fuel2_name != "")
    {
        ifuel2 = network_spec_index(fuel2_name);
        if (ifuel2 < 0)
            species_defined = false;
    }

    int ifuel3{-1};
    if (fuel3_name != "")
    {
        ifuel3 = network_spec_index(fuel3_name);
        if (ifuel3 < 0)
            species_defined = false;
    }

    int ifuel4{-1};
    if (fuel4_name != "")
    {
        ifuel4 = network_spec_index(fuel4_name);
        if (ifuel4 < 0)
            species_defined = false;
    }


    int iash1 = network_spec_index(ash1_name);
    if (iash1 < 0)
        species_defined = false;

    int iash2{-1};
    if (ash2_name != "")
    {
        iash2 = network_spec_index(ash2_name);
        if (iash2 < 0)
            species_defined = false;
    }

    int iash3{-1};
    if (ash3_name != "")
    {
        iash3 = network_spec_index(ash3_name);
        if (iash3 < 0)
            species_defined = false;
    }

    int iash4{-1};
    if (ash4_name != "")
    {
        iash4 = network_spec_index(ash4_name);
        if (iash4 < 0)
            species_defined = false;
    }

    if (!species_defined)
    {
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << " " << ifuel4 << " " << std::endl;
        std::cout << iash1 << " " << iash2 << " " << iash3 << " " << iash4 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    Real xn[NumSpec] = {0.0_rt};

    // set the composition of the underlying star

    Real xn_star[NumSpec] = {smallx};

    xn_star[iash1] = ash1_frac;
    if (ash2_name != "")
        xn_star[iash2] = ash2_frac;
    if (ash3_name != "")
        xn_star[iash3] = ash3_frac;
    if (ash4_name != "")
        xn_star[iash4] = ash4_frac;


    // the composition of the first layer above the star
    // by definition it is just made of fuel2

    Real xn_base[NumSpec] = {smallx};
    xn_base[ifuel2] = 1.0;

    // and the composition of the top layer, for which the fractions are specified in inputs

    Real xn_top_layer[NumSpec] = {smallx};
    xn_top_layer[ifuel1] = fuel1_frac;
    xn_top_layer[ifuel2] = fuel2_frac;
    xn_top_layer[ifuel3] = fuel3_frac;
    xn_top_layer[ifuel4] = fuel4_frac;


    // check if they sum to 1
    double sum{0.0};
    for (auto e : xn_star) {
        sum += e;
    }
    if (std::abs(sum) - 1.0_rt > NumSpec * smallx) {
        amrex::Error("ERROR: ash mass fractions don't sum to 1");
    }

    sum = 0.0;
    for (auto e : xn_base) {
        if (std::abs(sum) - 1.0_rt > NumSpec * smallx) {
            amrex::Error("ERROR: fuel mass fractions don't sum to 1");
        }   
    }

    sum = 0.0;
    for (auto e : xn_top_layer) {
        if (std::abs(sum) - 1.0_rt > NumSpec * smallx) {
            amrex::Error("ERROR: fuel mass fractions don't sum to 1");
        }   
    }

    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    if (nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xznl_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznr_hse;

    Array2D<Real, 0, NPTS_MODEL-1, 0, nvar-1> model_hse;

    // compute the coordinates of the new gridded function

    Real dCoord = (xmax - xmin) / static_cast<Real>(nx);

    for (int i = 0; i < nx; ++i) {
        xznl_hse(i) = xmin + (static_cast<Real>(i)) * dCoord;
        xzn_hse(i) = xmin + (static_cast<Real>(i) + 0.5_rt) * dCoord;
        xznr_hse(i) = xmin + (static_cast<Real>(i) + 1.0_rt) * dCoord;
    }

    // put the model onto our new uniform grid

    bool fluff = false;

    // determine the conditions at the base

    eos_t eos_state;
    eos_state.T = T_base;
    eos_state.rho = dens_base;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    eos(eos_input_rt, eos_state);

    // store the conditions at the base -- we'll use the entropy later
    // to constrain the isentropic layer
  
    Real pres_base = eos_state.p;
    Real entropy_base = eos_state.s;


    // set an initial temperature profile and composition

    double a = 0.0;

    for (int i = 0; i < nx; ++i) {

        // Composition
        // hyperbolic tangent transition between star and layer1
        // sigmoid between layer1 and layer2

        for (int n = 0; n < NumSpec; ++n) {

            // tanh part
            a = xn_star[n] + 
                0.5_rt * (xn_base[n] - xn_star[n]) *
                (1.0_rt + std::tanh((xzn_hse(i) - (xmin + H_star - delta) + delta) / delta)); 

            // add sigmoid
            model_hse(i, ispec+n) = xn_top_layer[n] + (a - xn_top_layer[n]) /
                pow(1 + pow(xzn_hse(i)/H_sig, alpha_sig), beta_sig);
        }

        // Temperature
        // hyperbolic tangent transition for the temperature

        model_hse(i, itemp) = T_star + 0.5_rt * (T_base - T_star) *
            (1.0_rt + std::tanh((xzn_hse(i) - (xmin + H_star - delta) + delta) / delta));

        // the density and pressure will be determined via HSE,
        // for now, set them to the base conditions
     
        model_hse(i, idens) = dens_base;
        model_hse(i, ipres) = pres_base;

    }


    // find the index of the base height
    int index_base = -1;

    if (index_base_from_temp == 1) {
        for (int i = 0; i < nx; ++i) {
            if (model_hse(i, itemp) > 0.9995 * T_base) {
                index_base = i+1;
                break;  
            }
        }
    }
    else {
        for (int i = 0; i < nx; ++i) {
            if (xzn_hse(i) >= xmin + H_star + delta) {
                index_base = i+1;
                break;
            }
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: base_height not found on grid");
    }

    std::cout << "index_base = " << index_base << "; x=" << xzn_hse(index_base) << std::endl;


    // find the index of the radiative zone height
    int index_rad = -1;

    if (do_isothermal_instead) {
        index_rad = nx+1; // temporarily put this outside of the grid
    } else {
    
        if (H_rad < H_star) {
            amrex::Error("ERROR: H_rad<H_star. Adiabatic (base) layer should be below radiative layer");
        }

        for (int i = index_base; i < nx; ++i) {
            if (xzn_hse(i) >= xmin + H_rad) {
                index_rad = i;
                break;
            }
        }

        if (index_rad == -1) {
            amrex::Error("ERROR: rad_height not found on grid");
        }

        std::cout << "index_rad = " << index_rad << "; x=" << xzn_hse(index_rad) << std::endl;
    }


    // make the base thermodynamics consistent for this base point -- that is
    // what we will integrate from!

    eos_state.rho = model_hse(index_base, idens);
    eos_state.T = model_hse(index_base, itemp);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_hse(index_base, ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model_hse(index_base, ipres) = eos_state.p;




    // HSE solve

    // the HSE state will be done putting creating an isentropic state until
    // reaching H_rad. Then, we set the radiative profile using the luminosity.

    // once the density goes below low_density_cutoff, we stop HSE

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;

    Real p_want;
    Real drho;
    Real dtemp;
    Real entropy;

    // integrate up until radiative zone
    for (int i = index_base+1; i < index_rad+1; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        // compute the gravitation acceleration at the lower edge
        Real g_zone;
        if (do_invsq_grav == 1) {
            g_zone = -C::Gconst * M_enclosed / std::pow(xznl_hse(i), 2);
        } else {
            g_zone = g_const;
        }

        // we've already set initial guesses for density, temperature, and  composition
        dens_zone = model_hse(i, idens);
        temp_zone = model_hse(i, itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, ispec+n);
        }

        // iteration loop
        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            // get the pressure we want from the HSE equation, just the
            // zone below the current.  Note, we are using an average of
            // the density of the two zones as an approximation of the
            // interface value -- this means that we need to iterate for
            // find the density and pressure that are consistent

            // furthermore, we need to get the entropy that we need,
            // which will come from adjusting the temperature in
            // addition to the density.

            // HSE differencing
        
            p_want = model_hse(i-1, ipres) +
            delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * g_zone;

            // now we have two functions to zero:
            //   A = p_want - p(rho,T)
            //   B = entropy_base - s(rho,T)
            // We use a two dimensional Taylor expansion and find the deltas
            // for both density and temperature

            // now we know the pressure and the entropy that we want, so we
            // need to find the temperature and density through a two
            // dimensional root find

            // (t, rho) -> (p, s)
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            entropy = eos_state.s;
            pres_zone = eos_state.p;

            Real A = p_want - pres_zone;
            Real B = entropy_base - entropy;

            // Solve F(rho,T) = [A(rho,T),B(rho,T)] = [0,0] with Newton-Raphson method
            // [drho,dT] = -inv(J_F) [A,B]
            // J_F is the jacobian matrix {{dA/drho,dA/dT},{dB/drho,dB/dT}}
            // Work out algebraically the expression below for the updates drho and dT

            // EOS quantities
            Real dpt = eos_state.dpdT; // dp/dT
            Real dpd = eos_state.dpdr; // dp/drho
            Real dst = eos_state.dsdT; // ds/dT
            Real dsd = eos_state.dsdr; // ds/drho

            dtemp = ((dsd / (dpd - 0.5_rt * delx * g_zone)) * A - B)/
                (dsd * dpt / (dpd - 0.5_rt * delx * g_zone) - dst);

            drho = (A - dpt * dtemp) / (dpd - 0.5_rt * delx * g_zone);

            // Update density and temperature in the zone
            // This assumes both values are going down
            dens_zone = amrex::max(0.9_rt * dens_zone,
                                    amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            temp_zone = amrex::max(0.9_rt * temp_zone,
                                    amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));


            // check if the density falls below our minimum cut-off --
            // if so, floor it
            if (dens_zone < low_density_cutoff) {
                dens_zone = low_density_cutoff;
                temp_zone = T_lo;
                converged_hse = true;
                fluff = true;

                // raise a warning because this happens before the radiative part
                std::cout << "Warning: reached cut-off density before radiative zone" << std::endl;

                break;
            }

            // Stop if converged
            if (std::abs(drho) < TOL * dens_zone &&
                std::abs(dtemp) < TOL*temp_zone) {
                converged_hse = true;
                break;
            }

            // If reached T_lo and we intend to do isothermal, set index_rad to previous zone
            // break and restart this iteration in the next part of the calculation
            if (do_isothermal_instead && temp_zone < T_lo) {
                temp_zone = T_lo;
                index_rad = i-1;
                converged_hse = true;
                std::cout << "Starting isothermal at x=" << xzn_hse(index_rad) << std::endl;
                break;
            }
        
        }

        if (! converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "isentropic layer" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << p_want << " " << entropy_base << " " << entropy << std::endl;
            std::cout << drho << " " << dtemp << std::endl;
            amrex::Error("Error: HSE non-convergence");
        }

        // call the EOS one more time for this zone and then go on to the next
        // (t, rho) -> (p)
     
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone
        model_hse(i, idens) = dens_zone;
        model_hse(i, itemp) = temp_zone;
        model_hse(i, ipres) = pres_zone;

        // to make this process converge faster, set the density in the
        // next zone to the density in this zone
        model_hse(i+1,idens) = dens_zone;

    }

    std::cout << "Reached end of adiabatic region at r=" << xzn_hse(index_rad) <<" T=" << model_hse(index_rad,itemp) << std::endl;


    // integrate up to the top of the model

    // Some extra vars needed here
    Real rad_D = 0.0;
    Real t4 = 0.0;

    for (int i = index_rad; i < nx; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        // compute the gravitation acceleration at the lower edge
        Real g_zone;
        if (do_invsq_grav == 1) {
            g_zone = -C::Gconst * M_enclosed / std::pow(xznl_hse(i), 2);
        } else {
            g_zone = g_const;
        }

        // we've already set initial guesses for density, temperature, and  composition
        dens_zone = model_hse(i, idens);
        temp_zone = model_hse(i, itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, ispec+n);
        }

        // iteration loop
        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            if (fluff) {
                dens_zone = low_density_cutoff;
                temp_zone = T_lo;

            } else {

                p_want = model_hse(i-1, ipres) +
                    delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * g_zone;

                if (do_isothermal_instead) {

                    temp_zone = T_lo;
                    
                } else {

                    // In optically thick medium, radiation energy density is transported
                    // according to a diffusion law:
                    // Flux = -D d(T^4)/dx
                    // where the diffusion coefficient is D = ac/(3*kappa*rho)
                    // kappa is not included in EOS, just use classical electron scattering
                    // kap = 0.2*(1+X), where X is the hydrogen fraction
                    // std::cout << "hydrogen1 index in the network is:" << Species::H1-1 << std::endl; 

                    // kap = 0.2*(1+model_hse(i, ispec + Species::H1-1));
                    // temp_zone = model_hse(i-1, itemp) - delx * 3.0_rt*kap*flux*
                    //     0.5_rt*(dens_zone + model_hse(i-1, idens)) / 
                    //     (4.0_rt * C::a_rad * C::c_light * pow(model_hse(i-1, itemp), 3));
                    // t4 = pow(model_hse(i-1, itemp), 4) - delx * 3.0_rt*kap*F_rad*
                    //     0.5_rt*(dens_zone + model_hse(i-1, idens)) / (C::a_rad * C::c_light);

                    // Average the diffusion coefficient over two zones
                    rad_D = C::a_rad * C::c_light /
                            (3 * 0.2 * 0.5 * (
                                (1 + model_hse(i, ispec+Species::H1-1)) * dens_zone +
                                (1 + model_hse(i-1, ispec+Species::H1-1)) * model_hse(i-1, idens)
                            ));

                    // Update T^4
                    t4 = pow(model_hse(i-1, itemp), 4) - delx*F_rad/rad_D;

                    // std::cout << pow(t4,0.25) << std::endl;
                    // break;

                    if (t4 < pow(T_lo, 4.0_rt)) {
                    // if (t4 < 0.0) {
                    // if (true) {
                        // amrex::Error("Error: T<0");
                        temp_zone = T_lo; // set to isothermal
                    } else {
                        temp_zone = pow(t4, 0.25_rt);
                    }
                
                    // std::cout << "T=" << temp_zone << std::endl;

                }

            }

            // (t, rho) -> (p)
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            entropy = eos_state.s;
            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * g_zone);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                    amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            if (std::abs(drho) < TOL * dens_zone) {
                converged_hse = true;
                break;
            }

            if (dens_zone < low_density_cutoff) {
                dens_zone = low_density_cutoff;
                temp_zone = T_lo;
                converged_hse = true;
                fluff = true;
                break;
            }
        }

        if (! converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "x = " << xzn_hse(i) << "(radiative layer)" << std::endl;
            std::cout << "rho,T = " << dens_zone << " " << temp_zone << std::endl;
            std::cout << "drho,dT = " << drho << " " << dtemp << std::endl;
            amrex::Error("Error: HSE non-convergence");
        }

        // call the EOS one more time for this zone and then go on to the next
        // (t, rho) -> (p)
     
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone
        model_hse(i, idens) = dens_zone;
        model_hse(i, itemp) = temp_zone;
        model_hse(i, ipres) = pres_zone;

        // to make this process converge faster, set the density in the
        // next zone to the density in this zone
        model_hse(i+1,idens) = dens_zone;

    }


    // integrate down -- using the temperature profile defined above
  
    for (int i = index_base-1; i >= 0; --i) {

        Real delx = xzn_hse(i+1) - xzn_hse(i);

        // compute the gravitation acceleration at the upper edge

        Real g_zone;
        if (do_invsq_grav == 1) {
            g_zone = -C::Gconst * M_enclosed / std::pow(xznr_hse(i), 2);
        } else {
            g_zone = g_const;
        }

        // we already set the temperature and composition profiles
        temp_zone = model_hse(i, itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, ispec+n);
        }

        // use our previous initial guess for density
     
        dens_zone = model_hse(i+1, idens);


        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            // get the pressure we want from the HSE equation, just the
            // zone below the current.  Note, we are using an average of
            // the density of the two zones as an approximation of the
            // interface value -- this means that we need to iterate for
            // find the density and pressure that are consistent

            // HSE differencing
        
            p_want = model_hse(i+1, ipres) -
                delx * 0.5_rt * (dens_zone + model_hse(i+1, idens)) * g_zone;

            // we will take the temperature already defined in model_hse
            // so we only need to zero:
            //   A = p_want - p(rho)

            // (t, rho) -> (p)
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            Real A = p_want - pres_zone;

            drho = A / (dpd + 0.5_rt * delx * g_zone);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                   amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            if (std::abs(drho) < TOL * dens_zone) {
                converged_hse = true;
                break;
            }

        }

        if (! converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "base isothermal layer" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << p_want << std::endl;
            std::cout << drho << std::endl;
            amrex::Error("Error: HSE non-convergence");
        }


        // call the EOS one more time for this zone and then go on to the next
        // (t, rho) -> (p)
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone
        model_hse(i, idens) = dens_zone;
        model_hse(i, itemp) = temp_zone;
        model_hse(i, ipres) = pres_zone;

    }

    
    // auto deltastr = num_to_unitstring(delta);
    auto dxstr = num_to_unitstring(dCoord);

    // std::string outfile = model_prefix + ".hse.tanh.delta_" + deltastr + ".dx_" + dxstr;
    // std::string outfile = model_prefix + "_" + dxstr + ".hse";
    std::string outfile = model_prefix;
    if (! do_isothermal_instead) {
        std::stringstream ss;
        ss << std::setprecision(4) << F_rad;
        outfile += "_Frad_" + ss.str();
    }
    outfile += "_" + dxstr + ".hse";

    std::string outfile2 = outfile + ".extras";

    std::ofstream of(outfile);
    std::ofstream of2(outfile2);

    of << "# npts = " << nx << std::endl;
    of << "# num of variables = " << nvar << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < nx; ++i) {
        of << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, idens) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, itemp) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, ipres) << " ";
        for (int n = 0; n < NumSpec; ++n) {
            of << std::setprecision(12) << std::setw(20) << model_hse(i, ispec+n) << " ";
        }
        of << std::endl;
    }

    of.close();

    std::cout << "Saved model to " << outfile << std::endl; 

    // some metadata
    of << "# generated by toy_atm_layered" << std::endl;
    //of << "# inputs file: " << params_file << std::endl;

    // extras file

    of2 << "# npts = " << nx << std::endl;
    of2 << "# num of variables = " << 2 << std::endl;
    of2 << "# entropy" << std::endl;
    of2 << "# c_s" << std::endl;

    for (int i = 0; i < nx; ++i) {
        eos_state.rho = model_hse(i, idens);
        eos_state.T = model_hse(i, itemp);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model_hse(i, ispec+n);
        }

        eos(eos_input_rt, eos_state);

        of2 << std::setprecision(12) << std::setw(20) << xzn_hse(i) << std::endl;
        of2 << std::setprecision(12) << std::setw(20) << eos_state.s << std::endl;
        of2 << std::setprecision(12) << std::setw(20) << eos_state.cs << std::endl;                
    }

    // compute the maximum HSE error
  
    Real max_hse_error = -1.e30;

    for (int i = 1; i < nx-1; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        // compute the gravitation acceleration at the lower edge
     
        Real g_zone;
        if (do_invsq_grav == 1) {
            g_zone = -C::Gconst * M_enclosed / std::pow(xznl_hse(i), 2);
        } else {
            g_zone = g_const;
        }

        Real dpdr = (model_hse(i, ipres) - model_hse(i-1, ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, idens) + model_hse(i-1, idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, idens) > low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
  
}
