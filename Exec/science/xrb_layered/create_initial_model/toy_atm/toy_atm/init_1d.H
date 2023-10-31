//  Adapted from original AMReX-Astro code
//  (github.com/AMReX-Astro/initial_models/tree/main/toy_atm/init_1d.H)
//  by Simon Guichandut (2023).

//  Create a 1-d hydrostatic, atmosphere with an isothermal region
//  (T_star) representing the NS, a hyperbolic tangent rise to a
//  peak temperature (T_base) representing the base of an accreted
//  layer at r=H_star, an isentropic profile representing the convective
//  region at ignition, then a radiative profile from r=H_rad. 

//  The composition transitions from pure ash, to pure helium, to a solar 
//  and CNO mixture, with initial hydrogen fraction X0 and CNO fraction
//  Zcno (inputs file).  The CNO is split between o14 and o15 with ratio
//  determined by their beta decay times (70.606/122.24). Because hot CNO burns
//  hydrogen at a constant rate:
//  X = X0(1-y/yd),
//  where yd is the depletion column depth, another input parameter.
//  
//  In the radiative layer,
//  dF/dy = -eps,
//  where eps is the constant hot CNO burning rate.

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
//         |        .           \    +--------- X0 (hydrogen)
//         |        .            \  /
//         |        .             \/
//         |        .             /\
//         |        .            /  \
//         |        .           /    +--------- 1-X0-Zcno (helium)
//         |        .          /     +--------- Zcno (o14/o15)
//    0.0  +--------+---------+-----------------> r
//         |<H_star>|    
//         |---<H_rad>----| 
//         |------<y_d>-----|     
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
    int ih1 = network_spec_index("hydrogen-1");
    int ihe4 = network_spec_index("helium-4");
    int io14 = network_spec_index("oxygen-14");
    int io15 = network_spec_index("oxygen-15");
    int iash = network_spec_index(problem_rp::ash_name);

    Real xn[NumSpec] = {0.0_rt};

    // set the composition of the underlying star

    Real xn_star[NumSpec] = {smallx};
    xn_star[iash] = 1.0;

    // the composition of the first layer above the star
    // First we'll set the cno to the correct relative abundance, and the rest will be helium
    // Later, we'll replace some of the helium with hydrogen according
    // to column depth, or pressure

    Real xn_base[NumSpec] = {smallx};
    xn_base[io15] = problem_rp::Zcno / (1.0 + 70.606/122.24);
    xn_base[io14] = problem_rp::Zcno - xn_base[io15];
    xn_base[ihe4] = 1.0 - problem_rp::Zcno;

    // check if they sum to 1
    double sum{0.0};
    for (auto e : xn_base) {
        sum += e;
    }
    if (std::abs(sum) - 1.0_rt > NumSpec * smallx) {
            amrex::Error("ERROR: fuel mass fractions don't sum to 1");
    }   



    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xznl_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznr_hse;

    Array2D<Real, 0, NPTS_MODEL-1, 0, nvar-1> model_hse;

    // compute the coordinates of the new gridded function

    Real dCoord = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    for (int i = 0; i < problem_rp::nx; ++i) {
        xznl_hse(i) = problem_rp::xmin + (static_cast<Real>(i)) * dCoord;
        xzn_hse(i) = problem_rp::xmin + (static_cast<Real>(i) + 0.5_rt) * dCoord;
        xznr_hse(i) = problem_rp::xmin + (static_cast<Real>(i) + 1.0_rt) * dCoord;
    }

    // put the model onto our new uniform grid

    bool fluff = false;

    // determine the conditions at the base

    eos_t eos_state;
    eos_state.T = problem_rp::T_base;
    eos_state.rho = problem_rp::dens_base;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    eos(eos_input_rt, eos_state);

    // store the conditions at the base -- we'll use the entropy later
    // to constrain the isentropic layer
  
    Real pres_base = eos_state.p;
    Real entropy_base = eos_state.s;


    // set an initial temperature profile and composition

    for (int i = 0; i < problem_rp::nx; ++i) {

        // Make hyperbolic tangent transition between star and base

        // Composition
        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i, ispec+n) = xn_star[n] + 
                0.5_rt * (xn_base[n] - xn_star[n]) *
                (1.0_rt + std::tanh((xzn_hse(i) - (problem_rp::xmin + problem_rp::H_star - problem_rp::delta) + problem_rp::delta) / problem_rp::delta)); 
        } 

        // Temperature
        model_hse(i, itemp) = problem_rp::T_star + 0.5_rt * (problem_rp::T_base - problem_rp::T_star) *
            (1.0_rt + std::tanh((xzn_hse(i) - (problem_rp::xmin + problem_rp::H_star - problem_rp::delta) + problem_rp::delta) / problem_rp::delta));


        // the density and pressure will be determined via HSE,
        // for now, set them to the base conditions
        model_hse(i, idens) = problem_rp::dens_base;
        model_hse(i, ipres) = pres_base;
    }


    // find the index of the base height
    int index_base = -1;

    if (problem_rp::index_base_from_temp == 1) {
        for (int i = 0; i < problem_rp::nx; ++i) {
            if (model_hse(i, itemp) > 0.9995 * problem_rp::T_base) {
                index_base = i+1;
                break;  
            }
        }
    } else {
        for (int i = 0; i < problem_rp::nx; ++i) {
            if (xzn_hse(i) >= problem_rp::xmin + problem_rp::H_star + problem_rp::delta) {
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

    if (problem_rp::H_rad < problem_rp::H_star) {
        amrex::Error("ERROR: H_rad<H_star. Adiabatic (base) layer should be below radiative layer");
    }

    for (int i = index_base; i < problem_rp::nx; ++i) {
        if (xzn_hse(i) >= problem_rp::xmin + problem_rp::H_rad) {
            index_rad = i;
            break;
        }
    }

    if (index_rad == -1) {
        amrex::Error("ERROR: rad_height not found on grid");
    }

    std::cout << "index_rad = " << index_rad << "; x=" << xzn_hse(index_rad) << std::endl;


    // make the base thermodynamics consistent for this base point -- that is
    // what we will integrate from!

    eos_state.rho = model_hse(index_base, idens);
    eos_state.T = model_hse(index_base, itemp);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_hse(index_base, ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model_hse(index_base, ipres) = eos_state.p;

    std::cout << "base column depth y8 = " << -1.*eos_state.p/problem_rp::g_const/1.e8 << std::endl;


    // HSE solve

    // the HSE state will be done putting creating an isentropic state until
    // reaching H_rad. Then, we set the radiative profile using the flux

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
            delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * problem_rp::g_const;

            // pressure also sets the hydrogen fraction
            // (+ P/g/yd because g is defined negative)
            xn[ih1] = problem_rp::X0 * (1 + p_want/problem_rp::g_const/problem_rp::yd);
            xn[ihe4] = 1.0 - problem_rp::Zcno - xn[ih1];
            if (xn[ih1]<0) {
                xn[ih1] = 0;
                xn[ihe4] = 1.0 - problem_rp::Zcno;
            }


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

            dtemp = ((dsd / (dpd - 0.5_rt * delx * problem_rp::g_const)) * A - B)/
                (dsd * dpt / (dpd - 0.5_rt * delx * problem_rp::g_const) - dst);

            drho = (A - dpt * dtemp) / (dpd - 0.5_rt * delx * problem_rp::g_const);

            // Update density and temperature in the zone
            // This assumes both values are going down
            dens_zone = amrex::max(0.9_rt * dens_zone,
                                    amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            temp_zone = amrex::max(0.9_rt * temp_zone,
                                    amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));


            // check if the density falls below our minimum cut-off --
            // if so, floor it
            if (dens_zone < problem_rp::low_density_cutoff) {
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::T_lo;
                converged_hse = true;
                fluff = true;

                // raise a warning because this happens before the radiative part
                std::cout << "\nWarning: reached cut-off density before radiative zone!!" << std::endl;

                break;
            }

            // Stop if converged
            if (std::abs(drho) < TOL * dens_zone &&
                std::abs(dtemp) < TOL*temp_zone) {
                converged_hse = true;
                break;
            }
        
        }

        if (! converged_hse) {
            std::cout << "\nError zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "x = " << xzn_hse(i) << "(isentropic layer)" << std::endl;
            std::cout << "\ndens_zone   temp_zone" << std::endl;
            std::cout << dens_zone << "   " << temp_zone << std::endl;
            std::cout << "\np_want   entropy_base   entropy" << std::endl;
            std::cout << p_want << " " << entropy_base << " " << entropy << std::endl;
            std::cout << "\nX(h1)   X(he4)   X(ash)" << std::endl;
            std::cout << xn[ih1] << " " << xn[ihe4] << " " << xn[iash] << std::endl;
            std::cout << "\ndrho   dT" << std::endl;
            std::cout << drho << " " << dtemp << std::endl;
            amrex::Error("Error: HSE non-convergence\n");
        }

        // call the EOS one more time for this zone and then go on to the next
        // (t, rho) -> (p)
     
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
            model_hse(i, ispec+n) = xn[n];
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

    std::cout << "Reached end of adiabatic region at r=" << xzn_hse(index_rad) <<" T8=" << model_hse(index_rad,itemp)/1.e8 << std::endl;
    std::cout << "radiative column depth y8 = " << -1.*pres_zone/problem_rp::g_const/1.e8 << std::endl;




    // Above is the radiative layer

    // Get the radiative flux at the edge of the adiabatic zone

    int i0 = index_rad;

    Real rad_D = C::a_rad * C::c_light /
        (3 * 0.2 * 0.5 * (
            (1 + model_hse(i0, ispec+Species::H1-1)) * model_hse(i0, idens) +
            (1 + model_hse(i0-1, ispec+Species::H1-1)) * model_hse(i0-1, idens)
        ));

    Real dT4_dr = (pow(model_hse(i0, itemp), 4) - pow(model_hse(i0-1, itemp), 4)) / 
                  (xzn_hse(i0) - xzn_hse(i0-1));

    Real F0 = problem_rp::F0;
    if (F0 == -1) {
        F0 = -rad_D * dT4_dr;
    }
    std::cout << "Flux at bottom of rad zone F22 = " << F0/1.e22 << std::endl;

    Array1D<Real, 0, NPTS_MODEL-1> rad_flux;
    for (int i = 0; i < problem_rp::nx; ++i) {
        rad_flux(i) = F0;
    }


    // Now integrate up to the top of the model
    Real t4 = 0.0;

    for (int i = index_rad; i < problem_rp::nx; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);


        // we've already set initial guesses for density, temperature, and composition
        dens_zone = model_hse(i, idens);
        temp_zone = model_hse(i, itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, ispec+n);
        }

        // iteration loop
        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            if (fluff) {
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::T_lo;

                // if here, assume we're at low column and X is constant
                xn[ih1] = problem_rp::X0;
                xn[ihe4] = 1.0 - problem_rp::X0 - problem_rp::Zcno;

            } else {

                p_want = model_hse(i-1, ipres) +
                    delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * problem_rp::g_const;

                // pressure also sets the hydrogen fraction
                xn[ih1] = problem_rp::X0 * (1 + p_want/problem_rp::g_const/problem_rp::yd);
                xn[ihe4] = 1.0 - problem_rp::Zcno - xn[ih1];
                if (xn[ih1] < 0) {
                    xn[ih1] = 0;
                    xn[ihe4] = 1.0 - problem_rp::Zcno;
                }


                // In optically thick medium, radiation energy density is transported
                // according to a diffusion law:
                // Flux = -D d(T^4)/dx
                // where the diffusion coefficient is D = ac/(3*kappa*rho)
                // kappa is not included in EOS, just use classical electron scattering
                // kap = 0.2*(1+X), where X is the hydrogen fraction
                // Flux increments from constant hot CNO burning, i.e.
                // dF/dy = -5.8e13 (Zcno/0.01) (or 0 if y>yd)

                // std::cout << "hydrogen1 index in the network is:" << Species::H1-1 << std::endl; 

                // Average the diffusion coefficient over two zones
                rad_D = C::a_rad * C::c_light /
                        (3 * 0.2 * 0.5 * (
                            (1 + xn[ih1]) * dens_zone +
                            (1 + model_hse(i-1, ispec+Species::H1-1)) * model_hse(i-1, idens)
                        ));

                // Update the flux
                if (xn[ih1] < smallx) {
                    rad_flux(i) = rad_flux(i-1);
                } else {
                    rad_flux(i) = rad_flux(i-1) + 0.5 * (dens_zone + model_hse(i-1, idens)) * delx *
                                                    5.8e13 * problem_rp::Zcno/0.01;
                }

                // Update T^4
                t4 = pow(model_hse(i-1, itemp), 4) - 0.5 * (rad_flux(i) + rad_flux(i-1)) / rad_D;

                // std::cout << pow(t4,0.25) << std::endl;
                // break;

                if (t4 < pow(problem_rp::T_lo, 4.0_rt)) {
                // if (t4 < 0.0) {
                // if (true) {
                    // amrex::Error("Error: T<0");
                    temp_zone = problem_rp::T_lo; // set to isothermal
                } else {
                    temp_zone = pow(t4, 0.25_rt);
                }
            
                // std::cout << "T=" << temp_zone << std::endl;
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

            drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * problem_rp::g_const);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                    amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            // if (!fluff) {
            //     std::cout << rad_flux(i) << std::endl;
            // }

            if (std::abs(drho) < TOL * dens_zone) {
                converged_hse = true;
                break;
            }

            if (dens_zone < problem_rp::low_density_cutoff) {
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::T_lo;
                converged_hse = true;
                fluff = true;
                break;
            }
        }

        if (! converged_hse) {
            std::cout << "\nError zone " << i << " did not converge in init_1d" << std::endl;
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
            model_hse(i, ispec+n) = xn[n];
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

        std::cout << "top column depth y = " << -1.*pres_zone/problem_rp::g_const << std::endl;


    // integrate down -- using the temperature profile defined above
  
    for (int i = index_base-1; i >= 0; --i) {

        Real delx = xzn_hse(i+1) - xzn_hse(i);

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
                delx * 0.5_rt * (dens_zone + model_hse(i+1, idens)) * problem_rp::g_const;

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

            drho = A / (dpd + 0.5_rt * delx * problem_rp::g_const);

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

    // std::string outfile = problem_rp::model_prefix + ".hse.tanh.delta_" + deltastr + ".dx_" + dxstr;
    std::string outfile = problem_rp::model_prefix + "_" + dxstr + ".hse";

    std::string outfile2 = outfile + ".extras";

    std::ofstream of(outfile);
    std::ofstream of2(outfile2);

    of << "# npts = " << problem_rp::nx << std::endl;
    of << "# num of variables = " << nvar << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < problem_rp::nx; ++i) {
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

    of2 << "# npts = " << problem_rp::nx << std::endl;
    of2 << "# num of variables = " << 2 << std::endl;
    of2 << "# entropy" << std::endl;
    of2 << "# c_s" << std::endl;

    for (int i = 0; i < problem_rp::nx; ++i) {
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

    for (int i = 1; i < problem_rp::nx-1; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        Real dpdr = (model_hse(i, ipres) - model_hse(i-1, ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, idens) + model_hse(i-1, idens)) * problem_rp::g_const;

        if (dpdr != 0.0_rt && model_hse(i+1, idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
  
}