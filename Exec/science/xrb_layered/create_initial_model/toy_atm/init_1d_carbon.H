// See the header of init_1d.H for the full explanation
// The difference here is that carbon is included in the helium layer
// following a power law Zc(y) = (y/y0)^del
// Then the helium fraction is Y = 1-X-Zc-Zcno
// Zcno is still just made up of o14/o15
//
//         X
//         ^
//         |  star  | reg 1 | reg 2  |  region 3
//    1.0  +--------*  -------- (he4)
//         |  ash   . /        \
//         |        . \         \    +--------- X0 (hydrogen)
//         |        .  \         \  /
//         |        .   \ (c12)   \/
//         |        .    \        /\
//         |        .     \      /  \
//         |        .      \    /    +--------- 1-X0-Zcno (helium)
//         |        .----------/--------------- Zcno (o14/o15)
//    0.0  +--------+---------+-----------------> r
//         |<H_star>|
//         |------<y_d>------|
//         |
//

#include <AMReX_Array.H>

#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <model_util.H>

using namespace amrex;

AMREX_INLINE void init_1d_carbon()
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
    int ic12 = network_spec_index("carbon-12");
    int io14 = network_spec_index("oxygen-14");
    int io15 = network_spec_index("oxygen-15");
    int iash = network_spec_index(problem_rp::ash_name);

    Real xn[NumSpec] = {0.0_rt};

    // set the composition of the underlying star

    Real xn_star[NumSpec] = {smallx};
    xn_star[iash] = 1.0;

    // the composition of the first layer above the star
    // First we'll set the cno to the correct relative abundance of oxygen,
    // the prescribed carbon power law, and the rest will be helium
    // Later, we'll replace some of the helium/carbon with hydrogen according
    // to column depth, or pressure

    Real xn_base[NumSpec] = {smallx};
    xn_base[io15] = problem_rp::Zcno / (1.0 + 70.621/122.266);
    xn_base[io14] = problem_rp::Zcno - xn_base[io15];
    xn_base[ic12] = pow(problem_rp::y1/problem_rp::carbon_y, problem_rp::carbon_del);
    xn_base[ihe4] = 1.0 - xn_base[ic12] - problem_rp::Zcno;

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
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    // We know the pressure (from y1) and the temperature
    Real pres_base = -problem_rp::g_const * problem_rp::y1;

    // However there is no (p,T) input option, so we need to iterate
    // to find the base density.
    // We know that it is roughly 1.3-1.5e6 (Initial_Model notebook)

    Real err_min = 1.0e99;
    Real n_best = 0;
    for (int n=0; n<200; ++n) {
        eos_state.rho = (1.3 + n*0.001)*1.e6;
        eos(eos_input_rt, eos_state);

        Real error = std::abs(eos_state.p - pres_base);
        if (error < err_min) {
            err_min = error;
            n_best = n;
        }
    }

    Real dens_base = (1.3 + n_best*0.001)*1.e6;

    std::cout << "Best fitting base density = " << dens_base << std::endl;

    eos_state.rho = dens_base;
    eos(eos_input_rt, eos_state);

    pres_base = eos_state.p;


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
        model_hse(i, idens) = dens_base;
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


    // make the base thermodynamics consistent for this base point -- that is
    // what we will integrate from!

    eos_state.rho = model_hse(index_base, idens);
    eos_state.T = model_hse(index_base, itemp);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_hse(index_base, ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model_hse(index_base, ipres) = eos_state.p;

    Real y_base = -1.*eos_state.p/problem_rp::g_const;
    std::cout << "Recovered base column depth y8 = " << y_base/1.e8 << std::endl;


    // HSE solve

    // the HSE state will be done by fixing the temperature according to the power-law
    // then iterating until the pressure from HSE and from EOS agree

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real y_zone;

    Real p_want;
    Real y_want;
    Real drho;
    Real dtemp;

    // Keep track of which power-law we are currently in
    // We know which one is first
    Real current_T0 = eos_state.T;
    Real current_y0 = y_base;
    Real current_del = problem_rp::del1;

    bool cross_y2 = false;
    bool cross_y3 = false;


    // Integrate up to the top of the model

    for (int i = index_base+1; i < problem_rp::nx; ++i) {

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

                // if here, assume we're at low column and X is constant, and no carbon
                xn[ih1] = problem_rp::X0;
                xn[ihe4] = 1.0 - problem_rp::X0 - problem_rp::Zcno;
                xn[ic12] = 0.0;

            } else {

                p_want = model_hse(i-1, ipres) +
                    delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * problem_rp::g_const;

                // The column depth determines both the temperature and hydrogen fraction
                y_want = -p_want/problem_rp::g_const;
                // std::cout << y_want << " " << current_del << std::endl;

                temp_zone = current_T0 * pow(y_want/current_y0, current_del);

                if (y_want > problem_rp::y2) {
                    xn[ic12] = pow(y_want/problem_rp::carbon_y, problem_rp::carbon_del);
                } else {
                    xn[ic12] = 0.0;
                }

                xn[ih1] = problem_rp::X0 * (1 - y_want/problem_rp::yd);
                xn[ihe4] = 1.0 - problem_rp::Zcno - xn[ic12] - xn[ih1];
                if (xn[ih1] < 0) {
                    xn[ih1] = 0;
                    xn[ihe4] = 1.0 - problem_rp::Zcno - xn[ic12];
                }

                // set to isothermal if hit threshold
                if (temp_zone < problem_rp::T_lo) {
                    temp_zone = problem_rp::T_lo;
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

            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * problem_rp::g_const);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                    amrex::min(dens_zone + drho, 1.1_rt * dens_zone));


            if (std::abs(drho) < TOL * dens_zone) {
                converged_hse = true;
                break;
            }

            if (dens_zone <= problem_rp::low_density_cutoff) {
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::T_lo;
                converged_hse = true;
                if (!fluff) {
                    std::cout << "Hit fluff at x=" << xzn_hse(i) << ", y=" << y_want << std::endl;
                }
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

        // Update the power-law
        y_zone = -1*pres_zone/problem_rp::g_const;

        if ( !cross_y2 and y_zone < problem_rp::y2 ) {
            current_T0 = temp_zone;
            current_y0 = y_zone;
            current_del = problem_rp::del2;
            cross_y2 = true;

            std::cout << "crossed y2 at x=" << xzn_hse(i) << std::endl;

        } else if ( !cross_y3 and y_zone < problem_rp::y3 ) {
            current_T0 = temp_zone;
            current_y0 = y_zone;
            current_del = problem_rp::del3;
            cross_y3 = true;

            std::cout << "crossed y3 at x=" << xzn_hse(i) << std::endl;
        }

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
