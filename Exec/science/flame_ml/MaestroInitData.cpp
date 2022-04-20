
#include <Maestro.H>
using namespace amrex;

// initializes data on a specific level
void Maestro::InitLevelData(const int lev, const Real time, const MFIter& mfi,
                            const Array4<Real> scal, const Array4<Real> vel) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::InitLevelData()", InitLevelData);

    const auto tileBox = mfi.tilebox();

    const auto vel_fuel_loc = vel_fuel;

    // initialize velocity
    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        for (auto n = 0; n < AMREX_SPACEDIM - 1; ++n) {
            vel(i, j, k, n) = 0.0;
        }
        vel(i, j, k, AMREX_SPACEDIM - 1) = vel_fuel_loc;
    });
    
    const auto s0_arr = s0_init.const_array();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
	int r = AMREX_SPACEDIM == 2 ? j : k;
	if (AMREX_SPACEDIM == 2) {
	    r = i < 64 ? max(j-i, 0) : max(j-127+i,0);
	}

        // set the scalars using s0
        scal(i, j, k, Rho) = s0_arr(lev, r, Rho);
        scal(i, j, k, RhoH) = s0_arr(lev, r, RhoH);
        scal(i, j, k, Temp) = s0_arr(lev, r, Temp);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) = s0_arr(lev, r, FirstSpec + comp);
        }
        // initialize pi to zero for now
        scal(i, j, k, Pi) = 0.0;
    });
    
    /*
    ///  COPIED FROM MaestroBaseState.CPP  ///
    // figure out the indices for different species
    const int ic12 = network_spec_index("carbon-12");
    const int io16 = network_spec_index("oxygen-16");
    const int img24 = network_spec_index("magnesium-24");

    if (ic12 < 0 || io16 < 0 || img24 < 0) {
        Abort("ERROR: species indices not defined");
    }

    // length of the domain
    const Real rlen = geom[lev].ProbHi(AMREX_SPACEDIM - 1) -
		      geom[lev].ProbLo(AMREX_SPACEDIM - 1);

    // width of the domain
    const Real rwid = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);

    // figure out the thermodynamics of the fuel and ash state
    RealVector xn_fuel(NumSpec, 0.0);
    RealVector xn_ash(NumSpec, 0.0);

    eos_t eos_state;

    // fuel
    xn_fuel[ic12] = xc12_fuel;
    xn_fuel[io16] = 1.0 - xc12_fuel;

    eos_state.rho = dens_fuel;
    eos_state.T = temp_fuel;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_fuel[comp];
    }

    eos(eos_input_rt, eos_state);

    Real p_ambient = eos_state.p;
    Real rhoh_fuel = dens_fuel * eos_state.h;

    // ash
    xn_ash[io16] = 1.0 - xc12_fuel;
    xn_ash[img24] = xc12_fuel;

    eos_state.rho = dens_fuel;  // initial guess
    eos_state.T = temp_ash;
    for (auto comp = 0; comp < NumSpec; ++comp) {
        eos_state.xn[comp] = xn_ash[comp];
    }
    eos_state.p = p_ambient;

    eos(eos_input_tp, eos_state);

    Real dens_ash = eos_state.rho;
    Real rhoh_ash = dens_ash * eos_state.h;

    const auto dx = geom[lev].CellSizeArray();

    ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
	Real x = (Real(i) + 0.5) * dx[0];
	Real y_offset = 64 * sin(x * M_PI / rwid);
	
	// height above the bottom of the domain
        Real rloc = geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
	    (Real(j) + 0.5 - y_offset) * dx[AMREX_SPACEDIM - 1];

        if (rloc <
            geom[lev].ProbLo(AMREX_SPACEDIM - 1) + interface_pos_frac * rlen) {
            // fuel
            scal(i, j, k, Rho) = dens_fuel;
            scal(i, j, k, RhoH) = rhoh_fuel;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) = dens_fuel * xn_fuel[comp];
            }

        } else {
            // ash
            scal(i, j, k, Rho) = dens_ash;
            scal(i, j, k, RhoH) = rhoh_ash;
            for (auto comp = 0; comp < NumSpec; ++comp) {
                scal(i, j, k, FirstSpec + comp) = dens_fuel * xn_ash[comp];
            }
        }

        // give the temperature a smooth profile
        scal(i, j, k, Temp) =
            temp_fuel +
            (temp_ash - temp_fuel) * 0.5 *
                (1.0 + tanh((rloc - (geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                                     interface_pos_frac * rlen)) /
                            (smooth_len_frac * rlen)));

        // give the carbon mass fraction a smooth profile too
        RealVector xn_smooth(NumSpec, 0.0);

        xn_smooth[ic12] =
            xn_fuel[ic12] +
            (xn_ash[ic12] - xn_fuel[ic12]) * 0.5 *
                (1.0 + tanh((rloc - (geom[lev].ProbLo(AMREX_SPACEDIM - 1) +
                                     interface_pos_frac * rlen)) /
                            (smooth_len_frac * rlen)));

        xn_smooth[io16] = xn_fuel[io16];
        xn_smooth[img24] = 1.0 - xn_smooth[ic12] - xn_smooth[io16];

        // get the new density and enthalpy
	eos_t eos_state_new;
        eos_state_new.rho = scal(i, j, k, Rho);
        eos_state_new.T = scal(i, j, k, Temp);
        eos_state_new.p = p_ambient;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            eos_state_new.xn[comp] = xn_smooth[comp];
        }

        // (T,p) --> rho, h
        eos(eos_input_tp, eos_state_new);

        scal(i, j, k, Rho) = eos_state_new.rho;
        scal(i, j, k, RhoH) = eos_state_new.rho * eos_state_new.h;
        for (auto comp = 0; comp < NumSpec; ++comp) {
            scal(i, j, k, FirstSpec + comp) =
                eos_state_new.rho * xn_smooth[comp];
        }
    });
    */
}

void Maestro::InitLevelDataSphr(const int lev, const Real time, MultiFab& scal,
                                MultiFab& vel) {
    Abort("InitLevelDataSphr not implemented");
}
