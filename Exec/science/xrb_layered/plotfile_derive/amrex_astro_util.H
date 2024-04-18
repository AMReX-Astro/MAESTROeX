#ifndef AMREX_ASTRO_UTIL_H
#define AMREX_ASTRO_UTIL_H

#include <iostream>
#include <regex>
#include <string>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

#include <network.H>

using namespace amrex;

///
/// Gets the variable ``varname`` from the ``job_info`` file and returns as a
/// string
///
inline
std::string GetVarFromJobInfo (const std::string pltfile, const std::string varname) {
    std::string filename = pltfile + "/job_info";
    std::regex re("\\s*" + varname + "\\s*=\\s*([\\w\\d\\.\\-\\+]+)");

    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
        std::stringstream buf;
        buf << jobfile.rdbuf();
        std::string file_contents = buf.str();

        if (std::regex_search(file_contents, m, re)) {
            return m[1];
        } else {
            Print() << "Unable to find " << varname << " in job_info file!" << std::endl;
        }
    } else {
        Print() << "Could not open job_info file!" << std::endl;
    }

    return "";
}

///
/// Get the center from the job info file and return as a Real Vector
///
inline
Vector<Real> GetCenter (const std::string pltfile) {
    auto center_str = GetVarFromJobInfo(pltfile, "center");

    // split string
    std::istringstream iss {center_str};
    Vector<Real> center;

    std::string s;
    while (std::getline(iss, s, ','))
        center.push_back(stod(s));

    return center;
}

///
/// return the radial coordinate of a zone from the center and the
/// volume of the zone
///
inline
std::pair<Real, Real> get_coord_info(const Array<Real, AMREX_SPACEDIM>& p,
                                     const Vector<Real>& center,
                                     const Array<Real, AMREX_SPACEDIM>& dx_level,
                                     const int coord, const bool sphr) {


    // compute radial coord and index

    Real r_zone{0.0_rt};
    Real vol = std::numeric_limits<Real>::quiet_NaN();

#if AMREX_SPACEDIM == 1
    // 1-d spherical geometry / spherical Sedov explosion

    AMREX_ASSERT(coord == 2);

    r_zone = p[0] - center[0];
    Real r_r = problo[0]+static_cast<Real>(i+1)*dx_level[0];
    Real r_l = problo[0]+static_cast<Real>(i)*dx_level[0];
    vol = (4.0_rt/3.0_rt) * M_PI * dx_level[0] *
        (r_r*r_r + r_l*r_r + r_l*r_l);

#elif AMREX_SPACEDIM == 2
    if (sphr) {
        // 2-d axisymmetric geometry / spherical Sedov explosion

        AMREX_ASSERT(coord == 1);

        // axisymmetric V = pi (r_r**2 - r_l**2) * dz
        //                = pi dr * dz * (r_r + r_l)
        //                = 2 pi r dr dz

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));
        vol = 2 * M_PI * p[0] * dx_level[0] * dx_level[1];

    } else {
        // 2-d Cartesian geometry / cylindrical Sedov explosion

        AMREX_ASSERT(coord == 0);

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));
        vol = dx_level[0] * dx_level[1];

    }

#else

    AMREX_ASSERT(coord == 0);

    vol = dx_level[0] * dx_level[1] * dx_level[2];

    if (sphr) {
        // 3-d Cartesian geometry / spherical Sedov explosion

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]) +
                           (p[2] - center[2]) * (p[2] - center[2]));

    } else {
        // 3-d Cartesian geometry / cylindrical Sedov explosion

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));

    }
#endif

    return {r_zone, vol};

}


///
/// return the index of the density variable by searching through
/// the list of variables in the plotfile.
///
inline int
get_dens_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "density");
    if (idx == var_names_pf.cend()) {
        // MAESTROeX uses "rho"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "rho");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the density component");
        }
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

///
/// return the index of the temperature variable by searching through
/// the list of variables in the plotfile.
///
inline int
get_temp_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "Temp");
    if (idx == var_names_pf.cend()) {
        // for MAESTROeX, we'll default to "tfromp"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "tfromp");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the temperature component");
        }
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

///
/// return the index of the pressure variable by searching through
/// the list of variables in the plotfile.
///
inline int
get_pres_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "pressure");
    if (idx == var_names_pf.cend()) {
        // for MAESTROeX, we'll default to "p0"
        idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "p0");
        if (idx == var_names_pf.cend()) {
            amrex::Error("Error: could not find the pressure component");
        }
    }
    return std::distance(var_names_pf.cbegin(), idx);
}


///
/// return the index of the first species component by searching through
/// the list of variables in the plotfile.
///
inline int
get_spec_index(const std::vector<std::string>& var_names_pf) {

    std::string first_spec_name = "X(" + short_spec_names_cxx[0] + ")";
    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), first_spec_name);
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find the first species");
    }
    int spec_comp = std::distance(var_names_pf.cbegin(), idx);

    // safety check -- make sure the species in the plotfile are identical to
    // those defined in the network we built this tool with.

    for (int n = 0; n < NumSpec; ++n) {
        std::string current_spec_name = "X(" + short_spec_names_cxx[n] + ")";
        if (current_spec_name != var_names_pf[spec_comp+n]) {
            std::cout << current_spec_name << std::endl;
            std::cout << var_names_pf[spec_comp+n] << std::endl;
            std::cout << "Make sure to compile with the same network as the plotfile" << std::endl;
            std::cout << "    make NETWORK_DIR=aprox13   (for example)" << std::endl;
            amrex::Error("Error: species don't match");
        }
    }
    return spec_comp;
}

#endif


