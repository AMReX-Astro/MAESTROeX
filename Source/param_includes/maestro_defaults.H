
// This file is automatically created by parse_maestro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

amrex::Real Maestro::stop_time = -1.0;
int         Maestro::max_step = 1;
amrex::Real Maestro::cflfac = 0.8;
int         Maestro::regrid_int = -1;
int         Maestro::plot_int = 10;
std::string Maestro::plot_base_name = "plt";
int         Maestro::init_iter = 4;
int         Maestro::init_divu_iter = 4;
int         Maestro::restart = -1;
int         Maestro::restart_into_finer = 0;
int         Maestro::do_initial_projection = 1;
int         Maestro::do_reflux = 1;
