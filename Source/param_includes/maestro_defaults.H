
// This file is automatically created by parse_maestro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

std::string Maestro::model_file = "model.hse";
int         Maestro::perturb_model = 0;
int         Maestro::print_init_hse_diag = 0;
amrex::Real Maestro::stop_time = -1.0;
int         Maestro::max_step = 1;
amrex::Real Maestro::cflfac = 0.8;
int         Maestro::regrid_int = -1;
int         Maestro::plot_int = 10;
std::string Maestro::plot_base_name = "plt";
int         Maestro::do_reflux = 1;
amrex::Real Maestro::small_temp = 5.e6;
amrex::Real Maestro::small_dens = 1.e-5;
int         Maestro::use_tfromp = 0;
int         Maestro::use_eos_e_instead_of_h = 0;
int         Maestro::use_pprime_in_tfromp = 0;
