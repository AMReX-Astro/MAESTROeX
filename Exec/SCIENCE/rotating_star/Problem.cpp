#include "Maestro.H"

#ifdef DO_PROBLEM_POST_TIMESTEP

using namespace amrex;

void
Maestro::problem_post_timestep() {

    static const std::string script_name = "post_to_slack.sh";

    if ( (plot_int > 0 && istep % plot_int == 0) ||
         (plot_deltat > 0 && std::fmod(t_new, plot_deltat) < dt) ||
         (istep == max_step) || (t_old >= stop_time) )
    {
        // wrote a plotfile
        std::string command = "bash " + script_name + ' ' + t_new + ' ' + istep + ' ' + plot_base_name;

        system(command.c_str());
    }

    if ( (small_plot_int > 0 && istep % small_plot_int == 0) ||
         (small_plot_deltat > 0 && std::fmod(t_new, small_plot_deltat) < dt) ||
         (istep == max_step)  || (t_old >= stop_time) )
    {
        // wrote a small plotfile
        std::string command = "bash " + script_name + ' ' + t_new + ' ' + istep + ' ' + small_plot_base_name;

        system(command.c_str());
    }
}

#endif
