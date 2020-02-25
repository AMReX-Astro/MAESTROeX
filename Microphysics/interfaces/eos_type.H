#ifndef _eos_type_H_
#define _eos_type_H_

#include <AMReX.H>
#include <network.H>

struct eos_t {
    amrex::Real rho;
    amrex::Real T;
    amrex::Real p;
    amrex::Real e;
    amrex::Real h;
    amrex::Real s;
    amrex::Real xn[NumSpec];
    amrex::Real aux[NumAux];

    amrex::Real dpdT;
    amrex::Real dpdr;
    amrex::Real dedT;
    amrex::Real dedr;
    amrex::Real dhdT;
    amrex::Real dhdr;
    amrex::Real dsdT;
    amrex::Real dsdr;
    amrex::Real dpde;
    amrex::Real dpdr_e;

    amrex::Real cv;
    amrex::Real cp;
    amrex::Real xne;
    amrex::Real xnp;
    amrex::Real eta;
    amrex::Real pele;
    amrex::Real ppos;
    amrex::Real mu;
    amrex::Real mu_e;
    amrex::Real y_e;
    amrex::Real gam1;
    amrex::Real cs;

    amrex::Real abar;
    amrex::Real zbar;

#ifdef EXTRA_THERMO
    amrex::Real dpdA;
    amrex::Real dpdZ;
    amrex::Real dedA;
    amrex::Real dedZ;
#endif

    amrex::Real conductivity;

};

enum eos_input_t {eos_input_rt = 0,
                  eos_input_rh,
                  eos_input_tp,
                  eos_input_rp,
                  eos_input_re,
                  eos_input_ps,
                  eos_input_ph,
                  eos_input_th};

// these are used to allow for a generic interface to the
// root finding

enum eos_iter {itemp = 1,
               idens = 2,
               iener = 3,
               ienth = 4,
               ientr = 5,
               ipres = 6};

#endif
