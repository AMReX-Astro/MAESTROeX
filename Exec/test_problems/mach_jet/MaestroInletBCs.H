#ifndef _MaestroInletBCs_H_
#define _MaestroInletBCs_H_

struct InletBCs
{
    amrex::Real INLET_RHO;
    amrex::Real INLET_RHOH;
    amrex::Real INLET_TEMP;
    amrex::Real INLET_CS;
};

#endif