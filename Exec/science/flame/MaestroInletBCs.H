#ifndef MaestroInletBCs_H
#define MaestroInletBCs_H

#include <AMReX_BLFort.H>

namespace InletBCs {

extern amrex::Real INLET_RHO;
extern amrex::Real INLET_RHOH;
extern amrex::Real INLET_TEMP;
extern RealVector INLET_RHOX;
extern amrex::Real INLET_VEL;
};  // namespace InletBCs

#endif
