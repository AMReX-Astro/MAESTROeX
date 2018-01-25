
#include <Maestro.H>

using namespace amrex;

void
Maestro::PutInPertForm() {

}

void
Maestro::ConvertRhoXToX(Vector<MultiFab>& scal,
                        bool flag)
{
    if (flag) {
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                MultiFab::Divide(scal[lev],scal[lev],Rho,comp,1,0);
            }
        }
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            for (int comp=FirstSpec; comp<FirstSpec+NumSpec; ++comp) {
                MultiFab::Multiply(scal[lev],scal[lev],Rho,comp,1,0);
            }
        }
    }

    // average down data
    AverageDown(scal,FirstSpec,NumSpec);
}

void
Maestro::ConvertRhoHToH(Vector<MultiFab>& scal,
                        bool flag)
{
    if (flag) {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Divide(scal[lev],scal[lev],Rho,RhoH,1,0);
        }
    }
    else {
        for (int lev=0; lev<=finest_level; ++lev) {
            MultiFab::Multiply(scal[lev],scal[lev],Rho,RhoH,1,0);
        }
    }

    // average down data
    AverageDown(scal,RhoH,1);
}
