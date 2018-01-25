
#include <Maestro.H>

using namespace amrex;

void
Maestro::PutInPertForm(Vector<MultiFab>& scal, 
		       const Vector<Real>& s0, 
		       const int& comp, 
		       bool flag) {
    // place 1d array onto a cartesian grid
    Vector<MultiFab> s0_cart(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	s0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
    }
    
    // s0 is not edge centered
    // note that bcs parameter is not used
    Put1dArrayOnCart(s0,s0_cart,0,0);

    if (flag) {
	for (int lev = 0; lev <= finest_level; ++lev) {
	    MultiFab::Subtract(scal[lev],s0_cart[lev],0,comp,1,0);
	}
    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
	    MultiFab::Add(scal[lev],s0_cart[lev],0,comp,1,0);
	}
    }

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
