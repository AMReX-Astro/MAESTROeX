#include <Radial.H>
#include <Maestro.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;


// ---------------------------------
// Write 1D radial diagnostics
// ---------------------------------
void test ()
{
    // exact-solution problem to test subroutines that compute
    // diagnostics dependent on velocity only;
    
    // make BoxArray and Geometry
    int n_cell = 128;
    int max_grid_size = 32;
    BoxArray ba;
    Geometry tgeom;
    {
	IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
	IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
	Box domain(dom_lo, dom_hi);

	ba.define(domain);
	ba.maxSize(max_grid_size);

	// define physical box, [-100, 100] in each direction
	RealBox real_box({AMREX_D_DECL(-100.0_rt,-100.0_rt,-100.0_rt)},
		         {AMREX_D_DECL( 100.0_rt, 100.0_rt, 100.0_rt)});

	tgeom.define(domain, real_box);
    }
    DistributionMapping dm(ba);

    // define density and pressure
    MultiFab rho0(ba, dm,       1, 0);
    MultiFab p0  (ba, dm,       1, 0);
    MultiFab rhoX(ba, dm, NumSpec, 0);
    rho0.setVal(1.);
    p0.setVal(1.);
    rhoX.setVal(0.);
    rhoX.setVal(1., 0);  // first species: rhoX = 1
    
    // define velocities
    MultiFab u_mf(ba, dm, AMREX_SPACEDIM, 0);
    MultiFab w0_mf(ba, dm, AMREX_SPACEDIM, 0);
    w0_mf.setVal(0.);

    const auto& center_p = center;

    const auto dx = pgeom[lev].CellSizeArray();
    const auto prob_lo = pgeom[lev].ProbLoArray();

    for ( MFIter mfi(vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

	// Get the index space of the valid region
	const Box& tileBox = mfi.tilebox();
	
	const Array4<Real> vel_arr = u_mf.array(mfi);
	
	AMREX_PARALLEL_FOR_3D(tileBox, i, j, k, {	
	    Real x = prob_lo[0] + (Real(i)+0.5) * dx[0] - center_p[0];
	    Real y = prob_lo[1] + (Real(j)+0.5) * dx[1] - center_p[1];
	    Real z = prob_lo[2] + (Real(k)+0.5) * dx[2] - center_p[2];
	    Real radius = std::sqrt(x*x + y*y + z*z);

	    if (radius < 50.0) {
		vel_arr(i,j,k) = radius;
	    } else {
		vel_arr(i,j,k) = 50.0*50.0/radius;
	    }
	});
    }        

    
    // put all variables into a single multifab
    MultiFab s0(ba, dm, 2 + NumSpec + 2*AMREX_SPACEDIM, 0);
    MultiFab::Copy(s0, rho0, 0, 0, 1, 0);
    MultiFab::Copy(s0, p0,   0, 1, 1, 0);
    MultiFab::Copy(s0, rhoX, 0,         2,        NumSpec, 0);
    MultiFab::Copy(s0, u_mf, 0, 2+NumSpec, AMREX_SPACEDIM, 0);
    MultiFab::Copy(s0, w0_mf,0, 2+NumSpec+AMREX_SPACEDIM, AMREX_SPACEDIM, 0);

    
    // variable names for plotfile
    Vector<std::string> varnames(2 + NumSpec + 2*AMREX_SPACEDIM);
    int cnt = 0;
    varnames[cnt++] = "rho";
    varnames[cnt++] = "p0";

    for (int i = 0; i < NumSpec; i++) {
	int len = 20;
	Vector<int> int_spec_names(len);
	//
	// This call return the actual length of each string in "len"
	//
	get_spec_names(int_spec_names.dataPtr(), &i, &len);
	auto* spec_name = new char[len + 1];
	for (int j = 0; j < len; j++) {
	    spec_name[j] = int_spec_names[j];
	}
	spec_name[len] = '\0';
	std::string spec_string = "rhoX(";
	spec_string += spec_name;
	spec_string += ')';

	varnames[cnt++] = spec_string;
	delete[] spec_name;
    }
    
    // add velocities
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	std::string x = "vel";
	x += (120 + i);
	varnames[cnt++] = x;
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	std::string x = "w0";
	x += (120 + i);
	varnames[cnt++] = x;
    }

    // write plotfile
    std::string testfilename = "test_plt";
    
    WriteSingleLevelPlotfile(testfilename,
			     s0, varnames,
			     tgeom, 0, 0);

    WriteJobInfo(testfilename);
}

