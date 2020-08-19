//
// Process a plotfile to produce diagnostics for rotational problems.
//
#include <fstream>
#include <iostream>

#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include <Postprocess.H>

using namespace amrex;

std::string inputs_name = "";


int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	Postprocess postproc;

	postproc.init();

	postproc.diag();
	
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}

