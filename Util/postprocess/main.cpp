//
// Process a plotfile to produce diagnostics for rotational problems.
//
#include <fstream>
#include <iostream>
#include <regex>

#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include <Maestro.H>
#include <Maestro_F.H>
#include <Radial.H>

using namespace amrex;

std::string inputs_name = "";
BaseStateGeometry base_geom;
int finest_level = 0;
Vector<BoxArray> grid;
Vector<DistributionMapping> dmap;
Vector<Geometry> pgeom;
GpuArray<Real,3> center;

// Helper subroutines
int GetdrdxFac (const std::string pltfile);
bool GetOctant (const std::string pltfile);


int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// plotfile names
	std::string iFile;

	// read in parameters from inputs file
	ParmParse pp;
	pp.query("infile", iFile);
	if (iFile.empty()) {
	    Print() << "WARNING: Plotfile was not specified!\n";
	    Print() << "Running exact solution test problem" << std::endl;
	    test();
	    iFile = "test_plt0000000";
	} else {
	    // initialize species
	    maestro_network_init();
	}

	// optional parameters
	std::string imFile;
	std::string deltat, numfiles;
	pp.query("modelfile", imFile); 
	pp.query("dt", deltat);
	pp.query("nfiles", numfiles);
	
	// read input grid 
	amrex::PlotFileData pltfile(iFile);
	finest_level = pltfile.finestLevel();
	const auto probLo = pltfile.probLo();
	const auto probHi = pltfile.probHi();
	
	grid.resize(finest_level+1);
	dmap.resize(finest_level+1);
	pgeom.resize(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    grid[lev] = pltfile.boxArray(lev);
	    dmap[lev] = pltfile.DistributionMap(lev);

	    Box domain = pltfile.probDomain(lev);
	    RealBox real_box(probLo, probHi);
	    pgeom[lev].define(domain,&real_box);
	}
		
	// read input MultiFabs
	Vector<MultiFab> rho_mf(finest_level+1);
	Vector<MultiFab> p0_mf(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    rho_mf[lev] = pltfile.get(lev, "rho");
	    p0_mf[lev] = pltfile.get(lev, "p0");
	}

	Vector<MultiFab> rhoX_mf(finest_level+1);
	Vector<MultiFab> u_mf(finest_level+1);
	Vector<MultiFab> w0_mf(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    rhoX_mf[lev].define(grid[lev], dmap[lev], NumSpec, 0);
	    u_mf[lev]   .define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
	    w0_mf[lev]  .define(grid[lev], dmap[lev], AMREX_SPACEDIM, 0);
	}

	for (int i = 0; i < NumSpec; ++i) {
	    int len = 20;
	    Vector<int> int_spec_names(len);
	    //
	    // This call return the actual length of each string in "len"
	    //
	    get_spec_names(int_spec_names.dataPtr(),&i,&len);
	    char* spec_name = new char[len+1];
	    for (int j = 0; j < len; j++)
		spec_name[j] = int_spec_names[j];
	    spec_name[len] = '\0';
	    std::string spec_string = "rhoX(";
	    spec_string += spec_name;
	    spec_string += ')';

	    for (int lev = 0; lev <= finest_level; ++lev) {
		MultiFab::Copy(rhoX_mf[lev], pltfile.get(lev,spec_string), 0, i, 1, 0);
	    }
	    
	    delete [] spec_name;
	}

	// velocities
	for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	    std::string x = "vel";
	    std::string w = "w0";
	    x += (120+i);
	    w += (120+i);

	    for (int lev = 0; lev <= finest_level; ++lev) {
		MultiFab::Copy(u_mf[lev], pltfile.get(lev,x), 0, i, 1, 0);
		MultiFab::Copy(w0_mf[lev], pltfile.get(lev,w), 0, i, 1, 0);
	    }
	}

	// setup geometry
	const auto domainBoxFine = pltfile.probDomain(finest_level);
	const auto dxFine = pltfile.cellSize(finest_level);
	
	// assume spherical
	int drdxfac = GetdrdxFac(iFile);
	bool octant = GetOctant(iFile);
	base_geom.max_radial_level = 0;
	base_geom.dr_fine = dxFine[0] / drdxfac;
	int domhi = domainBoxFine.bigEnd(0)+1;

	// assume constant dr
	if (!octant) {
	    base_geom.nr_irreg = int(round( (3*(domhi/2-0.5)*(domhi/2-0.5)-0.75)/2.0 ));
	} else {
	    base_geom.nr_irreg = int(round( (3*(domhi-0.5)*(domhi-0.5)-0.75)/2.0 ));
	}

	double lenx, leny, lenz, max_dist;
	if (octant) {
	    lenx = probHi[0] - probLo[0];
	    leny = probHi[1] - probLo[1];
	    lenz = probHi[2] - probLo[2];
	} else {
	    lenx = 0.5*(probHi[0] - probLo[0]);
	    leny = 0.5*(probHi[1] - probLo[1]);
	    lenz = 0.5*(probHi[2] - probLo[2]);
	}
	max_dist = sqrt(lenx*lenx + leny*leny + lenz*lenz);
	base_geom.nr_fine = int(max_dist / base_geom.dr_fine) + 1;
	
	// Print() << "Check: " << NumSpec << "; " << base_geom.dr_fine << ", " << base_geom.nr_fine << std::endl;
	base_geom.Init(base_geom.max_radial_level, base_geom.nr_fine, base_geom.dr_fine, base_geom.nr_irreg, pgeom, finest_level, center);

	// Radial states
	BaseState<Real> rho0(base_geom.max_radial_level+1,base_geom.nr_fine);
	BaseState<Real> p0(base_geom.max_radial_level+1,base_geom.nr_fine);

	Average(rho_mf, rho0, 0);
	Average(p0_mf, p0, 0);

	
	// Write radial output file
	Print() << "Writing radial diag file" << std::endl;
	WriteRadialFile(iFile, rho0, p0, u_mf, w0_mf);

	// Write 2D slice output file
	if (deltat.empty() || numfiles.empty()) {
	    Print() << "Writing 2D slice file for single plot file" << std::endl;
	    Write2dSliceFile (iFile, u_mf, w0_mf);
	} else {
	    Print() << "Writing 2D slice file for " << stoi(numfiles)+1 << " plot files" << std::endl;
	    Write2dSliceFile (iFile, u_mf, w0_mf, stoi(deltat), stoi(numfiles));
	}

	// Write diag file for initial model if specified
	if (!imFile.empty()) {
	    Print() << "Writing diag file for initial model" << std::endl;
	    WriteModelDiagFile(imFile);
	}
	
	
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}

///
/// Gets the variable ``varname`` from the ``job_info`` file and returns as a
/// string
///
std::string GetVarFromJobInfo (const std::string pltfile, const std::string& varname) {
    std::string filename = pltfile + "/job_info";
    std::regex re("(?:[ \\t]*)" + varname + "\\s*=\\s*(.*)\\s*\\n");
    //std::regex re(varname + " = ([^ ]*)\\n");
    
    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
	std::stringstream buf;
	buf << jobfile.rdbuf();
	std::string file_contents = buf.str();

	if (std::regex_search(file_contents, m, re)) {
	    return m[1];
	} else {
	    Print() << "Unable to find " << varname << " in job_info file!" << std::endl;
	}
    } else {
	Print() << "Could not open job_info file!" << std::endl;
    }

    return "";
}

// Get drdxfac from the job info file
int GetdrdxFac (const std::string pltfile) {
    auto drdxfac_str = GetVarFromJobInfo(pltfile, "maestro.drdxfac");
    // Print() << "drdxfac_str = " << drdxfac_str << std::endl;

    if (drdxfac_str == "") {
	return 5;
    }
    
    // retrieve first number
    std::istringstream iss {drdxfac_str};
    int drdxfac;
    std::string s;
    if (std::getline(iss, s, ' '))
	drdxfac = stoi(s);
    
    return drdxfac;
}

// Get octant from the job info file
bool GetOctant (const std::string pltfile) {
    auto octant_str = GetVarFromJobInfo(pltfile, "maestro.octant");
    
    if (octant_str == "true")
	return true;
    else
	// Print() << "octant_str = " << octant_str << std::endl;
	return false;
}
