#include <BaseState.H>
#include <BaseStateGeometry.H>

extern BaseStateGeometry base_geom;
extern int finest_level;
extern amrex::GpuArray<amrex::Real,3> center;
extern amrex::Vector<amrex::BoxArray> grid;
extern amrex::Vector<amrex::DistributionMapping> dmap;
extern amrex::Vector<amrex::Geometry> pgeom;


// Read parameters from file
std::string GetVarFromJobInfo (const std::string pltfile, const std::string& varname);


// write radial plotfile to disk
void WriteRadialFile (const std::string& plotfilename,
		      const BaseState<amrex::Real>& rho0_in,
		      const BaseState<amrex::Real>& p0_in, 
		      const amrex::Vector<amrex::MultiFab>& u_in,
		      const amrex::Vector<amrex::MultiFab>& w0_in);

// write additional diagnostics from model file
void WriteModelDiagFile (const std::string& plotfilename);

// read radial file from disk
void ReadBaseCCFile (const std::string& plotfilename,
		     BaseState<amrex::Real>& r_in,
		     BaseState<amrex::Real>& rho0_in,
		     BaseState<amrex::Real>& p0_in);

// read input model file from disk
void ReadModelFile (const std::string& plotfilename,
		     BaseState<amrex::Real>& r_in,
		     BaseState<amrex::Real>& rho0_in,
		     BaseState<amrex::Real>& p0_in);

// postprocessing subroutines
void MakeRadialNFreq (const std::string& plotfilename,
		      const BaseState<amrex::Real>& p0_s,
		      const BaseState<amrex::Real>& rho0_s, 
		      BaseState<amrex::Real>& freq0,
		      BaseState<amrex::Real>& entropy_s, BaseState<amrex::Real>& grav_s);

void MakeRadialGrav (const BaseState<amrex::Real>& rho0_in,
		     BaseState<amrex::Real>& grav,
		     const BaseStateArray<amrex::Real>& r_cc_loc);

void MakeConvectionVel (const amrex::Vector<amrex::MultiFab>& velr,
			BaseState<amrex::Real>& vel_conv);

void MakeRadialRotationRatio (const BaseState<amrex::Real>& s0_in,
			      const amrex::Vector<amrex::MultiFab>& omega,
			      BaseState<amrex::Real>& ratio_omega);

void MakeVelrc (const amrex::Vector<amrex::MultiFab>& vel,
		const amrex::Vector<amrex::MultiFab>& w0rcart,
		amrex::Vector<amrex::MultiFab>& rad_vel,
		amrex::Vector<amrex::MultiFab>& circ_vel);

void MakeRotationRate (const std::string& plotfilename,
		       const amrex::Vector<amrex::MultiFab>& vel,
		       const amrex::Vector<amrex::MultiFab>& w0cart,
		       amrex::Vector<amrex::MultiFab>& omega);

void MakeLatShear (const amrex::Vector<amrex::MultiFab>& omega_in,
		   BaseState<amrex::Real>& shear,
		   const BaseStateArray<amrex::Real>& r_cc_loc);

// ----------------------------
// in Slice.cpp
// ----------------------------

// write 2d slice plotfile to disk
void Write2dSliceFile (const std::string& plotfilename,
		       const amrex::Vector<amrex::MultiFab>& u_in,
		       const amrex::Vector<amrex::MultiFab>& w0_in, 
		       const int deltat = 0, const int nfiles = 0);

void MakeMeridionalCirculation (const amrex::Vector<amrex::MultiFab>& vel,
				const amrex::Vector<amrex::MultiFab>& w0rcart,
				amrex::Vector<amrex::MultiFab>& radii_vel,
				const int max_grid_size,
				amrex::Geometry& geomFine);


// ----------------------------
// in Average.cpp
// ----------------------------

// compute radially-averaged values
void Average (const amrex::Vector<amrex::MultiFab>& phi,
	      BaseState<amrex::Real>& phibar,
	      int comp);

// compute averaged value onto r-theta (x-z) plane
void Average2d (const amrex::Vector<amrex::MultiFab>& phi,
		amrex::Vector<amrex::MultiFab>& phibar,
		int barcomp, 
		const amrex::Vector<amrex::Box>& bardomain,
		int comp);
