#include <BaseState.H>
#include <BaseStateGeometry.H>

extern BaseStateGeometry base_geom;
extern int finest_level;
extern amrex::GpuArray<amrex::Real,3> center;
//extern amrex::Vector<BoxArray> grid;
//extern amrex::Vector<DistributionMapping> dmap;
extern amrex::Vector<amrex::Geometry> pgeom;


// Read parameters from file
std::string GetVarFromJobInfo (const std::string pltfile, const std::string& varname);

// compute radially-averaged values
void Average (const amrex::Vector<amrex::MultiFab>& phi,
	      BaseState<amrex::Real>& phibar,
	      int comp);


// write radial plotfile to disk
void WriteRadialFile (const std::string& plotfilename,
		      const BaseState<amrex::Real>& rho0_in,
		      const BaseState<amrex::Real>& p0_in, 
		      const amrex::Vector<amrex::MultiFab>& u_in,
		      const amrex::Vector<amrex::MultiFab>& w0_in);

// postprocessing subroutines
void MakeRadialNFreq (const std::string& plotfilename,
		      const BaseState<amrex::Real>& p0_s,
		      const BaseState<amrex::Real>& rho0_s, 
		      BaseState<amrex::Real>& freq0);

void MakeConvectionVel (const amrex::Vector<amrex::MultiFab>& velr,
			BaseState<amrex::Real>& vel_conv);

void MakeVelrc (const amrex::Vector<amrex::MultiFab>& vel,
		const amrex::Vector<amrex::MultiFab>& w0rcart,
		amrex::Vector<amrex::MultiFab>& rad_vel,
		amrex::Vector<amrex::MultiFab>& circ_vel);
