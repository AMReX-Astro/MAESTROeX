#include <BaseState.H>
#include <BaseStateGeometry.H>

extern BaseStateGeometry base_geom;
extern int finest_level;

// compute radially-averaged values
void Average (const amrex::Vector<amrex::Geometry>& geom,
	      const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& center,
	      const amrex::Vector<amrex::MultiFab>& phi,
	      BaseState<amrex::Real>& phibar,
	      int comp);


// write radial plotfile to disk
void WriteRadialFile (const std::string& plotfilename,
		      const BaseState<amrex::Real>& rho0_in,
		      const BaseState<amrex::Real>& rhoh0_in,
		      const BaseState<amrex::Real>& p0_in,
		      const amrex::Vector<amrex::MultiFab>& u_in,
		      const amrex::Vector<amrex::MultiFab>& w0_in);

// postprocessing subroutines
void MakeRadialNFreq (const BaseState<amrex::Real>& p0_s,
		      const BaseState<amrex::Real>& rho0_s, 
		      BaseState<amrex::Real>& freq0);


