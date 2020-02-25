#include <eos_data.H>

// Declare extern variables

bool EOSData::initialized;
AMREX_GPU_MANAGED amrex::Real EOSData::mintemp;
AMREX_GPU_MANAGED amrex::Real EOSData::maxtemp;
AMREX_GPU_MANAGED amrex::Real EOSData::mindens;
AMREX_GPU_MANAGED amrex::Real EOSData::maxdens;
AMREX_GPU_MANAGED amrex::Real EOSData::minx;
AMREX_GPU_MANAGED amrex::Real EOSData::maxx;
AMREX_GPU_MANAGED amrex::Real EOSData::minye;
AMREX_GPU_MANAGED amrex::Real EOSData::maxye;
AMREX_GPU_MANAGED amrex::Real EOSData::mine;
AMREX_GPU_MANAGED amrex::Real EOSData::maxe;
AMREX_GPU_MANAGED amrex::Real EOSData::minp;
AMREX_GPU_MANAGED amrex::Real EOSData::maxp;
AMREX_GPU_MANAGED amrex::Real EOSData::mins;
AMREX_GPU_MANAGED amrex::Real EOSData::maxs;
AMREX_GPU_MANAGED amrex::Real EOSData::minh;
AMREX_GPU_MANAGED amrex::Real EOSData::maxh;
