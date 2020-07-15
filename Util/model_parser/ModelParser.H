#ifndef _MODEL_PARSER_H_
#define _MODEL_PARSER_H_

#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <BaseState.H>
#include <MaestroUtil.H>
#include <network_properties.H>
#include <state_indices.H>

/// Define Real vector types for CUDA-compatability. If `AMREX_USE_CUDA`, then
/// this will be stored in CUDA managed memory.
#ifdef AMREX_USE_CUDA
typedef amrex::Gpu::ManagedVector<amrex::Real> RealVector;
typedef amrex::Gpu::ManagedVector<int> IntVector;
#else
typedef amrex::Vector<amrex::Real> RealVector;
typedef amrex::Vector<int> IntVector;
#endif

class ModelParser {
   public:
    ModelParser() { model_initialized = false; };

    void ReadFile(const std::string& model_file);

    amrex::Real Interpolate(const amrex::Real r, const int ivar,
                            bool interpolate_top = false);

    // // arrays for storing the model data
    amrex::Vector<RealVector> model_state;
    RealVector model_r;

    int npts_model;

    bool model_initialized;

    // integer keys for indexing the model_state array
    static constexpr int idens_model = 0;
    static constexpr int itemp_model = 1;
    static constexpr int ipres_model = 2;
    static constexpr int ispec_model = 3;
    static constexpr int nvars_model = 3 + NumSpec;
};

#endif