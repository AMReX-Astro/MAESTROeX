#ifndef _MODEL_PARSER_H_
#define _MODEL_PARSER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <state_indices.H>
#include <network_properties.H>
#include <AMReX.H>

/// Define Real vector types for CUDA-compatability. If `AMREX_USE_CUDA`, then
/// this will be stored in CUDA managed memory.
#ifdef AMREX_USE_CUDA
typedef amrex::Gpu::ManagedVector<amrex::Real> RealVector;
typedef amrex::Gpu::ManagedVector<int> IntVector;
#else
typedef amrex::Vector< amrex::Real > RealVector;
typedef amrex::Vector< int > IntVector;
#endif

class ModelParser
{
public:

    ModelParser() {model_initialized = false};

    void ReadFile(const string& model_file);

    amrex::Real Interpolate(const amrex::Real r, const int ivar, 
                            bool interpolate_top=false);

    // arrays for storing the model data
    Vector<RealVector> model_state;
    RealVector model_r;

    static int npts_model;

    static bool model_initialized;

private:
    // integer keys for indexing the model_state array
    constexpr int idens_model = 0;
    constexpr int itemp_model = 1;
    constexpr int ipres_model = 2;
    constexpr int ispec_model = 3;
    constexpr int nvars_model = 3 + NumSpec;
}

#endif