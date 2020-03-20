#ifndef BaseState_H_
#define BaseState_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>

/// Define Real vector types for CUDA-compatability. If `AMREX_USE_CUDA`, then
/// this will be stored in CUDA managed memory.
#ifdef AMREX_USE_CUDA
typedef amrex::Gpu::ManagedVector<amrex::Real> RealVector;
typedef amrex::Gpu::ManagedVector<int> IntVector;
#else
typedef amrex::Vector< amrex::Real > RealVector;
typedef amrex::Vector< int > IntVector;
#endif

template <class T>
class BaseState
{
public:

    /*
      public member functions
    */

    /// constructor
    BaseState();

    BaseState(const int num_levs);

    BaseState(const int num_levs, const int length, const int ncomp=1);

    /// copy constructor. This makes a deep copy of the src.
    BaseState(const BaseState<T>& src);

    /// destructor
    virtual ~BaseState ();

    /// allocate memory for the BaseState
    void define(const int num_levs, const int length, const int ncomp=1);

    /// number of levels in base 
    int nLevels() { return data.size(); };
    int len() { return len; };
    int nComp() { return nvar; };

    /// returns a references to the nth component at position i and level lev.
    T& operator() (const int lev, const int i, const int n=0);

    amrex::Gpu::ManagedVector<T> operator() (const int lev, const int i, const int ncomp, const int start_comp);

    /// same but for const BaseStates
    const T& operator() (const int lev, const int i, const int n=0);

private:

    amrex::Vector< amrex::Gpu::ManagedVector<T> > data;

    static int nvar;
    static int len;
};

template <class T>
T&
BaseState<T>::operator() (const int lev, const int i, const int n) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    return data[lev][nvar*i + n];
}

template <class T>
amrex::Gpu::ManagedVector<T>
BaseState<T>::operator() (const int lev, const int i, const int ncomp, const int start_comp) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(ncomp < this->nvar && ncomp > 0);
    AMREX_ASSERT(start_comp < this-> nvar && start_comp >= 0);
    AMREX_ASSERT(ncomp+start_comp < this->nvar);

    amrex::Gpu::ManagedVector<T> vec(ncomp);

    for (auto comp = 0; comp < ncomp; ++comp) {
        vec[comp] = data[lev][nvar*i + start_comp+comp];
    }
}

#endif