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

    /// set to some scalar value 
    void setVal(const T& val);

    /// number of levels in base 
    int nLevels() { return data.size(); };
    int len() { return len; };
    int nComp() { return nvar; };

    /// returns a references to the nth component at position i and level lev.
    T& operator() (const int lev, const int i, const int n=0);

    amrex::Gpu::ManagedVector<T> operator() (const int lev, const int i, const int ncomp, const int start_comp);

    /// same but for const BaseStates
    const T& operator() (const int lev, const int i, const int n=0);

    /// scalar addition to the whole base state 
    friend BaseState<T>& operator+ (const T& val, const BaseState<T>& p);
    BaseState<T>& operator+= (const T& val);

    /// element-wise addition 
    friend BaseState<T>& operator+ (BaseState<T>& lhs, const BaseState<T>& rhs);
    BaseState<T>& operator+= (const BaseState<T>& rhs);

    /// scalar subtraction from the whole base state 
    friend BaseState<T>& operator- (const T& val, const BaseState<T>& p);
    BaseState<T>& operator-= (const T& val);

    /// element-wise subtraction
    friend BaseState<T>& operator- (BaseState<T>& lhs, const BaseState<T>& rhs);
    BaseState<T>& operator-= (const BaseState<T>& rhs);

    /// scalar multiplication of the whole base state 
    friend BaseState<T>& operator* (const T& val, const BaseState<T>& p);
    BaseState<T>& operator*= (const T& val);

    /// element-wise multiplication 
    friend BaseState<T>& operator* (BaseState<T>& lhs, const BaseState<T>& rhs);
    BaseState<T>& operator*= (const BaseState<T>& rhs);

    /// scalar division of the whole base state 
    friend BaseState<T>& operator/ (const T& val, const BaseState<T>& p);
    BaseState<T>& operator/= (const T& val);

    /// element-wise division 
    friend BaseState<T>& operator/ (BaseState<T>& lhs, const BaseState<T>& rhs);
    BaseState<T>& operator/= (const BaseState<T>& rhs);

    /// comparison operator
    friend bool operator==(BaseState<T>& lhs, const BaseState<T>& rhs);



protected:

    amrex::Vector< amrex::Gpu::ManagedVector<T> > data;

    static int nvar;
    static int len;
};


#endif