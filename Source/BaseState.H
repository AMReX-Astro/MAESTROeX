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
    BaseState() {};

    BaseState(const int num_levs);

    BaseState(const int num_levs, const int length, const int ncomp=1);

    /// copy constructor. This makes a deep copy of the src.
    BaseState(const BaseState<T>& src);

    /// destructor
    virtual ~BaseState () {};

    /// allocate memory for the BaseState
    void define(const int num_levs, const int length, const int ncomp=1);

    /// set to some scalar value 
    void setVal(const T& val);

    /// set the comp'th component to some scalar value
    void setVal(const int comp, const T& val);

    /// number of levels in base 
    int nLevels() const { return data.size(); };
    int length() const { return len; };
    int nComp() const { return nvar; };

    /// returns a references to the nth component at position i and level lev.
    T& operator() (const int lev, const int i, const int n=0);

    T* AMREX_RESTRICT operator() (const int lev, const int i, const int ncomp, const int start_comp);

    T* ptr (const int lev, const int i, const int n=0);

    amrex::Gpu::ManagedVector<T>& operator[](const int lev) { return data[lev]; };

    const amrex::Gpu::ManagedVector<T>& operator[](const int lev) const { return data[lev]; };

    /// scalar addition to the whole base state 
    template <class U>
    friend BaseState<U> operator+ (const U val, const BaseState<U>& p);
    BaseState<T>& operator+= (const T val);

    /// element-wise addition 
    template <class U>
    friend BaseState<U> operator+ (const BaseState<U>& lhs, const BaseState<U>& rhs);
    BaseState<T>& operator+= (const BaseState<T>& rhs);

    /// scalar subtraction from the whole base state 
    template <class U>
    friend BaseState<U> operator- (const U val, const BaseState<U>& p);
    BaseState<T>& operator-= (const T val);

    /// element-wise subtraction
    template <class U>
    friend BaseState<U> operator- (const BaseState<U>& lhs, const BaseState<U>& rhs);
    BaseState<T>& operator-= (const BaseState<T>& rhs);

    /// scalar multiplication of the whole base state 
    template <class U>
    friend BaseState<U> operator* (const U val, const BaseState<U>& p);
    BaseState<T>& operator*= (const T val);

    /// element-wise multiplication 
    template <class U>
    friend BaseState<U> operator* (const BaseState<U>& lhs, const BaseState<U>& rhs);
    BaseState<T>& operator*= (const BaseState<T>& rhs);

    /// scalar division of the whole base state 
    template <class U>
    friend BaseState<U> operator/ (const U val, const BaseState<U>& p);
    BaseState<T>& operator/= (const T val);

    /// element-wise division 
    template <class U>
    friend BaseState<U> operator/ (const BaseState<U>& lhs, const BaseState<U>& rhs);
    BaseState<T>& operator/= (const BaseState<T>& rhs);

    /// comparison operator
    template <class U>
    friend bool operator==(const BaseState<U>& lhs, const BaseState<U>& rhs);

    /// comparison operator
    template <class U>
    friend bool operator!=(const BaseState<U>& lhs, const BaseState<U>& rhs);

protected:

    amrex::Vector< amrex::Gpu::ManagedVector<T> > data;

    int nvar;
    int len;
};

template <class T>
BaseState<T>::BaseState(const int num_levs) 
{
    data.resize(num_levs);
}

template <class T>
BaseState<T>::BaseState(const int num_levs, const int length, const int ncomp) 
: len(length), nvar(ncomp)
{
    data.resize(num_levs);

    for (auto l = 0; l < num_levs; ++l) {
        data[l].resize(len*nvar);
    }
}

template <class T>
BaseState<T>::BaseState(const BaseState<T>& src) 
: len(src.length()), nvar(src.nComp())
{
    data.resize(src.data.size());

    for (auto l = 0; l < src.data.size(); ++l) {
        data[l].resize(len*nvar);
        for (auto i = 0; i < len*nvar; ++i) {
            data[l][i] = src[l][i];
        }
    }
}

template <class T>
void
BaseState<T>::define(const int num_levs, const int length, const int ncomp)
{
    data.resize(num_levs);
    len = length;
    nvar = ncomp;

    for (auto l = 0; l < num_levs; ++l) {
        data[l].resize(len*nvar);
    }
}

template <class T>
void
BaseState<T>::setVal(const T& val)
{
    for (auto l = 0; l < data.size(); ++l) {
        AMREX_PARALLEL_FOR_1D(nvar*len, i, {
            data[l][i] = val;
        });
    }
}

template <class T>
void
BaseState<T>::setVal(const int comp, const T& val)
{
    for (auto l = 0; l < data.size(); ++l) {
        AMREX_PARALLEL_FOR_1D(len, i, {
            data[l][comp + i*nvar] = val;
        });
    }
}

template <class T>
T&
BaseState<T>::operator() (const int lev, const int i, const int n) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    return data[lev][nvar*i + n];
}

template <class T>
T* AMREX_RESTRICT
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

    return vec.dataPtr();
}

template <class T>
T*
BaseState<T>::ptr (const int lev, const int i, const int n) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    T* AMREX_RESTRICT p = data[lev].dataPtr();

    return p + n + nvar*i;
}


template <class T>
BaseState<T>
operator+ (const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s += val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator+= (const T val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] += val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator+ (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s += rhs;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator+= (const BaseState<T>& rhs) {
    AMREX_ASSERT(data.size() == rhs.data.size());
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] += rhs[l][i];
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator- (const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s -= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator-= (const T val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] -= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator- (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s -= rhs;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator-= (const BaseState<T>& rhs) {
    AMREX_ASSERT(data.size() == rhs.data.size());
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] -= rhs[l][i];
        }
    }
    return *this;
}


template <class T>
BaseState<T>
operator* (const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s *= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator*= (const T val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] *= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator* (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s *= rhs;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator*= (const BaseState<T>& rhs) {
    AMREX_ASSERT(data.size() == rhs.data.size());
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] *= rhs[l][i];
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator/ (const T val, const BaseState<T>& p) {
    // make a deep copy
    BaseState<T> s(p);
    s /= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/= (const T val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] /= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>
operator/ (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    BaseState<T> s(lhs);
    s /= rhs;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/= (const BaseState<T>& rhs) {
    AMREX_ASSERT(data.size() == rhs.data.size());
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*len; ++i) {
            data[l][i] /= rhs[l][i];
        }
    }
    return *this;
}

template <class T>
bool
operator== (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.data.size() == rhs.data.size());
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    bool equiv = true;
    for (auto l = 0; l < lhs.data.size(); ++l) {
        for (auto i = 0; i < lhs.nvar*lhs.len; ++i) {
            if (lhs[l][i] != rhs[l][i]) {
                equiv = false;
                break;
            }
            if (!equiv) break;
        }
    }
    return equiv;
}

template <class T>
bool
operator!= (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.data.size() == rhs.data.size());
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    bool equiv = false;
    for (auto l = 0; l < lhs.data.size(); ++l) {
        for (auto i = 0; i < lhs.nvar*lhs.len; ++i) {
            if (lhs[l][i] != rhs[l][i]) {
                return true;
            }
        }
    }
    return equiv;
}
#endif