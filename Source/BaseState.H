#ifndef BaseState_H_
#define BaseState_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

template <class T>
class BaseState
{
public:

    /*
      public member functions
    */

    BaseState() {};

    BaseState(const int num_levs, const int length=1, const int ncomp=1);
    
    /// create a BaseState object from an existing Fortran-friendly base state
    BaseState(const amrex::Gpu::ManagedVector<T>& src, const int num_levs, const int length, const int ncomp=1);

    BaseState(const amrex::Vector<T>& src, const int num_levs, const int length, const int ncomp=1);

    /// copy constructor. This makes a deep copy of the src.
    BaseState(const BaseState<T>& src);

    // /// destructor
    // ~BaseState () { };

    /// allocate memory for the BaseState
    void define(const int num_levs, const int length=1, const int ncomp=1);

    void resize(const int num_levs, const int length=1, const int ncomp=1);

    /// set to some scalar value 
    void setVal(const T& val);

    /// set the comp'th component to some scalar value
    void setVal(const int comp, const T& val);

    /// number of levels in base 
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
    int nLevels() const { return nlev; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int length() const { return len; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int nComp() const { return nvar; };

    /// deep copy 
    void copy(const BaseState<T>& src);
    void copy(const amrex::Gpu::ManagedVector<T>& src);
    void copy(const amrex::Vector<T>& src);

    void toVector(amrex::Vector<T>& vec);
    void toVector(amrex::Gpu::ManagedVector<T>& vec);

    /// swap the data with src
    void swap(BaseState<T>& src);

    /// returns a reference to the nth component at position i and level lev.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (const int lev, const int i=0, const int n=0) const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* operator() (const int lev, const int i, const int ncomp, const int start_comp);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* ptr (const int lev, const int i=0, const int n=0);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator[](const int i) { return dptr[i]; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const T& operator[](const int i) const { return dptr[i]; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T* dataPtr () const noexcept { return this->dptr; }

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
    BaseState<T> operator/ (const T val);
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

    amrex::Gpu::ManagedVector<T> base_data;

    T*  dptr;

    int nlev;
    int nvar;
    int len;    
};

template <class T>
BaseState<T>::BaseState(const int num_levs, const int length, const int ncomp) 
: nlev(num_levs), len(length), nvar(ncomp)
{
    base_data.resize(nlev*len*nvar);
    base_data.shrink_to_fit();
    dptr = base_data.dataPtr();
}

template <class T>
BaseState<T>::BaseState(const amrex::Gpu::ManagedVector<T>& src, const int num_levs, const int length, const int ncomp)
: nlev(num_levs), len(length), nvar(ncomp)
{
    base_data.resize(nlev*len*nvar);
    base_data.shrink_to_fit();
    dptr = base_data.dataPtr();
    // need to switch from Fortran array-ordering to C++ ordering
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto n = 0; n < nvar; ++n) {
                base_data[n + nvar*(r + len*l)] = src[l + nlev*(r + len*n)];
            }
        }
    }
}

template <class T>
BaseState<T>::BaseState(const amrex::Vector<T>& src, const int num_levs, const int length, const int ncomp)
: nlev(num_levs), len(length), nvar(ncomp)
{
    base_data.resize(nlev*len*nvar);
    base_data.shrink_to_fit();
    dptr = base_data.dataPtr();
    // need to switch from Fortran array-ordering to C++ ordering
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto n = 0; n < nvar; ++n) {
                base_data[n + nvar*(r + len*l)] = src[l + nlev*(r + len*n)];
            }
        }
    }
}

template <class T>
BaseState<T>::BaseState(const BaseState<T>& src) 
: nlev(src.nLevels()), len(src.length()), nvar(src.nComp())
{
    base_data.resize(nlev*len*nvar);
    base_data.shrink_to_fit();
    dptr = base_data.dataPtr();
    for (auto i = 0; i < nlev*len*nvar; i++) {
        base_data[i] = src[i];
    }
}

template <class T>
void
BaseState<T>::define(const int num_levs, const int length, const int ncomp)
{
    nlev = num_levs;
    len = length;
    nvar = ncomp;

    base_data.resize(nlev*len*nvar);
    base_data.shrink_to_fit();
    dptr = base_data.dataPtr();
}

template <class T>
void
BaseState<T>::resize(const int num_levs, const int length, const int ncomp)
{
    if (base_data.size() > 0) {
        base_data.resize(0);
    }
    this->define(num_levs, length, ncomp);
}

template <class T>
void
BaseState<T>::setVal(const T& val)
{
    T * AMREX_RESTRICT ptr = dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] = val;
    });

    amrex::Gpu::synchronize();
}

template <class T>
void
BaseState<T>::setVal(const int comp, const T& val)
{
    T * AMREX_RESTRICT ptr = dptr;
    const auto ncomp = nvar;
    AMREX_PARALLEL_FOR_1D(len*nlev, i, {
        ptr[comp + i*ncomp] = val;
    });

    amrex::Gpu::synchronize();
}

template <class T>
void
BaseState<T>::copy(const BaseState<T>& src)
{
    AMREX_ASSERT(nlev == src.nlev);
    AMREX_ASSERT(nvar == src.nvar);
    AMREX_ASSERT(len == src.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT src_ptr = src.dptr;
    AMREX_PARALLEL_FOR_1D(len*nlev*nvar, i, {
        ptr[i] = src_ptr[i];
    });

    amrex::Gpu::synchronize();
}

template <class T>
void
BaseState<T>::copy(const amrex::Gpu::ManagedVector<T>& src)
{
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar*(r + len*l)] = src[l + nlev*(r + len*comp)];
            } 
        }
    }
}

template <class T>
void
BaseState<T>::copy(const amrex::Vector<T>& src)
{
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                base_data[comp + nvar*(r + len*l)] = src[l + nlev*(r + len*comp)];
            } 
        }
    }
}

template <class T>
void
BaseState<T>::toVector(amrex::Vector<T>& vec)
{
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                vec[l + nlev*(r + len*comp)] = base_data[comp+nvar*(r + len*l)];
            }
        }
    }
}

template <class T>
void
BaseState<T>::toVector(amrex::Gpu::ManagedVector<T>& vec)
{
    for (auto l = 0; l < nlev; ++l) {
        for (auto r = 0; r < len; ++r) {
            for (auto comp = 0; comp < nvar; ++comp) {
                vec[l + nlev*(r + len*comp)] = base_data[comp+nvar*(r + len*l)];
            }
        }
    }
}

template <class T>
void
BaseState<T>::swap(BaseState<T>& src)
{
    AMREX_ASSERT(nlev == src.nlev);
    AMREX_ASSERT(nvar == src.nvar);
    AMREX_ASSERT(len == src.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT src_ptr = src.dptr;
    AMREX_PARALLEL_FOR_1D(len*nlev*nvar, i, {
        T temp = ptr[i];
        ptr[i] = src_ptr[i];
        src_ptr[i] = temp;
    });

    amrex::Gpu::synchronize();
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T&
BaseState<T>::operator() (const int lev, const int i, const int n) const {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT((i < this->len && i >= 0) || (lev == 0 && i < this->nlev));
    AMREX_ASSERT(n < this->nvar && n >= 0);

    return dptr[(lev*len + i)*nvar + n];
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T* 
BaseState<T>::operator() (const int lev, const int i, const int ncomp, const int start_comp) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(ncomp < this->nvar && ncomp > 0);
    AMREX_ASSERT(start_comp < this-> nvar && start_comp >= 0);
    AMREX_ASSERT(ncomp+start_comp < this->nvar);

    amrex::Gpu::ManagedVector<T> vec(ncomp);

    for (auto comp = 0; comp < ncomp; ++comp) {
        vec[comp] = dptr[(lev*len + i)*nvar + start_comp+comp];
    }

    return vec.dataPtr();
}

template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T*
BaseState<T>::ptr (const int lev, const int i, const int n) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    T* p = dptr[(lev*len + i)*nvar + n];

    return p;
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
    T * AMREX_RESTRICT ptr = dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] += val;
    });
    amrex::Gpu::synchronize();
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
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT rhs_ptr = rhs.dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] += rhs_ptr[i];
    });
    amrex::Gpu::synchronize();
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
    T * AMREX_RESTRICT ptr = dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] -= val;
    });
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
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT rhs_ptr = rhs.dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] -= rhs_ptr[i];
    });
    amrex::Gpu::synchronize();
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
    T * AMREX_RESTRICT ptr = dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] *= val;
    });
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
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT rhs_ptr = rhs.dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] *= rhs_ptr[i];
    });
    amrex::Gpu::synchronize();
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
BaseState<T>
BaseState<T>::operator/ (const T val) {
    // make a deep copy
    BaseState<T> s(*this);
    s /= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/= (const T val) {
    T * AMREX_RESTRICT ptr = dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] /= val;
    });
    amrex::Gpu::synchronize();
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
    AMREX_ASSERT(nlev == rhs.nlev);
    AMREX_ASSERT(nvar == rhs.nvar);
    AMREX_ASSERT(len == rhs.len);

    T * AMREX_RESTRICT ptr = dptr;
    T * AMREX_RESTRICT rhs_ptr = rhs.dptr;
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] /= rhs_ptr[i];
    });
    amrex::Gpu::synchronize();
    return *this;
}

template <class T>
bool
operator== (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.nlev == rhs.nlev);
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    for (auto i = 0; i < lhs.nvar*lhs.len*lhs.nlev; ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

template <class T>
bool
operator!= (const BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.nlev == rhs.nlev);
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    for (auto i = 0; i < lhs.nvar*lhs.len*lhs.nlev; ++i) {
        if (lhs[i] != rhs[i]) {
            return true;
        }
    }
    return false;
}
#endif