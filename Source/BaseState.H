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

    /// constructor
    // constexpr BaseState() : dptr(nullptr) {};

    BaseState(const int num_levs, const int length, const int ncomp=1);

    /// copy constructor. This makes a deep copy of the src.
    BaseState(const BaseState<T>& src);

    // /// destructor
    // ~BaseState () { };

    // /// allocate memory for the BaseState
    // inline
    // void define(const int num_levs, const int length, const int ncomp=1);

    /// set to some scalar value 
    void setVal(const T val);

    /// set the comp'th component to some scalar value
    void setVal(const int comp, const T val);

    /// number of levels in base 
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE 
    int nLevels() const { return nlev; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int length() const { return len; };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int nComp() const { return nvar; };

    /// returns a reference to the nth component at position i and level lev.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (int lev, int i, int n=0) const;

    AMREX_GPU_HOST_DEVICE inline
    T* operator() (const int lev, const int i, const int ncomp, const int start_comp);

    AMREX_GPU_HOST_DEVICE inline
    T* ptr (const int lev, const int i, const int n=0);

    AMREX_GPU_HOST_DEVICE inline
    T& operator[](const int i) { return dptr[i]; };

    AMREX_GPU_HOST_DEVICE inline
    const T& operator[](const int i) const { return dptr[i]; };

//     /// scalar addition to the whole base state 
//     template <class U>
//     friend BaseState<U> operator+ (const U val, const BaseState<U>& p);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator+= (const T val);

//     /// element-wise addition 
//     template <class U>
//     friend BaseState<U> operator+ (const BaseState<U>& lhs, const BaseState<U>& rhs);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator+= (const BaseState<T>& rhs);

//     /// scalar subtraction from the whole base state 
//     template <class U>
//     friend BaseState<U> operator- (const U val, const BaseState<U>& p);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator-= (const T val);

//     /// element-wise subtraction
//     template <class U>
//     friend BaseState<U> operator- (const BaseState<U>& lhs, const BaseState<U>& rhs);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator-= (const BaseState<T>& rhs);

//     /// scalar multiplication of the whole base state 
//     template <class U>
//     friend BaseState<U> operator* (const U val, const BaseState<U>& p);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator*= (const T val);

//     /// element-wise multiplication 
//     template <class U>
//     friend BaseState<U> operator* (const BaseState<U>& lhs, const BaseState<U>& rhs);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator*= (const BaseState<T>& rhs);

//     /// scalar division of the whole base state 
//     template <class U>
//     friend BaseState<U> operator/ (const U val, const BaseState<U>& p);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator/= (const T val);

//     /// element-wise division 
//     template <class U>
//     friend BaseState<U> operator/ (const BaseState<U>& lhs, const BaseState<U>& rhs);
//     AMREX_GPU_HOST_DEVICE inline
//     BaseState<T>& operator/= (const BaseState<T>& rhs);

    /// comparison operator
    template <class U>
    friend bool operator==(const BaseState<U>& lhs, const BaseState<U>& rhs);

//     /// comparison operator
//     template <class U>
//     friend bool operator!=(const BaseState<U>& lhs, const BaseState<U>& rhs);

    

// // protected:

    // amrex::GpuArray<T,nlev*len*nvar> base_data;

    amrex::Gpu::ManagedVector<T> base_data;

    T*  dptr;
//     amrex::Real * AMREX_RESTRICT real_ptr;

    int nlev;
    int nvar;
    int len;


    // amrex::MFIter MakeMFIter () {
    //     amrex::Vector< amrex::MultiFab> m_dummy_mf(1);
    //     m_dummy_mf[0].define(grids[0], dmap[0], 1, 0);
    //     return amrex::MFIter(m_dummy_mf[0], false);
    // }

    
};



// template <class T, int nlev, int len, int nvar>
// BaseState<T,nlev,len,nvar>::BaseState()
// {

// }

template <class T>
BaseState<T>::BaseState(const int num_levs, const int length, const int ncomp) 
: nlev(num_levs), len(length), nvar(ncomp)
{
    base_data.resize(nlev*len*nvar);
    dptr = base_data.dataPtr();
}

// template <class T>
// BaseState<T>::BaseState(const BaseState<T>& src) 
// : nlev(src.nLevels()), len(src.length()), nvar(src.nComp())
// {
//     base_data.resize(nlev*len*nvar);
//     dptr = base_data.dataPtr();
//     for (auto i = 0; i < nlev*len*nvar; ++i) {
//         dptr[i] = src[i];
//     }
// }

// template <class T>
// inline
// void
// BaseState<T>::define(const int num_levs, const int length, const int ncomp)
// {
//     nlev = num_levs;
//     len = length;
//     nvar = ncomp;

//     base_data.resize(nlev*len*nvar);
//     dptr = base_data.dataPtr();
// }

template <class T>
void
BaseState<T>::setVal(const T val)
{
    T * AMREX_RESTRICT ptr = base_data.dataPtr();
    AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
        ptr[i] = val;
    });

    amrex::Gpu::synchronize();

    // AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
    //     dptr[i] = val;
    // });

    // T* AMREX_RESTRICT p = this->base_data.data();
    // const int l = this->nvar*this->len*this->nlev;
    // amrex::ParallelFor(l,
    // [=] AMREX_GPU_DEVICE (long idx)
    // {
    //     p[idx] = val;
    // });
    amrex::Print() << "print data here " << base_data[0] << std::endl;

    // amrex::Print() << "print pointer " << p[0] << std::endl;
    // // amrex::Print() << "print dptr " << dptr[0] << std::endl;

    // AMREX_PARALLEL_FOR_1D(nvar*len*nlev, i, {
    //     real_ptr[i] = val;
    // });

    // amrex::Print() << "print real data here " << real_base[0] << std::endl;

    // amrex::Print() << "print real pointer " << real_ptr[0] << std::endl;



    // amrex::Gpu::ManagedVector<T> testvec(nvar*len*nlev, 0.0);

    // T* tp = testvec.dataPtr();
    // AMREX_PARALLEL_FOR_1D(l, i, {
    //     tp[i] = val;
    // });

    // amrex::Print() << "print tp here " << testvec[0] << std::endl;


    // for (auto i = 0; i < nvar*len*nlev; ++i) {
    //     dptr[i] = val;
    // }

    // amrex::Print() << "print data again here " << base_data[0] << std::endl;

    // amrex::Print() << "val = " << val << std::endl;
}

template <class T>
void
BaseState<T>::setVal(const int comp, const T val)
{
    T * AMREX_RESTRICT ptr = base_data.dataPtr();
    AMREX_PARALLEL_FOR_1D(len*nlev, i, {
        ptr[comp + i*nvar] = val;
    });
}

template <class T>
AMREX_GPU_HOST_DEVICE inline
T&
BaseState<T>::operator() (int lev, int i, int n) const {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    return dptr[(lev*len + i)*nvar + n];
}

template <class T>
AMREX_GPU_HOST_DEVICE inline
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
AMREX_GPU_HOST_DEVICE inline
T*
BaseState<T>::ptr (const int lev, const int i, const int n) {
    AMREX_ASSERT(lev >= 0);
    AMREX_ASSERT(i < this->len && i >= 0);
    AMREX_ASSERT(n < this->nvar && n >= 0);

    T* p = dptr[(lev*len + i)*nvar + n];

    return p;// + n + nvar*i;
}


// template <class T>
// inline
// BaseState<T>
// operator+ (const T val, const BaseState<T>& p) {
//     // make a deep copy
//     BaseState<T> s(p);
//     s += val;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator+= (const T val) {
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] += val;
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator+ (const BaseState<T>& lhs, const BaseState<T>& rhs) {
//     // make a deep copy
//     BaseState<T> s(lhs);
//     s += rhs;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator+= (const BaseState<T>& rhs) {
//     AMREX_ASSERT(nlev == rhs.nlev);
//     AMREX_ASSERT(nvar == rhs.nvar);
//     AMREX_ASSERT(len == rhs.len);

//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] += rhs[i];
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator- (const T val, const BaseState<T>& p) {
//     // make a deep copy
//     BaseState<T> s(p);
//     s -= val;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator-= (const T val) {
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] -= val;
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator- (const BaseState<T>& lhs, const BaseState<T>& rhs) {
//     // make a deep copy
//     BaseState<T> s(lhs);
//     s -= rhs;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator-= (const BaseState<T>& rhs) {
//     AMREX_ASSERT(nlev == rhs.nlev);
//     AMREX_ASSERT(nvar == rhs.nvar);
//     AMREX_ASSERT(len == rhs.len);

//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] += rhs[i];
//     }
//     return *this;
// }


// template <class T>
// inline
// BaseState<T>
// operator* (const T val, const BaseState<T>& p) {
//     // make a deep copy
//     BaseState<T> s(p);
//     s *= val;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator*= (const T val) {
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] *= val;
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator* (const BaseState<T>& lhs, const BaseState<T>& rhs) {
//     // make a deep copy
//     BaseState<T> s(lhs);
//     s *= rhs;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator*= (const BaseState<T>& rhs) {
//     AMREX_ASSERT(nlev == rhs.nlev);
//     AMREX_ASSERT(nvar == rhs.nvar);
//     AMREX_ASSERT(len == rhs.len);
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] += rhs[i];
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator/ (const T val, const BaseState<T>& p) {
//     // make a deep copy
//     BaseState<T> s(p);
//     s /= val;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator/= (const T val) {
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] /= val;
//     }
//     return *this;
// }

// template <class T>
// inline
// BaseState<T>
// operator/ (const BaseState<T>& lhs, const BaseState<T>& rhs) {
//     // make a deep copy
//     BaseState<T> s(lhs);
//     s /= rhs;
//     return s;
// }

// template <class T>
// AMREX_GPU_HOST_DEVICE inline
// BaseState<T>&
// BaseState<T>::operator/= (const BaseState<T>& rhs) {
//     AMREX_ASSERT(nlev == rhs.nlev);
//     AMREX_ASSERT(nvar == rhs.nvar);
//     AMREX_ASSERT(len == rhs.len);
//     for (auto i = 0; i < nvar*len*nlev; ++i) {
//         dptr[i] /= rhs[i];
//     }
//     return *this;
// }

template <class T>
inline
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

// template <class T>
// inline
// bool
// operator!= (const BaseState<T>& lhs, const BaseState<T>& rhs) {
//     AMREX_ASSERT(lhs.nlev == rhs.nlev);
//     AMREX_ASSERT(lhs.nvar == rhs.nvar);
//     AMREX_ASSERT(lhs.len == rhs.len);

//     for (auto i = 0; i < lhs.nvar*lhs.len*lhs.nlev; ++i) {
//         if (lhs[i] != rhs[i]) {
//             return true;
//         }
//     }
//     return false;
// }




#endif