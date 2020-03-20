
#include <BaseState.H>

using namespace amrex;

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
: len(src.len()), nvar(src.nComp())
{
    data.resize(src.size());

    for (auto l = 0; l < src.size(); ++l) {
        data[l].resize(len*nvar);
        for (auto i = 0; i < len*nvar; ++i) {
            data[l][i] = src[l][i];
        }
    }
}

template <class T>
void
BaseState<T>::define(const int num_levs, const int length, const int ncomp)
: len(length), nvar(ncomp)
{
    data.resize(num_levs);

    for (auto l = 0; l < num_levs; ++l) {
        data[l].resize(len*nvar);
    }
}

template <class T>
void
BaseState<T>::setVal(const T& val)
{
    for (auto l = 0; l < data.size(); ++l) {
        AMREX_PARALLEL_FOR_1D(nvar*ncomp, i, {
            data[l][i] = val;
        })
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

    return vec;
}


template <class T>
BaseState<T>&
BaseState<T>::operator+ (const T& val, const BaseState<T>& p) {
    // make a deep copy
    auto s = BaseState<T>(p);
    s += val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator+= (const T& val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] += val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator+ (BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    auto s = BaseState<T>(lhs);
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
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] += rhs[l][i];
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator- (const T& val, const BaseState<T>& p) {
    // make a deep copy
    auto s = BaseState<T>(p);
    s -= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator-= (const T& val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] -= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator- (BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    auto s = BaseState<T>(lhs);
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
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] -= rhs[l][i];
        }
    }
    return *this;
}


template <class T>
BaseState<T>&
BaseState<T>::operator* (const T& val, const BaseState<T>& p) {
    // make a deep copy
    auto s = BaseState<T>(p);
    s *= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator*= (const T& val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] *= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator* (BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    auto s = BaseState<T>(lhs);
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
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] *= rhs[l][i];
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/ (const T& val, const BaseState<T>& p) {
    // make a deep copy
    auto s = BaseState<T>(p);
    s /= val;
    return s;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/= (const T& val) {
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] /= val;
        }
    }
    return *this;
}

template <class T>
BaseState<T>&
BaseState<T>::operator/ (BaseState<T>& lhs, const BaseState<T>& rhs) {
    // make a deep copy
    auto s = BaseState<T>(lhs);
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
        for (auto i = 0; i < nvar*ncomp; ++i) {
            data[l][i] /= rhs[l][i];
        }
    }
    return *this;
}

template <class T>
bool
BaseState<T>::operator== (BaseState<T>& lhs, const BaseState<T>& rhs) {
    AMREX_ASSERT(lhs.data.size() == rhs.data.size());
    AMREX_ASSERT(lhs.nvar == rhs.nvar);
    AMREX_ASSERT(lhs.len == rhs.len);

    bool equiv = true;
    for (auto l = 0; l < data.size(); ++l) {
        for (auto i = 0; i < nvar*ncomp; ++i) {
            if (lhs[l][i] != rhs[l][i]) {
                equiv = false;
                break;
            }
            if (!equiv) break;
        }
    }
    return equiv;
}