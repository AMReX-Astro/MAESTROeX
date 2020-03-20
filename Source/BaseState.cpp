
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
T&
BaseState<T>::operator() (const int lev, const int i, const int n=0)
{
    return data[lev][i*nvar + n];
}



// template <class T>
// BaseState<T>::~BaseState() {
// }

