#ifndef PhysBCFunctMaestro_H_
#define PhysBCFunctMaestro_H_

#include <AMReX_PhysBCFunct.H>

extern "C"
{
    typedef void (*BndryFuncDefaultMaestro)(Real* data, AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                     const int* dom_lo, const int* dom_hi,
                                     const Real* dx, const Real* grd_lo,
                                     const Real* time, const int* bc,
                                     const int icomp);
    typedef void (*BndryFunc3DDefaultMaestro)(Real* data, const int* lo, const int* hi,
                                       const int* dom_lo, const int* dom_hi,
                                       const Real* dx, const Real* grd_lo,
                                       const Real* time, const int* bc,
                                       const int icomp);
}

//! This version calls function working on array
class BndryFuncArrayMaestro : amrex::BndryFuncArray
{
public:
    BndryFuncArray () noexcept {}
    BndryFuncArray (BndryFuncDefaultMaestro inFunc) noexcept : m_func(inFunc) {}
    BndryFuncArray (BndryFunc3DDefaultMaestro inFunc) noexcept : m_func3D(inFunc) {}

    void operator() (Box const& bx, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int orig_comp);

    bool RunOnGPU () const noexcept { return m_run_on_gpu; }
    void setRunOnGPU (bool b) noexcept { m_run_on_gpu = b; }

protected:
    BndryFuncDefaultMaestro   m_func   = nullptr;
    BndryFunc3DDefaultMaestro m_func3D = nullptr;
    bool m_run_on_gpu = false;
};

void
BndryFuncArrayMaestro::operator () (Box const& /*bx*/, FArrayBox& dest,
                                    const int dcomp, const int numcomp,
                                    Geometry const& geom, const Real time,
                                    const Vector<BCRec>& bcr, const int bcomp, // BCRec for this box
                                    const int /*orig_comp*/)
{
    BL_ASSERT(m_func != nullptr || m_func3D != nullptr);

    const int*  lo = dest.loVect();
    const int*  hi = dest.hiVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    const int* dom_hi = domain.hiVect();
    const Real* dx = geom.CellSize();

    Real grd_lo[AMREX_SPACEDIM];
    const Real* problo = geom.ProbLo();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        grd_lo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }

    static_assert(sizeof(BCRec) == 2*AMREX_SPACEDIM*sizeof(int),
                  "Let us know if this assertion fails");

    for (int icomp = 0; icomp < numcomp; ++icomp)
    {
        Real* data = dest.dataPtr(dcomp+icomp);
        if (m_func != nullptr) {
            m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),
                   dom_lo,dom_hi,
                   dx,grd_lo,&time,bcr[bcomp+icomp].vect(),dcomp+icomp);
        } else {
            m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
                     AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                     AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcr[bcomp+icomp].vect());
        }
    }
}

using PhysBCFunctMaestro = amrex::PhysBCFunct<amrex::BndryFuncArrayMaestro>;

#endif
