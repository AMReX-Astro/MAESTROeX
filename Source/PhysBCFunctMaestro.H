#ifndef PhysBCFunctMaestro_H_
#define PhysBCFunctMaestro_H_

#include <AMReX_PhysBCFunct.H>

extern "C" {
typedef void (*BndryFuncDefaultMaestro)(
    const amrex::Array4<amrex::Real>& scal, const amrex::Box& bx,
    const amrex::Box& domainBox, const amrex::Real* dx, const amrex::BCRec bcs,
    const amrex::Real* gridlo, const int comp);
}

/// This version calls function working on array
class BndryFuncArrayMaestro : public amrex::BndryFuncArray {
   public:
    BndryFuncArrayMaestro() noexcept {}
    BndryFuncArrayMaestro(BndryFuncDefaultMaestro inFunc) noexcept
        : m_func(inFunc) {}

    void operator()(amrex::Box const& bx, amrex::FArrayBox& dest,
                    const int dcomp, const int numcomp,
                    amrex::Geometry const& geom, [[maybe_unused]] const amrex::Real time,
                    const amrex::Vector<amrex::BCRec>& bcr,
                    const int bcomp,  // BCRec for this box
                    [[maybe_unused]] const int orig_comp) {
        BL_ASSERT(m_func != nullptr);

        const int* lo = dest.loVect();
        const amrex::Box& domain = geom.Domain();
        const int* dom_lo = domain.loVect();
        const amrex::Real* dx = geom.CellSize();

        amrex::Real grd_lo[AMREX_SPACEDIM];
        const amrex::Real* problo = geom.ProbLo();
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            grd_lo[i] = problo[i] + dx[i] * (lo[i] - dom_lo[i]);
        }

        static_assert(sizeof(amrex::BCRec) == 2 * AMREX_SPACEDIM * sizeof(int),
                      "Let us know if this assertion fails");

        for (int icomp = 0; icomp < numcomp; ++icomp) {
            m_func(dest.array(dcomp + icomp), bx, domain, dx,
                   bcr[bcomp + icomp], grd_lo, dcomp + icomp);
        }
    }

   protected:
    BndryFuncDefaultMaestro m_func = nullptr;
};

using PhysBCFunctMaestro = amrex::PhysBCFunct<BndryFuncArrayMaestro>;

#endif
