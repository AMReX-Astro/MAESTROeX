#ifndef PhysBCFunctMaestro_H_
#define PhysBCFunctMaestro_H_

#include <AMReX_PhysBCFunct.H>

extern "C"
{
typedef void (*BndryFuncDefaultMaestro)(amrex::Real* data,
                                        AMREX_ARLIM_P (lo), AMREX_ARLIM_P (hi),
                                        const int* dom_lo, const int* dom_hi,
                                        const amrex::Real* dx, const amrex::Real* grd_lo,
                                        const amrex::Real* time, const int* bc,
                                        const int icomp);
typedef void (*BndryFunc3DDefaultMaestro)(amrex::Real* data, const int* lo, const int* hi,
                                          const int* dom_lo, const int* dom_hi,
                                          const amrex::Real* dx, const amrex::Real* grd_lo,
                                          const amrex::Real* time, const int* bc,
                                          const int icomp);
}

/// This version calls function working on array
class BndryFuncArrayMaestro
    : public amrex::BndryFuncArray
{
public:
    BndryFuncArrayMaestro () noexcept {}
	BndryFuncArrayMaestro (BndryFuncDefaultMaestro inFunc) noexcept: m_func(inFunc) {
	}
	BndryFuncArrayMaestro (BndryFunc3DDefaultMaestro inFunc) noexcept: m_func3D(inFunc) {
	}

    void operator () (amrex::Box const& bx, amrex::FArrayBox& dest,
                        const int dcomp, const int numcomp,
                        amrex::Geometry const& geom, const amrex::Real time,
                        const amrex::Vector<amrex::BCRec>&bcr,
                        const int bcomp,  // BCRec for this box
                        const int orig_comp)
    {
    	BL_ASSERT(m_func != nullptr || m_func3D != nullptr);

    	const int*  lo = dest.loVect();
    	const int*  hi = dest.hiVect();
    	const amrex::Box& domain = geom.Domain();
    	const int* dom_lo = domain.loVect();
    	const int* dom_hi = domain.hiVect();
    	const amrex::Real* dx = geom.CellSize();

    	amrex::Real grd_lo[AMREX_SPACEDIM];
    	const amrex::Real* problo = geom.ProbLo();
    	for (int i = 0; i < AMREX_SPACEDIM; i++)
    	{
    		grd_lo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    	}

    	static_assert(sizeof(amrex::BCRec) == 2*AMREX_SPACEDIM*sizeof(int),
    	              "Let us know if this assertion fails");

    	for (int icomp = 0; icomp < numcomp; ++icomp)
    	{
    		amrex::Real* data = dest.dataPtr(dcomp+icomp);
    		if (m_func != nullptr) {
    			m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),
    			       dom_lo,dom_hi,
    			       dx,grd_lo,&time,bcr[bcomp+icomp].vect(),dcomp+icomp+1);
                        // pass in +1 for the component since fortran uses 1-based component index
    		} else {
    			m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
    			         AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
    			         AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),&time,bcr[bcomp+icomp].vect(),
    			         dcomp+icomp+1);
                        // pass in +1 for the component since fortran uses 1-based component index
    		}
    	}
    }

protected:
	BndryFuncDefaultMaestro m_func     = nullptr;
	BndryFunc3DDefaultMaestro m_func3D = nullptr;
};

using PhysBCFunctMaestro = amrex::PhysBCFunct<BndryFuncArrayMaestro>;

#endif
