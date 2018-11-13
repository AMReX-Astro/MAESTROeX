
#include <PhysBCFunctMaestro.H>

using namespace amrex;

#if 0

PhysBCFunctMaestro::PhysBCFunctMaestro (const Geometry& geom,
                                        const Vector<BCRec>& bcr,
                                        const BndryFunctBase& func)
    : m_geom(geom), m_bcr(bcr), m_bc_func(func.clone())
{
}

void
PhysBCFunctMaestro::define (const Geometry& geom,
                            const Vector<BCRec>& bcr,
                            const BndryFunctBase& func)
{
    m_geom = geom;
    m_bcr = bcr;
    m_bc_func.reset(func.clone());
}

void
PhysBCFunctMaestro::FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time)
{
    // timer for profiling
    BL_PROFILE_VAR("PhysBCFunctMaestro::FillBoundary",PhysBC_FillBoundary);

    if (mf.nGrow() == 0) return;

    if (m_geom.isAllPeriodic()) return;

    const Box&     domain      = m_geom.Domain();
    const int*     dlo         = domain.loVect();
    const int*     dhi         = domain.hiVect();
    const Real*    dx          = m_geom.CellSize();
    const RealBox& prob_domain = m_geom.ProbDomain();
    const Real*    problo      = prob_domain.lo();

    // create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomain.grow(i, mf.nGrow());
        }
    }

    for (int comp=dcomp; comp<dcomp+ncomp; ++comp)
    {

        int bccomp = comp-dcomp;

#ifdef _OPENMP
#pragma omp parallel
#endif
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& dest = mf[mfi];
            const Box& bx = dest.box();

            // if there are cells not in the valid + periodic grown box
            // we need to fill them here
            if (!gdomain.contains(bx)) {

                const int* fablo = bx.loVect();
                const int* fabhi = bx.hiVect();

                Real xlo[AMREX_SPACEDIM];
                for (int i = 0; i < AMREX_SPACEDIM; i++) {
                    xlo[i] = problo[i] + dx[i]*(fablo[i]-dlo[i]);
                }

                (*m_bc_func)(dest.dataPtr(comp), fablo, fabhi, dlo, dhi,
                             dx, xlo, &time, m_bcr[bccomp].vect());
            }
        }
    }
}

#endif
