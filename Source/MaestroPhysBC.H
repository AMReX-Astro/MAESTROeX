#ifndef MaestroPhysBC_H_
#define MaestroPhysBC_H_

#include <AMReX_PhysBCFunct.H>

using namespace amrex;

class MaestroPhysBC
    : public PhysBCFunctBase
{
public:
    MaestroPhysBC () {}
    virtual ~MaestroPhysBC () {}
    MaestroPhysBC (const Geometry& geom, const Array<BCRec>& bcr);
    void define (const Geometry& geom, const Array<BCRec>& bcr);
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time) override;
    using PhysBCFunctBase::FillBoundary;
private:
    Geometry     m_geom;
    Array<BCRec> m_bcr;
};

#endif
