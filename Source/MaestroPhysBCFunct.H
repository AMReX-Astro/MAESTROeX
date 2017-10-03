#ifndef MaestroPhysBCFunct_H_
#define MaestroPhysBCFunct_H_

#include <AMReX_PhysBCFunct.H>

using namespace amrex;

class MaestroPhysBCFunct
    : public PhysBCFunctBase
{
public:
    MaestroPhysBCFunct () {}
    virtual ~MaestroPhysBCFunct () {}
    MaestroPhysBCFunct (const Geometry& geom, const Array<BCRec>& bcr, const BndryFunctBase& func);
    void define (const Geometry& geom, const Array<BCRec>& bcr, const BndryFunctBase& func);
    virtual void FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time) override;
    using PhysBCFunctBase::FillBoundary;
private:
    Geometry     m_geom;
    Array<BCRec> m_bcr;
    std::unique_ptr<BndryFunctBase> m_bc_func;
};

#endif
