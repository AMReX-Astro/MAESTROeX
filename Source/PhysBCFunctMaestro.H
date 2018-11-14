#ifndef PhysBCFunctMaestro_H_
#define PhysBCFunctMaestro_H_

#include <AMReX_PhysBCFunct.H>

using PhysBCFunctMaestro = amrex::PhysBCFunct<amrex::BndryFuncArray>;

#if 0
class PhysBCFunctMaestro
    : public amrex::PhysBCFunctBase
{
public:
    PhysBCFunctMaestro () {}
    virtual ~PhysBCFunctMaestro () {}

    PhysBCFunctMaestro (const amrex::Geometry& geom,
                        const amrex::Vector<amrex::BCRec>& bcr,
                        const amrex::BndryFunctBase& func);

    void define (const amrex::Geometry& geom,
                 const amrex::Vector<amrex::BCRec>& bcr,
                 const amrex::BndryFunctBase& func);

    virtual void FillBoundary (amrex::MultiFab& mf,
                               int dcomp,
                               int ncomp,
                               amrex::Real time) override;

    using PhysBCFunctBase::FillBoundary;

private:
    amrex::Geometry                        m_geom;
    amrex::Vector<amrex::BCRec>            m_bcr;
    std::unique_ptr<amrex::BndryFunctBase> m_bc_func;
};
#endif

#endif
