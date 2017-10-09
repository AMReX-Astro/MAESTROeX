#ifndef MaestroPhysBCFunct_H_
#define MaestroPhysBCFunct_H_

#include <AMReX_PhysBCFunct.H>

class MaestroPhysBCFunct
    : public amrex::PhysBCFunctBase
{
public:
    MaestroPhysBCFunct () {}
    virtual ~MaestroPhysBCFunct () {}

    MaestroPhysBCFunct (const amrex::Geometry& geom,
                        const amrex::Vector<amrex::BCRec>& bcr,
                        const amrex::BndryFunctBase& func);

    void define (const amrex::Geometry& geom,
                 const amrex::Vector<amrex::BCRec>& bcr,
                 const amrex::BndryFunctBase& func);

    virtual void FillBoundary (amrex::MultiFab& mf,
                               int dcomp,
                               int ncomp,
                               Real time) override;

    using PhysBCFunctBase::FillBoundary;

private:
    amrex::Geometry                        m_geom;
    amrex::Vector<amrex::BCRec>            m_bcr;
    std::unique_ptr<amrex::BndryFunctBase> m_bc_func;
};

#endif
