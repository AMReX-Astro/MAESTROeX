
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void
Maestro::ScalarFill(const Array4<Real>& scal, const Box& bx, 
                    const Box& domainBox, const Real* dx, 
                    const BCRec bcs, 
                    const Real* gridlo, const int comp)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ScalarFill()", ScalarFill);

    fab_filcc(bx, scal, 1, domainBox, dx, gridlo, &bcs);

    FillExtBC(scal, bx, domainBox, bcs, comp, false);
}


void
Maestro::VelFill(const Array4<Real>& vel, 
                 const Box& bx, const Box& domainBox, 
                 const Real* dx, 
                 const BCRec bcs, const Real * gridlo)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelFill()", VelFill);

    fab_filcc(bx, vel, 1, domainBox, dx, gridlo, &bcs);

    FillExtBC(vel, bx, domainBox, bcs, 0, true);
}

void
Maestro::FillExtBC(const Array4<Real>& q, const Box& bx, 
                   const Box& domainBox, const BCRec bcs, 
                   const int comp, const bool is_vel)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::FillExtBC()", FillExtBC);

    // do nothing if there are no exterior boundaries
    bool found_ext_boundary = false;
    for (auto i = 0; i < AMREX_SPACEDIM; ++i) {
        if (bcs.lo(i) == BCType::ext_dir || bcs.hi(i) == BCType::ext_dir) {
            found_ext_boundary = true;
            break;
        }
    }
    if (!found_ext_boundary) return;

    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    const auto lo = bx.loVect3d();
    const auto hi = bx.hiVect3d();

    const Box qbx(q);
    const auto q_lo = qbx.loVect3d();
    const auto q_hi = qbx.hiVect3d();

//     const auto is = max(q_lo[0], domlo[0]);
//     const auto ie = min(q_hi[0], domhi[0]);

//     const auto js = max(q_lo[1], domlo[1]);
//     const auto je = min(q_hi[1], domhi[1]); 

// #if AMREX_SPACEDIM == 3
//     const auto ks = max(q_lo[2], domlo[2]);
//     const auto ke = min(q_hi[2], domhi[2]); 
// #endif

    if (lo[0] < domlo[0]) {
        auto imin = lo[0];
        auto imax = domlo[0];

        if (bcs.lo(0) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (i < imax) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[0,0] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }

    if (hi[0] > domhi[0]) {
        auto imin = domhi[0]+1;
        auto imax = hi[0];

        if (bcs.hi(0) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (i > imin) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[0,1] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }

    if (lo[1] < domlo[1]) {
        auto jmin = lo[1];
        auto jmax = domlo[1];

        if (bcs.lo(1) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (j < jmax) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[1,0] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }

    if (hi[1] > domhi[1]) {
        auto jmin = domhi[1]+1;
        auto jmax = hi[1];

        if (bcs.hi(1) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (j > jmin) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[1,1] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }

#if AMREX_SPACEDIM == 3

    if (lo[2] < domlo[2]) {
        auto kmin = lo[2];
        auto kmax = domlo[2];

        if (bcs.lo(2) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (k < kmax) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[2,0] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }

    if (hi[2] > domhi[2]) {
        auto kmin = domhi[2]+1;
        auto kmax = hi[2];

        if (bcs.hi(2) == BCType::ext_dir) {
            if (comp == Pi || is_vel) {
                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    if (k > kmin) {
                        q(i,j,k) = 0.0;
                    }
                })
            } else {
                Abort("MaestroBCFill bc[2,1] - must supply Dirichlet boundary conditions for scalar");
            } 
        }
    }
#endif
}