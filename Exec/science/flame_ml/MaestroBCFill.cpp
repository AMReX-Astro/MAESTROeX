
#include <Maestro.H>
#include <Maestro_F.H>

#include <MaestroInletBCs.H>

using namespace amrex;
using namespace InletBCs;

void Maestro::ScalarFill(const Array4<Real>& scal, const Box& bx,
                         const Box& domainBox, const Real* dx, const BCRec bcs,
                         const Real* gridlo, const int comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ScalarFill()", ScalarFill);

    fab_filcc(bx, scal, 1, domainBox, dx, gridlo, &bcs);

    FillExtBC(scal, bx, domainBox, dx, bcs, comp, false);
}

void Maestro::VelFill(const Array4<Real>& vel, const Box& bx,
                      const Box& domainBox, const Real* dx, const BCRec bcs,
                      const Real* gridlo, const int comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::VelFill()", VelFill);

    fab_filcc(bx, vel, 1, domainBox, dx, gridlo, &bcs);

    FillExtBC(vel, bx, domainBox, dx, bcs, comp, true);
}

void Maestro::FillExtBC(const Array4<Real>& q, const Box& bx,
                        const Box& domainBox, const Real* dx, const BCRec bcs,
                        const int comp, const bool is_vel) {
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

    // get parameters for EXT_DIR bcs
    const Real INLET_RHO_l = INLET_RHO;
    const Real INLET_RHOH_l = INLET_RHOH;
    const Real INLET_TEMP_l = INLET_TEMP;
    RealVector INLET_RHOX_l(NumSpec);
    for (int i = 0; i < NumSpec; ++i) INLET_RHOX_l[i] = INLET_RHOX[i];
    const Real INLET_VEL_l = INLET_VEL;

    const auto domlo = domainBox.loVect3d();
    const auto domhi = domainBox.hiVect3d();

    const auto lo = bx.loVect3d();
    const auto hi = bx.hiVect3d();

    if (lo[0] < domlo[0]) {
        auto imax = domlo[0];

        if (bcs.lo(0) == BCType::ext_dir) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (i < imax) {
                    q(i, j, k) = 1.0;
                }
            });
            Gpu::synchronize();
        }
    }

    if (hi[0] > domhi[0]) {
        auto imin = domhi[0];

        if (bcs.hi(0) == BCType::ext_dir) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (i > imin) {
                    q(i, j, k) = 1.0;
                }
            });
            Gpu::synchronize();
        }
    }

    if (lo[1] < domlo[1]) {
        auto jmax = domlo[1];

        if (bcs.lo(1) == BCType::ext_dir) {
            if (is_vel) {
                if (comp == 0) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (j < jmax) {
                            q(i, j, k) = 0.0;
                        }
                    });
                } else if (comp == 1) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (j < jmax) {
                            q(i, j, k) = INLET_VEL_l;
                        }
                    });
                }
                Gpu::synchronize();
            } else {
                if (INLET_VEL != 0.0) {
                    if (comp == Rho) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (j < jmax) {
                                            q(i, j, k) = INLET_RHO_l;
                                        }
                                    });
                    } else if (comp == RhoH) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (j < jmax) {
                                            q(i, j, k) = INLET_RHOH_l;
                                        }
                                    });
                    } else if (comp >= FirstSpec &&
                               comp < FirstSpec + NumSpec) {
                        Real* AMREX_RESTRICT INLET_RHOX_p =
                            INLET_RHOX.dataPtr();
                        ParallelFor(
                            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                if (j < jmax) {
                                    q(i, j, k) = INLET_RHOX_p[comp - FirstSpec];
                                }
                            });
                    } else if (comp == Temp) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (j < jmax) {
                                            q(i, j, k) = INLET_TEMP_l;
                                        }
                                    });
                    } else {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (j < jmax) {
                                            q(i, j, k) = 0.0;
                                        }
                                    });
                    }
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (j < jmax) {
                            q(i, j, k) = 0.0;
                        }
                    });
                }
                Gpu::synchronize();
            }  // end if is_vel
        }
    }

    if (hi[1] > domhi[1]) {
        auto jmin = domhi[1];

        if (bcs.hi(1) == BCType::ext_dir) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (j > jmin) {
                    q(i, j, k) = 1.0;
                }
            });
        }
        Gpu::synchronize();
    }

#if AMREX_SPACEDIM == 3

    if (lo[2] < domlo[2]) {
        auto kmax = domlo[2];

        if (bcs.lo(2) == BCType::ext_dir) {
            if (is_vel) {
                if (comp == 0 || comp == 1) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (k < kmax) {
                            q(i, j, k) = 0.0;
                        }
                    });
                } else if (comp == 2) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (k < kmax) {
                            q(i, j, k) = INLET_VEL_l;
                        }
                    });
                }
            } else {
                if (INLET_VEL != 0.0) {
                    if (comp == Rho) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (k < kmax) {
                                            q(i, j, k) = INLET_RHO_l;
                                        }
                                    });
                    } else if (comp == RhoH) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (k < kmax) {
                                            q(i, j, k) = INLET_RHOH_l;
                                        }
                                    });
                    } else if (comp >= FirstSpec &&
                               comp < FirstSpec + NumSpec) {
                        Real* AMREX_RESTRICT INLET_RHOX_p =
                            INLET_RHOX.dataPtr();
                        ParallelFor(
                            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                if (k < kmax) {
                                    q(i, j, k) = INLET_RHOX_p[comp - FirstSpec];
                                }
                            });
                    } else if (comp == Temp) {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (k < kmax) {
                                            q(i, j, k) = INLET_TEMP_l;
                                        }
                                    });
                    } else {
                        ParallelFor(bx,
                                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        if (k < kmax) {
                                            q(i, j, k) = 0.0;
                                        }
                                    });
                    }
                } else {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        if (k < kmax) {
                            q(i, j, k) = 0.0;
                        }
                    });
                }
            }  // end if is_vel
            Gpu::synchronize();
        }
    }

    if (hi[2] > domhi[2]) {
        auto kmin = domhi[2];

        if (bcs.hi(2) == BCType::ext_dir) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k > kmin) {
                    q(i, j, k) = 0.0;
                }
            });
        }
        Gpu::synchronize();
    }
#endif
}
