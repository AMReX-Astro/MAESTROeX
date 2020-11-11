
#include <Maestro.H>

using namespace amrex;

void Maestro::Slopex(const Box& bx, Array4<Real> const s,
                     Array4<Real> const slx, const Box& domainBox,
                     const Vector<BCRec>& bcs, int ncomp, int bc_start_comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Slopex()", Slopex);

    /////////////
    // x-dir
    /////////////

    int ilo = domainBox.loVect()[0];
    int ihi = domainBox.hiVect()[0];

    // create lo and hi vectors
    IntVector bclo(ncomp);
    IntVector bchi(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        bclo[i] = bcs[bc_start_comp + i].lo()[0];
        bchi[i] = bcs[bc_start_comp + i].hi()[0];
    }

    int* AMREX_RESTRICT bclo_p = bclo.dataPtr();
    int* AMREX_RESTRICT bchi_p = bchi.dataPtr();

    if (slope_order == 0) {
        // 1st order

        ParallelFor(bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        slx(i, j, k, n) = 0.0;
                    });

    } else if (slope_order == 2) {
        // 2nd order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                Real del = 0.5 * (s(i + 1, j, k, n) - s(i - 1, j, k, n));
                Real dpls = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                Real dmin = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                Real slim =
                    amrex::min(amrex::Math::abs(dpls), amrex::Math::abs(dmin));
                slim = dpls * dmin > 0.0 ? slim : 0.0;
                Real sflag = amrex::Math::copysign(1.0, del);
                slx(i, j, k, n) =
                    sflag * amrex::min(slim, amrex::Math::abs(del));

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (i == ilo - 1) {
                        slx(i, j, k, n) = 0.0;
                    } else if (i == ilo) {
                        del = (s(i + 1, j, k, n) + 3.0 * s(i, j, k, n) -
                               4 * s(i - 1, j, k, n)) /
                              3.0;
                        dpls = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        slx(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (i == ihi + 1) {
                        slx(i, j, k, n) = 0.0;
                    } else if (i == ihi) {
                        del = -(s(i - 1, j, k, n) + 3.0 * s(i, j, k, n) -
                                4 * s(i + 1, j, k, n)) /
                              3.0;
                        dpls = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                        dmin = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        slx(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }
            });

    } else if (slope_order == 4) {
        // 4th order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                // left
                Real dcen = 0.5 * (s(i, j, k, n) - s(i - 2, j, k, n));
                Real dmin = 2.0 * (s(i - 1, j, k, n) - s(i - 2, j, k, n));
                Real dpls = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                Real dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                Real dflag = amrex::Math::copysign(1.0, dcen);
                Real dxl = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // right
                dcen = 0.5 * (s(i + 2, j, k, n) - s(i, j, k, n));
                dmin = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                dpls = 2.0 * (s(i + 2, j, k, n) - s(i + 1, j, k, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);
                Real dxr = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // center
                dcen = 0.5 * (s(i + 1, j, k, n) - s(i - 1, j, k, n));
                dmin = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                dpls = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);

                Real ds = 4.0 / 3.0 * dcen - (dxr + dxl) / 6.0;
                slx(i, j, k, n) =
                    dflag * amrex::min(amrex::Math::abs(ds), dlim);

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (i == ilo - 1) {
                        slx(i, j, k, n) = 0.0;
                    } else if (i == ilo) {
                        Real del = -16.0 / 15.0 * s(i - 1, j, k, n) +
                                   0.5 * s(i, j, k, n) +
                                   2.0 / 3.0 * s(i + 1, j, k, n) -
                                   0.1 * s(i + 2, j, k, n);
                        dmin = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                        dpls = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? dlim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        slx(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    } else if (i == ilo + 1) {
                        // Recalculate the slope at lo(1)+1 using the revised dxl
                        Real del = -16.0 / 15.0 * s(i - 2, j, k, n) +
                                   0.5 * s(i - 1, j, k, n) +
                                   2.0 / 3.0 * s(i, j, k, n) -
                                   0.1 * s(i + 1, j, k, n);
                        dmin = 2.0 * (s(i - 1, j, k, n) - s(i - 2, j, k, n));
                        dpls = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? dlim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        dxl = sflag * amrex::min(slim, amrex::Math::abs(del));

                        ds = 4.0 / 3.0 * dcen - (dxr + dxl) / 6.0;
                        slx(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (i == ihi + 1) {
                        slx(i, j, k, n) = 0.0;
                    } else if (i == ihi) {
                        Real del = -(-16.0 / 15.0 * s(i + 1, j, k, n) +
                                     0.5 * s(i, j, k, n) +
                                     2.0 / 3.0 * s(i - 1, j, k, n) -
                                     0.1 * s(i - 2, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i - 1, j, k, n));
                        dpls = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? dlim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        slx(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    } else if (i == ihi - 1) {
                        // Recalculate the slope at hi(1)-1 using the revised dxr
                        Real del = -(-16.0 / 15.0 * s(i + 2, j, k, n) +
                                     0.5 * s(i + 1, j, k, n) +
                                     2.0 / 3.0 * s(i, j, k, n) -
                                     0.1 * s(i - 1, j, k, n));
                        dmin = 2.0 * (s(i + 1, j, k, n) - s(i, j, k, n));
                        dpls = 2.0 * (s(i + 2, j, k, n) - s(i + 1, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? dlim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        dxr = sflag * amrex::min(slim, amrex::Math::abs(del));

                        ds = 4.0 / 3.0 * dcen - (dxl + dxr) / 6.0;
                        slx(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }
            });
    }
}

void Maestro::Slopey(const Box& bx, Array4<Real> const s,
                     Array4<Real> const sly, const Box& domainBox,
                     const Vector<BCRec>& bcs, int ncomp, int bc_start_comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Slopey()", Slopey);

    /////////////
    // y-dir
    /////////////

    int jlo = domainBox.loVect()[1];
    int jhi = domainBox.hiVect()[1];

    // create lo and hi vectors
    IntVector bclo(ncomp);
    IntVector bchi(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        bclo[i] = bcs[bc_start_comp + i].lo()[1];
        bchi[i] = bcs[bc_start_comp + i].hi()[1];
    }

    int* AMREX_RESTRICT bclo_p = bclo.dataPtr();
    int* AMREX_RESTRICT bchi_p = bchi.dataPtr();

    if (slope_order == 0) {
        // 1st order

        ParallelFor(bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        sly(i, j, k, n) = 0.0;
                    });

    } else if (slope_order == 2) {
        // 2nd order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                Real del = 0.5 * (s(i, j + 1, k, n) - s(i, j - 1, k, n));
                Real dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                Real dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                Real slim =
                    amrex::min(amrex::Math::abs(dpls), amrex::Math::abs(dmin));
                slim = dpls * dmin > 0.0 ? slim : 0.0;
                Real sflag = amrex::Math::copysign(1.0, del);
                sly(i, j, k, n) =
                    sflag * amrex::min(slim, amrex::Math::abs(del));

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (j == jlo - 1) {
                        sly(i, j, k, n) = 0.0;
                    } else if (j == jlo) {
                        del = (s(i, j + 1, k, n) + 3.0 * s(i, j, k, n) -
                               4.0 * s(i, j - 1, k, n)) /
                              3.0;
                        dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        sly(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (j == jhi + 1) {
                        sly(i, j, k, n) = 0.0;
                    } else if (j == jhi) {
                        del = -(s(i, j - 1, k, n) + 3.0 * s(i, j, k, n) -
                                4.0 * s(i, j + 1, k, n)) /
                              3.0;
                        dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        sly(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }
            });

    } else if (slope_order == 4) {
        // 4th order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                // left
                Real dcen = 0.5 * (s(i, j, k, n) - s(i, j - 2, k, n));
                Real dmin = 2.0 * (s(i, j - 1, k, n) - s(i, j - 2, k, n));
                Real dpls = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                Real dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                Real dflag = amrex::Math::copysign(1.0, dcen);
                Real dyl = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // right
                dcen = 0.5 * (s(i, j + 2, k, n) - s(i, j, k, n));
                dmin = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                dpls = 2.0 * (s(i, j + 2, k, n) - s(i, j + 1, k, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);
                Real dyr = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // center
                dcen = 0.5 * (s(i, j + 1, k, n) - s(i, j - 1, k, n));
                dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);

                Real ds = 4.0 / 3.0 * dcen - (dyr + dyl) / 6.0;
                sly(i, j, k, n) =
                    dflag * amrex::min(amrex::Math::abs(ds), dlim);

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (j == jlo - 1) {
                        sly(i, j, k, n) = 0.0;
                    } else if (j == jlo) {
                        Real del = -16.0 / 15.0 * s(i, j - 1, k, n) +
                                   0.5 * s(i, j, k, n) +
                                   2.0 / 3.0 * s(i, j + 1, k, n) -
                                   0.1 * s(i, j + 2, k, n);
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                        dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        sly(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    } else if (j == jlo + 1) {
                        // Recalculate the slope at lo(2)+1 using the revised dyl
                        Real del = -16.0 / 15.0 * s(i, j - 2, k, n) +
                                   0.5 * s(i, j - 1, k, n) +
                                   2.0 / 3.0 * s(i, j, k, n) -
                                   0.1 * s(i, j + 1, k, n);
                        dmin = 2.0 * (s(i, j - 1, k, n) - s(i, j - 2, k, n));
                        dpls = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        dyl = sflag * amrex::min(slim, amrex::Math::abs(del));
                        ds = 4.0 / 3.0 * dcen - (dyr + dyl) / 6.0;
                        sly(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (j == jhi + 1) {
                        sly(i, j, k, n) = 0.0;
                    } else if (j == jhi) {
                        Real del = -(-16.0 / 15.0 * s(i, j + 1, k, n) +
                                     0.5 * s(i, j, k, n) +
                                     2.0 / 3.0 * s(i, j - 1, k, n) -
                                     0.1 * s(i, j - 2, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j - 1, k, n));
                        dpls = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        sly(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    } else if (j == jhi - 1) {
                        // Recalculate the slope at lo(2)+1 using the revised dyr
                        Real del = -(-16.0 / 15.0 * s(i, j + 2, k, n) +
                                     0.5 * s(i, j + 1, k, n) +
                                     2.0 / 3.0 * s(i, j, k, n) -
                                     0.1 * s(i, j - 1, k, n));
                        dmin = 2.0 * (s(i, j + 1, k, n) - s(i, j, k, n));
                        dpls = 2.0 * (s(i, j + 2, k, n) - s(i, j + 1, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        dyr = sflag * amrex::min(slim, amrex::Math::abs(del));
                        ds = 4.0 / 3.0 * dcen - (dyl + dyr) / 6.0;
                        sly(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }
            });
    }
}

#if (AMREX_SPACEDIM == 3)
void Maestro::Slopez(const Box& bx, Array4<Real> const s,
                     Array4<Real> const slz, const Box& domainBox,
                     const Vector<BCRec>& bcs, int ncomp, int bc_start_comp) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Slopez()", Slopez);

    /////////////
    // z-dir
    /////////////

    int klo = domainBox.loVect()[2];
    int khi = domainBox.hiVect()[2];

    // create lo and hi vectors
    IntVector bclo(ncomp);
    IntVector bchi(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        bclo[i] = bcs[bc_start_comp + i].lo()[2];
        bchi[i] = bcs[bc_start_comp + i].hi()[2];
    }

    int* AMREX_RESTRICT bclo_p = bclo.dataPtr();
    int* AMREX_RESTRICT bchi_p = bchi.dataPtr();

    if (slope_order == 0) {
        // 1st order

        ParallelFor(bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        slz(i, j, k, n) = 0.0;
                    });

    } else if (slope_order == 2) {
        // 2nd order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                Real del = 0.5 * (s(i, j, k + 1, n) - s(i, j, k - 1, n));
                Real dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                Real dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                Real slim =
                    amrex::min(amrex::Math::abs(dpls), amrex::Math::abs(dmin));
                slim = dpls * dmin > 0.0 ? slim : 0.0;
                Real sflag = amrex::Math::copysign(1.0, del);
                slz(i, j, k, n) =
                    sflag * amrex::min(slim, amrex::Math::abs(del));

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (k == klo - 1) {
                        slz(i, j, k, n) = 0.0;
                    } else if (k == klo) {
                        del = (s(i, j, k + 1, n) + 3.0 * s(i, j, k, n) -
                               4.0 * s(i, j, k - 1, n)) /
                              3.0;
                        dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        slz(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (k == khi + 1) {
                        slz(i, j, k, n) = 0.0;
                    } else if (k == khi) {
                        del = -(s(i, j, k - 1, n) + 3.0 * s(i, j, k, n) -
                                4.0 * s(i, j, k + 1, n)) /
                              3.0;
                        dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                        slim = amrex::min(amrex::Math::abs(dpls),
                                          amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        sflag = amrex::Math::copysign(1.0, del);
                        slz(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));
                    }
                }
            });

    } else if (slope_order == 4) {
        // 4th order

        ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                // left
                Real dcen = 0.5 * (s(i, j, k, n) - s(i, j, k - 2, n));
                Real dmin = 2.0 * (s(i, j, k - 1, n) - s(i, j, k - 2, n));
                Real dpls = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                Real dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                Real dflag = amrex::Math::copysign(1.0, dcen);
                Real dzl = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // right
                dcen = 0.5 * (s(i, j, k + 2, n) - s(i, j, k, n));
                dmin = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                dpls = 2.0 * (s(i, j, k + 2, n) - s(i, j, k + 1, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);
                Real dzr = dflag * amrex::min(dlim, amrex::Math::abs(dcen));

                // center
                dcen = 0.5 * (s(i, j, k + 1, n) - s(i, j, k - 1, n));
                dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                dlim =
                    amrex::min(amrex::Math::abs(dmin), amrex::Math::abs(dpls));
                dlim = dpls * dmin > 0.0 ? dlim : 0.0;
                dflag = amrex::Math::copysign(1.0, dcen);

                Real ds = 4.0 / 3.0 * dcen - (dzr + dzl) / 6.0;
                slz(i, j, k, n) =
                    dflag * amrex::min(amrex::Math::abs(ds), dlim);

                if (bclo_p[n] == EXT_DIR || bclo_p[n] == HOEXTRAP) {
                    if (k == klo - 1) {
                        slz(i, j, k, n) = 0.0;
                    } else if (k == klo) {
                        Real del = -16.0 / 15.0 * s(i, j, k - 1, n) +
                                   0.5 * s(i, j, k, n) +
                                   2.0 / 3.0 * s(i, j, k + 1, n) -
                                   0.1 * s(i, j, k + 2, n);
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                        dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        slz(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));

                    } else if (k == klo + 1) {
                        // Recalculate the slope at lo(2)+1 using the revised dzl
                        dmin = 2.0 * (s(i, j, k - 1, n) - s(i, j, k - 2, n));
                        dpls = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        dzl = slz(i, j, k - 1, n);
                        ds = 4.0 / 3.0 * dcen - (dzr + dzl) / 6.0;
                        slz(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }

                if (bchi_p[n] == EXT_DIR || bchi_p[n] == HOEXTRAP) {
                    if (k == khi + 1) {
                        slz(i, j, k, n) = 0.0;
                    } else if (k == khi) {
                        Real del = -(-16.0 / 15.0 * s(i, j, k + 1, n) +
                                     0.5 * s(i, j, k, n) +
                                     2.0 / 3.0 * s(i, j, k - 1, n) -
                                     0.1 * s(i, j, k - 2, n));
                        dmin = 2.0 * (s(i, j, k, n) - s(i, j, k - 1, n));
                        dpls = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        Real sflag = amrex::Math::copysign(1.0, del);
                        slz(i, j, k, n) =
                            sflag * amrex::min(slim, amrex::Math::abs(del));

                    } else if (k == khi - 1) {
                        // Recalculate the slope at lo(3)+1 using the revised dzr
                        dmin = 2.0 * (s(i, j, k + 1, n) - s(i, j, k, n));
                        dpls = 2.0 * (s(i, j, k + 2, n) - s(i, j, k + 1, n));
                        Real slim = amrex::min(amrex::Math::abs(dpls),
                                               amrex::Math::abs(dmin));
                        slim = dpls * dmin > 0.0 ? slim : 0.0;
                        dzr = slz(i, j, k + 1, n);
                        ds = 4.0 / 3.0 * dcen - (dzl + dzr) / 6.0;
                        slz(i, j, k, n) =
                            dflag * amrex::min(amrex::Math::abs(ds), dlim);
                    }
                }
            });
    }
}
#endif