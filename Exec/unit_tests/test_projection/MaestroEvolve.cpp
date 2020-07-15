
#include <Maestro.H>
#include <Problem_F.H>

using namespace amrex;

// advance solution to final time
void Maestro::Evolve() {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::Evolve()", Evolve);

    Print() << "Calling Evolve()" << std::endl;

    int project_type;
    get_project_type(&project_type);

    // -------------------------------------------------------------------------
    //  allocate arrays
    // -------------------------------------------------------------------------

    Print() << "...allocate arrays" << std::endl;

    Vector<MultiFab> umid(finest_level + 1);
    Vector<MultiFab> gphi(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_old(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_mid(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > umac_new(finest_level + 1);
    Vector<std::array<MultiFab, AMREX_SPACEDIM> > gphi_mac(finest_level + 1);
    Vector<MultiFab> utemp(finest_level + 1);

    if (project_type == 1) {
        // HG projection.  Velocities are cell-centered
        for (int lev = 0; lev <= finest_level; ++lev) {
            umid[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
            gphi[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);

            umid[lev].setVal(0.);
            gphi[lev].setVal(0.);
        }
    } else {
        // MAC projection.  Velocities are nodal in respective dimension
        for (int lev = 0; lev <= finest_level; ++lev) {
            AMREX_D_TERM(
                umac_old[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_old[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_old[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                umac_mid[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_mid[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_mid[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                umac_new[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , umac_new[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , umac_new[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););
            AMREX_D_TERM(
                gphi_mac[lev][0].define(convert(grids[lev], nodal_flag_x),
                                        dmap[lev], 1, ng_s);
                , gphi_mac[lev][1].define(convert(grids[lev], nodal_flag_y),
                                          dmap[lev], 1, ng_s);
                , gphi_mac[lev][2].define(convert(grids[lev], nodal_flag_z),
                                          dmap[lev], 1, ng_s););

            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                umac_old[lev][d].setVal(0.);
                umac_mid[lev][d].setVal(0.);
                umac_new[lev][d].setVal(0.);
                gphi_mac[lev][d].setVal(0.);
            }
            utemp[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, ng_s);
            utemp[lev].setVal(0.);
        }
    }

    // -------------------------------------------------------------------------
    //  initialize velocity field
    // -------------------------------------------------------------------------

    Print() << "...initialize velocity field" << std::endl;

    if (project_type == 1) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& vel = uold[lev];
            const Real* dx = geom[lev].CellSize();

            // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(vel, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                init_vel(ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_3D(vel[mfi]),
                         ZFILL(dx));
            }
        }

        AverageDown(uold, 0, AMREX_SPACEDIM);
        FillPatch(t_old, uold, uold, uold, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

    } else {
        // need to initialize the mac velocity
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_old[lev][0];
            MultiFab& vmac_mf = umac_old[lev][1];
            MultiFab& wmac_mf = umac_old[lev][2];
            const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                init_mac_vel(ARLIM_3D(lo), ARLIM_3D(hi),
                             BL_TO_FORTRAN_3D(umac_mf[mfi]),
                             BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                             BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                             ZFILL(dx));
            }
        }

        if (finest_level == 0) {
            // fill periodic ghost cells
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    umac_old[lev][d].FillBoundary(geom[lev].periodicity());
                }
            }
            // fill ghost cells behind physical boundaries
            FillUmacGhost(umac_old);
        } else {
            // edge_restriction for velocities
            AverageDownFaces(umac_old);
            // fill level n ghost cells using interpolation from level n-1 data
            FillPatchUedge(umac_old);
        }
    }

    BaseState<Real> dummy;

    if (project_type == 1) {
        // write uold
        WritePlotFile(0, t_new, dt, dummy, dummy, dummy, dummy, uold, uold,
                      uold);

        // copy the velocity field over to the intermediate state, umid
        for (int lev = 0; lev <= finest_level; ++lev)
            MultiFab::Copy(umid[lev], uold[lev], 0, 0, AMREX_SPACEDIM, ng_s);
    } else {
        // convect MAC field to cell-centered

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_old[lev][0];
            MultiFab& vmac_mf = umac_old[lev][1];
            MultiFab& wmac_mf = umac_old[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
                                  BL_TO_FORTRAN_3D(umac_mf[mfi]),
                                  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                                  BL_TO_FORTRAN_3D(utemp_mf[mfi]));
            }
        }
        // write umac_old
        WritePlotFile(0, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // copy the velocity field over to the intermediate state, umid
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                MultiFab::Copy(umac_mid[lev][comp], umac_old[lev][comp], 0, 0,
                               1, ng_s);
        }
    }

    //--------------------------------------------------------------------------
    // 'pollute' the velocity field by adding the gradient of a scalar
    //--------------------------------------------------------------------------

    Print() << "...pollute velocity field" << std::endl;

    if (project_type == 1) {
        // add_grad_scalar

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& gphi_mf = gphi[lev];
            MultiFab& umid_mf = umid[lev];
            const Real* dx = geom[lev].CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(gphi_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                add_grad_scalar(ARLIM_3D(lo), ARLIM_3D(hi),
                                BL_TO_FORTRAN_FAB(gphi_mf[mfi]),
                                BL_TO_FORTRAN_FAB(umid_mf[mfi]),
                                phys_bc.dataPtr(), ZFILL(dx));
            }
        }

        // fill ghosts and boundary
        AverageDown(umid, 0, AMREX_SPACEDIM);
        FillPatch(t_old, umid, umid, umid, 0, 0, AMREX_SPACEDIM, 0, bcs_u);

        // write umid
        WritePlotFile(1, t_new, dt, dummy, dummy, dummy, dummy, umid, uold,
                      uold);
        // write gphi
        WritePlotFile(2, t_new, dt, dummy, dummy, dummy, dummy, gphi, uold,
                      uold);

        // copy the velocity field over to the final state, unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(unew[lev], umid[lev], 0, 0, AMREX_SPACEDIM, ng_s);

            // swap pointers so NodalProj works properly
            std::swap(uold[lev], unew[lev]);
        }

    } else {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& gphix_mac_mf = gphi_mac[lev][0];
            MultiFab& gphiy_mac_mf = gphi_mac[lev][1];
#if (AMREX_SPACEDIM == 3)
            MultiFab& gphiz_mac_mf = gphi_mac[lev][2];
#endif
            MultiFab& umac_mid_mf = umac_mid[lev][0];
            MultiFab& vmac_mid_mf = umac_mid[lev][1];
#if (AMREX_SPACEDIM == 3)
            MultiFab& wmac_mid_mf = umac_mid[lev][2];
#endif
            const Real* dx = geom[lev].CellSize();
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(gphix_mac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                // NOTE: I have no idea what 'box_phys_bc' is in MAESTROeX, so am
                // just going to pass in bcs_u[0].data, as is used for velpred
                // and mkutrans

                add_grad_scalar_mac(ARLIM_3D(lo), ARLIM_3D(hi),
                                    BL_TO_FORTRAN_3D(gphix_mac_mf[mfi]),
                                    BL_TO_FORTRAN_3D(gphiy_mac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                    BL_TO_FORTRAN_3D(gphiz_mac_mf[mfi]),
#endif
                                    BL_TO_FORTRAN_3D(umac_mid_mf[mfi]),
                                    BL_TO_FORTRAN_3D(vmac_mid_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                    BL_TO_FORTRAN_3D(wmac_mid_mf[mfi]),
#endif
                                    phys_bc.dataPtr(), bcs_u[0].data(),
                                    ZFILL(dx));
            }
        }

        // average down and do boundaries
        if (finest_level == 0) {
            // fill periodic ghost cells
            for (int lev = 0; lev <= finest_level; ++lev) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    umac_mid[lev][d].FillBoundary(geom[lev].periodicity());
                }
            }
            // fill ghost cells behind physical boundaries
            FillUmacGhost(umac_mid);
        } else {
            // edge_restriction for velocities
            AverageDownFaces(umac_mid);
            // fill level n ghost cells using interpolation from level n-1 data
            FillPatchUedge(umac_mid);
        }

        // write umac_mid
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_mid[lev][0];
            MultiFab& vmac_mf = umac_mid[lev][1];
            MultiFab& wmac_mf = umac_mid[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
                                  BL_TO_FORTRAN_3D(umac_mf[mfi]),
                                  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                                  BL_TO_FORTRAN_3D(utemp_mf[mfi]));
            }
        }

        WritePlotFile(1, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // write gphi_mac
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = gphi_mac[lev][0];
            MultiFab& vmac_mf = gphi_mac[lev][1];
            MultiFab& wmac_mf = gphi_mac[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
                                  BL_TO_FORTRAN_3D(umac_mf[mfi]),
                                  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                                  BL_TO_FORTRAN_3D(utemp_mf[mfi]));
            }
        }

        WritePlotFile(2, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);

        // copy the velocity field over to the final state, unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                MultiFab::Copy(umac_new[lev][comp], umac_mid[lev][comp], 0, 0,
                               1, umac_new[lev][comp].nGrow());
        }
    }

    //--------------------------------------------------------------------------
    // project out the divergent portion of the velocity field
    //--------------------------------------------------------------------------

    Print() << "...projection" << std::endl;

    if (project_type == 1) {
        // hgprojection -- here pi is nodal and u is cell-centered

        for (int lev = 0; lev <= finest_level; ++lev) {
            // build the density used in the projection -- we are just doing
            // constant density, so set it to 1
            sold[lev].setVal(0.);
            snew[lev].setVal(0.);
            sold[lev].setVal(1., Rho, 1, 1);
            snew[lev].setVal(1., Rho, 1, 1);
            pi[lev].setVal(0.);
            gpi[lev].setVal(0.);
            rhcc_for_nodalproj[lev].setVal(0.);
        }

        // build the coefficient in the divergence.  We are doing
        // divergence-free (incompressible), so set beta0 = 1
        beta0_old.setVal(1.0);
        beta0_new.setVal(1.0);

        t_new = t_old + 1.;

        // NodalProj is going to operate on unew, where we've temporarily stored the initial data. Let's instead copy this initial data to umid.
        for (int lev = 0; lev <= finest_level; ++lev)
            std::swap(umid[lev], unew[lev]);

        // hgproject
        NodalProj(initial_projection_comp, rhcc_for_nodalproj);

        // swap pointers to restore initial data to uold and new data to unew
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(uold[lev], umid[lev]);
            std::swap(umid[lev], unew[lev]);
        }

        // write unew
        WritePlotFile(3, t_new, dt, dummy, dummy, dummy, dummy, unew, uold,
                      uold);

        // I think now can compare to uold and see if it's the same?
        {
            int lev = finest_level;
            Real norm = 0.;
            umid[lev].setVal(0.);

            MultiFab::Copy(umid[lev], uold[lev], 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Subtract(umid[lev], unew[lev], 0, 0, AMREX_SPACEDIM, 0);

            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp)
                norm += umid[lev].norm2(comp) / uold[lev].norm2(comp);

            Print() << "\nRelative error = " << norm << std::endl;

            norm = 0.;
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                norm = umid[lev].norm0(comp);

                Print() << "\nAbsolute error, dim " << comp << " = " << norm
                        << std::endl;
            }
        }

    } else {
        // mac projection -- here pi is cell-centered and u is MAC

        Vector<MultiFab> macpi(finest_level + 1);
        Vector<MultiFab> macrhs(finest_level + 1);

        beta0_old.setVal(1.0);
        beta0_new.setVal(1.0);

        for (int lev = 0; lev <= finest_level; ++lev) {
            // cell-centered MultiFabs
            macpi[lev].define(grids[lev], dmap[lev], 1, 1);
            macpi[lev].setVal(0.);
            macrhs[lev].define(grids[lev], dmap[lev], 1, 1);
            macrhs[lev].setVal(0.);

            sold[lev].setVal(0.);
            snew[lev].setVal(0.);
            sold[lev].setVal(1., Rho, 1, 1);
            snew[lev].setVal(1., Rho, 1, 1);
        }

        // macproject
        auto is_predictor = 0;
        MacProj(umac_new, macpi, macrhs, beta0_old, is_predictor);

        // I think now can compare to umac_old and see if it's the same?
        {
            int lev = finest_level;
            Real norm = 0.;
            for (auto comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                umac_mid[lev][comp].setVal(0.);
                MultiFab::Copy(umac_mid[lev][comp], umac_old[lev][comp], 0, 0,
                               1, 0);
                MultiFab::Subtract(umac_mid[lev][comp], umac_new[lev][comp], 0,
                                   0, 1, 0);

                norm +=
                    umac_mid[lev][comp].norm2() / umac_old[lev][comp].norm2();
            }
            Print() << "\nRelative error = " << norm << std::endl;
        }

        // write umac_new
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& umac_mf = umac_new[lev][0];
            MultiFab& vmac_mf = umac_new[lev][1];
            MultiFab& wmac_mf = umac_new[lev][2];
            MultiFab& utemp_mf = utemp[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(umac_mf, true); mfi.isValid(); ++mfi) {
                const Box& tilebox = mfi.tilebox();
                const int* lo = tilebox.loVect();
                const int* hi = tilebox.hiVect();

                convert_MAC_to_cc(ARLIM_3D(lo), ARLIM_3D(hi),
                                  BL_TO_FORTRAN_3D(umac_mf[mfi]),
                                  BL_TO_FORTRAN_3D(vmac_mf[mfi]),
#if (AMREX_SPACEDIM == 3)
                                  BL_TO_FORTRAN_3D(wmac_mf[mfi]),
#endif
                                  BL_TO_FORTRAN_3D(utemp_mf[mfi]));
            }
        }

        WritePlotFile(3, t_new, dt, dummy, dummy, dummy, dummy, utemp, uold,
                      uold);
    }
}
