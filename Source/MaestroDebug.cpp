
#include <Maestro.H>
#include <Maestro_F.H>

using namespace amrex;

void Maestro::PrintBase(const RealVector& base, const bool is_cell_centered) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintBase()", PrintBase);

    const int max_lev = base_geom.max_radial_level + 1;

    for (auto lev = 0; lev <= base_geom.finest_radial_level; ++lev) {
        for (auto i = 0; i <= base_geom.numdisjointchunks(lev, 0); ++i) {
            auto lo = base_geom.r_start_coord(lev, i);
            auto hi = is_cell_centered ? base_geom.r_end_coord(lev, i)
                                       : base_geom.r_end_coord(lev, i) + 1;
            for (auto r = lo; r <= hi; ++r) {
                Print() << "base lev, r " << lev << ", " << r << ", "
                        << base[lev + max_lev * r] << std::endl;
            }
        }
    }
}

void Maestro::PrintBase(const BaseState<Real>& base_s,
                        const bool is_cell_centered) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintBase()", PrintBase);

    const auto base = base_s.const_array();

    for (auto lev = 0; lev <= base_geom.finest_radial_level; ++lev) {
        for (auto i = 0; i <= base_geom.numdisjointchunks(lev, 0); ++i) {
            auto lo = base_geom.r_start_coord(lev, i);
            auto hi = is_cell_centered ? base_geom.r_end_coord(lev, i)
                                       : base_geom.r_end_coord(lev, i) + 1;
            for (auto r = lo; r <= hi; ++r) {
                std::cout << std::setprecision(16) << "base lev, r " << lev
                          << ", " << r << ", " << base(lev, r) << std::endl;
            }
        }
    }
}

// print out the contents of a Vector of MultiFabs
void Maestro::PrintMF(const Vector<MultiFab>& MF) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintMF()", PrintMF);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const MultiFab& MF_mf = MF[lev];

        const BoxArray& ba = MF_mf.boxArray();
        const DistributionMapping& dm = MF_mf.DistributionMap();
        const int myProc = ParallelDescriptor::MyProc();

        for (int i = 0; i < ba.size(); ++i) {
            if (dm[i] == myProc) {
                // we want all processors to write, not just the IOProcessor
                std::cout << "Grid #" << i << std::endl;
                std::cout << "Processor #" << myProc << std::endl;

                const Box& validBox = ba[i];

                auto lo = validBox.loVect3d();
                auto hi = validBox.hiVect3d();

                std::cout << "Level " << lev << std::endl;
                std::cout << "valid box ";
                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    std::cout << "(" << lo[n] << ", " << hi[n] << ")  ";
                }
                std::cout << std::endl;

                const Array4<const Real> MF_arr = MF_mf.array(i);

                for (auto comp = 0; comp < MF_mf.nComp(); ++comp) {
                    for (auto k = lo[2]; k <= hi[2]; ++k) {
                        for (auto j = lo[1]; j <= hi[1]; ++j) {
                            for (auto ii = lo[0]; ii <= hi[0]; ++ii) {
#if (AMREX_SPACEDIM == 2)
                                std::cout << "lev, i, j, comp" << lev << " "
                                          << ii << " " << j << " " << comp
                                          << " " << MF_arr(ii, j, k, comp)
                                          << std::endl;
#else
                                std::cout << "lev, i, j, k, comp" << lev << " "
                                          << ii << " " << j << " " << k << " "
                                          << comp << " "
                                          << MF_arr(ii, j, k, comp)
                                          << std::endl;
#endif
                            }
                        }
                    }
                }
            }
            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

void Maestro::PrintBA(const Vector<BoxArray>& ba) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintBA()", PrintBA);

    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int i = 0; i < ba[lev].size(); ++i) {
            Print() << "Grid #" << i << std::endl;

            const Box& validBox = ba[lev][i];
            auto lo = validBox.loVect3d();
            auto hi = validBox.hiVect3d();

            Print() << "Level " << lev << std::endl;
            Print() << "valid box ";
            for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                Print() << "(" << lo[n] << ", " << hi[n] << ")  ";
            }
        }
    }
}

void Maestro::PrintEdge(
    const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& EDGE, int dir) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PrintEdge()", PrintEdge);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // get references to the MultiFabs at level lev
        const MultiFab& EDGE_mf = EDGE[lev][dir];

        const BoxArray& ba = EDGE_mf.boxArray();
        const DistributionMapping& dm = EDGE_mf.DistributionMap();
        const int myProc = ParallelDescriptor::MyProc();

        for (int i = 0; i < ba.size(); ++i) {
            if (dm[i] == myProc) {
                // we want all processors to write, not just the IOProcessor
                std::cout << "Grid #" << i << std::endl;
                std::cout << "Processor #" << myProc << std::endl;

                // EDGE BASED
                const Box& validBox = ba[i];
                auto lo = validBox.loVect3d();
                auto hi = validBox.hiVect3d();

                std::cout << "Level " << lev << std::endl;
                std::cout << "valid box ";
                for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
                    std::cout << "(" << lo[n] << ", " << hi[n] << ")  ";
                }
                std::cout << std::endl;

                const Array4<const Real> EDGE_arr = EDGE_mf.array(i);

                for (auto comp = 0; comp < EDGE_mf.nComp(); ++comp) {
                    for (auto k = lo[2]; k <= hi[2]; ++k) {
                        for (auto j = lo[1]; j <= hi[1]; ++j) {
                            for (auto ii = lo[0]; ii <= hi[0]; ++ii) {
#if (AMREX_SPACEDIM == 2)
                                std::cout << "lev, i, j, comp" << lev << " "
                                          << ii << " " << j << " " << comp
                                          << " " << EDGE_arr(ii, j, k, comp)
                                          << std::endl;
#else
                                std::cout << "lev, i, j, k, comp" << lev << " "
                                          << ii << " " << j << " " << k << " "
                                          << comp << " "
                                          << EDGE_arr(ii, j, k, comp)
                                          << std::endl;
#endif
                            }
                        }
                    }
                }
            }
            // add this barrier so only one grid gets printed out at a time
            ParallelDescriptor::Barrier();
        }
    }
}

// utility to write out a multilevel multifab to a plotfile
void Maestro::WriteMF(const Vector<MultiFab>& mf, const std::string& name) {
    int nComp = mf[0].nComp();

    Vector<std::string> varnames;
    varnames.resize(nComp);
    for (int i = 0; i < nComp; ++i) {
        varnames[i] = "MultiFab_" + std::to_string(i);
    }

    // temporary MultiFab to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level + 1);

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] = new MultiFab((mf[i]).boxArray(),
                                       (mf[i]).DistributionMap(), nComp, 0);
    }

    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->ParallelCopy((mf[i]), 0, 0, nComp);
    }

    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf.push_back(plot_mf_data[i]);
    }

    Vector<int> step_array;
    step_array.resize(maxLevel() + 1, 0);

    WriteMultiLevelPlotfile(name, finest_level + 1, plot_mf, varnames, Geom(),
                            0., step_array, refRatio());
}
