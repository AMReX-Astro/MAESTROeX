#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

#include <extern_parameters.H>

#include <network.H>
#include <eos.H>
#include <conductivity.H>

#include <fundamental_constants.H>

#include <amrex_astro_util.H>

using namespace amrex;

inline int
get_vx_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "velx");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find velx component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

inline int
get_vy_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "vely");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find vely component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

inline int
get_dT_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "tpert");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find tpert component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

inline int
get_Hnuc_index(const std::vector<std::string>& var_names_pf) {

    auto idx = std::find(var_names_pf.cbegin(), var_names_pf.cend(), "Hnuc");
    if (idx == var_names_pf.cend()) {
        amrex::Error("Error: could not find Hnuc component");
    }
    return std::distance(var_names_pf.cbegin(), idx);
}

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile(diag_rp::plotfile);

    if (pltfile.empty()) {
        std::cout << "no plotfile specified" << std::endl;
        std::cout << "use: diag.plotfile=plt00000 (for example)" << std::endl;
        amrex::Error("no plotfile");
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    std::string outfile = pltfile + "/derived";
    std::cout << outfile << std::endl;

    PlotFileData pf(pltfile);

    const int ndims = pf.spaceDim();
    AMREX_ALWAYS_ASSERT(ndims <= AMREX_SPACEDIM);

    const int nlevs = pf.finestLevel() + 1;

    Vector<std::string> varnames;
    varnames = pf.varNames();

    // find variable indices
    // We want:
    // density, temperature, pressure, species
    // vertical velocity, temperature perturbation
    // we will assume here that the species are contiguous, so we will find
    // the index of the first species

    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = get_dens_index(var_names_pf);
    int temp_comp = get_temp_index(var_names_pf);
    int pres_comp = get_pres_index(var_names_pf);
    int spec_comp = get_spec_index(var_names_pf);
    int dT_comp = get_dT_index(var_names_pf);

    int vy_comp = get_vy_index(var_names_pf);
    //int vx_comp = get_vx_index(var_names_pf);

    //int Hnuc_comp = get_Hnuc_index(var_names_pf);

    // create the variable names we will derive and store in the output
    // file

    Vector<std::string> gvarnames;

    gvarnames.push_back("del");
    gvarnames.push_back("del_ad");
    gvarnames.push_back("del_ledoux");
    gvarnames.push_back("vy2");
    gvarnames.push_back("Fconv_dT");
    gvarnames.push_back("Fconv_T");
    gvarnames.push_back("Fconv_mlt");
    gvarnames.push_back("Fconv_mlt_v");
    gvarnames.push_back("Fconv_mlt_vabs");
    gvarnames.push_back("v2_mlt");
    gvarnames.push_back("Fkin");
    gvarnames.push_back("Frad");
    gvarnames.push_back("Fh1");
    gvarnames.push_back("Fc12");
    gvarnames.push_back("Hp");
    gvarnames.push_back("cp");
    //gvarnames.push_back("Da");
    //gvarnames.push_back("Ri");

    // interpret the boundary conditions

    BCRec bcr_default;
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    IntVect ng(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim < ndims) {
            bcr_default.setLo(idim, BCType::hoextrapcc);
            bcr_default.setHi(idim, BCType::hoextrapcc);
        } else {
            bcr_default.setLo(idim, BCType::int_dir);
            bcr_default.setHi(idim, BCType::int_dir);
            is_periodic[idim] = 1;
            ng[idim] = 0;
        }
    }

    // we need the variables constructed with ghost cells

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {

        // output MultiFab

        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), gvarnames.size(), 0);

        Vector<BCRec> bcr{bcr_default};
        auto is_per = is_periodic;

        Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                         pf.coordSys(), is_per);
        geom.push_back(vargeom);

        PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
            (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

        // fill the pressure and temperature mfs with ghost cells
        // we also need all of the species

        MultiFab temp_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab dens_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab pres_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab species_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), NumSpec, ng);
        //MultiFab vx_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab vy_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        MultiFab dT_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
        //MultiFab Hnuc_mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);

        if (ilev == 0) {

            // temperature
            {
                MultiFab smf = pf.get(ilev, var_names_pf[temp_comp]);
                FillPatchSingleLevel(temp_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // density
            {
                MultiFab smf = pf.get(ilev, var_names_pf[dens_comp]);
                FillPatchSingleLevel(dens_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // pressure
            {
                MultiFab smf = pf.get(ilev, var_names_pf[pres_comp]);
                FillPatchSingleLevel(pres_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // species
            {
                for (int n = 0; n < NumSpec; ++n) {
                    MultiFab smf = pf.get(ilev, var_names_pf[spec_comp+n]);
                    FillPatchSingleLevel(species_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                         0, n, 1, vargeom, physbcf, 0);
                }
            }

            // v
            //{
            //    MultiFab smf = pf.get(ilev, var_names_pf[vx_comp]);
            //    FillPatchSingleLevel(vx_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
            //                         0, 0, 1, vargeom, physbcf, 0);
            //}

            {
                MultiFab smf = pf.get(ilev, var_names_pf[vy_comp]);
                FillPatchSingleLevel(vy_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // dT
            {
                MultiFab smf = pf.get(ilev, var_names_pf[dT_comp]);
                FillPatchSingleLevel(dT_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            }

            // Hnuc
            //{
            //    MultiFab smf = pf.get(ilev, var_names_pf[Hnuc_comp]);
            //    FillPatchSingleLevel(Hnuc_mf, ng, Real(0.0), {&smf}, {Real(0.0)},
            //                         0, 0, 1, vargeom, physbcf, 0);
            //}



        } else {
            auto* mapper = (Interpolater*)(&cell_cons_interp);

            IntVect ratio(pf.refRatio(ilev-1));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }

            Geometry cgeom(pf.probDomain(ilev-1), RealBox(pf.probLo(),pf.probHi()),
                           pf.coordSys(), is_per);
            PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cphysbcf
                (cgeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

            // temperature
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[temp_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[temp_comp]);
                FillPatchTwoLevels(temp_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // density
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[dens_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[dens_comp]);
                FillPatchTwoLevels(dens_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // pressure
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[pres_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[pres_comp]);
                FillPatchTwoLevels(pres_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // species
            {
                for (int n = 0; n < NumSpec; ++n) {
                    MultiFab cmf = pf.get(ilev-1, var_names_pf[spec_comp+n]);
                    MultiFab fmf = pf.get(ilev  , var_names_pf[spec_comp+n]);
                    FillPatchTwoLevels(species_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                       {&fmf}, {Real(0.0)}, 0, n, 1, cgeom, vargeom,
                                       cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
                }
            }

            // v
            //{
            //    MultiFab cmf = pf.get(ilev-1, var_names_pf[vx_comp]);
            //    MultiFab fmf = pf.get(ilev  , var_names_pf[vx_comp]);
            //    FillPatchTwoLevels(vx_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
            //                       {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
            //                       cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            //}
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[vy_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[vy_comp]);
                FillPatchTwoLevels(vy_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }


            // dT
            {
                MultiFab cmf = pf.get(ilev-1, var_names_pf[dT_comp]);
                MultiFab fmf = pf.get(ilev  , var_names_pf[dT_comp]);
                FillPatchTwoLevels(dT_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            // Hnuc
            //{
            //    MultiFab cmf = pf.get(ilev-1, var_names_pf[Hnuc_comp]);
            //    MultiFab fmf = pf.get(ilev  , var_names_pf[Hnuc_comp]);
            //    FillPatchTwoLevels(Hnuc_mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
            //                       {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
            //                       cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            //}

        }

        auto const& dx = pf.cellSize(ilev);

        const MultiFab& lev_data_mf = pf.get(ilev);

        for (MFIter mfi(temp_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();

            // output storage
            auto const& ga = gmf[ilev].array(mfi);

            // temperature, pressure and composition with ghost cells
            auto const& T = temp_mf.const_array(mfi);
            auto const& Rho = dens_mf.const_array(mfi);
            auto const& P = pres_mf.const_array(mfi);
            auto const& X = species_mf.const_array(mfi);
            //auto const& Hnuc = Hnuc_mf.const_array(mfi);
            //auto const& Vx = vx_mf.const_array(mfi);

            // all of the data without ghost cells
            const auto& fab = lev_data_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {

                Real dT_dy = 0.0;
                //Real drho_dy = 0.0;
                //Real dvx_dy = 0.0;
                Real del = 0.0;

                // y is the vertical
                dT_dy = (T(i,j+1,k) - T(i,j-1,k)) / (2.0*dx[1]);
                //drho_dy = (Rho(i,j+1,k) - Rho(i,j-1,k)) / (2.0*dx[1]);
                //dvx_dy = (Vx(i,j+1,k) - Vx(i,j-1,k)) / (2.0*dx[1]);
                Real dp = P(i,j+1,k) - P(i,j-1,k);
                if (dp != 0.0) {
                    del = (T(i,j+1,k) - T(i,j-1,k)) / dp * (P(i,j,k) / T(i,j,k));
                    }

                Real pres = fab(i,j,k,pres_comp);
                Real rho  = fab(i,j,k,dens_comp);
                Real temp = fab(i,j,k,temp_comp);
                //Real velx   = fab(i,j,k,vx_comp);
                Real vely   = fab(i,j,k,vy_comp);
                Real delT   = fab(i,j,k,dT_comp);
                //Real Hnuc   = fab(i,j,k,Hnuc_comp);

                // Make EOS
                eos_t eos_state;
                eos_state.rho = rho;
                eos_state.T = temp;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                }
                eos(eos_input_rt, eos_state);
                conductivity(eos_state);

                // Derive from EOS
                Real cp = eos_state.cp;
                Real chi_T = eos_state.dpdT * eos_state.T/eos_state.p;
                Real chi_rho = eos_state.dpdr * eos_state.rho/eos_state.p;
                Real delta = chi_T/chi_rho; // dlnd/dlnT = T/d dd/dT = T/d (dP/dT)/(dP/dd) = T/d chi_T/chi_d
                Real del_ad = eos_state.p * chi_T / (eos_state.rho * eos_state.T * cp * chi_rho);



                // Stuff for del_ledoux
                // del_ledoux = del_ad + B, where B is the composition term
                // We calculate it like MESA, Paxton+ 2013 Equation 8
                // but we do a centered difference

                Real lnP_plus{0.0};  // pressure "above"
                Real lnP_minus{0.0};  // pressure "below"

                Real lnPalt_plus{0.0};  // pressure with "above" species
                Real lnPalt_minus{0.0};  // pressure with "below" species

                // Assuming non-spherical, ndims=2
                // y is the vertical

                lnP_plus = std::log(P(i,j+1,k));
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j+1,k,n);
                }
                eos(eos_input_rt, eos_state);
                lnPalt_plus = std::log(eos_state.p);

                lnP_minus = std::log(P(i,j-1,k));
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = X(i,j-1,k,n);
                }
                eos(eos_input_rt, eos_state);
                lnPalt_minus = std::log(eos_state.p);

                Real denom = lnP_plus - lnP_minus;
                Real B{0.0};
                if (denom != 0.0) {
                    B = -1 / chi_T * (lnPalt_plus - lnPalt_minus) / denom;
                }
                Real del_ledoux = del_ad + B;


                // Restore eos
                eos_state.rho = rho;
                eos_state.T = temp;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = fab(i,j,k,spec_comp+n);
                }
                eos(eos_input_rt, eos_state);
                conductivity(eos_state);

                // Convective grads
                // varnames del,del_ad,del_ledoux,del-del_ad,del-del_ledoux
                ga(i,j,k,0) = del;
                ga(i,j,k,1) = del_ad;
                ga(i,j,k,2) = del_ledoux;
                //ga(i,j,k,3) = del-del_ad;
                //ga(i,j,k,4) = del-del_ledoux;

                // Velocities : "vx2","vy2","v2","vy3","vyabs3"
                //ga(i,j,k,5) = pow(velx,2);
                // ga(i,j,k,6) = pow(vely,2);
                // ga(i,j,k,7) = pow(velx,2) + pow(vely,2);
                // ga(i,j,k,8) = pow(vely,3);
                // ga(i,j,k,9) = pow(std::abs(vely),2);

                ga(i,j,k,3) = pow(vely,2);
                //ga(i,j,k,6) = pow(vely,3);
                //ga(i,j,k,7) = pow(std::abs(vely),2);

                // Gravity/Scale height
                // auto grav = GetVarFromJobInfo(pltfile, "maestro.grav_const"); // really slow!!
                // Real g = std::abs(std::stod(grav));
                // std::cout << g << std::endl;
                Real g = 1.29e14;
                Real Hp = pres/(rho*g);


                // Fluxes

                // Convective heat flux : "Fconv_dT"
                ga(i,j,k,4) = rho * cp * vely * delT;

                // Alternate version that should be equivalent : "Fconv_T"
                // because <vdT> = <v(T-T0)> = <vT> - <v>T0 = <vT> (the vertical velocity averages to zero)
                ga(i,j,k,5) = rho * cp * vely * temp;

                // MLT heat flux in the efficient regime : "Fconv_mlt"
                ga(i,j,k,6) = 0.5 * rho * cp * temp * vely * (del-del_ad);

                // MLT flux as a function of velocity : "Fconv_mlt_v3" and "Fconv_mlt_vabs3"
                ga(i,j,k,7) = 4 * rho * cp * temp * pow(vely,3) / (delta * g * Hp);
                ga(i,j,k,8) = 4 * rho * cp * temp * pow(std::abs(vely),3) / (delta * g * Hp);

                // MLT velocity squared : "v2_mlt"
                ga(i,j,k,9) = g * delta * (del-del_ad) * Hp/8;

                // Kinetic flux : "Fkin"
                ga(i,j,k,10) = 0.5 * rho * pow(vely,3);

                // Radiative flux  : "Frad"
                // conductivity is k = 4*a*c*T^3/(kap*rho)
                // see Microphysics/conductivity/stellar/actual_conductivity.H
                ga(i,j,k,11) = -eos_state.conductivity * dT_dy;

                // Hydrogen flux : "Fh1"
                ga(i,j,k,12) = rho * vely * fab(i,j,k,spec_comp+0); // this is rho*v*X, not rho*v*dX

                // Hydrogen flux : "Fc12"
                ga(i,j,k,13) = rho * vely * fab(i,j,k,spec_comp+2);

                // Other (scale height, cp, Damkholer number, Richardson number)
                ga(i,j,k,14) = Hp;
                ga(i,j,k,15) = cp;
                //ga(i,j,k,20) = Hp * Hnuc / (vely*cp*temp);
                //ga(i,j,k,21) = g * drho_dy / (rho * pow(dvx_dy, 2));

            });
        }
    }

    Vector<int> level_steps;
    Vector<IntVect> ref_ratio;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        level_steps.push_back(pf.levelStep(ilev));
        if (ilev < pf.finestLevel()) {
            ref_ratio.push_back(IntVect(pf.refRatio(ilev)));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ref_ratio[ilev][idim] = 1;
            }
        }
    }

    WriteMultiLevelPlotfile(outfile, nlevs, GetVecOfConstPtrs(gmf), gvarnames,
                            geom, pf.time(), level_steps, ref_ratio);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv);

    // initialize the runtime parameters

    init_extern_parameters();

    // initialize C++ Microphysics

    eos_init(diag_rp::small_temp, diag_rp::small_dens);
    network_init();

    main_main();
    amrex::Finalize();
}
