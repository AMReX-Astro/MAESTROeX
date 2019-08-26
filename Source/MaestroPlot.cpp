#include <Maestro.H>
#include <Maestro_F.H>
#include <MaestroPlot.H>
#include <AMReX_buildInfo.H>
#include <iterator>     // std::istream_iterator

using namespace amrex;

// write a small plotfile to disk
void Maestro::WriteSmallPlotFile (const int step,
                                  const Real t_in,
                                  const Real dt_in,
                                  const RealVector& rho0_in,
                                  const RealVector& rhoh0_in,
                                  const RealVector& p0_in,
                                  const RealVector& gamma1bar_in,
                                  const Vector<MultiFab>& u_in,
                                  Vector<MultiFab>& s_in,
                                  const Vector<MultiFab>& S_cc_in)
{
	WritePlotFile (step, t_in, dt_in, rho0_in, rhoh0_in, p0_in,
	               gamma1bar_in, u_in, s_in, S_cc_in, true);
}

// write plotfile to disk
void
Maestro::WritePlotFile (const int step,
                        const Real t_in,
                        const Real dt_in,
                        const RealVector& rho0_in,
                        const RealVector& rhoh0_in,
                        const RealVector& p0_in,
                        const RealVector& gamma1bar_in,
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& s_in,
                        const Vector<MultiFab>& S_cc_in,
                        const bool is_small)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WritePlotFile()",WritePlotFile);

	// wallclock time
	const Real strt_total = ParallelDescriptor::second();

	std::string plotfilename;

	if (!is_small) {
		plotfilename = plot_base_name;
	} else {
		plotfilename = small_plot_base_name;
	}

	if (step == plotInitData) {
		if (plotfilename.back() == '_') {
			plotfilename += "InitData";
		} else {
			plotfilename += +"_InitData";
		}

	}
	else if (step == plotInitProj) {
		if (plotfilename.back() == '_') {
			plotfilename += "after_InitProj";
		} else {
			plotfilename += +"_after_InitProj";
		}
	}
	else if (step == plotDivuIter) {
		if (plotfilename.back() == '_') {
			plotfilename += "after_DivuIter";
		} else {
			plotfilename += +"_after_DivuIter";
		}
	}
	else {
		PlotFileName(step, &plotfilename);
	}

	// convert rho0 to multi-D MultiFab
	Vector<MultiFab> rho0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		rho0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(rho0_in,rho0_cart,0,0);

	// convert rhoh0 to multi-D MultiFab
	Vector<MultiFab> rhoh0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		rhoh0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(rhoh0_in,rhoh0_cart,0,0);

	// convert p0 to multi-D MultiFab
	Vector<MultiFab> p0_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		p0_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(p0_in,p0_cart,0,0);

	// convert gamma1bar to multi-D MultiFab
	Vector<MultiFab> gamma1bar_cart(finest_level+1);
	for (int lev=0; lev<=finest_level; ++lev) {
		gamma1bar_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	}
	Put1dArrayOnCart(gamma1bar_in,gamma1bar_cart,0,0);

	int nPlot = 0;
	const auto& varnames = PlotFileVarNames(&nPlot);

	const auto& mf = PlotFileMF(nPlot,t_in,dt_in,rho0_cart,rhoh0_cart,p0_cart,
	                            gamma1bar_cart,u_in,s_in,p0_in,gamma1bar_in,
	                            S_cc_in);

	// WriteMultiLevelPlotfile expects an array of step numbers
	Vector<int> step_array;
	step_array.resize(maxLevel()+1, step);

	if (!is_small) {
		WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
		                        Geom(), t_in, step_array, refRatio());
	} else {
		int nSmallPlot = 0;
		const auto& small_plot_varnames = SmallPlotFileVarNames(&nSmallPlot,
		                                                        varnames);

		const auto& small_mf = SmallPlotFileMF(nPlot, nSmallPlot, mf, varnames,
		                                       small_plot_varnames);

		WriteMultiLevelPlotfile(plotfilename, finest_level+1, small_mf,
		                        small_plot_varnames, Geom(), t_in, step_array,
		                        refRatio());

		for (int i = 0; i <= finest_level; ++i)
			delete small_mf[i];
	}

	WriteJobInfo(plotfilename);

	VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

	// write out the cell-centered base state
	if (ParallelDescriptor::IOProcessor()) {

	  for (int lev=0; lev<=max_radial_level; ++lev) {

	    std::ofstream BaseCCFile;
	    BaseCCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
	    std::string BaseCCFileName(plotfilename + "/BaseCC_");
	    std::string levStr = std::to_string(lev);
	    BaseCCFileName.append(levStr);
	    BaseCCFile.open(BaseCCFileName.c_str(), std::ofstream::out   |
			    std::ofstream::trunc |
			    std::ofstream::binary);
	    if( !BaseCCFile.good()) {
	      amrex::FileOpenFailed(BaseCCFileName);
	    }

	    BaseCCFile.precision(17);

	    BaseCCFile << "r_cc  rho0  rhoh0  p0  gamma1bar \n";

	    int nr = nr_fine / pow(2,(max_radial_level-lev));

	    for (int i=0; i<nr; ++i) {
	      BaseCCFile << r_cc_loc[lev+(max_radial_level+1)*i] << " "
			 << rho0_in[lev+(max_radial_level+1)*i] << " "
			 << rhoh0_in[lev+(max_radial_level+1)*i] << " "
			 << p0_in[lev+(max_radial_level+1)*i] << " "
			 << gamma1bar_in[lev+(max_radial_level+1)*i] << "\n";
	    }
	  }
	}

	// write out the face-centered base state
	if (ParallelDescriptor::IOProcessor()) {

	  for (int lev=0; lev<=max_radial_level; ++lev) {

	    std::ofstream BaseFCFile;
	    BaseFCFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
	    std::string BaseFCFileName(plotfilename + "/BaseFC_");
	    std::string levStr = std::to_string(lev);
	    BaseFCFileName.append(levStr);
	    BaseFCFile.open(BaseFCFileName.c_str(), std::ofstream::out   |
			    std::ofstream::trunc |
			    std::ofstream::binary);
	    if( !BaseFCFile.good()) {
	      amrex::FileOpenFailed(BaseFCFileName);
	    }

	    BaseFCFile.precision(17);

	    BaseFCFile << "r_edge  w0 \n";

	    int nr = nr_fine / pow(2,(max_radial_level-lev));

	    for (int i=0; i<nr+1; ++i) {
	      BaseFCFile << r_edge_loc[lev+(max_radial_level+1)*i] << " "
			 << w0[lev+(max_radial_level+1)*i] << "\n";
	    }
	  }
	}

	// wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;

	// print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total,ParallelDescriptor::IOProcessorNumber());
	if (maestro_verbose > 0) {
		Print() << "Time to write plotfile: " << end_total << '\n';
	}

	for (int i = 0; i <= finest_level; ++i) {
		delete mf[i];
	}

}


// get plotfile name
void
Maestro::PlotFileName (const int lev, std::string* plotfilename)
{
	*plotfilename = Concatenate(*plotfilename, lev, 7);
}

// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::PlotFileMF (const int nPlot,
                     const Real t_in,
                     const Real dt_in,
                     const Vector<MultiFab>& rho0_cart,
                     const Vector<MultiFab>& rhoh0_cart,
                     const Vector<MultiFab>& p0_cart,
                     const Vector<MultiFab>& gamma1bar_cart,
                     const Vector<MultiFab>& u_in,
                     Vector<MultiFab>& s_in,
                     const RealVector& p0_in,
                     const RealVector& gamma1bar_in,
                     const Vector<MultiFab>& S_cc_in)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileMF()",PlotFileMF);

	// MultiFab to hold plotfile data
	Vector<const MultiFab*> plot_mf;

	// temporary MultiFab to hold plotfile data
	Vector<MultiFab*> plot_mf_data(finest_level+1);

	// temporary MultiFab for calculations
	Vector<MultiFab> tempmf(finest_level+1);
	Vector<MultiFab> tempmf_scalar1(finest_level+1);
	Vector<MultiFab> tempmf_scalar2(finest_level+1);
	RealVector tempbar_plot ((max_radial_level+1)*nr_fine);
	tempbar_plot.shrink_to_fit();
	std::fill(tempbar_plot.begin(), tempbar_plot.end(), 0.);

	int dest_comp = 0;

	// build temporary MultiFab to hold plotfile data
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot,0);
		tempmf[i].define(grids[i],dmap[i],AMREX_SPACEDIM,0);

		tempmf_scalar1[i].define(grids[i],dmap[i],1,0);
		tempmf_scalar2[i].define(grids[i],dmap[i],1,0);
	}

	// velocity
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(u_in[i],0,dest_comp,AMREX_SPACEDIM);
	}
	dest_comp += AMREX_SPACEDIM;

	// magvel
	MakeMagvel(u_in, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	++dest_comp;

	// momentum = magvel * rho
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		MultiFab::Multiply(*plot_mf_data[i], s_in[i], Rho, dest_comp, 1, 0);
	}
	++dest_comp;

	// vorticity
	MakeVorticity(u_in, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	++dest_comp;

	// rho
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Rho,dest_comp,1);
	}
	++dest_comp;

	// rhoh
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
	}
	++dest_comp;

	// h
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
		MultiFab::Divide(*plot_mf_data[i], s_in[i], Rho, dest_comp, 1, 0);
	}
	++dest_comp;

	// rhoX
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],FirstSpec,dest_comp,NumSpec);
	}
	dest_comp += NumSpec;

	if (plot_spec) {
		// X
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(s_in[i],FirstSpec,dest_comp,NumSpec);
			for (int comp=0; comp<NumSpec; ++comp) {
				MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
			}
		}
		dest_comp += NumSpec;

		// abar
		MakeAbar(s_in, tempmf);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;

	}

	Vector<MultiFab> stemp             (finest_level+1);
	Vector<MultiFab> rho_Hext          (finest_level+1);
	Vector<MultiFab> rho_omegadot      (finest_level+1);
	Vector<MultiFab> rho_Hnuc          (finest_level+1);

	for (int lev=0; lev<=finest_level; ++lev) {
		stemp             [lev].define(grids[lev], dmap[lev],   Nscal, 0);
		rho_Hext          [lev].define(grids[lev], dmap[lev],       1, 0);
		rho_omegadot      [lev].define(grids[lev], dmap[lev], NumSpec, 0);
		rho_Hnuc          [lev].define(grids[lev], dmap[lev],       1, 0);
	}

	if (dt_in < small_dt) {
		React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, small_dt);
	} else {
		React(s_in, stemp, rho_Hext, rho_omegadot, rho_Hnuc, p0_in, dt_in*0.5);
	}

	if (plot_spec || plot_omegadot) {
		// omegadot
		if (plot_omegadot) {
			for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i]->copy(rho_omegadot[i],0,dest_comp,NumSpec);
				for (int comp=0; comp<NumSpec; ++comp) {
					MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
				}
			}
			dest_comp += NumSpec;
		}
	}

	if (plot_Hext) {
		// Hext
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(rho_Hext[i],0,dest_comp,1);
			MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp,1,0);
		}
		++dest_comp;
	}

	if (plot_Hnuc) {
		// Hnuc
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(rho_Hnuc[i],0,dest_comp,1);
			MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp,1,0);
		}
		++dest_comp;
	}

	if (plot_eta) {
		// eta_rho
		Put1dArrayOnCart(etarho_cc,tempmf,1,0,bcs_u,0,1);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;
	}

	// compute tfromp
	TfromRhoP(s_in,p0_in);
	// tfromp
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
	}
	++dest_comp;

	// compute tfromh
	TfromRhoH(s_in,p0_in);
	for (int i = 0; i <= finest_level; ++i) {
		// tfromh
		plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
	}
	++dest_comp;

	// deltap
	// compute & copy tfromp
	PfromRhoH(s_in,s_in,tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		// tfromh
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],p0_cart[i],0,dest_comp,1,0);
	}
	++dest_comp;

	// deltaT
	// compute & copy tfromp
	TfromRhoP(s_in,p0_in);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
	}
	// compute tfromh
	TfromRhoH(s_in,p0_in);
	// compute deltaT = (tfromp - tfromh) / tfromh
	for (int i = 0; i <= finest_level; ++i) {
		MultiFab::Subtract(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
		MultiFab::Divide(*plot_mf_data[i],s_in[i],Temp,dest_comp,1,0);
	}
	++dest_comp;

	// restore tfromp if necessary
	if (use_tfromp) {
		TfromRhoP(s_in,p0_in);
	}

	// pi
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
	}
	++dest_comp;

	// pioverp0
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
		MultiFab::Divide(*plot_mf_data[i], p0_cart[i], 0, dest_comp, 1, 0);
	}
	++dest_comp;

	// p0pluspi
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Pi,dest_comp,1);
		MultiFab::Add(*plot_mf_data[i], p0_cart[i], 0, dest_comp, 1, 0);
	}
	++dest_comp;

	if (plot_gpi) {
		// gpi
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(gpi[i],0,dest_comp,AMREX_SPACEDIM);
		}
		dest_comp += AMREX_SPACEDIM;
	}

	// rhopert
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],Rho,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],rho0_cart[i],0,dest_comp,1,0);
	}
	++dest_comp;

	// rhohpert
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(s_in[i],RhoH,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],rhoh0_cart[i],0,dest_comp,1,0);
	}
	++dest_comp;

	// tpert
	{
		Average(s_in, tempbar_plot, Temp);
		Put1dArrayOnCart(tempbar_plot,tempmf,0,0,bcs_f,0);

		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(s_in[i],Temp,dest_comp,1);
			MultiFab::Subtract(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
		}
	}
	++dest_comp;

    if (plot_base_state) {
        // rho0, rhoh0, h0 and p0
        for (int i = 0; i <= finest_level; ++i) {
        	plot_mf_data[i]->copy( rho0_cart[i],0,dest_comp,1);
        	plot_mf_data[i]->copy(rhoh0_cart[i],0,dest_comp+1,1);
        	plot_mf_data[i]->copy(rhoh0_cart[i],0,dest_comp+2,1);

        	// we have to use protected_divide here to guard against division by zero
        	// in the case that there are zeros rho0
        	MultiFab& plot_mf_data_mf = *plot_mf_data[i];
#ifdef _OPENMP
#pragma omp parallel
#endif
        	for ( MFIter mfi(plot_mf_data_mf, true); mfi.isValid(); ++mfi ) {
                plot_mf_data_mf[mfi].protected_divide(plot_mf_data_mf[mfi], dest_comp, dest_comp+2);
        	}

        	plot_mf_data[i]->copy(p0_cart[i],0,dest_comp+3,1);
        }
        dest_comp += 4;
    }

	Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
	Vector<MultiFab> w0r_cart(finest_level+1);

	if (spherical == 1) {

		for (int lev=0; lev<=finest_level; ++lev) {
			// w0mac will contain an edge-centered w0 on a Cartesian grid,
			// for use in computing divergences.
			AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1); ,
			             w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1); ,
			             w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1); );
			for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
				w0mac[lev][idim].setVal(0.);
			}

			// w0r_cart is w0 but onto a Cartesian grid in cell-centered as
			// a scalar.  Since w0 is the radial expansion velocity, w0r_cart
			// is the radial w0 in a zone
			w0r_cart[lev].define(grids[lev], dmap[lev], 1, 0);
			w0r_cart[lev].setVal(0.);
		}

		if (evolve_base_state == 1) {
			MakeW0mac(w0mac);
			Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0,1);
		}

		// Mach number
		MachfromRhoHSphr(s_in,u_in,p0_in,w0r_cart,tempmf);
	} else {
		// Mach number
		MachfromRhoH(s_in,u_in,p0_in,tempmf);
	}

	// MachNumber
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	++dest_comp;

	// deltagamma
	MakeDeltaGamma(s_in, p0_in, p0_cart, gamma1bar_in, gamma1bar_cart, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	++dest_comp;

	// entropy
	MakeEntropy(s_in, tempmf);
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	++dest_comp;

	// entropypert = (entropy - entropybar) / entropybar
	{
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}

		Average(tempmf, tempbar_plot, 0);
		Put1dArrayOnCart(tempbar_plot,tempmf,0,0,bcs_f,0);

		for (int i = 0; i <= finest_level; ++i) {
			MultiFab::Subtract(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
			MultiFab::Divide(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
		}
	}
	++dest_comp;

	if (plot_pidivu) {
		// pidivu
		MakePiDivu(u_in, s_in, tempmf);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;
	}

    // processor number of each tile
    if (plot_processors) {
        for (int i = 0; i <= finest_level; ++i) {
            (*plot_mf_data[i]).setVal(ParallelDescriptor::MyProc());
    	}
        ++dest_comp;
    }

	if (plot_ad_excess) {
		// ad_excess
		MakeAdExcess(s_in, tempmf);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;
	}

	// S
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(S_cc_in[i],0,dest_comp,1);
	}
	++dest_comp;

	// soundspeed
	if (plot_cs) {
		CsfromRhoH(s_in, p0_in, p0_cart, tempmf);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;
	}

  // gravitational_acceleration
	if (plot_grav) {
		MakeGrav(rho0_new, tempmf);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		}
		++dest_comp;
	}

        if (plot_base_state) {
            // w0
            Put1dArrayOnCart(w0,tempmf,1,1,bcs_u,0,1);
            for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,AMREX_SPACEDIM);
            }
            dest_comp += AMREX_SPACEDIM;

            // divw0
            MakeDivw0(w0mac, tempmf);
            for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
            }
            dest_comp++;
        }

	// thermal
	Vector<MultiFab> Tcoeff            (finest_level+1);
	Vector<MultiFab> hcoeff            (finest_level+1);
	Vector<MultiFab> Xkcoeff           (finest_level+1);
	Vector<MultiFab> pcoeff            (finest_level+1);

	for (int lev=0; lev<=finest_level; ++lev) {
		Tcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
		hcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
		Xkcoeff           [lev].define(grids[lev], dmap[lev], NumSpec, 1);
		pcoeff            [lev].define(grids[lev], dmap[lev],       1, 1);
	}

	if (use_thermal_diffusion) {
		MakeThermalCoeffs(s_in,Tcoeff,hcoeff,Xkcoeff,pcoeff);
		MakeExplicitThermal(tempmf,s_in,Tcoeff,hcoeff,Xkcoeff,pcoeff,p0_in,0);
	} else {
		for (int lev=0; lev<=finest_level; ++lev) {
			Tcoeff[lev].setVal(0.);
			tempmf[lev].setVal(0.);
		}
	}
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
	}
	dest_comp++;

	// conductivity
	for (int i = 0; i <= finest_level; ++i) {
		tempmf[i].setVal(0.);
		plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
		MultiFab::Subtract(*plot_mf_data[i],Tcoeff[i],0,dest_comp,1,0);
	}
	dest_comp++;

	// radial and circular velocities
	if (spherical == 1) {
		MakeVelrc(u_in, w0r_cart, tempmf, tempmf_scalar1);
		for (int i = 0; i <= finest_level; ++i) {
			plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
			plot_mf_data[i]->copy(tempmf_scalar1[i],0,dest_comp+1,1);
		}
		dest_comp += 2;
	}

	if (do_sponge) {
		init_sponge(rho0_old.dataPtr());
		MakeSponge(tempmf);

		if (plot_sponge_fdamp) {
			// compute f_damp assuming sponge=1/(1+dt*kappa*fdamp)
			// therefore fdamp = (1/sponge-1)/(dt*kappa)
			for (int i = 0; i <= finest_level; ++i) {
				// scalar1 = 1
				tempmf_scalar1[i].setVal(1.);
				// scalar2 = dt * kappa
				tempmf_scalar2[i].setVal(dt * sponge_kappa);
				// plot_mf = 1
				plot_mf_data[i]->copy(tempmf_scalar1[i],0,dest_comp,1);
				// plot_mf = 1/sponge
				MultiFab::Divide(*plot_mf_data[i],tempmf[i],0,dest_comp,1,0);
				// plot_mf = 1/sponge - 1
				MultiFab::Subtract(*plot_mf_data[i],tempmf_scalar1[i],0,dest_comp,1,0);
				// plot_mf = (1/sponge-1)/(dt*kappa)
				MultiFab::Divide(*plot_mf_data[i],tempmf_scalar2[i],0,dest_comp,1,0);
			}
		} else {
			for (int i = 0; i <= finest_level; ++i) {
				plot_mf_data[i]->copy(tempmf[i],0,dest_comp,1);
			}
		}
		dest_comp++;
	}

	// add plot_mf_data[i] to plot_mf
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf.push_back(plot_mf_data[i]);
		// delete [] plot_mf_data[i];
	}

	return plot_mf;

}


// this takes the multifab of all variables and extracts those
// required for the small plot file
Vector<const MultiFab*>
Maestro::SmallPlotFileMF(const int nPlot, const int nSmallPlot,
                         Vector<const MultiFab*> mf,
                         const Vector<std::string> varnames,
                         const Vector<std::string> small_plot_varnames)
{

	// timer for profiling
	BL_PROFILE_VAR("Maestro::SmallPlotFileMF()",SmallPlotFileMF);

	// MultiFab to hold plotfile data
	Vector<const MultiFab*> plot_mf;

	// temporary MultiFabs to hold plotfile data
	Vector<MultiFab*> plot_mf_data(finest_level+1);

	int dest_comp = 0;

	// build temporary MultiFab to hold plotfile data
	for (int i = 0; i <= finest_level; ++i) {
		plot_mf_data[i] = new MultiFab(mf[i]->boxArray(),
		                               mf[i]->DistributionMap(),nSmallPlot,0);

	}

	for (auto it=small_plot_varnames.begin();
	     it!=small_plot_varnames.end(); ++it) {

		for (auto n = 0; n < nPlot; n++) {
			if (*it == varnames[n]) {
				for (int i = 0; i <= finest_level; ++i)
					plot_mf_data[i]->copy(*(mf[i]), n, dest_comp, 1);
				++dest_comp;
				break;
			}
		}
	}

	// add plot_mf_data[i] to plot_mf
	for (int i = 0; i <= finest_level; ++i)
		plot_mf.push_back(plot_mf_data[i]);

	return plot_mf;

}

// set plotfile variable names
Vector<std::string>
Maestro::PlotFileVarNames (int * nPlot) const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

	// velocities (AMREX_SPACEDIM)
	// magvel, momentum
	// rho, rhoh, h, rhoX, tfromp, tfromh, deltap, deltaT Pi (Nscal+4 -- the extra 4 are h, tfromh, deltap and deltaT)
	// rho' and rhoh' and t' (3)
	// pioverp0, p0pluspi (2)
	// MachNumber, deltagamma, entropy, entropypert, S
	// thermal, conductivity

	(*nPlot) = AMREX_SPACEDIM + Nscal + 19;

	if (plot_spec) (*nPlot) += NumSpec + 1; // X + 1 (abar)
	if (plot_spec || plot_omegadot) (*nPlot) += NumSpec; // omegadot

	if (plot_Hext) (*nPlot)++;
	if (plot_Hnuc) (*nPlot)++;
	if (plot_eta) (*nPlot)++;
	if (plot_gpi) (*nPlot) += AMREX_SPACEDIM;
	// rho0, rhoh0, h0, p0, w0, divw0 (5+AMREX_SPACEDIM)
        if (plot_base_state) (*nPlot) += AMREX_SPACEDIM + 5;
	if (plot_cs) (*nPlot)++;
  if (plot_grav) (*nPlot)++;
	if (plot_ad_excess) (*nPlot)++;
	if (plot_pidivu) (*nPlot)++;
    if (plot_processors) (*nPlot)++;
	if (spherical == 1) (*nPlot) += 2; // radial_velocity, circ_velocity
	if (do_sponge) (*nPlot)++;

	Vector<std::string> names(*nPlot);

	int cnt = 0;

	// add velocities
	for (int i=0; i<AMREX_SPACEDIM; ++i) {
		std::string x = "vel";
		x += (120+i);
		names[cnt++] = x;
	}

	names[cnt++] = "magvel";
	names[cnt++] = "momentum";

	names[cnt++] = "vort";

	// density and enthalpy
	names[cnt++] = "rho";
	names[cnt++] = "rhoh";
	names[cnt++] = "h";

	for (int i = 0; i < NumSpec; i++) {
		int len = 20;
		Vector<int> int_spec_names(len);
		//
		// This call return the actual length of each string in "len"
		//
		get_spec_names(int_spec_names.dataPtr(),&i,&len);
		char* spec_name = new char[len+1];
		for (int j = 0; j < len; j++)
			spec_name[j] = int_spec_names[j];
		spec_name[len] = '\0';
		std::string spec_string = "rhoX(";
		spec_string += spec_name;
		spec_string += ')';

		names[cnt++] = spec_string;
		delete [] spec_name;
	}

	if (plot_spec) {
		for (int i = 0; i < NumSpec; i++) {
			int len = 20;
			Vector<int> int_spec_names(len);
			//
			// This call return the actual length of each string in "len"
			//
			get_spec_names(int_spec_names.dataPtr(),&i,&len);
			char* spec_name = new char[len+1];
			for (int j = 0; j < len; j++) {
				spec_name[j] = int_spec_names[j];
			}
			spec_name[len] = '\0';
			std::string spec_string = "X(";
			spec_string += spec_name;
			spec_string += ')';

			names[cnt++] = spec_string;

			delete [] spec_name;
		}

		names[cnt++] = "abar";
	}

	if (plot_spec || plot_omegadot) {
		for (int i = 0; i < NumSpec; i++) {
			int len = 20;
			Vector<int> int_spec_names(len);
			//
			// This call return the actual length of each string in "len"
			//
			get_spec_names(int_spec_names.dataPtr(),&i,&len);
			char* spec_name = new char[len+1];
			for (int j = 0; j < len; j++) {
				spec_name[j] = int_spec_names[j];
			}
			spec_name[len] = '\0';
			std::string spec_string = "omegadot(";
			spec_string += spec_name;
			spec_string += ')';

			names[cnt++] = spec_string;

			delete [] spec_name;
		}
	}

	if (plot_Hext) names[cnt++] = "Hext";
	if (plot_Hnuc) names[cnt++] = "Hnuc";
	if (plot_eta) names[cnt++] = "eta_rho";

	names[cnt++] = "tfromp";
	names[cnt++] = "tfromh";
	names[cnt++] = "deltap";
	names[cnt++] = "deltaT";
	names[cnt++] = "Pi";
	names[cnt++] = "pioverp0";
	names[cnt++] = "p0pluspi";

	if (plot_gpi) {
		// add gpi
		for (int i=0; i<AMREX_SPACEDIM; ++i) {
			std::string x = "gpi";
			x += (120+i);
			names[cnt++] = x;
		}
	}

	names[cnt++] = "rhopert";
	names[cnt++] = "rhohpert";
	names[cnt++] = "tpert";
        if (plot_base_state) {
            names[cnt++] = "rho0";
            names[cnt++] = "rhoh0";
            names[cnt++] = "h0";
            names[cnt++] = "p0";
        }
	names[cnt++] = "MachNumber";
	names[cnt++] = "deltagamma";
	names[cnt++] = "entropy";
	names[cnt++] = "entropypert";
	if (plot_pidivu) names[cnt++] = "pi_divu";
    if (plot_processors) names[cnt++] = "processor_number";
	if (plot_ad_excess) names[cnt++] = "ad_excess";
	names[cnt++] = "S";

	if (plot_cs) names[cnt++] = "soundspeed";

  if (plot_grav) names[cnt++] = "maggrav";

        if (plot_base_state) {
            // w0 and divw0
            for (int i=0; i<AMREX_SPACEDIM; ++i) {
		std::string x = "w0";
		x += (120+i);
		names[cnt++] = x;
            }
            names[cnt++] = "divw0";
        }

	names[cnt++] = "thermal";
	names[cnt++] = "conductivity";

	if (spherical == 1) {
		names[cnt++] = "radial_velocity";
		names[cnt++] = "circ_velocity";
	}

	if (do_sponge) {
		if (plot_sponge_fdamp) {
			names[cnt++] = "sponge_fdamp";
		} else {
			names[cnt++] = "sponge";
		}
	}

	return names;

}

// set plotfile variable names
Vector<std::string>
Maestro::SmallPlotFileVarNames (int * nPlot, Vector<std::string> varnames) const
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::SmallPlotFileVarNames()",SmallPlotFileVarNames);

    Vector<std::string> names(*nPlot);

    ParmParse pp("maestro");

    int nPltVars = pp.countval("small_plot_vars");

    if (nPltVars > 0) { // small_plot_vars defined in inputs file

        std::string nm;

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("small_plot_vars", nm, i);

            if (nm == "ALL")
                return varnames;
            else if (nm == "NONE") {
                names.clear();
                return names;
            } else {
                // test to see if it's a valid varname by iterating over
                // varnames
                auto found_name = false;
                for (auto it=varnames.begin(); it!=varnames.end(); ++it) {
                    if (nm == *it) {
                        names.push_back(nm);
                        found_name = true;
                        break;
                    }
                }

                if (!found_name)
                    Print() << "Small plot file variable " << nm << " is invalid\n";
            }
        }
    } else {
        // use default value of small_plot_vars which is a string that needs to be split
        std::stringstream sstream(small_plot_vars);
        std::string nm;

        while (sstream >> nm) {
            // test to see if it's a valid varname by iterating over
            // varnames
            auto found_name = false;
            for (auto it=varnames.begin(); it!=varnames.end(); ++it) {
                if (nm == *it) {
                    names.push_back(nm);
                    found_name = true;
                    break;
                }
            }

			if (!found_name)
				Print() << "Small plot file variable " << nm << " is invalid\n";
		}
	}

    names.shrink_to_fit();
    *nPlot = names.size();

    return names;

}

void
Maestro::WriteJobInfo (const std::string& dir) const
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::WriteJobInfo()",WriteJobInfo);

	if (ParallelDescriptor::IOProcessor())
	{
		// job_info file with details about the run
		std::ofstream jobInfoFile;
		std::string FullPathJobInfoFile = dir;

		std::string PrettyLine = std::string(78, '=') + "\n";
		std::string OtherLine = std::string(78, '-') + "\n";
		std::string SkipSpace = std::string(8, ' ');

		FullPathJobInfoFile += "/job_info";
		jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

		// job information
		jobInfoFile << PrettyLine;
		jobInfoFile << " MAESTROeX Job Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "job name: " << job_name << "\n\n";
		jobInfoFile << "inputs file: " << inputs_name << "\n\n";

		jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
		jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";

		jobInfoFile << "tile size: ";
		for (int d=0; d<AMREX_SPACEDIM; ++d) {
			jobInfoFile << FabArrayBase::mfiter_tile_size[d] << " ";
		}
		jobInfoFile << "\n";
#endif
		jobInfoFile << "\n";
		jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
		        getCPUTime()/3600.0;

		jobInfoFile << "\n\n";

		// plotfile information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Plotfile Information\n";
		jobInfoFile << PrettyLine;

		time_t now = time(0);

		// Convert now to tm struct for local timezone
		tm* localtm = localtime(&now);
		jobInfoFile << "output data / time: " << asctime(localtm);

		char currentDir[FILENAME_MAX];
		if (getcwd(currentDir, FILENAME_MAX)) {
			jobInfoFile << "output dir:         " << currentDir << "\n";
		}

		jobInfoFile << "\n\n";

		// build information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Build Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
		jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
		jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
		jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
		jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
		jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
		jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

		jobInfoFile << "\n";

		jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
		jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

		jobInfoFile << "\n";

		for (int n = 1; n <= buildInfoGetNumModules(); n++) {
			jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
		}

		const char* githash1 = buildInfoGetGitHash(1);
		const char* githash2 = buildInfoGetGitHash(2);
		const char* githash3 = buildInfoGetGitHash(3);
		if (strlen(githash1) > 0) {
			jobInfoFile << "MAESTROeX git describe: " << githash1 << "\n";
		}
		if (strlen(githash2) > 0) {
			jobInfoFile << "AMReX git describe: " << githash2 << "\n";
		}
		if (strlen(githash3) > 0) {
			jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
		}

		const char* buildgithash = buildInfoGetBuildGitHash();
		const char* buildgitname = buildInfoGetBuildGitName();
		if (strlen(buildgithash) > 0) {
			jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
		}

		jobInfoFile << "\n\n";

		// grid information
		jobInfoFile << PrettyLine;
		jobInfoFile << " Grid Information\n";
		jobInfoFile << PrettyLine;

		for (int i = 0; i <= finest_level; i++)
		{
			jobInfoFile << " level: " << i << "\n";
			jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
			jobInfoFile << "   maximum zones   = ";
			for (int n = 0; n < BL_SPACEDIM; n++)
			{
				jobInfoFile << geom[i].Domain().length(n) << " ";
			}
			jobInfoFile << "\n\n";
		}

		jobInfoFile << " Boundary conditions\n";
		Vector<int> lo_bc_out(BL_SPACEDIM), hi_bc_out(BL_SPACEDIM);
		ParmParse pp("maestro");
		pp.getarr("lo_bc",lo_bc_out,0,BL_SPACEDIM);
		pp.getarr("hi_bc",hi_bc_out,0,BL_SPACEDIM);


// these names correspond to the integer flags setup in the
// Castro_setup.cpp
		const char* names_bc[] =
		{ "interior", "inflow", "outflow",
		  "symmetry", "slipwall", "noslipwall" };


		jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
		jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
		if (BL_SPACEDIM >= 2) {
			jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
			jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
		}
		if (BL_SPACEDIM == 3) {
			jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
			jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
		}

		jobInfoFile << "\n\n";


// species info
		Real Aion = 0.0;
		Real Zion = 0.0;

		int mlen = 20;

		jobInfoFile << PrettyLine;
		jobInfoFile << " Species Information\n";
		jobInfoFile << PrettyLine;

		jobInfoFile <<
		        std::setw(6) << "index" << SkipSpace <<
		        std::setw(mlen+1) << "name" << SkipSpace <<
		        std::setw(7) << "A" << SkipSpace <<
		        std::setw(7) << "Z" << "\n";
		jobInfoFile << OtherLine;

		for (int i = 0; i < NumSpec; i++)
		{

			int len = mlen;
			Vector<int> int_spec_names(len);
			//
			// This call return the actual length of each string in "len"
			//
			get_spec_names(int_spec_names.dataPtr(),&i,&len);
			char* spec_name = new char[len+1];
			for (int j = 0; j < len; j++)
				spec_name[j] = int_spec_names[j];
			spec_name[len] = '\0';

			// get A and Z
			get_spec_az(&i, &Aion, &Zion);

			jobInfoFile <<
			        std::setw(6) << i << SkipSpace <<
			        std::setw(mlen+1) << std::setfill(' ') << spec_name << SkipSpace <<
			        std::setw(7) << Aion << SkipSpace <<
			        std::setw(7) << Zion << "\n";
			delete [] spec_name;
		}
		jobInfoFile << "\n\n";

		// runtime parameters
		jobInfoFile << PrettyLine;
		jobInfoFile << " Inputs File Parameters\n";
		jobInfoFile << PrettyLine;

		ParmParse::dumpTable(jobInfoFile, true);

		jobInfoFile.close();

		// now the external parameters
		const int jobinfo_file_length = FullPathJobInfoFile.length();
		Vector<int> jobinfo_file_name(jobinfo_file_length);

		for (int i = 0; i < jobinfo_file_length; i++)
			jobinfo_file_name[i] = FullPathJobInfoFile[i];

		runtime_pretty_print(jobinfo_file_name.dataPtr(), &jobinfo_file_length);
	}
}

void
Maestro::MakeMagvel (const Vector<MultiFab>& vel,
                     Vector<MultiFab>& magvel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeMagvel()",MakeMagvel);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);

#if (AMREX_SPACEDIM == 3)
	if (spherical == 1) {
		for (int lev=0; lev<=finest_level; ++lev) {
			w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);
			w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);
			w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1);
		}
		MakeW0mac(w0mac);
	}
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& magvel_mf = magvel[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_magvel(AMREX_INT_ANYD(tileBox.loVect()),
                            AMREX_INT_ANYD(tileBox.hiVect()),
                            lev,
				            BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
				            w0.dataPtr(),
				            BL_TO_FORTRAN_ANYD(magvel_mf[mfi]));
			}

		} else {

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				MultiFab& w0macx_mf = w0mac[lev][0];
				MultiFab& w0macy_mf = w0mac[lev][1];
				MultiFab& w0macz_mf = w0mac[lev][2];

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_magvel_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                 AMREX_INT_ANYD(tileBox.hiVect()),
				                 BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
				                 BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
				                 BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
				                 BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
				                 BL_TO_FORTRAN_ANYD(magvel_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(magvel,0,1);
	FillPatch(t_old,magvel,magvel,magvel,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}


void
Maestro::MakeVelrc (const Vector<MultiFab>& vel,
                    const Vector<MultiFab>& w0rcart,
                    Vector<MultiFab>& rad_vel,
                    Vector<MultiFab>& circ_vel)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVelrc()",MakeVelrc);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& radvel_mf = rad_vel[lev];
		MultiFab& circvel_mf = circ_vel[lev];
		const MultiFab& w0rcart_mf = w0rcart[lev];
		const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
			make_velrc(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
			           BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
			           BL_TO_FORTRAN_ANYD(w0rcart_mf[mfi]),
			           BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
			           BL_TO_FORTRAN_ANYD(radvel_mf[mfi]),
			           BL_TO_FORTRAN_ANYD(circvel_mf[mfi]));
		}
	}

	// average down and fill ghost cells
	AverageDown(rad_vel,0,1);
	FillPatch(t_old,rad_vel,rad_vel,rad_vel,0,0,1,0,bcs_f);
	AverageDown(circ_vel,0,1);
	FillPatch(t_old,circ_vel,circ_vel,circ_vel,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}


void
Maestro::MakeAdExcess (const Vector<MultiFab>& state,
                       Vector<MultiFab>& ad_excess)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeAdExcess()",MakeAdExcess);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		MultiFab& ad_excess_mf = ad_excess[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
		if (spherical == 0) {
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_ad_excess(AMREX_INT_ANYD(tileBox.loVect()),
                               AMREX_INT_ANYD(tileBox.hiVect()),
				               BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
				               BL_TO_FORTRAN_ANYD(ad_excess_mf[mfi]));
			}

		} else {

			const MultiFab& normal_mf = normal[lev];

#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_ad_excess_sphr(AMREX_INT_ANYD(tileBox.loVect()),
				                    AMREX_INT_ANYD(tileBox.hiVect()),
				                    BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
				                    BL_TO_FORTRAN_ANYD(normal_mf[mfi]),
				                    BL_TO_FORTRAN_ANYD(ad_excess_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(ad_excess,0,1);
	FillPatch(t_old,ad_excess,ad_excess,ad_excess,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}


void
Maestro::MakeGrav (const RealVector& rho0,
                   Vector<MultiFab>& grav)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeGrav()",MakeGrav);

  RealVector grav_cell( (max_radial_level+1)*nr_fine );
  grav_cell.shrink_to_fit();

  make_grav_cell(grav_cell.dataPtr(),
                 rho0.dataPtr(),
                 r_cc_loc.dataPtr(),
                 r_edge_loc.dataPtr());

  Put1dArrayOnCart(grav_cell,grav,0,0,bcs_f,0);

	// average down and fill ghost cells
	AverageDown(grav,0,1);
	FillPatch(t_old,grav,grav,grav,0,0,1,0,bcs_f);
}


void
Maestro::MakeVorticity (const Vector<MultiFab>& vel,
                        Vector<MultiFab>& vorticity)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeVorticity()",MakeVorticity);

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		MultiFab& vorticity_mf = vorticity[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(vel_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.

            // NOTE: do not offload to the gpu as make_vorticity_3d contains nested
            // functions
			make_vorticity(ARLIM_3D(tileBox.loVect()),
                           ARLIM_3D(tileBox.hiVect()),
			               BL_TO_FORTRAN_3D(vel_mf[mfi]), ZFILL(dx),
			               BL_TO_FORTRAN_3D(vorticity_mf[mfi]), phys_bc.dataPtr());
		}


	}

	// average down and fill ghost cells
	AverageDown(vorticity,0,1);
	FillPatch(t_old,vorticity,vorticity,vorticity,0,0,1,0,bcs_f);
}

void
Maestro::MakeDeltaGamma (const Vector<MultiFab>& state,
                         const RealVector& p0,
                         const Vector<MultiFab>& p0_cart,
                         const RealVector& gamma1bar,
                         const Vector<MultiFab>& gamma1bar_cart,
                         Vector<MultiFab>& deltagamma)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeDeltaGamma()",MakeDeltaGamma);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		MultiFab& deltagamma_mf = deltagamma[lev];

		if (spherical == 0) {

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_deltagamma(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
                                lev,
				                BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
				                p0.dataPtr(), gamma1bar.dataPtr(),
				                BL_TO_FORTRAN_ANYD(deltagamma_mf[mfi]));
			}

		} else {

			const MultiFab& p0cart_mf = p0_cart[lev];
			const MultiFab& gamma1barcart_mf = gamma1bar_cart[lev];

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_deltagamma_sphr(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
				                     BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
				                     BL_TO_FORTRAN_ANYD(p0cart_mf[mfi]), BL_TO_FORTRAN_ANYD(gamma1barcart_mf[mfi]),
				                     BL_TO_FORTRAN_ANYD(deltagamma_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(deltagamma,0,1);
	FillPatch(t_old,deltagamma,deltagamma,deltagamma,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MakeEntropy (const Vector<MultiFab>& state,
                      Vector<MultiFab>& entropy)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeEntropy()",MakeEntropy);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		MultiFab& entropy_mf = entropy[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(state_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
			make_entropy(AMREX_INT_ANYD(tileBox.loVect()),
			             AMREX_INT_ANYD(tileBox.hiVect()),lev,
			             BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
			             BL_TO_FORTRAN_ANYD(entropy_mf[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(entropy,0,1);
	FillPatch(t_old,entropy,entropy,entropy,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MakeDivw0 (const Vector<std::array<MultiFab, AMREX_SPACEDIM> >& w0mac,
                    Vector<MultiFab>& divw0)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakeDivw0()",MakeDivw0);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		MultiFab& divw0_mf = divw0[lev];

		if (spherical == 0) {

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(divw0_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();
				const Real* dx = geom[lev].CellSize();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_divw0(AMREX_INT_ANYD(tileBox.loVect()),
                           AMREX_INT_ANYD(tileBox.hiVect()),lev,
				           w0.dataPtr(), AMREX_REAL_ANYD(dx),
				           BL_TO_FORTRAN_ANYD(divw0_mf[mfi]));
			}

		} else {

			const MultiFab& w0macx_mf = w0mac[lev][0];
			const MultiFab& w0macy_mf = w0mac[lev][1];
			const MultiFab& w0macz_mf = w0mac[lev][2];

			// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
			for ( MFIter mfi(divw0_mf, true); mfi.isValid(); ++mfi ) {

				// Get the index space of the valid region
				const Box& tileBox = mfi.tilebox();
				const Real* dx = geom[lev].CellSize();

				// call fortran subroutine
				// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
				// lo/hi coordinates (including ghost cells), and/or the # of components
				// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
				make_divw0_sphr(AMREX_INT_ANYD(tileBox.loVect()),
                                AMREX_INT_ANYD(tileBox.hiVect()),
				                BL_TO_FORTRAN_ANYD(w0macx_mf[mfi]),
				                BL_TO_FORTRAN_ANYD(w0macy_mf[mfi]),
				                BL_TO_FORTRAN_ANYD(w0macz_mf[mfi]),
                                AMREX_REAL_ANYD(dx),
				                BL_TO_FORTRAN_ANYD(divw0_mf[mfi]));
			}
		}
	}

	// average down and fill ghost cells
	AverageDown(divw0,0,1);
	FillPatch(t_old,divw0,divw0,divw0,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MakePiDivu (const Vector<MultiFab>& vel,
                     const Vector<MultiFab>& state,
                     Vector<MultiFab>& pidivu)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakePiDivu()",MakePiDivu);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& vel_mf = vel[lev];
		const MultiFab& state_mf = state[lev];
		MultiFab& pidivu_mf = pidivu[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(pidivu_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
			make_pidivu(AMREX_INT_ANYD(tileBox.loVect()), AMREX_INT_ANYD(tileBox.hiVect()),
			            BL_TO_FORTRAN_ANYD(vel_mf[mfi]),
                        AMREX_REAL_ANYD(dx),
			            BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
			            BL_TO_FORTRAN_ANYD(pidivu_mf[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(pidivu,0,1);
	FillPatch(t_old,pidivu,pidivu,pidivu,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}

void
Maestro::MakeAbar (const Vector<MultiFab>& state,
                   Vector<MultiFab>& abar)
{
	// timer for profiling
	BL_PROFILE_VAR("Maestro::MakePiDivu()",MakeAbar);

#ifdef AMREX_USE_CUDA
    auto not_launched = Gpu::notInLaunchRegion();
    // turn on GPU
    if (not_launched) Gpu::setLaunchRegion(true);
#endif

	for (int lev=0; lev<=finest_level; ++lev) {

		// get references to the MultiFabs at level lev
		const MultiFab& state_mf = state[lev];
		MultiFab& abar_mf = abar[lev];

		// Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
		for ( MFIter mfi(abar_mf, true); mfi.isValid(); ++mfi ) {

			// Get the index space of the valid region
			const Box& tileBox = mfi.tilebox();
			const Real* dx = geom[lev].CellSize();

			// call fortran subroutine
			// use macros in AMReX_ArrayLim.H to pass in each FAB's data,
			// lo/hi coordinates (including ghost cells), and/or the # of components
			// We will also pass "validBox", which specifies the "valid" region.
#pragma gpu box(tileBox)
			make_abar(AMREX_INT_ANYD(tileBox.loVect()),
                      AMREX_INT_ANYD(tileBox.hiVect()),
			          BL_TO_FORTRAN_ANYD(state_mf[mfi]), state_mf.nComp(),
			          BL_TO_FORTRAN_ANYD(abar_mf[mfi]));
		}

	}

	// average down and fill ghost cells
	AverageDown(abar,0,1);
	FillPatch(t_old,abar,abar,abar,0,0,1,0,bcs_f);

#ifdef AMREX_USE_CUDA
    // turn off GPU
    if (not_launched) Gpu::setLaunchRegion(false);
#endif
}
