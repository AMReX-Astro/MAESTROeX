
#include <Maestro.H>
#include <AMReX_buildInfo.H>

using namespace amrex;

// write plotfile to disk
void
Maestro::WriteDiagFile (const int step,
                        const Real t_in,
                        const Vector<Real>& rho0_in,
                        const Vector<Real>& p0_in,
                        const Vector<MultiFab>& u_in,
                        Vector<MultiFab>& s_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::WriteDiagFile()",WriteDiagFile);

    // w0mac will contain an edge-centered w0 on a Cartesian grid,   
    // for use in computing divergences.
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
    
    // w0r_cart is w0 but onto a Cartesian grid in cell-centered as
    // a scalar.  Since w0 is the radial expansion velocity, w0r_cart
    // is the radial w0 in a zone
    Vector<MultiFab> w0r_cart(finest_level+1);

    if (spherical == 1) {

	for (int lev=0; lev<=finest_level; ++lev) {
	    AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);,
			 w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);,
			 w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1););

	    for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
		w0mac[lev][idim].setVal(0.);
	    }

	    w0r_cart[lev].define(grids[lev], dmap[lev], 1, 0);
	    w0r_cart[lev].setVal(0.);
	}

	// put w0 on Cartesian edges as a vector  
	MakeW0mac(w0mac);

	// put w0 in Cartesian cell-centers as a scalar (the radial 
	// expansion velocity)
	Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_u,0);

    } else {
	Abort("ERROR: WriteDiagFile() not supported for non-spherical geometry");
    }

    // initialize diagnosis variables
    Real T_max=0.0, T_center=0.0; 
    Vector<Real> coord_Tmax(AMREX_SPACEDIM,0.0);
    Vector<Real> vel_Tmax(AMREX_SPACEDIM,0.0);
    Real Mach_max=0.0;

    for (int lev=0; lev<=finest_level; ++lev) {
	
	// get references to the MultiFabs at level lev
	const MultiFab& sin_mf = s_in[lev];
	const MultiFab& uin_mf = u_in[lev];
        const MultiFab& w0macx_mf = w0mac[lev][0];
	const MultiFab& w0macy_mf = w0mac[lev][1];
	const MultiFab& w0macz_mf = w0mac[lev][2];
	const MultiFab& w0rcart_mf = w0r_cart[lev];
	const MultiFab& normal_mf = normal[lev];
	const Real* dx = geom[lev].CellSize();

        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
        for ( MFIter mfi(sin_mf); mfi.isValid(); ++mfi ) {

	    // Get the index space of the valid region
            const Box& validBox = mfi.validbox();

	    // call fortran subroutine
            // use macros in AMReX_ArrayLim.H to pass in each FAB's data, 
            // lo/hi coordinates (including ghost cells), and/or the # of components
	    diag_sphr(&lev, ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
		      BL_TO_FORTRAN_FAB(sin_mf[mfi]),
		      rho0_in.dataPtr(), p0_in.dataPtr(), 
		      BL_TO_FORTRAN_3D(uin_mf[mfi]), 
		      BL_TO_FORTRAN_3D(w0macx_mf[mfi]), 
		      BL_TO_FORTRAN_3D(w0macy_mf[mfi]),
		      BL_TO_FORTRAN_3D(w0macz_mf[mfi]), 
		      BL_TO_FORTRAN_3D(w0rcart_mf[mfi]),
		      dx, 
		      BL_TO_FORTRAN_3D(normal_mf[mfi]), 
		      &T_max, coord_Tmax.dataPtr(), vel_Tmax.dataPtr(), 
		      &T_center, &Mach_max);
	}
    }

    // sum T_center over all processors
    ParallelDescriptor::ReduceRealSum(T_center);

    // find the largest Mach number over all processors
    ParallelDescriptor::ReduceRealMax(Mach_max);

    // for T_max, we want to know where the hot spot is, so we do a
    // gather on the temperature and find the index corresponding to
    // the maxiumum.  We then pack the coordinates and velocities
    // into a local array and gather that to the I/O processor and
    // pick the values corresponding to the maximum.
    int nprocs = ParallelDescriptor::NProcs();
    int ioproc = ParallelDescriptor::IOProcessorNumber();

    Vector<Real> Tmax_global(nprocs);
    
    if (nprocs == 1) {
	Tmax_global[0] = T_max;
    } else {
	ParallelDescriptor::Gather(&T_max, 1, &Tmax_global[0], 1, ioproc);
    }

    // determine index of max global T
    int index_max = 0;
    Real T_max_data = 0.0;
    for (int ip=0; ip<nprocs; ++ip) {
	if (Tmax_global[ip] > T_max_data) {
	    T_max_data = Tmax_global[ip];
	    index_max = ip;
	}
    }
    
    // T_max_coords will contain both the coordinate information and
    // the velocity information, so there are 2*dm values on each
    // proc
    Vector<Real> T_max_coords(2*AMREX_SPACEDIM);
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
	T_max_coords[i] = coord_Tmax[i];
	T_max_coords[i+3] = vel_Tmax[i];
    }

    Vector<Real> Tmax_coords_global(2*AMREX_SPACEDIM*nprocs);

    if (nprocs == 1) {
	for (int i=0; i<2*AMREX_SPACEDIM; ++i) {
	    Tmax_coords_global[i] = T_max_coords[i];
	}
    } else {
	ParallelDescriptor::Gather(&T_max_coords[0], 6, &Tmax_coords_global[0], 6, ioproc);
    }

    // initialize global variables
    Vector<Real> coord_Tmax_data(AMREX_SPACEDIM);
    Vector<Real> vel_Tmax_data(AMREX_SPACEDIM);
    Real Rloc_Tmax, vr_Tmax;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
	coord_Tmax_data[i] = Tmax_coords_global[2*AMREX_SPACEDIM*index_max+i];
	vel_Tmax_data[i] = Tmax_coords_global[2*AMREX_SPACEDIM*index_max+i+3];
    }

    //
    // reduce the current level's data with the global data
    //
    if (ParallelDescriptor::IOProcessor()) {

	// if T_max_data is the new max, then copy the location as well
	if ( T_max_data > T_max ) {
	    T_max = T_max_data;
	    
	    for (int i=0; i<AMREX_SPACEDIM; ++i) {
		coord_Tmax[i] = coord_Tmax_data[i];
		vel_Tmax[i] = vel_Tmax_data[i];
	    }

	    // compute center of domain
	    Vector<Real> center(3, 0.0);
	    const Real* probLo = geom[0].ProbLo();
	    const Real* probHi = geom[0].ProbHi();

	    for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
		center[idim] = 0.5*(*(probLo+idim) + *(probHi+idim));
	    }

	    // compute the radius of the bubble from the center
	    Rloc_Tmax = sqrt( (coord_Tmax[0] - center[0])*(coord_Tmax[0] - center[0]) + 
			      (coord_Tmax[1] - center[1])*(coord_Tmax[1] - center[1]) + 
			      (coord_Tmax[2] - center[2])*(coord_Tmax[2] - center[2]) );

	    // use the coordinates of the hot spot and the velocity components
            // to compute the radial velocity at the hotspot
	    vr_Tmax = ((coord_Tmax[0] - center[0])/Rloc_Tmax)*vel_Tmax[0] + 
		      ((coord_Tmax[1] - center[1])/Rloc_Tmax)*vel_Tmax[1] + 
		      ((coord_Tmax[2] - center[2])/Rloc_Tmax)*vel_Tmax[2];
	    
	    T_center /= 8.0;  // ncenter = 8
	}

    }


    // write out diagnosis data
    if (ParallelDescriptor::IOProcessor()) {
	const std::string& diagfilename = "diag_temp.out";
	std::ofstream diagfile;
	
	if (step == 0) {
	    // create file after initialization
	    diagfile.open(diagfilename, std::ofstream::out | 
			       std::ofstream::trunc | std::ofstream::binary);

	    // write variable names
	    diagfile << "time \t max{T} \t x(max{T}) \t y(max{T}) \t z(max{T}) \t";
	    diagfile << "vx(max{T}) \t vy(max{T}) \t vz(max{T}) \t";
	    diagfile << "R(max{T}) \t vr(max{T}) \t T_center" << endl;

	} else {
	    // append to file 
	    diagfile.open(diagfilename, std::ofstream::out | 
			       std::ofstream::app | std::ofstream::binary);
	}

	diagfile.precision(15);
	diagfile << t_in << "\t" << T_max <<"\t"; 
	diagfile << coord_Tmax[0] << "\t" << coord_Tmax[1] << "\t" << coord_Tmax[2] << "\t";
	diagfile << vel_Tmax[0] << "\t" << vel_Tmax[1] << "\t" << vel_Tmax[2] << "\t"; 
	diagfile << Rloc_Tmax << "\t" << vr_Tmax << "\t" << T_center << endl;

	// close file
	diagfile.close();
    }

}


// put together a vector of multifabs for writing
Vector<const MultiFab*>
Maestro::DiagFileMF (const Vector<MultiFab>& p0_cart,
                     const Vector<MultiFab>& rho0_cart,
                     const Vector<MultiFab>& u_in,
                           Vector<MultiFab>& s_in,
                     const Vector<Real>& p0_in)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::DiagFileMF()",DiagFileMF);

    // velocities (AMREX_SPACEDIM)
    // rho, rhoh, rhoX, tfromp, tfromh, Pi (Nscal+1)
    // X (NumSpec)
    // rho0, p0 (2)
    // deltaT (1)
    // w0 (AMREX_SPACEDIM)
    int nPlot = 2*AMREX_SPACEDIM + Nscal + NumSpec + 4;

    // MultiFab to hold plotfile data
    Vector<const MultiFab*> plot_mf;

    // temporary MultiFab to hold plotfile data
    Vector<MultiFab*> plot_mf_data(finest_level+1);

    // temporary MultiFab for calculations
    Vector<MultiFab> tempmf(finest_level+1);

    int dest_comp = 0;

    // build temporary MultiFab to hold plotfile data
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i] = new MultiFab((s_in[i]).boxArray(),(s_in[i]).DistributionMap(),nPlot         ,0);
        tempmf[i].define(grids[i],dmap[i],AMREX_SPACEDIM,0);
    }

    // velocity
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((u_in[i]),0,dest_comp,AMREX_SPACEDIM);
    }
    dest_comp += AMREX_SPACEDIM;

    // rho
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),Rho,dest_comp,1);
    }
    ++dest_comp;

    // rhoh
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),RhoH,dest_comp,1);
    }
    ++dest_comp;

    // rhoX
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),FirstSpec,dest_comp,NumSpec);
    }
    dest_comp += NumSpec;

    // X
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),FirstSpec,dest_comp,NumSpec);
        for (int comp=0; comp<NumSpec; ++comp) {
            MultiFab::Divide(*plot_mf_data[i],s_in[i],Rho,dest_comp+comp,1,0);
        }
    }
    dest_comp += NumSpec;

    // compute tfromp
    TfromRhoP(s_in,p0_in);

    // tfromp
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
    }
    ++dest_comp;

    // compute tfromh
    TfromRhoH(s_in,p0_in);

    for (int i = 0; i <= finest_level; ++i) {
        // tfromh
        plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
    }
    ++dest_comp;

    // deltaT
    // compute & copy tfromp
    TfromRhoP(s_in,p0_in);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((s_in[i]),Temp,dest_comp,1);
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
        plot_mf_data[i]->copy((s_in[i]),Pi,dest_comp,1);
    }
    ++dest_comp;

    // rho0 and p0
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((rho0_cart[i]),0,dest_comp  ,1);
        plot_mf_data[i]->copy((  p0_cart[i]),0,dest_comp+1,1);
    }
    dest_comp += 2;

    if (spherical == 1) {
	Vector<std::array< MultiFab, AMREX_SPACEDIM > > w0mac(finest_level+1);
	Vector<MultiFab> w0r_cart(finest_level+1);

	for (int lev=0; lev<=finest_level; ++lev) {
	    // w0mac will contain an edge-centered w0 on a Cartesian grid,
	    // for use in computing divergences.
	    AMREX_D_TERM(w0mac[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 1);,
			 w0mac[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 1);,
			 w0mac[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 1););
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
	    Put1dArrayOnCart(w0,w0r_cart,1,0,bcs_f,0);
	}
    }   // spherical

    // w0
    Put1dArrayOnCart(w0,tempmf,1,1,bcs_u,0);
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf_data[i]->copy((tempmf[i]),0,dest_comp,AMREX_SPACEDIM);
    }
    dest_comp += AMREX_SPACEDIM;

    // add plot_mf_data[i] to plot_mf
    for (int i = 0; i <= finest_level; ++i) {
        plot_mf.push_back(plot_mf_data[i]);
    }

    return plot_mf;

}

// set plotfile variable names
Vector<std::string>
Maestro::DiagFileVarNames () const
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::PlotFileVarNames()",PlotFileVarNames);

    // velocities (AMREX_SPACEDIM)
    // rho, rhoh, rhoX, tfromp, tfromh, Pi (Nscal+1)
    // X (NumSpec)
    // rho0, p0 (2)
    // deltaT (1)
    // w0 (AMREX_SPACEDIM)
    int nPlot = 2*AMREX_SPACEDIM + Nscal + NumSpec + 4;
    Vector<std::string> names(nPlot);

    int cnt = 0;

    // add velocities
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        names[cnt++] = x;
    }

    // density and enthalpy
    names[cnt++] = "rho";
    names[cnt++] = "rhoh";

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
    }

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
    }

    names[cnt++] = "tfromp";
    names[cnt++] = "tfromh";
    names[cnt++] = "deltaT";
    names[cnt++] = "Pi";

    names[cnt++] = "rho0";
    names[cnt++] = "p0";

    // add w0
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "w0";
        x += (120+i);
        names[cnt++] = x;
    }

    return names;
}
