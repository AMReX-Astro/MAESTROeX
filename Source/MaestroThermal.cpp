
#include <Maestro.H>
#include <AMReX_VisMF.H>

using namespace amrex;

////////////////////////////////////////////////////////////////////////////
// Compute the quantity: thermal = del dot kappa grad T
// 
//  if temp_diffusion_formulation = 1, then we compute this directly.
//  if temp_diffusion_formulation = 2, then we compute the algebraically
//     equivalent form with grad h - grad X_k - grad p_0 formulation
///////////////////////////////////////////////////////////////////////////
void 
Maestro::MakeExplicitThermal(Vector<MultiFab>& thermal, 
			     const Vector<MultiFab>& scal, 
			     const Vector<MultiFab>& Tcoeff, 
			     const Vector<MultiFab>& hcoeff, 
			     const Vector<MultiFab>& Xkcoeff, 
			     const Vector<MultiFab>& pcoeff, 
			     const Vector<Real>& p0, 
			     int temp_formulation) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeExplicitThermal()",MakeExplicitThermal);

    for (int lev = 0; lev <= finest_level; ++lev) {
	thermal[lev].setVal(0.);
    }

    // Variable to store values to be acted upon by (div B grad) operator
    Vector<MultiFab> phi(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
	phi[lev].define(grids[lev], dmap[lev], 1, 0);
	phi[lev].setVal(0.);
    }

    // 
    // Compute thermal = div B grad phi using MLABecLaplacian class
    //  
    LPInfo info;

    // turn off multigrid coarsening since no actual solve is performed
    info.setMaxCoarseningLevel(0);

    MLABecLaplacian mlabec(geom, grids, dmap, info);

    // order of stencil
    int stencil_order = 2;
    mlabec.setMaxOrder(stencil_order);

    for (int lev = 0; lev <= finest_level; ++lev) {
	mlabec.setLevelBC(lev, &thermal[lev]);
    }

    if (temp_formulation == 1) 
    {
	// compute div Tcoeff grad T
	// alpha is set to zero 
	mlabec.setScalars(0.0, -1.0);

	// set value of phi
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Copy(phi[lev],scal[lev],Temp,0,1,0);
	}

	ApplyThermal(mlabec, thermal, Tcoeff, phi, Temp); 
    }
    else { // if temp_formulation == 2

	// temporary residual variable
	Vector<MultiFab> resid(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    resid[lev].define(grids[lev], dmap[lev],1, 0);
	    resid[lev].setVal(0.);
	}

	// 1. Compute div hcoeff grad h
	mlabec.setScalars(0.0, -1.0);

	// set value of phi
	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Copy(phi[lev],scal[lev],RhoH,0,1,0);
	    MultiFab::Divide(phi[lev],scal[lev],Rho,0,1,0);
	}

	ApplyThermal(mlabec, resid, hcoeff, phi, RhoH); 

	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Add(thermal[lev],resid[lev],0,0,1,0);
	}

	// 2. Compute -div Xkcoeff grad Xk
	mlabec.setScalars(0.0, 1.0);

	// temporary Xk coeff variable
	Vector<MultiFab> Xicoeff(finest_level+1);
	for (int lev = 0; lev <= finest_level; ++lev) {
	    Xicoeff[lev].define(grids[lev], dmap[lev], 1, 1);
	}

	for (int nspec=FirstSpec; nspec<FirstSpec+NumSpec; ++nspec) {
	    // set value of phi
	    for (int lev=0; lev<=finest_level; ++lev) {
		MultiFab::Copy(phi[lev],scal[lev],nspec,0,1,0);
		MultiFab::Divide(phi[lev],scal[lev],Rho,0,1,0);
		MultiFab::Copy(Xicoeff[lev],Xkcoeff[lev],nspec-FirstSpec,0,1,1);
	    }

	    ApplyThermal(mlabec, resid, Xicoeff, phi, nspec); 

	    for (int lev=0; lev<=finest_level; ++lev) {
		MultiFab::Add(thermal[lev],resid[lev],0,0,1,0);
	    }
	}

	// 3. Compute -div pcoeff grad p0
	mlabec.setScalars(0.0, 1.0);
	
	// set value of phi
	Put1dArrayOnCart(p0, phi, 0, 0, bcs_s, Temp);

	ApplyThermal(mlabec, resid, pcoeff, phi, Temp); 

	for (int lev=0; lev<=finest_level; ++lev) {
	    MultiFab::Add(thermal[lev],resid[lev],0,0,1,0);
	}
    } // end if

}

// Use apply() to construct the form of the conduction term. 
// apply() forms the generic quantity:
//
//   (alpha * A - beta * div B grad) phi = RHS
void Maestro::ApplyThermal(MLABecLaplacian& mlabec, 
			   Vector<MultiFab>& thermalout, 
			   const Vector<MultiFab>& coeff,
			   Vector<MultiFab>& phi, 
			   int comp) 
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ApplyThermal()",ApplyThermal);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
	if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else {
	    // lo-side BCs
            if (bcs_s[comp].lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (bcs_s[comp].lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

	    // hi-side BCs
            if (bcs_s[comp].hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (bcs_s[comp].hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
	}
    }

    mlabec.setDomainBC(mlmg_lobc,mlmg_hibc);
    
    // coefficients for solver
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > face_bcoef(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        AMREX_D_TERM(face_bcoef[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);,
                     face_bcoef[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);,
                     face_bcoef[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0););
    }

    // average face-centered B coefficients
    PutDataOnFaces(coeff, face_bcoef, 1);

    // set coefficient matrix
    for (int lev = 0; lev <= finest_level; ++lev) {
	mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef[lev]));
    }

    // build an MLMG solver
    MLMG thermal_mlmg(mlabec);

    thermal_mlmg.apply(GetVecOfPtrs(thermalout), GetVecOfPtrs(phi));
}


////////////////////////////////////////////////////////////////////////////
// create the coefficients for grad{T}, grad{h}, grad{X_k}, and grad{p_0}
// for the thermal diffusion term in the enthalpy equation.  
////////////////////////////////////////////////////////////////////////////
void 
Maestro::MakeThermalCoeffs(const Vector<MultiFab>& scal,
			   Vector<MultiFab>& Tcoeff,
			   Vector<MultiFab>& hcoeff, 
			   Vector<MultiFab>& Xkcoeff, 
			   Vector<MultiFab>& pcoeff) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeThermalCoeffs()",MakeThermalCoeffs);

    for (int lev = 0; lev <= finest_level; ++lev) 
    {
	// get references to the MultiFabs at level lev
	const MultiFab& scal_mf = scal[lev];
	MultiFab& Tcoeff_mf  = Tcoeff[lev];
	MultiFab& hcoeff_mf  = hcoeff[lev];
	MultiFab& Xkcoeff_mf = Xkcoeff[lev];
	MultiFab& pcoeff_mf  = pcoeff[lev];

	Print() << "... Level " << lev << " create thermal coeffs:" << endl;

	// loop over boxes
	for ( MFIter mfi(scal_mf); mfi.isValid(); ++mfi) {

	    // Get the index space of valid region
	    const Box& validBox = mfi.validbox();

	    // call fortran subroutine
	    make_thermal_coeffs(ARLIM_3D(validBox.loVect()),ARLIM_3D(validBox.hiVect()), 
				BL_TO_FORTRAN_FAB(scal_mf[mfi]), 
				BL_TO_FORTRAN_3D(Tcoeff_mf[mfi]), 
				BL_TO_FORTRAN_3D(hcoeff_mf[mfi]), 
				BL_TO_FORTRAN_FAB(Xkcoeff_mf[mfi]),
				BL_TO_FORTRAN_3D(pcoeff_mf[mfi]));

	}
    }
    
}

////////////////////////////////////////////////////////////////////////////
// ThermalConduct implements thermal diffusion in the enthalpy equation.
// This is an implicit solve, using the multigrid solver.  This updates
// the enthalpy only. 
////////////////////////////////////////////////////////////////////////////
void
Maestro::ThermalConduct (const Vector<MultiFab>& s1,
			 Vector<MultiFab>& s2,
			 const Vector<MultiFab>& hcoeff1,
			 const Vector<MultiFab>& Xkcoeff1,
			 const Vector<MultiFab>& pcoeff1,
			 const Vector<MultiFab>& hcoeff2, 
			 const Vector<MultiFab>& Xkcoeff2, 
			 const Vector<MultiFab>& pcoeff2,
			 const Vector<Real>& p0_old,
			 const Vector<Real>& po_new)
{
    // timer for profiling
    BL_PROFILE_VAR("Maestro::ThermalConduct()",ThermalConduct);

    // Dummy coefficient matrix, holds all zeros
    Vector<MultiFab> Dcoeff(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
	Dcoeff[lev].define(grids[lev], dmap[lev], 1, 1);
	Dcoeff[lev].setVal(0.);
    }

    // solverrhs will hold solver RHS = (rho h)^1  + 
    //           dt/2 div . ( hcoeff1 grad h^1) -
    //           dt/2 sum_k div . (Xkcoeff2 grad X_k^2 + Xkcoeff1 grad X_k^1) -
    //           dt/2 div . ( pcoeff2 grad p_0^new + pcoeff1 grad p_0^old)
    Vector<MultiFab> solverrhs(finest_level+1);
    Vector<MultiFab> resid(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        solverrhs[lev].define(grids[lev], dmap[lev], 1, 0);
	resid[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    // compute RHS = (rho h)^1
    for (int lev = 0; lev <= finest_level; ++lev) {
	MultiFab::Copy(solverrhs[lev],s1[lev],RhoH,0,1,0);
    }

    // compute resid = div(hcoeff1 grad h^1) - sum_k div(Xkcoeff1 grad Xk^1) - div(pcoeff1 grad p0_old)
    MakeExplicitThermal(resid,s1,Dcoeff,hcoeff1,Xkcoeff1,pcoeff1,p0_old,2);

    // RHS = solverrhs + dt/2 * resid1
    for (int lev = 0; lev <= finest_level; ++lev) {
	MultiFab::LinComb(solverrhs[lev],1.0,solverrhs[lev],0,dt/2.0,resid[lev],0,0,1,0);
    }

    // compute resid = 0 - sum_k div(Xkcoeff2 grad Xk^2) - div(pcoeff2 grad p0_new)
    MakeExplicitThermal(resid,s2,Dcoeff,Dcoeff,Xkcoeff2,pcoeff2,p0_new,2);

    // RHS = solverrhs + dt/2 * resid2
    for (int lev = 0; lev <= finest_level; ++lev) {
	MultiFab::LinComb(solverrhs[lev],1.0,solverrhs[lev],0,dt/2.0,resid[lev],0,0,1,0);
    }

    // LHS coefficients for solver
    Vector<MultiFab> acoef(finest_level+1);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > face_bcoef(finest_level+1);
    Vector<MultiFab> phi(finest_level+1);
    for (int lev=0; lev<=finest_level; ++lev) {
        acoef[lev].define(grids[lev], dmap[lev], 1, 0);
        AMREX_D_TERM(face_bcoef[lev][0].define(convert(grids[lev],nodal_flag_x), dmap[lev], 1, 0);,
                     face_bcoef[lev][1].define(convert(grids[lev],nodal_flag_y), dmap[lev], 1, 0);,
                     face_bcoef[lev][2].define(convert(grids[lev],nodal_flag_z), dmap[lev], 1, 0););
	phi[lev].define(grids[lev], dmap[lev], 1, 0);
    }

    // set cell-centered A coefficient to zero
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Copy(acoef[lev],s2[lev],Rho,0,1,0);
    }

    // average face-centered Bcoefficients
    PutDataOnFaces(hcoeff2, face_bcoef, 1);

    // initialize value of phi to h^(2) as a guess
    for (int lev=0; lev<=finest_level; ++lev) {
	MultiFab::Copy(phi[lev],s2[lev],RhoH,0,1,0);
	MultiFab::Divide(phi[lev],s2[lev],Rho,0,1,0);
    }

    // 
    // Set up implicit solve using MLABecLaplacian class
    //
    LPInfo info;
    MLABecLaplacian mlabec(geom, grids, dmap, info);

    // order of stencil
    int linop_maxorder = 2;
    mlabec.setMaxOrder(linop_maxorder);

    // set boundaries for mlabec using enthalpy bc's
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
	if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else {
	    // lo-side BCs
            if (bcs_s[RhoH].lo(idim) == Outflow) {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            } else if (bcs_s[RhoH].lo(idim) == Inflow) {
                mlmg_lobc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }

	    // hi-side BCs
            if (bcs_s[RhoH].hi(idim) == Outflow) {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            } else if (bcs_s[RhoH].hi(idim) == Inflow) {
                mlmg_hibc[idim] = LinOpBCType::inflow;
            } else {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
            }
	}
    }

    mlabec.setDomainBC(mlmg_lobc,mlmg_hibc);

    for (int lev = 0; lev <= finest_level; ++lev) {
	mlabec.setLevelBC(lev, &phi[lev]);
    }

    mlabec.setScalars(1.0, dt/2.0);

    for (int lev = 0; lev <= finest_level; ++lev) {
	mlabec.setACoeffs(lev, acoef[lev]);
	mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef[lev]));
    }

    // solve A - dt/2 * div B grad phi = RHS

    // build an MLMG solver
    MLMG thermal_mlmg(mlabec);

    // set solver parameters
    thermal_mlmg.setVerbose(mg_verbose);
    thermal_mlmg.setCGVerbose(cg_verbose);

    // tolerance parameters taken from original MAESTRO fortran code
    Real thermal_tol_abs = -1.e0;
    for (int lev = 0; lev <= finest_level; ++lev) {
	thermal_tol_abs = std::max(thermal_tol_abs, phi[lev].norm0());
    }
    const Real solver_tol_abs = eps_mac*thermal_tol_abs;
    const Real solver_tol_rel = eps_mac;

    // solve for phi
    thermal_mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(solverrhs), solver_tol_rel, solver_tol_abs);
    
    // load new rho*h into s2
    for (int lev = 0; lev <= finest_level; ++lev) {
	MultiFab::Copy(s2[lev],phi[lev],0,RhoH,1,0);
	MultiFab::Multiply(s2[lev],s2[lev],Rho,RhoH,1,0);
    }

    // average fine data onto coarser cells
    AverageDown(s2,RhoH,1);

    // fill ghost cells
    FillPatch(t_old,s2,s2,s2,RhoH,RhoH,1,RhoH,bcs_s);
}
