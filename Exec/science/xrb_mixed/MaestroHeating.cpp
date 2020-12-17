
#include <Maestro.H>

using namespace amrex;

// compute heating term, rho_Hext
void Maestro::MakeHeating(Vector<MultiFab>& rho_Hext,
                          const Vector<MultiFab>& scal) {
    // timer for profiling
    BL_PROFILE_VAR("Maestro::MakeHeating()", MakeHeating);

    static int ilh1, ihe4, ilc12, iln14, ilo16;
    static bool firstCall = true;
    
    if (firstCall) {
	ilh1 = network_spec_index("hydrogen-1");
	ihe4 = network_spec_index("helium-4");
	ilc12 = network_spec_index("carbon-12");
	iln14 = network_spec_index("nitrogen-14");
	ilo16 = network_spec_index("oxygen-16");
	
	firstCall = false;
    }
    
    for (int lev = 0; lev <= finest_level; ++lev) {
        // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(rho_Hext[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // Get the index space of the valid region
            const Box& tileBox = mfi.tilebox();
	    
            const Array4<const Real> scal_arr = scal[lev].array(mfi);
            const Array4<Real> rho_Hext_arr = rho_Hext[lev].array(mfi);

            ParallelFor(tileBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
		Real rho = scal_arr(i, j, k, Rho);
		Real temp;
		// if (drive_initial_convection)
		//    temp = tempbar_init(j)
		// else
		    temp = scal_arr(i,j,k,Temp);
		// endif

		Real T6 = temp / 1.e6;
		Real T613 = std::pow(T6, 1.0/3.0);

		// total CNO abundance
		Real X_CNO = (scal_arr(i,j,k,FirstSpec+ilc12) + 
			      scal_arr(i,j,k,FirstSpec+iln14) + 
			      scal_arr(i,j,k,FirstSpec+ilo16)) / rho;

		// H abundance
		Real X_1 = scal_arr(i,j,k,FirstSpec+ilh1) / rho;

		// CNO heating from Kippenhahn & Weigert, Eq. 18.65
		Real g14 = 1.0 + 2.7e-3*T613 - 7.78e-3*T613*T613 - 1.49e-4*T6;

		// The temperature-insensitive Hot CNO cycle energy generation
		// this factor is from Wallace & Woosley (1981) ApJS 45
		Real factor = 5.86e15;
		Real eps_CNO = factor*X_CNO;
    
		// The triple alpha reaction energy generation
		// taken from Arnett's book pg 225
		// needs more accurate screening factor (f3a); just setting it to unity now
		const Real eps_0 = 3.9e11;
		const Real f3a = 1.0;
		Real t8i = 1.e8/temp;
		Real Y = scal_arr(i,j,k,FirstSpec+ihe4)/rho;		
		eps_CNO += eps_0*f3a*rho*rho*Y*Y*Y*std::exp(-42.94*t8i)*t8i*t8i*t8i;

		rho_Hext_arr(i, j, k) = rho * eps_CNO;
            });
        }
    }

    // average down
    AverageDown(rho_Hext, 0, 1);
}
