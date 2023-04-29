# 23.03

  * Fixed some clang-tidy issues and compiler warnings (#356, #359,
    #360, #364)

  * Updated xrb_mixed so it works again (#330, #362)

  * Fixed some memory access issues (#337, #338, #349, #352, #353)

# 22.07

  * No changes since last release

# 22.06

  * Bug fix for running code on GPUs

  * Fix precision bug for checkpoint files

# 22.05

  * No changes since last release

# 22.03

  * No changes since last release

# 22.02

  * Fix bug in restart functionality for flame problem

  * Sync up with changes to Microphysics

# 22.01

  * Remove all Fortran Microphysics

  * Remove old Fortran constants module, sponge, inlet BC, extern
    parameter support, and files that are no longer used

  * Add problem namespace to the problem runtime parameters

  * Convert embedded probins to C++ inputs for all problems

  * Change networks for test problems: xrb_mixed (xrb_simple ->
    rprox), wdconvect (ignition_chamulak -> ignition_simple)

  * Remove Travis references (Travis is no longer used)

# 21.12

  * Remove Fortran files that are no longer used

  * Remove Make.cuda_rules

# 21.11

  * No changes since last release

# 21.10

  * Bug fix for NSE problems

# 21.09

  * No changes since last release

# 21.08

  * Pressure correction for impenetrable top boundary for planar problems

# 21.07

  * Fix sphinx4 latex macro rendering

# 21.06

  * Change blocking factor to 8

# 21.05

  * Fix labels for tfromp and tfromh in output files

  * Remove reference to rtol_temp and atol_temp to sync with Microphysics

# 21.04

  * Add a yt VR of ECSN radial velocity

# 21.03

  * Bug fixes for simplified SDC algorithm

  * GPU bug fixes

# 21.02

  * Bug fix for exact base state algorithm

  * Sync up with latest updates of Microphysics

  * Converted all test problem routines to C++

# 21.01

  * Converted all remaining routines to C++

  * Fixed exact base state algorithm and combined MaestroAdvanceIrreg with MaestroAdvanceAvg subroutine

  * Converted all kernel launches to use amrex::ParallelFor for easier debugging

  * C++17 is now the minimum supported C++ standard.

# 20.12

  * Set default in SDC algorithm to use split projection for both planar and spherical problems

  * Bug fix: spherical SDC with split projection

  * Added test problem: imposed_external_heating

  * Fixed create_release action and release action python script

# 20.11

  * Update cpp linter action

# 20.10

  * Added input parameter to toggle centrifugal force

  * Enabled control for toggling AMReX fancy bottom solver

  * Updated test problem: fully convective star

# 20.09

  * Bug fix for SDC multilevel problems

  * Fixed non-deterministic issue with GPU runs

  * Parameters header files are now created at build time

  * Updated Urca test problem

# 20.08

  * Bug fix: tagging issues are resolved when finest_level < max_radial_level

  * Fixed dpdt implementation and SDC algorithm for thermal diffusion

  * Various bug fixes

  * Implemented clang-format

  * Added Github actions: version release, checkout submodule development branches

# 20.07

  * The master branch has been renamed the main branch. If you have an
    existing clone of MAESTROeX, then do the following to update for this
    change. First, do `git checkout master` if you're not already on the
    old master branch. Then do `git pull`. This will gather the updates
    to the repo, but will fail with the message `Your configuration specifies
    to merge with the ref 'refs/heads/master' from the remote, but no such ref
    was fetched.` Then you can simply do `git checkout main` and your local
    repo should automatically switch to that branch and track updates from
    the upstream repo on GitHub. If you like, you can then delete the old
    master branch with `git branch -D master`.

  * Fixed unit test: test_basestate

# 20.06

  * Completed converting base state vectors to C indexing

  * SDC bug fixes: use_tfromp=T

  * Fixed probin issue when compiling with PGI

  * Added AMReX as a submodule

  * Switch docs build to a Github action instead of Travis

# 20.05

  * Created new base state and base state geometry classes

  * Started to convert base state vectors (based on Fortran indexing) to new base state class that uses C indexing

  * Offloaded more routines to GPU: InitData, MakeThermalCoeffs

# 20.04a

  * Include latest release from Microphysics submodule

# 20.04

  * Offloaded more routines to GPU: MakePsi, MakeEtarho, FillUmacGhost, BCFill, Sponge, Dt, BaseState, EOS routines, Average

  * Converted code to use Microphysics directly

  * Fixed GPU and OMP reduction calls

# 20.03

  * Offloaded more routines to GPU: MakeUtrans, MacProj, NodalProj, Makew0, MakeBeta0, AdvectBase, EnforceHSE

  * Fixed non-deterministic OMP+MPI behavior

  * Update regression tests to cover wider range of MPI + OMP + multilevel

  * Merged all parameters into a single module

  * Add test problem: fully_convective_star

  * Various bug fixes

# 20.02

  * Offloaded more routines to GPU: Fill3dData, MakeEdgeScal, VelPred, PPM

  * Various bugfixes

  * Fix vorticity to use correct domain boundaries

  * Add travis clang static analysis

  * Update analysis scripts for Urca fields

# 20.01

  * Various initialization and checkpoint-restart fixes

  * Converted all .f90 files to .F90 files

  * Various tiling and GPU updates and bugfixes

  * Offloading more routines to GPU, including nodal and cell-centered solvers

  * Options to apply cutoff densities for heating and burning to a bound range of densiites

  * Work on the rotating_star problem - models, heating profiles, initialization

  * Rename subdirectories within Exec/ to lower-case to sync with Castro

# 19.12

  * Add option to switch to SDC algorithm (USE_SDC)

  * Continue to port hydro subroutines to GPU

  * Set default Microphysics integrator to VODE90

  * Add Zenodo DOI for JOSS publication

  * GPU-specific changes to Makefile

  * Bugfix: w0 is computed even when base cutoff density is outside domain

  * Update test: rotating_star

# 19.11

  * Port 3d hydro subroutine to GPU

  * Add enthalpy_pred_type options: predict_T_then_rhohprime,
    predict_T_then_h, predict_Tprime_then_h

  * Updates to documentation for JOSS paper, including new logo

  * Add scripts for plotting and comparing plotfiles

  * Bugfixes: boundary condition for thermal diffusion, fill w0_cart
    after restart/regrid

  * Update tests: rotating_star, code_comp, multilevel reacting_bubble


# 19.10

  * Move base state from 1D array to cartesian grid in many subroutines

  * Start porting hydro subroutines to GPU (ongoing)

  * Bug fix: fill w0_cart based on edge-centered w0, not cell-centered

  * Add option to output Nonaka plot (USE_NONAKA_PLOT)

  * Add diagnostics for planar geometries

  * Updates to documentation and JOSS paper


# 19.09

  * Minor changes to some problem setups and bugfixes

  * Remove 1D support (with compile-time error)

  * Significant updates to sphinx documentation, including doxygen

  * JOSS paper added; ready for submission


# 19.08

  * Option to plot gravity

  * Added checkpoint at simulation time level option

  * Fixed a ppm_type=2 index bug

  * GPU ports: FirstDtm EstDt

  * Various GPU optimizations

  * Boundary condition bugs in make_edge_scal


# 19.07

  * Add BSD-3 license

  * Minor bugfixes at initialization

  * Add diagnostic outputs of radially averaged profiles

  * Add tests: test_stability, code_comp (code comparison project for Exeter)


# 19.06

  * Implement multiple tagging criteria for multiple levels of AMR

  * Specify how often to print diagnostics using sum_interval and
    sum_per

  * Plot files now include base states output

  * Bugfix to new temporal algorithm in MaestroAdvanceAvg.cpp

  * Many thermodynamics and source term subroutines offloaded to GPU:
    burner_loop, make_thermal_coeffs, mk_sponge, rhoh_vs_t,
    make_gamma, mkrhohforce, make_S, update_scal, update_vel,
    compute_dt, make_pi_cc

  * Add cuda managed attributes to meth_params.F90 using
    parse_maestro_params.py

  * Add RealVector typedef so that Real Vectors can be CUDA managed

  * Move put_1d_array_on_sphr routine calls to C++ to avoid calls on the GPU

  * If USE_CUDA=TRUE, then set USE_GPU_PRAGMA=TRUE and add CUDA to
    define

  * Remove loop over level in 1d divu_cart and gradp0

  * Offloaded plotting routines in make_plot_variables.F90 to GPU


# 19.05

  * Travis documentation updates

  * Plotfile and small_plotfile logic fixed

  * Additional variables in plotfiles

  * Some progress toward the do_smallscale option

  * Bugfix to parsing strings with spaces in parse_maestro_params.py

  * Added mach_jet problem and inflow routine

  * Added fill_ext_bc routines framework to local copies of
    bc_fill_nd.F90 to support Dirichlet boundary conditions

  * Custom BndryFuncArrayMaestro that allows the bc fill routines to
    be component-aware

  * MLMG boundary conditions fixes

  * Started porting flame, xrb_mixed, sub_chandra problems

  * Make ppm_trace_foces true by default

  * Fixed refluxing bug, and implemented reflux_type (0=none,
    1=average-down, 2=reflux); default is 1 for now

  * Force the user to supply anelastic_cutoff_density and
    base_cutoff_density

  * Set burning_cutoff_density to base_cutoff_density if not specified


# 19.04

  * New temporal algorithm is now implemented correctly by separating
    full velocity state into w0 and perturbational velocity during
    each projection solve

  * Additional diagnostics and timing outputs

  * New rotating star model files

  * Fix test_react unit tests


# 19.03

  * Tagging changes for multilevel

  * Planar multilevel now works

  * Spherical multilevel works in limited capacity (still debugging)

  * New input files and initial model files

  * Fix unit tests: test_advect, test_diffusion, test_projection


# 19.02

  * Various planar multilevel bug fixes, and general multilevel bug
    fixes

  * Switch to averaging down edge fluxes instead of using the
    refluxing routines

  * New rotating star examples

  * New unit tests: test_advect, test_eos, test_project, test_react,
    test_diffusion, test_basestate


# 19.01

  * Add more plotfile variables

  * Modification to dp0/dt treatment in new temporal algorithm

  * Embed probin file at the bottom of inputs files

  * Multilevel planar and spherical work/bug fixes

  * Sphinx-based user's guide on website enabled


# 18.12

  * Miscellaneous bug fixes

  * Work on evolve_bases_state for irregular base state

  * Simplified boundary conditions class.

  * Minor changes to regression tests

  * Update interface to StarKiller Microphysics conductivity routine

    * Put conductivity in eos_type. Revise interface to conductivity routine

    * Update actual_conductivity interface in constant conductivity


# 18.11

  * Implemented OpenMP in 1D routines; tiling OpenMP in 3D routines

  * Added problem-specific runtime parameters to each problem
    directory

  * Eliminated dependency on FBoxLib; now only requires AMReX

  * Switch probin files to all use the namelist 'probin', not
    'probdata'

  * fix mask bug in reactions

  * corner ghost cell fix in fill_umac_ghost

  * various minor changes to wdconvect test problems

  * fixed a bug with use_exact_base_state where edge states weren't
    being computed properly

  * Irregular base state developments

  * Update burn_type to be compatible with Microphysics developments

  * Debugging utility for writing out a vector of MultiFabs.

  * Increased diagnostics in Diag

  * Simplify conductivity interface to mimic the eos interface.


# 18.09

  * New test problems: rt, double_bubble, shear jet, test_convect

  * Some multilevel bug fixes (multilevel still WIP)

  * More progress on the new use_exact_base_state algorithm (still
    WIP)


# 18.08

  * Load balancing improvements

  * Spherical variable-base state implemented (still testing and
    debugging)

  * Added files to help with AMR debugging with constant base-state
    spherical case

  * Implemented correct corner coupling for 3d conservative quantity
    advection

  * Implemented dpdt term for variable-base state spherical MAC
    projection


# 18.07

  * Initial commit of development/master branches.

  * Planar and spherical work for single level grids.

  * wdconvect and reacting bubble tests implemented.

  * AMR and variable spherical base state spacing are WIP.
