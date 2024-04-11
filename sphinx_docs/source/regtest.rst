*********************************************
Regression Testing and Continuous Integration
*********************************************

Regression Testing
==================

Nightly regression tests are run on MAESTROeX using the AMReX
regression testing framework:
https://github.com/amrex-codes/regression_testing/ .  We use the
nightly approach rather than continuous integration on github because
we need to ensure that we test the changes to AMReX and Microphysics
together with changes in MAESTROeX.

The basic flow of the regression testing framework is as follows:

* check out the ``development`` branch of MAESTROeX, AMReX, and Microphysics

* for each test we defined:

    * ``make realclean`` in the build directory

    * build the executable with any test-specific options (MPI, network, ...)

    * run the test using the inputs file listed in the test definition

    * compare zone-by-zone the results of the test to the stored
      benchmark.  Any differences, no matter how small, result in an
      error and test failure

* generate the test results webpage.

Tests are run both at Stony Brook at LBNL and the output is presented here:

  * LBNL: https://ccse.lbl.gov/pub/RegressionTesting/MAESTROeX/

  * Stony Brook: http://groot.astro.sunysb.edu/MAESTROeX/test-suite/gfortran/


If you want to run your own version of the test suite, you can follow
the quickstart in the ``regression_testing`` README.  The inputs file used
at Stony Brook can be used as a starting point and is found here:

  https://github.com/amrex-astro/actual_test_files


Continuous Integration
======================

We use Github Actions to run integration tests on the code and to build and deploy the documentation.

Currently, Github Actions runs the `clang static analyzer
<https://clang-analyzer.llvm.org/>`_, which finds potential bugs in
the code.  This is run on pull requests to the MAESTROeX github
repo, and are run weekly on the development branch.
