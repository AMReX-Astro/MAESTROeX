[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Build Status](https://travis-ci.com/AMReX-Astro/MAESTROeX.svg?branch=development)](https://travis-ci.com/AMReX-Astro/MAESTROeX)

[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)

# MAESTROeX
*a C++/F90 low Mach number stellar hydrodynamics code*

MAESTROeX solves the equations of low Mach number hydrodynamics for
stratified atmospheres/full spherical stars with a general equation of state,
and nuclear reaction networks in an adaptive-grid finite-volume framework. It
includes reactions and thermal diffusion and can be used on anything
from a single core to 100,000s of processor cores with MPI + OpenMP.

A description of the algorithm and links to the algorithm papers can be found here:

http://amrex-astro.github.io/MAESTROeX/


## Getting started

- To stay up-to-date with MAESTROeX, you will want to periodically pull changes
from the repository by typing `git pull`.

- To get things running, you will need to get a copy of Microphysics and AMReX.
Both are available on github and can be obtained via:

```
git clone https://github.com/starkiller-astro/Microphysics.git
git clone https://github.com/AMReX-Codes/amrex.git
```

- You will then need to setup your shell environment to tell MAESTROeX where to
find AMReX and Microphysics. Define the `AMREX_HOME` environment variable to point
to the `amrex/` directory, and `MICROPHYSICS_HOME` environment variable to point
to the `Microphysics/` directory. For example, if your shell is Bash:

```
export AMREX_HOME="/path/to/amrex/"
export MICROPHYSICS_HOME='/path/to/Microphysics" 
```

Note: you must specify the full path to the directories. 
Do not use `∼` to refer to your home directory; the scripts used by 
the build system will not be able to process this.

- Change directory to correspond to the problem that you want to run. Each problem lives under 
one of three sub-directories of `MAESTROeX/Exec/`: `SCIENCE/`, `TEST_PROBLEMS/`, or `UNIT_TESTS/`.
Then build the executable and run it by specifying an input file. 

    * For example, to run the standard reacting_bubble problem in 2-D:

    ```
    cd MAESTROeX/Exec/TEST_PROBLEMS/reacting_bubble/
    make DIM=2
    ./Maestro2d.gnu.ex inputs_2d_C
    ```

- The plotfiles (named `pltXXXXXXX`) are in BoxLib/AMReX format and can be visualized 
using Amrvis, VisIt, and yt.


For more detailed instructions on how to run the code and available test problems, 
refer to MAESTROeX User's Guide: 

https://amrex-astro.github.io/MAESTROeX/docs/getting_started.html


## Regression and unit testing

MAESTROeX needs to be tested in tandem with the AMReX and Microphysics
repo updates.  This testing is done on local machines using the AMReX
regression test framework
(https://github.com/AMReX-Codes/regression_testing).  Tests are run
nightly and reported here:

https://ccse.lbl.gov/pub/RegressionTesting/MAESTROeX/
http://groot.astro.sunysb.edu/MAESTROeX/test-suite/gfortran/

A number of small unit tests exist in `Exec/UNIT_TESTS` for testing
physics solvers independently.


## Development Model:

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `master`.

  * Contributions are welcomed from anyone.  *Any contributions that
    have the potential to change answers should be done via pull
    requests.*   A pull request should be generated from your fork of
    MAESTROeX and target the `development` branch.  (If you mistakenly
    target `master`, we can change it for you.)

    If there are a number of small commits making up the PR, we may
    wish to squash commits upon merge to have a clean history.
    *Please ensure that your PR title and first post are descriptive,
    since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `master`, in coordination with `AMReX`,
    `MAESTROeX`, and `Microphysics`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `master` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of MAESTROeX.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 10 non-merge commits to `MAESTROeX/Source/` or
    `MAESTROeX/sphinx_docs/` or one of the problems that is not your
    own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * listed in the author list on the Zenodo DOI for the project
    (as given in the .zenodo.json file)

  * invited to co-author general code papers / proceedings describing
    MAESTROeX, its performance, etc.  (Note: science papers will always
    be left to the science leads to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.


## Getting help

Join the mailing list to ask for help or stay up-to-date:

https://groups.google.com/forum/#!forum/maestro-help
