[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![Build Status](https://travis-ci.com/AMReX-Astro/MAESTROeX.svg?branch=development)](https://travis-ci.com/AMReX-Astro/MAESTROeX)

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

Tests are run nightly and reported here:

https://ccse.lbl.gov/pub/RegressionTesting/MAESTROeX/


## Getting help

Join the mailing list to ask for help or stay up-to-date:

https://groups.google.com/forum/#!forum/maestro-help