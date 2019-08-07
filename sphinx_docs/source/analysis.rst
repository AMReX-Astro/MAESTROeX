.. _sec:analysis:

********
Analysis
********


The Postprocessing Routines
===========================

The BoxLib/Tools/Postprocessing/F_Src/ directory contains a large
number of Fortran-based analysis routines for BoxLib datasets. Many
of these can be used with both MAESTRO and the compressible
astrophysics code, CASTRO.

To compile any of the individual routines, edit the GNUmakefile
add uncomment the line beginning with ‘programs +=’ containing
the routine you want to build.

General Analysis Routines
-------------------------

The following routines are generally applicable for any BoxLib-based
plotfile. Typing the executable names without any arguments will
provide usage examples.

-  faverage.f90

   Laterally average each of the variables in a plotfile (works for
   both 2-d and 3-d). This is written with MAESTRO plane-parallel geometry plotfiles in mind, and the averaging is done
   over the coordinate direction(s) perpendicular to gravity.

-  fboxinfo.f90

   Print out some basic information about the number of boxes on each
   refinement level and (optionally) the bounds of each of the boxes.

-  fcompare.f90

   | Compare two plotfiles, zone-by-zone to machine precision, and report
     the L2-norm of the error (both absolute and relative) for each
     variable. This assumes that the grids are identical.
   | With the optional –zone_info var argument, where var
     is the name of a variable, it will also report the full state
     for the zone where var has the largest error.

   This is used by in the regression test suite in
   Parallel/util/regtests/.

-  fextract.f90

   Extract a 1-d line through a dataset (1-, 2-, or 3-d). This works
   with both uniformly-gridded or AMR datasets. For multi-dimensional
   datasets, the coordinate direction to extract along can be specified.
   The line is always taken through the center of the domain. Either
   a single variable or all variables, along with the coordinate
   information, are output to a file.

-  fextrema.f90

   Report the min and max of each variable (or only a single variable)
   in one or more plotfiles.

-  fsnapshot2d.f90, fsnapshot3d.f90

   Create an image (PPM file) of a single variable in a plotfile. For
   3-d, the slice plane through the center of the domain is specified.
   Separate routines exist for 2-d and 3-d datasets.

-  ftime.f90

   For each plotfile, simply print the simulation time.

-  fvarnames.f90

   Simply print out the list of variables stored in a plotfile.

Data Processing Example
-----------------------

The routine fspeciesmass2d.f90 in
F_src/tutorial serves as a well-commented example of how
to work with MAESTRO plotfile data. This routine simply computes
the total mass of a particular species on the entire domain for a 2-d
dataset. It is written to understand a multilevel (AMR) dataset, and
only considers the finest-available data at any physical location in
the computational domain.

fspeciesmass2d.f90 should provide a good starting point for
writing a new analysis routine for BoxLib data.

.. _analysis:sec:particles:

Particle routines
-----------------

The parseparticles.py routine in the python/ subdirectory
can read in MAESTRO particle files containing particle histories
(usually named timestamp_??). See the discussion in
§ \ `[arch:sec:particles] <#arch:sec:particles>`__ for details on initializing particles in
MAESTRO. The driver test_parseparticles.py shows shows how to
use this module to plot particle histories. Additional documentation
is available from the module itself. In the python environment,
type:

.. code:: python

    import parseparticles
    help(parseparticles)

to get information on the classes and functions provided by the
parseparticles module.

As a concrete example, running reacting_bubble with particles enabled
will seed particles in the initial hotspot. To plot the results,
first set your PYTHONPATH environment variable to point to the
AmrPostprocessing/python/ directory, for example:

::

    export PYTHONPATH="/home/username/development/AmrPostprocessing/python"

This will allow python to see the parseparticles.py routine.
For the reacting_bubble problem, the plotparticles.py routine shows
how to plot the particle histories and make an animation of the
particles colored by the ash mass fraction. This script is run as:

::

    ./plotparticles.py timestamp_*

Note, these python routines require the NumPy and matplotlib packages.
On a Fedora Linux system, the necessary routines can be installed via:

::

    yum install python-matplotlib lyx-fonts stix-fonts
