.. _sec:analysis:

********
Analysis
********


The Postprocessing Routines
===========================

The ``amrex/Tools/Plotfile/`` directory contains a large
number of C++-based analysis routines for AMReX datasets. Many
of these can be used with both MAESTROeX and the compressible
astrophysics code, CASTRO.

To compile these routines, run ``make`` in the directory.

General Analysis Routines
-------------------------

The following routines are generally applicable for any AMReX-based
plotfile. Typing the executable names without any arguments will
provide usage examples.

-  ``fboxinfo.cpp``

   | Print out some basic information about the number of boxes on each
     refinement level and (optionally) the bounds of each of the boxes.

-  ``fcompare.cpp``

   | Compare two plotfiles, zone-by-zone to machine precision, and report
     the L2-norm of the error (both absolute and relative) for each
     variable. This assumes that the grids are identical.
   | With the optional –zone_info var argument, where var
     is the name of a variable, it will also report the full state
     for the zone where var has the largest error.

-  ``fextract.cpp``

   | Extract a 1-d line through a dataset (1-, 2-, or 3-d). This works
     with both uniformly-gridded or AMR datasets. For multi-dimensional
     datasets, the coordinate direction to extract along can be specified.
     The line is always taken through the center of the domain. Either
     a single variable or all variables, along with the coordinate
     information, are output to a file.

-  ``fextrema.cpp``

   | Report the min and max of each variable (or only a single variable)
     in one or more plotfiles.

-  ``fsnapshot.cpp``

   | Create an image (PPM file) of a single variable in a plotfile. For
     3-d, the slice plane through the center of the domain is specified.
     Separate routines exist for 2-d and 3-d datasets.

-  ``ftime.f90``

   | For each plotfile, simply print the simulation time.

-  ``fvarnames.f90``

   | Simply print out the list of variables stored in a plotfile.

Horizontal averaging
--------------------

The ``amrex/Tools/Postprocessing/`` directory contains a few more complex
analysis routines for working with AMReX-based plotfiles.

Of particular use for MAESTROeX data is the ``HorizontalAvg.cpp`` routine.
This will laterally average each of the variables in a plotfile (works for
both 2-d and 3-d). It is written with MAESTROeX plane-parallel geometry
plotfiles in mind, and the averaging is done over the coordinate direction(s)
perpendicular to gravity.
