.. _sec:runtime_parameters:

******************
Runtime Parameters
******************


Introduction to Runtime Parameters
==================================

The runtime parameters are defined in the MAESTRO/_parameters
files, and any other problem-specific parameter files. These
parameters are then made available to the code through the
probin_module.

Any runtime parameters defined by the microphysics (either in the
MAESTRO/Microphysics/ source or the external Microphysics/
source are also parsed at build time and read in via the same
namelist in probin.f90. These microphysics runtime parameters
can be accessed via extern_probin_module.

Parameter definitions take the form of:

::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the priority is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.

The following tables list all of the general
MAESTRO runtime parameters.
These tables are generated automatically using the
rp.py script in MAESTRO/docs/runtime_parameters/ by parsing
the MAESTRO/_parameters file. The problem-specific parameters
are not shown here.


Parameters by Namespace
=======================

.. toctree::

   runtime_parameters
