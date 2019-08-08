.. _sec:runtime_parameters:

******************
Runtime Parameters
******************


Introduction to Runtime Parameters
==================================

MAESTROeX has 2 sets of runtime parameters—those controlled by
C++ and those controlled by Fortran. The C++ parameters are set
in the inputs file and managed by the AMReX ``ParmParse``
class. For MAESTROeX-specific parameters, we list the runtime
parameters in a file ``Source/param/_cpp_parameters`` and generate the
C++ code and headers using the ``Source/param/mk_params.sh`` script—note
this script needs to be run every time the ``_cpp_parameters``
file is updated.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files ``_parameters`` located in the different
directories that hold the microphysics code. At compile time, a
script in the AMReX bulid system, ``findparams.py``, locates all
of the ``_parameters`` files that are needed for the given choice
of network, integrator, and EOS, and assembles all of the runtime
parameters into a module named ``extern_probin_module`` (using the
``write_probin.py`` script). The parameters are set in the ``&probin`` namelist
in the problem inputs file.

C++ parameter format
--------------------

The C parameters take the form of::

    # comment describing the parameter
    name   type   default   need in Fortran?   ifdef    fortran name    fortran type

Here,

  * ``name`` is the name of the parameter that will be looked for
    in the inputs file.

    The name can actually take the form of ``(a, b)``, where ``a`` is
    the name to be used in the inputs file where the parameter is set
    and ``b`` is the name used within the MAESTROeX C++ class.  It is not
    recommended to name new parameters with this functionality—this
    was implemented for backwards compatibility.


  * ``type`` is one of int, Real, or string

  * ``default`` is the default value of the parameter.

The next columns are optional, but you need to fill in all of the
information up to and including any of the optional columns you need
(e.g., if you are going to provide the fortran name, you also need to
provide ``need in Fortran?`` and ``ifdef``.

  * ``need in Fortran?`` is ``y`` if the runtime parameter should be
    made available in Fortran (through ``meth_params_module``).

  * ``ifdef`` provides the name of a preprocessor name that should
    wrap this parameter definition—it will only be compiled in if that
    name is defined to the preprocessor.

  * ``fortran name`` is the name that the parameter should use in
    Fortran—by default it will be the same as ``name``.

  * ``fortran type`` is the data type of the parameter in Fortran—by
    default it will be the Fortran-equivalent to ``type``.

Finally, any comment (starting with ``#``) immediately before the
parameter definition will be used to generate the documentation
describing the parameters.

Microphysics/extern parameter format
------------------------------------

The microphysics/extern parameter definitions take the form of::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the ``priority`` is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.


Parameters by Namespace
=======================

The following tables list all of the general
MAESTROeX runtime parameters.
These tables are generated automatically using the
``rp.py`` script in ``MAESTROeX/sphinx_docs`` by parsing
the ``MAESTROeX/Source/param/_cpp_parameters`` file. The problem-specific parameters
are not shown here.

.. toctree::

   runtime_parameters
