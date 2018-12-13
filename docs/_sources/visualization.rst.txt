*************
Visualization
*************

Plotfiles
=========

MAESTRO outputs plotfile specifically for visualization and analysis.
Table \ `[vis:table:plotfile] <#vis:table:plotfile>`__ lists the quantities stored in a plotfile.
Not all of these may be present, dependent on what options were used in
creating the plotfile.

By default, plotfiles store double precision data, but if
single_prec_plotfiles = T is set, then the data is
converted to single precision before outputting—this is done to
reduce file sizes.

.. raw:: latex

   \footnotesize

.. table:: Plotfile quantities

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | plotfile variable     | description           | runtime parameter     |
   | name                  |                       | controlling output    |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | :math:`\ut`           | –                     |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | x_vel                 |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | y_vel                 | :math:`\vt`           | –                     |
   +-----------------------+-----------------------+-----------------------+
   | z_vel                 | :math:`\wt`           | 3-d runs only         |
   +-----------------------+-----------------------+-----------------------+
   | density               | :math:`\rho`          | –                     |
   +-----------------------+-----------------------+-----------------------+
   | rhoh                  | :math:`(\rho h)`      | use_tfromp = F or     |
   |                       |                       | (use_tfromp = T and   |
   |                       |                       | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | h                     | :math:`(\rho h)/\rho` | use_tfromp = F or     |
   |                       |                       | (use_tfromp = T and   |
   |                       |                       | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | X(*)                  | :math:`(\rho X_k)/\rh | plot_spec             |
   |                       | o`                    |                       |
   +-----------------------+-----------------------+-----------------------+
   | tracer                | tracers               | plot_trac             |
   +-----------------------+-----------------------+-----------------------+
   | w0_x                  | :math:`w_0 \er \cdot  | plot_base             |
   |                       | \ex`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | w0_y                  | :math:`w_0 \er \cdot  | plot_base             |
   |                       | \ey`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | w0_z                  | :math:`w_0 \er \cdot  | plot_base             |
   |                       | \ez`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | divw0                 | :math:`\nabla \cdot w | plot_base             |
   |                       | _0`                   |                       |
   +-----------------------+-----------------------+-----------------------+
   | rho0                  | :math:`\rho_0`        | plot_base             |
   +-----------------------+-----------------------+-----------------------+
   | rhoh0                 | :math:`(\rho h)_0`    | plot_base             |
   +-----------------------+-----------------------+-----------------------+
   | h0                    | :math:`(\rho h)_0/\rh | plot_base             |
   |                       | o_0`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | p0                    | :math:`p_0`           | plot_base             |
   +-----------------------+-----------------------+-----------------------+
   | radial_velocity       | :math:`\Ubt \cdot \er | spherical == 1        |
   |                       |  + w_0`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | circum_velocity       | :math:`| \Ubt - (\Ubt | spherical == 1        |
   |                       |  \cdot \er) \er |`    |                       |
   +-----------------------+-----------------------+-----------------------+
   | magvel                | :math:`| \Ubt + w_0 \ | –                     |
   |                       | er |`                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | momentum              | :math:`\rho | \Ubt +  | –                     |
   |                       | w_0 \er |`            |                       |
   +-----------------------+-----------------------+-----------------------+
   | vort                  | :math:`| \nabla \time | –                     |
   |                       | s \Ubt |`             |                       |
   +-----------------------+-----------------------+-----------------------+
   | S                     | :math:`S`             | –                     |
   +-----------------------+-----------------------+-----------------------+
   | rhopert               | :math:`\rho - \rho_0` | –                     |
   +-----------------------+-----------------------+-----------------------+
   | rhohpert              | :math:`(\rho h) - (\r | use_tfromp = F or     |
   |                       | ho h)_0`              | (use_tfromp = T and   |
   |                       |                       | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | tfromp                | :math:`T(\rho, p_0, X | –                     |
   |                       | _k)`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | tfromh                | :math:`T(\rho, h, X_k | use_tfromp = F or     |
   |                       | )`                    | (use_tfromp = T and   |
   |                       |                       | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | deltaT                | :math:`[T(\rho, h, X_ | use_tfromp = F or     |
   |                       | k) - T(\rho, p_0, X_k | (use_tfromp = T and   |
   |                       | )]/T(\rho, h, X_k)`   | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | deltap                | :math:`|p(\rho,h,X_k) | use_tfromp = F or     |
   |                       |  - p_0|/p_0`          | (use_tfromp = T and   |
   |                       |                       | plot_h_with_use_tfrom |
   |                       |                       | p                     |
   |                       |                       | = T)                  |
   +-----------------------+-----------------------+-----------------------+
   | tpert                 | :math:`T(\rho,h,X_k)  | –                     |
   |                       | - \overline{T}`       |                       |
   |                       | if use_tfromp = F;    |                       |
   |                       | :math:`T(\rho,p_0,X_k |                       |
   |                       | ) - \overline{T}`     |                       |
   |                       | otherwise             |                       |
   +-----------------------+-----------------------+-----------------------+
   | Machnumber            | :math:`| \Ubt + w_0 \ | –                     |
   |                       | er | / c(\rho,p_0,X_k |                       |
   |                       | )`                    |                       |
   +-----------------------+-----------------------+-----------------------+
   | soundspeed            | :math:`c(\rho,p_0,X_k | plot_cs               |
   |                       | )`                    |                       |
   +-----------------------+-----------------------+-----------------------+
   | deltagamma            | :math:`\Gamma_1(\rho, | –                     |
   |                       | p_0,X_k) - \gammabar` |                       |
   +-----------------------+-----------------------+-----------------------+
   | entropy               | :math:`s(\rho,p_0,X_k | –                     |
   |                       | )`                    |                       |
   +-----------------------+-----------------------+-----------------------+
   | entropypert           | :math:`[s(\rho,p_0,X_ | –                     |
   |                       | k) - \overline{s}]/\o |                       |
   |                       | verline{s}`           |                       |
   +-----------------------+-----------------------+-----------------------+
   | sponge or             | :math:`1/(1 + \Delta  | –                     |
   | sponge_fdamp          | t \kappa f_\mathrm{da |                       |
   |                       | mp})`                 |                       |
   |                       | by default;           |                       |
   |                       | :math:`f_\mathrm{damp |                       |
   |                       | }`                    |                       |
   |                       | if plot_sponge_fdamp  |                       |
   |                       | = T                   |                       |
   +-----------------------+-----------------------+-----------------------+
   | pi                    | :math:`\pi`           | –                     |
   +-----------------------+-----------------------+-----------------------+
   | gpi_x                 | :math:`\nabla \pi \cd | plot_gpi              |
   |                       | ot \ex`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | gpi_y                 | :math:`\nabla \pi \cd | plot_gpi              |
   |                       | ot \ey`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | gpi_z                 | :math:`\nabla \pi \cd | plot_gpi              |
   |                       | ot \ez`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | pioverp0              | :math:`\pi / p_0`     | plot_base             |
   +-----------------------+-----------------------+-----------------------+
   | p0pluspi              | :math:`p_0 + \pi`     | plot_base             |
   +-----------------------+-----------------------+-----------------------+
   | omegadot(*)           | :math:`\dot{\omega}_k | plot_omegadot         |
   |                       | `                     |                       |
   +-----------------------+-----------------------+-----------------------+
   | enucdot               | :math:`(\rho \Hnuc)/\ | plot_Hnuc             |
   |                       | rho`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | Hext                  | :math:`(\rho \Hext)/\ | plot_Hext             |
   |                       | rho`                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | eta_rho               | :math:`\etarho`       | plot_eta              |
   +-----------------------+-----------------------+-----------------------+
   | thermal               | :math:`\nabla \cdot \ | use_thermal_diffusion |
   |                       | kappa \nabla T`       |                       |
   +-----------------------+-----------------------+-----------------------+
   | conductivity          | :math:`\kappa(\rho, T | use_thermal_diffusion |
   |                       | , X_k)`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | ad_excess             | :math:`\nabla - \nabl | plot_ad_excess        |
   |                       | a_\mathrm{ad}`        |                       |
   +-----------------------+-----------------------+-----------------------+
   | particle_count        | number of particles   | use_particles         |
   |                       | in a cell             |                       |
   +-----------------------+-----------------------+-----------------------+
   | processor_number      | processor number      | plot_processors       |
   |                       | containing the cell’s |                       |
   |                       | data                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | pi_divu               | :math:`\pi \nabla \cd | plot_pidivu           |
   |                       | ot \tilde{\Ub}`       |                       |
   |                       | (a measure of energy  |                       |
   |                       | conservation)         |                       |
   +-----------------------+-----------------------+-----------------------+

.. _vis:sec:miniplotfile:

Mini vs. regular plotfiles
--------------------------

MAESTRO can manage two independent sets of plotfiles. This allows you to
output the default plotfiles, which contain a lot of variables, sparsely,
and output mini-plotfiles much more frequently. A mini-plotfile is controlled
by the analogous runtime parameters as the main plotfiles:

-  mini_plot_int is the interval in steps between successive plotfiles

-  mini_plot_deltat is the interval in time between successive plotfiles

-  mini_plot_base_name is the base name that prefixes the plotfiles. The
   default is miniplt

To set the fields that are stored in the mini plotfile, a second set
of runtime parameters is used: mini_plot_var1,
mini_plot_var2, :math:`\ldots`, mini_plot_var9. These can be set to
any of the following:

-  "density"

-  "species": this gets all of the mass fractions

-  the name of an individual species in the network (like "helium-4")

-  "radvel": this gets both the radial and circumferential velocity

-  "velocity": all three componets

-  "temperature": this is either :math:`T(\rho,p_0)` or :math:`T(\rho,h)`, depending
   on the value of use_tfromp

-  "enuc": the nuclear energy generation rate

-  "mach": the Mach number

the AMReX file format
---------------------

MAESTRO stores the plotfile data in a hierarchical directory format,
with each level’s data stored in a separate subdirectory. Some meta-data
is stored in the top-level to help interpret the structure. The basic format is:

-  plt00000/

   -  Header

   -  job_info

   -  Level_00/

      -  Cell_D_00000

      -  …

      -  Cell_H

   -  model_cc_00000

   -  model_ec_00000

In the main directory, Header contains the information required
to interpret the data stored on disk. We describe this below. The
job_info file is a plaintext file that contains a lot of
information about the run (where it was run, when it was run, compiler
options, runtime options, etc.). It is not needed to interpret the
data. Finally, the model_cc_00000 and model_ec_00000
contain the MAESTRO basestate information for the cell-centered and
edge-centered basestate quantities respectively. These files are not
typically used for visualization.

Header
~~~~~~

The main Header is written by fabio_ml_multifab_write_d by
processor 0. The information contained is the following:

    | NavierStokes-V1.1
    | *number of variables*
    | *variable 1 name*
    | *variable 2 name*
    | :math:`\vdots`
    | *last variable name*
    | *number of dimensions*
    | *simulation time*
    | *number of levels* (0-based)
    | *physical domain minimum coordinate* (dm numbers)
    | *physical domain maximum coordinate* (dm numbers)
    | *jump in refinement between level 0 and 1*
    | :math:`\vdots`
    | *jump in refinement between level :math:`n-2` and :math:`n-1`*
    |  

Visualizing with Amrvis
=======================

| Amrvis is a tool developed together with AMReX to visualize datasets
  from codes built around the AMReX library. You can download the
  Amrvis source from:
| https://ccse.lbl.gov/Downloads/downloadAmrvis.html
| Amrvis exists in the C++ AMReX framework, so the build system is
  slightly different. A different executable is needed for 2- vs. 3-d
  datasets. Edit the GNUmakefile and set the compilers (probably
  g++ and gfortran) and the dimensionality, and turn off any
  of the volume rendering options. You will need to have the Motif library
  installed on your system (or a replacement, such as lesstif.

Once the code is built, you visualize a dataset as:

::

    amrvis3d.Linux.g++.gfortran.ex pltfile

where pltfile is the name of the plotfile directory. Different
variables can be selected from the drop down menu at the top. Middle
and right clicking in 3-d select the slice planes, and shift + middle
or right will extract 1-d lines through the data. In 2-d, middle and
right clicking alone extract 1-d lines.

If Amrvis cannot find the Palette file, then the plots will be
in grayscale. To fix this, copy the amrvis.defaults and
Palette files to your home directory and edit amrvis.defaults so that
the palette line points to the Palette file, e.g.:

::

    palette               /home/username/Palette

Visualizing with VisIt
======================

.. _sec:vis:python:

Python visualization scripts
============================

AmrPostprocessing/python provides some simple commandline
tools for doing visualizations of AMReX plotfiles (note: a subset
of these are distributed directly with AMReX in amrex/Tools/Py_util/). The main drivers
are written in python and use a set of Fortran routines, compiled with
f2py to interface with the plotfile data. To use the routine,
you will need to have matplotlib and f2py installed. On a
machine running Fedora linux, you can install these packages via

::

    yum install python-matplotlib f2py

The library required by the python routines can be built by typing
’make’ in that directory. If successful, you should find
a library fsnapshot.so.

The path to fsnapshot.so should be included in your PYTHONPATH
environment variable. This can be done by adding:

::

    export PYTHONPATH="${PYTHONPATH}:/home/user/AmrPostprocessing/python}

to your .bashrc.

It is recommended that you use matplotlib version 1.2.0 or
higher. If the fonts look strange in the output files, you can try
installing the lyx-fonts package and deleting your
.matplotlib directory, and trying again.

plotsinglevar.py
----------------

plotsinglevar.py does visualizations of 2-d AMReX plotfiles,
and slices through 3-d AMReX plotfiles. A simple plot can be made
via:

::

    plotsinglevar.py --log -o test.png plt00000/ tfromp

This will make a plot of “tfromp” from the plotfile plt00000 with log scaling,
and store the output in test.png. See Figure \ `[fig:python] <#fig:python>`__.
If you don’t do ‘-o’, then a default output filename consisting of the
plotfile name + component will be used.

.. raw:: latex

   \centering

.. figure:: \visfigpath/plt00000_tfromp
   :alt: [fig:python] Plot of reacting_bubble done with the python
   script plotsinglevar.py.

   [fig:python] Plot of reacting_bubble done with the python
   script plotsinglevar.py.

If you list 2 different variables after the plotfile name, then they
will be plotted side-by-side in a single figure. For example,

::

    plotsinglevar.py plt00000/ tfromp enucdot

produces the output shown in figure \ `[fig:python_two] <#fig:python_two>`__.

.. raw:: latex

   \centering

.. figure:: \visfigpath/plt00000_tfromp_enucdot
   :alt: [fig:python_two] Plot of reacting_bubble done with the
   python script plotsinglevar.py showing 2 variables plotted
   from a single plotfile.

   [fig:python_two] Plot of reacting_bubble done with the
   python script plotsinglevar.py showing 2 variables plotted
   from a single plotfile.

Additional options include ‘-m’ to specify the minimum data
value, ‘-M’ to specify the maximum data value, and ‘–eps’
to make an EPS plot instead of PNG. Running the script with no parameters
will give the full list available options.

Limited 3-d support is available. When run as with a plotfile name
and variable, it will plot slices (:math:`x`-:math:`y`, :math:`x`-:math:`z`, and :math:`y`-:math:`z`)
through the center of the domain. The option ‘–origin’
will put the slices through the origin.

contourcompare.py
-----------------

contourcompare.py takes two or three plotfiles and a single variable as arguments
and plots contours of the datasets on the same set of axes. This is
form comparisons of different runs. Running the script with no parameters
will give the full list available options.

For example:

::

    contourcompare.py tfromp plt00000 other_plt00000

will make a contour plot of the variable tfromp from the data in
plt00000 and other_plt00000 shown on the same axes.

runtimevis.py
-------------

The runtimevis.py script is designed to be run from a submission
script to produce plots from plotfiles as they are produced. This is
accomplished by hooking it into the process scripts described in
Chapter \ `[ch:managingjobs] <#ch:managingjobs>`__.

The script itself reads in an inputs file, vis.in, that
describes the variables to plot. From 1 to 6 variables can be
plotting from a plotfile. The script does its best to organize them
in columns and rows to maximize the plot area. The image is always
output at 1280\ :math:`\times`\ 720 pixels, corresponding to 720p HD resolution.
For each variable, a block of the form:

::

    [varname]
    min = 1
    max = 2
    log = 1

is supplied. If min or max are omitted, then the data
limits are computed automatically. If log is omitted, then no
log is taken of the data before plotting. The script is then run as:

::

    runtimevis.py plt00000

[sec:vis:yt] Visualizing with yt
================================

ytis a Python package for analyzing and visualizing simulation data,
and understand that AMReX data from MAESTRO and CASTRO (along
with many other simulation codes). For more
information, see the ythomepage at http://yt-project.org/ and
:raw-latex:`\cite{yt}`.

Some sample scripts that use yt with MAESTRO data are contained in
MAESTRO/Util/yt/.

Installing yt
-------------

The easiest way to obtain ytis through the use of an installation script:

.. code:: bash

    > wget http://hg.yt-project.org/yt/raw/yt-3.0/doc/install_script.sh
    > bash install_script.sh

By default, this ytinstall script will download and install, in its
own isolated environment, all the secondary utilities needed to get
ytup and running. Note that this currently includes installing
*hdf5, zlib, bzip2, libpng, freetype, python, numpy, matplotlib,
mercurial, ipython, h5py, Cython, Forthon* as well as a yt
*mercurial* bundle of changes. You can turn off the automatic
installation of any of these particular packages by setting the
appropriate INST\_\* variable to zero in the install script; you
may have to then change some of \*\_DIR variables to point to
your own particular installation of that package. It is usually best
to just let ytinstall its own stuff, which ensures things are
working properly.

After the install script has finished, and assuming you let yt
install its own packages, you’ll need to *prepend* some
environment variables with ytlocations so that your system finds
those first and stops looking. ytprovides a simple script to
do this, which will be announced upon successful completion.

Working with yt
---------------

The ytinstallation provides both an interactive, *iPython*-like,
interface or the ability to import ytmodules for use in a batch
script. The interactive interface should be in your $PATH if
you’ve followed the instructions in the previous section; to start it,
simply type iyt in a terminal.

.. code:: python

    > iyt

    Welcome to yt!


    In [1]: 

Codes like Enzo use what are called *parameter files* to
describe the general information—number of levels, domain
dimensions, time, etc.—of a a dataset. ytlikes to grab this
information before working on any specific data; this is accomplished
via the convenience method load:

.. code:: python

    In [1]: pf = load("plt00166")

Note that some older versions of ytneeded an inputs file
in the same directory as the plotfile, but as of yt3.0, all the
necessary metadata is obtained from the job_info file inside
the plotfile directory.
The load method returns an instance of the StaticOutput class. One
of the easiest ways of handling plots is via a PlotCollection
object

.. code:: python

    In [2]: center = (pf.domain_right_edge + pf.domain_left_edge) / 2.0
    In [3]: pc = PlotCollection(pf,center)

By default, the PlotCollection constructor places the center of
the plot to be at the location of peak density. Here we have
calculated the center of the domain by accessing the lower and upper
domain boundaries via the numpy arrays
pf.domain_left_edge and pf.domain_right_edge,
respectively. Note that up until this point, ythas not actually
loaded the AMR dataset.

Now we would like to take some slices of tfromp in the dataset and
generate some 2-d plots. To do this, we will use the
PlotCollection’s add_slice method:

.. code:: python

    In [4]: pc.add_slice("tfromp",0)
    In [5]: pc.add_slice("tfromp",1)
    In [6]: pc.add_slice("tfromp",2)

The first call to add_slice builds an
AMRHierarchy object associated with pf. The
AMRHierarchy object contains information about the actual dataset,
such as its layout in the simulation domain or on disk. Building the
hierarchy is expensive, but once it is built the data it contains can
be accessed via attributes and dictionary lookup. In other words, the
subsequent add_slice operations are relatively
cheap. The first parameter to the add_slice
method is obviously the variable we want; the second (optional)
parameter specifies an coordinate axis orthogonal to the slice being
made—0 for x, 1 for y, 2 for z. Now we want to save the plots from
the PlotCollection; this is done with the
save method, which takes an optional basename for the generated files:

.. code:: python

    In [7]: pc.save("my_cool_images")

This generates the following image files:

.. code:: python

    Out[7]: 
    ['my_cool_images_Slice_x_tfromp.png',
     'my_cool_images_Slice_y_tfromp.png',
     'my_cool_images_Slice_z_tfromp.png']

Figure `[fig:yt_slice] <#fig:yt_slice>`__ shows an example of one of the slice images.
Note that this was a quick and dirty generation of the image—there
is a lot of space around the figure, which can be removed with various
options to the ytmethods. Also, the user can specify derived
variables, log-scaling, annotations, etc. For more information see the
documentation at http://yt-project.org/docs/dev/  .

[fig:yt_slice]
|Example slice through 3-d dataset with yt.|

When writing a script to use in batch mode, one has to manually import
the import modules needed to work with yt. As such, all scripts must
import from the yt.mods module, which is essentially a
convenience module that sets up the appropriate ytnamespace.
For completeness, below is a script containing our example above.

.. code:: python

    # load our namespace
    from yt.mods import *

    # the plotfile I'm interested in
    fn = "plt00166"

    # load it into a StaticOutput object
    pf = load(fn)

    # find the center of the domain
    center = (pf.domain_right_edge + pf.domain_left_edge) / 2.0

    # associate a PlotCollection with the pf
    pc = PlotCollection(pf,center)

    # add some slices of tfromp
    pc.add_slice("tfromp",0)
    pc.add_slice("tfromp",1)
    pc.add_slice("tfromp",2)

    # save our plots to a files with basename
    # "my_cool_images"
    pc.save("my_cool_images")

2-d datasets
------------

To visualize 2-d data, you can use the SlicePlot function,
picking the normal direction to be "z":

::

    from yt.mods import *
     pf = load("plt00085")
    s = SlicePlot(pf, "z", "tfromp")
    s.save("tfromp.eps")

This generates the figure shown in Figure \ `[fig:yt2d] <#fig:yt2d>`__.

.. raw:: latex

   \centering

.. figure:: \visfigpath/tfromp
   :alt: [fig:yt2d] Example slice through 3-d dataset with yt.

   [fig:yt2d] Example slice through 3-d dataset with yt.

Volume Rendering
----------------

.. |Example slice through 3-d dataset with yt.| image:: \visfigpath/my_cool_images_Slice_x_tfromp_trim

