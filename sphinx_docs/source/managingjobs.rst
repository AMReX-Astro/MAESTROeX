*************
Managing Jobs
*************

General Info
============

All plotfile directories have a job_info file which lists as
host of parameters about the simulation, including:

-  A descriptive name for the simulation (the job_name runtime
   parameter

-  The number of MPI tasks and OpenMP threads

-  The total amount of CPU-hours used up to this point

-  The data and time of the plotfile creation and the directory it was written to.

-  The build date, machine, and directory

-  The MAESTRO, AMReX, and other relevant git hashes for the source

-  The directories used in the build

-  The compilers used and compilation flags

-  The number of levels and boxes for the grid

-  The properties of the species carried

-  The tolerances for the MG solves

-  Any warnings from the initialization procedure (note: these are not currently stored on restart).

-  The value of all runtime parameters (even those that were not explicitly set in the inputs file), along with an indicator showing if the default value
   was overridden.

This file makes it easy to understand how to recreate the run that
produced the plotfile.

Linux boxes
===========

gfortran
--------

gfortran is probably the best-supported compiler for MAESTRO. Here are some
version-specific notes:

-  *gfortran 4.8.x*: This typically works well, but sometimes we get
   an error allocating memory in cluster_f.f90. This
   is a compiler bug (affecting atleast 4.8.2 and 4.8.3):

   The code runs without any problem if it is compiled with -O2
   -ftree-vectorize -fno-range-check (our default) but with
   cluster_f.f90 compiled with -O2 -ftree-vectorize
   -fno-range-check -fno-tree-pre. The “fno-tree-pre” option
   turns off “ftree-pre” that is turned on by “O2”

   GCC manual says,

       -ftree-pre
       Perform partial redundancy elimination (PRE) on trees. This flag is enabled by default at -O2 and -O3.

   gfortran 4.8.5 appears to work without issues

-  *gfortran 5.1.1*: These compilers have no known issues.

-  *gfortran 5.3.x*: These compilers have no known issues.

-  *gfortran 6.2*: These compilers work without any known issues.

   gfortran 6.2.1 is used for nightly regression testing.

PGI compilers
-------------

The AMReX floating point exception trapping is disabled with PGI
compilers earlier than version 16, due to problems with PGI using the
system header files. From version 16 onward, things should work.

There are no known issues with PGI 16.5 compilers—these are used
for nightly regression testing.

Working at OLCF (ORNL)
======================

Titan Compilers
---------------

The preferred compilers on Titan are the Cray compilers.
Cray 8.4.0 works well on titan/OLCF with MPI/OpenMP.

Monitoring Allocations
----------------------

The showusage and showusage -f commands give an
overview of the usage.

Automatic Restarting and Archiving of Data
------------------------------------------

The submission script titan.run and shell script
process.titan in Util/job_scripts/titan/
are designed to allow you to run MAESTRO with minimal interaction,
while being assured that the data is archived to HPSS on the OLCF
machines.

To use the scripts, first create a directory in HPSS that has the same
name as the directory on lustre you are running in (just the directory
name, not the full path). E.g. if you are running in a directory
call wdconvect_run, then do:

::

    hsi
    mkdir wdconvect_run

(Note: if the hsi command prompts you for your password, you will need
to talk to the OLCF help desk to ask for password-less access to
HPSS).

The script process.titan is called from titan.run and will
run in the background and continually wait until checkpoint or
plotfiles are created (actually, it always leaves the most recent one
alone, since data may still be written to it, so it waits until there
are more than 1 in the directory).

Then the script will use htar to archive the plotfiles and
checkpoints to HPSS. If the htar command was successful, then
the plotfiles are copied into a plotfile/ subdirectory. This is
actually important, since you don’t want to try archiving the data a
second time and overwriting the stored copy, especially if a purge
took place.

Additionally, if the ftime executable is in your path (
ftime.f90 lives in amrex/Tools/Postprocessing/F_src/), then
the script will create a file called ftime.out that lists the
name of the plotfile and the corresponding simulation time.

Finally, right when the job is submitted, the script will tar up all
of the diagnostic files created by diag.f90 and ftime.out
and archive them on HPSS. The .tar file is given a name that
contains the date-string to allow multiple archives to co-exist.

The titan.run submission script has code in it that will look at
the most recently generated checkpoint files, make sure that they were
written out correctly (it looks to see if there is a Header file,
since that is the last thing written), and automatically set the
–restart flag on the run line to restart from the most recent
checkpoint file. This allows you to job-chain a bunch of submission
and have them wait until the previous job finished and then
automatically queue up:

::

    qsub -W depend=afterany:<JOB-ID>  <QSUB SCRIPT>

where is the id number of the job that must complete
before the new submission runs and QSUB SCRIPT is the submission
script (e.g. titan.run). This way you can queue up a bunch of
runs and literally leave things alone and it will restart from the
right place automatically and store the data as it is generated.

When process.titan is running, it creates a lockfile (called
process.pid) that ensures that only one instance of the script
is running at any one time. Sometimes if the machine crashes, the
process.pid file will be left behind, in which case, the script
aborts. Just delete that if you know the script is not running.

The chainsub.sh script can be used to automatically launch a
number of jobs depending on a single, currently queued (or running)
job.

Profiling and Debugging on GPUs
-------------------------------

To get an idea of how code performs on Titan’s GPUs, there are a few tools
available. We’ll overview a few here.

Score-P with CUBE and vampir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Score-P is a profiling and tracing tool that can be used to instrument
code to generate data for other tools to analyze, such as CUBE and
vampir. These tools have been developed to analyze performance of HPC
codes that run on several nodes, not specifically for analyzing GPU usage.
Still, they do support some GPU analysis. In the next section, we’ll discuss
NVIDA’s tools specifically for analyzing GPU usage.
At the time of writing, Score-P usage is fairly well documented on OLCF’s
website here: https://www.olcf.ornl.gov/kb_articles/software-scorep/.
We’ll review the essentials here, but please see the link for more details.

To instrument a code with Score-P you must re-compile. First, your
desired modules will need to be loaded. Please note that *order is
important* — you need to load modules needed for compilation before loading
Score-P. The Score-P module is designed to detect the loaded
environment and will configure itself based on that. These tools have been
tested with the PGI 16.3.0 compilers, and we will use them in our examples.
One possible set of module loads is

.. code:: bash

    $ module load PrgEnv-pgi
    $ module swap pgi/15.7.0 pgi/16.3.0
    $ module load cudatoolkit
    $ module load scorep/3.0-rc1

In the above we’ve loaded version 3.0, release candidate 1, which added some
support for analyzing OpenACC code. The next step is to compile. You
essentially preface your normal compile (and link) line with the Score-P
executable and options. As an example using the Fortran wrapper required on
Titan, we have

.. code:: bash

    $ scorep --cuda --openacc -v ftn gpuprogram.f90

One way to achieve this in MAESTRO is to modify the appropriate make file. For
PGI, this would be $AMREX_HOME/Tools/F_mk/comps/Linux_pgi.mak. If
this proves useful, it may be worth it to build Score-P into the build
infrastructure.

Once compiled, we are ready to generate profiling and tracing data. Among those
that develop these tools, note that they draw a distinction between profiling
and tracing. Profiling generates a timing (or perhaps other metric) summary of
the entire program’s execution while tracing produces a timeline of the
execution. Score-P’s analysis is configured with environment variables.
Some of the key configuration variables used in testing include

.. code:: bash

    export SCOREP_ENABLE_PROFILING=yes
    export SCOREP_ENABLE_TRACING=yes
    export SCOREP_EXPERIMENT_DIRECTORY=profile-trace
    export SCOREP_CUDA_ENABLE=yes,kernel_counter,flushatexit
    export SCOREP_CUDA_BUFFER=200M
    export SCOREP_TOTAL_MEMORY=1G
    export SCOREP_OPENACC_ENABLE=yes

For a full listing and definition of configuration variables, execute

.. code:: bash

    $ scorep-info config-vars --full 

Except for very simple codes, you will never want to enable both tracing and
profiling. The overhead is too high, and the code will likely crash or be
excessively slow. Typically, it’s best to profile first and then trace. The
profiling data can be used to help configure tracing (as we’ll see shortly).

Once the configuration is set, simply run the code as you normally would.
Experience suggests you will need to load the same modules that were loaded for
compilation when executing. If analysis is being done through a batch script,
note that you cannot do a simple module load ... in the script. First you
need to do source /opt/modules/default/init/bash in the script, and then
module loads will work as usual.

After executing, analysis data will be stored in the specified
SCOREP_EXPERIMENT_DIRECTORY. With profiling, you
will see a file like profile.cubex. This can be opened with cube
(module load cube).

As mentioned, the profiling data can be used to get recommended settings for
tracing. Running

.. code:: bash

    $ scorep-score -r profile.cubex

will yield output showing estimated sizes for e.g. SCOREP_TOTAL_MEMORY.
It also list functions that are called many times. If you don’t care about
them and they’re slowing Score-P down (or making an outrageously large output
file), you can configure Score-P to ignore them in its analysis. To filter a set
of functions, you need to provide a filter file, for example

.. code:: bash

    $ export SCOREP_FILTERING_FILE=scorep.filter

where

.. code:: bash

    $ cat scorep.filter
    SCOREP_REGION_NAMES_BEGIN
     EXCLUDE
       matmul_sub
       matvec_sub
    SCOREP_REGION_NAMES_END

This would tell Score-P not to trace the routines matmul_sub and
matvec_sub. See the OLCF KnowledgeBase article and/or Score-P’s
help for more, but this doesn’t seem to be the best-documented aspect of the
program.

Running with tracing enabled will generate a traces.otf2 file that can be
inspected with vampir (module load vampir)

nvprof and nvvp
~~~~~~~~~~~~~~~

NVIDIA provides tools for specifically analyzing how your code utilizes their
GPUs. Score-P is a fully-featured profiler with some CUDA and OpenACC
support. It can be useful for providing context for GPU execution and it allows
you to, for example, see line numbers for OpenACC directives that are executed.
nvprof will only analyze GPU execution, but in exchange you get much more
detail than is available with Score-P. nvvp is NVIDIA’s visual
profiler. It can be used to read data generated by nvprof. Most useful
is the guided analysis it will perform, which analyzes your code’s GPU
performance for bottlenecks and suggests ways to improve performance. Both are
provided when you load the cudatoolkit module.

With nvprof, no instrumentation is necessary. Instead, you compile
normally and then run nvprof on the executable. As before, be sure when
executing to load the modules used at compile-time. Executing nvprof on
Titan’s compute nodes requires some unexpected options having to do with how
aprun and nvprof interact.

To get a basic overview printed to the terminal on Titan’s compute node, execute

.. code:: bash

    $ aprun -b nvprof --profile-child-processes ./gpuprogram.exe arg1 arg2... 

To generate tracing data for nvvp, execute

.. code:: bash

    $ aprun -b nvprof --profile-child-processes -o nvprof.timeline.out%p 
      ./gpuprogram.exe arg1 arg2... 

nvvp can then be used to read nvprof.timeline.out%p, where the
%p will be replaced with the process ID. You *must* include %p in
the output file’s name or the code will crash, even if you’re not running a
multi-process code.

To generate profile-like metric data for nvvp, execute

.. code:: bash

    $ aprun -b nvprof --profile-child-processes --analysis-metrics 
      -o nvprof.metrics.out%p ./gpuprogram.exe arg1 arg2... 

This is the output needed for nvvp’s guided analysis.

Target Metrics
~~~~~~~~~~~~~~

The output from profilers may be difficult to makes sense of. The purpose of
this section is to note different metrics and reasonable targets for them.
Note that these may be specific to the k20x hardware in Titan.

-  Threads per block: 256-512. Note that if your code requires many
   registers per thread, then this will limit the number of threads per block.

-  Occupancy: 60% is a reasonable target. We have had success with codes
   even achieving only 23% occupancy.

One very useful tool for determining target metrics and what is limiting your
performance is a spreadsheet developed by NVIDIA to calculate occupancy. Every
installation of the CUDA Toolkit should have this occupancy calculator in a
tools subdirectory. At time of writing, the calculator is also available at
this link:
http://developer.download.nvidia.com/compute/cuda/CUDA_Occupancy_calculator.xls.
The document is actually more than a simple calculator. It contains quite a
bit of interesting insight into optimizing a GPU code. More on occupancy can be
found here:
http://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html#occupancy.

Batch Submission of yt Visualization Scripts
--------------------------------------------

Rhea—preferred method
~~~~~~~~~~~~~~~~~~~~~

 this section needs to be updated. See the titan section

The best way to do visualization is to use rhea, the OLCF vis machine.
You need to build yt via the install_script.sh script *on
rhea*. It also must be on a *Lustre filesystem*, so it is seen by
the compute node. It is best to build it in your $PROJWORK directory,
since that has a longer time between purges. Once the installation is
complete, it will tell you what script to source to define the
necessary paths.

The scripts in MAESTRO/Util/job_scripts/rhea/ will handle the
visualization. On rhea, the job script gives you access to the
compute node, and then you can run serial jobs or a parallel job with
mpirun. The process-rhea.run script will request the
resources and run the parallel-yt-rhea script.
parallel-yt-rhea will launch the visualization process (defined
via the variables at the top of the script) on all the plotfiles
found to match the prefix defined in the script. Several serial
jobs are run together, with the script using lock files to keep track
of how many processors are in use. When processors free, the next
file in the list is processed, and so on, until there are no more
plotfiles left to process. If a .png file matching the
plotfile prefix is found, then that plotfile is skipped.

Note: the line in parallel-yt-rhea that sources the yt activate script may need to be modified to point to the
correct yt installation path.

Titan
~~~~~

You can run yt python scripts in the titan batch queues to do your
visualization. You need to install yt and all its dependencies
manually somewhere on a *Lustre filesystem*—this ensures that the
compute nodes can see it. A good choice is the project space, since
that has a longer purge window. The following procedure will setup
the development version of yt (from source)

-  create a directory in your $PROJWORK directory named yt/

-  in yt/, down load the yt install script:

   ::

       wget https://bitbucket.org/yt_analysis/yt/raw/yt/doc/install_script.sh

-  edit the script to use Conda to get the necessary dependencies
   and to build yt from source. This is accomplished by setting the
   following variables near the top of the script:

   ::

       INST_CONDA=1
       INST_YT_SOURCE=1

-  run the script:

   ::

       source install_script.sh

When the script is done, you will have a new python installation in a sub-directory
called yt-conda/ and the script will tell you how to modify your path
in your .bashrc

**Important: make sure that you are not loading any other
python environments in your .bashrc, e.g., via modules.**

To test thing out, start up the python shell, and try doing
import yt. If there are no errors, then you are good.

The python code vol.py and submission script
yt-vis-titan.run in MAESTRO/Util/job_scripts/titan/vis/ show
how to do a simple visualization in the batch queue using yt. Note
that vol.py is executable, and that we run it via aprun to
ensure that it is launched on the compute node.

The scripts vis-titan.run and parallel-yt-new in that
same directorywill manage the yt jobs by calling a
python script for each file that matches a pattern.
Note that the actual visualization command itself is launched by
parallel-yt-new, again using the aprun command. But
aprun can only launch a single job at a time, so this means
we cannot easily do (trivally) parallel visualization on a node. For
this reason, running on rhea is preferred.

Remote VisIt Visualization on Lens
----------------------------------

*Note: this information may be out-of-date. It is recommended that
yt be used instead.*

For large data files, visualization with VisIt should be done with
a client running on your local machine interfacing with VisIt running
on the remote machine. For the lens machine at NCCS, the proper setup
is described below.

First, on lens, in your .bashrc, add:

::

    export PATH="/sw/analysis-x64/visit/bin/":$PATH

(you would think that you could just add module load visit but this
does not seem to work with VisIt.

On your local machine, launch VisIt. Note: this procedure seems to
work with VisIt 2.4.2, but not VisIt 2.5.0 for some reason.

-  First we setup a new host

   -  From the menu, select *options :math:`\rightarrow` host profiles*

   -  Create a new host by clicking on the *New Host* button.

   -  Enter the *Host nickname* as lens

   -  Enter the *Remote host name* as lens.ccs.ornl.gov

   -  Enter the *Path to Visit installation* as /sw/analysis-x64/visit (not sure if this is needed)

   -  Make sure that your *username* is correct

   -  Check *Tunnel data connections through SSH*

-  Setup the *Launch Profiles*

   -  Click on the *Launch Profiles* tab

   -  Click on the *New Profile* button

   -  Enter the *Profile name* as parallel

   -  Click on the *Parallel* tab

   -  Check *Launch parallel engine*

   -  Select the *Parallel launch method* as qsub/mpirun

   -  Set the *Partition / Pool / Queue* to computation

   -  Change the *Default number of processors* to 8

   -  Set the *Default number of nodes* to 2

   -  Set the *Default Bank / Account* to AST006

   -  Set the *Default Time Limit* to 00:30:00

-  Click on *Apply* and *Post*

-  Save your changes by selecting *Options :math:`\rightarrow` Save Settings*

To do remote visualization, select *File :math:`\rightarrow` Open*.
From the drop down list at the top, select lens. You will be
prompted for your password. After that, you can navigate to the
directory on lens with the data.

To make a movie (output sequence of images):

-  save a view in VisIt you like as a session file (File :math:`\rightarrow` Save session).

-  On lens, create a file called files.visit which lists all
   of the files you want to visualize, one per line, with /Header
   after the filename. This can be done simply as:

   ::

         ls -1 | grep -v processed | awk '{print $1"/Header"}' > files.visit
         

   (note: the processed bit is for when you used the script above to
   automatically archive the data).

-  Edit the session file, searching for the name of the plotfile you
   originally visualized, and change it to read files.visit. Make
   sure that the path is correct. This may appear multiple times.

-  Relaunch VisIt locally and restore the session (File :math:`\rightarrow` Restore session). It will render the first image. Then reopen (File :math:`\rightarrow` ReOpen file). After this is done, the buttons that allow you to move through the files should become active (black).

-  Resave the session file

-  To generate the frames, you have 2 options:

   #. File :math:`\rightarrow` Save movie. Pick *New simple movie*,
      then set the format to *PNG* and add this to the output box by
      clicking the right arrow, then in the very last screen, select:
      Later, tell me the command to run.

      VisIt will pop up a box showing the command to run. Trying to
      get the currently running session of VisIt to generate the frames
      seems problamatic. Note: you will probably want to edit out the
      -v x.x.x argument in the commandline to not have it force
      to use a specific version.

   #. If the session file successfully does the remote visualization
      as desired, you can run the movie via the commandline with something like:

      ::

             visit -movie -format png -geometry 1080x1080 -output subchandra_cutoff3_ \
                 -start 0 -end 44 -sessionfile subchandra_radvel.session
             

Working at NERSC
================

edison compilers
----------------

The default compilers on edison are the Intel compilers, but
PGI and Cray also work well

-  Intel 15.0.1 works well on edison/NERSC with MPI/OpenMP

-  Intel 16.0.2 works fine.

-  Cray 8.4.x has worked in the past, but it has not been used
   at NERSC in a while.

Note: in order to compile, you will need to ensure that both the
python and python_base modules are loaded (via the
module command).

Running Jobs
------------

edison is configured with 24 cores per node split between two Intel
IvyBridge 12-core processors. Each processor connects to 1/2 of the
node’s memory and is called a NUMA node, so there are 2 NUMA nodes per
edison node. Best performance is seen when running with 6 or 12 threads.

Note: edison switched to SLURM as the batch system. Your job is submitted
using the sbatch command. Options to sbatch are specified at the
top of your submission script with #SBATCH as the prefix. These options
can be found on the sbatch manpage. For instance,

::

    #SBATCH -N 2
    #SBATCH -J myjob
    #SBATCH -A repo-name
    #SBATCH -p regular
    #SBATCH -t 12:00:00

will request 2 nodes (-N), under the account repo-name (-J),
in the regular queue, and for a 12-hour window -t.

If you are using OpenMP, then your script should set OMP_NUM_THREADS, e.g.,

::

    export OMP_NUM_THREADS=12

By default, SLURM will change directory into the submission directory. The
job is launched from your script using srun, e.g.:

::

    srun -n 48 ./main.Linux.Cray.mpi.exe inputs_3d

to run 48 MPI tasks (across the 2 nodes), or

::

    export OMP_NUM_THREADS=6
    srun -n 8 -c 6 ./main.Linux.Cray.mpi.omp.exe inputs_3d

to use 8 MPI tasks each with 6 threads.

The scripts in Util/job_scripts/edison/ provides some examples.

To chain jobs, such that one queues up after the previous job finished,
use the chainslurm.sh script in that same directory. You can view
the job dependency using:

::

    squeue -l -j job-id

where job-id is the number of the job.

Jobs are submitted with sbatch. A job can be canceled using
scancel, and the status can be checked using squeue -u
*username*.

.. _automatic-restarting-and-archiving-of-data-1:

Automatic Restarting and Archiving of Data
------------------------------------------

The same set of submission scripts described for titan are available
for edison at NERSC in Util/job_scripts/edison/. In particular,
the job submission script will set the restart command line parameters
to restart from the most recent checkpoint file in the output directory.

Note: NERSC does not allow for the process script that archives
to HPSS to run in the main job submission script. Instead, a separate
job needs to be run in the “xfer” queue. The script edison.xfer.slurm
in Util/job_scripts/edison/ shows how this works.

Jobs in the xfer queue start up quickly. The best approach is
to start one as you start your main job (or make it dependent on the
main job). The sample process.xrb script will wait for output
and then archive it as it is produced, using the techniques described
for titan above.

To check the status of a job in the xfer queue, use:

::

    squeue -u username -M all

Batch visualization using yt
----------------------------

yt can be built using the install_script.sh. It has been
tested using the build of yt from source and dependencies via conda,
by setting:

::

    INST_CONDA=1
    INST_YT_SOURCE=1

in the install_script.sh. Once these are set, run:

::

    source install_script.sh

Note: installation was done in the home directory.

This way of building yt installs it’s own python and support
libraries in a directory, yt-conda. **Important:** You need
to make sure that your start-up files (typically .bashrc.ext at
NERSC) don’t module load python or any python libraries, as this
will interfere with the conda installation. The install script will
direct you to add the install location to your path.

The scripts parallel-yt and process-edison.slurm in
Util/job_scripts/edison show how to invoke yt to loop over a
series of plotfiles and do visualization. A number of tasks are run
at once on the node, each operating on a separate file. The
parallel-yt script then calls vol.py to do the volume
rendering with yt. Note: it is important that srun be used to
launch the yt script to ensure that it is run on the compute node.

A simple test-yt.slurm script shows how to just call the
yt python script directly, using one node and 24 threads, again
using srun to execute on the compute node.

If you want to keep up with the development version of yt, then you
can update the source in yt-conda/bin/src/yt-hg, using:

::

    hg pull
    hg update yt

and then rebuild it via:

::

    python setup.py develop --user

Using the AmrPostprocesing python plotting scripts on hopper
------------------------------------------------------------

To build the fsnapshot.so library, you need to do:

::

    module load gcc

f2py is already in the path, so the library should then build without issue.

Then edit your .bashrc.ext file to set the PYTHONPATH to
the python_plotfile directory, e.g.:

::

    export PYTHONPATH="/global/homes/z/zingale/AmrPostprocessing/python"

and set the PATH to that directory,

::

    export PATH="/global/homes/z/zingale/AmrPostprocessing/python:$PATH"

To run the script, you need to do:

::

    module load matplotlib
    module load python

Remote visualization on hopper
------------------------------

VisIt is already configured to work with hopper. If the host does not appear
in your local version of visit, copy the host_nersc_hopper.xml file
from the .visit/allhosts/ directory under the system’s VisIt install path
to your :math:`\mathtt{\sim}`/.visit/hosts/ directory.

Working at NCSA (Blue Waters)
=============================

Overview
--------

Blue Waters consists of 22,640 Cray XE6 compute nodes and 4,224
Cray XK7 compute nodes.

Each XE node has two AMD Interlagos model 6276 compute units, each of
which has 16 integer cores (thus, a single node has a total of 32 integer
cores). Two integer cores share a multithreaded, 256-bit wide floating
point unit (FPU). If both integer cores have their own thread, each has access
to 128-bit floating point processing, whereas if only one thread is
assigned the process can access all 256 bits. In one major science
application on Blue Waters it was found that having an OpenMP thread for
each integer core gave the best performance, but when starting a new
application it’s best to experiment. One OpenMP thread per FPU may
be better in some cases.

Each compute unit is divided into two NUMA nodes. Cores in
the same NUMA region share a pool of L3 cache. For the same science
application as before it was found that the best performance was achieved
by assigning an MPI task to each NUMA node. Thus, each physical node
has four MPI tasks.

The XK nodes consist of one AMD Interlagos model 6276 compute unit
and an NVIDIA GK110 “Kepler” GPU accelerator (Tesla K20X). The
GPU is configured with 14 streaming multiprocessor units (SMXs), each
of which has 192 single-precision or 64 double-precision CUDA cores. Thus
there are a total of 2688 SP CUDA cores or 896 DP CUDA cores.

For more details, please see
https://bluewaters.ncsa.illinois.edu/user-guide

BW Compilers
------------

The Cray compilers are the default on blue waters, and version
8.3.3 works well with MAESTRO.

.. _monitoring-allocations-1:

Monitoring Allocations
----------------------

The usage command will list the current user’s usage and
usage -P *project* will
list the usage for all users in a project allocation named “project”.
