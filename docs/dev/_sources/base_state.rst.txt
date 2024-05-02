**********
Base State
**********

``BaseState`` and ``BaseStateArray``
====================================

In AMReX, grid data is stored in FABs. Data in individual cells can then be accessed using an ``Array4``, e.g.

.. code-block:: cpp

   MultiFab mf; // collection of FABs
   const int ncomp = mf.nComp();

   // iterate over FABs
   for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();
       const Array4<Real>& arr = mf.array(mfi);

       AMREX_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
       {
           arr(i,j,k,n) += 1;
       });
   }

In MAESTROeX, we use a similar technique to store and access the base state. Instead of a FAB, the base state is stored in a ``BaseState`` object. The ``BaseState`` class is templated, so the data stored in a ``BaseState`` object can be anything that acts like a scalar number (i.e. has routines defining basic arithmetic and logical operations).

The ``BaseState`` object can have 1, 2 or 3 dimensions, allowing it represent data on multiple levels, a single level, and/or with multiple components.

Instead of accessing data in individual cells using and ``Array4``, for a ``BaseState`` object we use a ``BaseStateArray``. For example

.. code-block:: cpp

   // create a BaseState object containing Reals
   // with 2 levels, 25 cells per level and 3 components
   const int nlev = 2;
   const int ncells = 25;
   const int ncomp = 3;

   BaseState<Real> base_state(nlev, ncells, ncomp);
   BaseStateArray<Real> base_arr = base_state.array();

   for (auto l = 0; l < nlev; ++l) {
       ParallelFor(ncells, [=] AMREX_GPU_DEVICE (int r)
       {
           base_arr(l,r,0) = 0.0;
           base_arr(l,r,1) = 1.0;
           base_arr(l,r,2) = 2.0;
       });
   }
   Gpu::synchronize();

Note here that we **must** call ``Gpu::synchronize()`` after the GPU kernel. For FABS, the ``MFIter`` loop calls ``Gpu::synchronize()`` at the end to ensure that all GPU streams are complete before moving on. However, as we often iterate over the base state outside of an ``MFIter`` loop, it's necessary that we call it explicitly here.


``BaseState`` arithmetic
========================

If performing simple arithmetic with ``BaseState`` s, then this can be done straightforwardly using the built in member functions. ``BaseState`` s can be added, subtracted, multiplied and divided by scalars and element-wise by other ``BaseState`` s.

.. code-block:: cpp

   const int nlevs = 5;
   const int ncells = 20;

   // if we don't define the number of components, it defaults to 1
   BaseState<Real> base_state(nlevs, ncells);
   BaseState<Real> other_state(base_state); // copy constructor performs a deep copy

   base_state.setVal(0.0); // set all cells to 0.0
   other_state.setVal(1.0); // set all cells to 1.0

   base_state += 0.5;

   BaseState<Real> another_state = base_state * (other_state - 2.5);

``BaseStateGeometry``
=====================

In AMReX, information about the grid geometry and routines associated with it are stored in the ``Geometry`` class. In MAESTROeX, information and routines concerning the geometry of the base state are stored in the ``BaseStateGeometry`` class.

Member variables of the ``BaseStateGeometry`` class are normally only ever read by routines of other classes. Consequently, ``BaseState`` objects storing the actual data are private member variables (so can only be accessed directly by ``BaseStateGeometry`` methods), however they have corresponding public ``BaseStateArray`` variables which can be accessed by other classes.

The base state geometry of a problem is initialized in ``Maestro::Setup()`` by ``BaseStateGeometry::Init()``. After it has been initialized, its member variables can be used. For example, a typical iteration over a ``BaseState`` would look like

.. code-block:: cpp

    BaseState<Real> base_data(base_geom.max_radial_level+1, base_geom.nr_fine);
    auto base_arr = base_data.array();

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        ParallelFor(base_geom.nr_fine, [=] AMREX_GPU_DEVICE (int r)
        {
            base_arr(l,r) = 0.0;
        });
        Gpu::synchronize();
    }


For multilevel problems, the iteration over the base state looks a bit different, as on each level we must make sure that we only consider base state cells covered by grids. The discrete blocks of base state cells covered by grids on each level are defined by ``base_geom.numdisjointchunks``, and the start and end points of each block/chunk by the variables ``base_geom.r_start_coord`` and ``base_geom.r_end_coord``.

.. code-block:: cpp

    BaseState<Real> base_data(base_geom.max_radial_level+1, base_geom.nr_fine);
    auto base_arr = base_data.array();

    for (auto l = 0; l <= base_geom.max_radial_level; ++l) {
        for (auto i = 1; i <= base_geom.numdisjointchunks(l); ++i) {
            const auto lo = base_geom.r_start_coord(l,i);
            const auto hi = base_geom.r_end_coord(l,i);

            ParallelFor(hi-lo+1, [=] AMREX_GPU_DEVICE (int j)
            {
                auto r = j + lo;
                base_arr(l,r) = 0.0;
            });
            Gpu::synchronize();
        }
    }