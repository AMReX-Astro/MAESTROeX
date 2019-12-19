---
title: 'MAESTROeX: A Massively Parallel Low Mach Number Astrophysical Solver'
tags:
  - C++
  - Fortran90
  - convection
  - hydrodynamics
  - nuclear reactions
  - nucleosynthesis
  - abundances
  - supernovae
authors:
  - name: Duoming Fan
    orcid: 0000-0002-3246-4315
    affiliation: 1
  - name: Andrew Nonaka
    orcid: 0000-0003-1791-0265
    affiliation: 1
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Donald Willcox
    orcid: 0000-0003-2300-5165
    affiliation: 1
  - name: Alice Harpole
    orcid: 0000-0002-1530-781X
    affiliation: 2
  - name: Michael Zingale
    orcid: 0000-0001-8401-030X
    affiliation: 2
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Department of Physics and Astronomy, Stony Brook University
    index: 2
date: 30 August 2019
bibliography: paper.bib
aas-journal: Astrophysical Journal
aas-doi: 10.3847/1538-4357/ab4f75
---

# Summary
``MAESTROeX`` is a massively parallel, finite-volume C++/F90 solver for low Mach number
astrophysical flows.  Our code utilizes a low Mach number equation set allowing for more
efficient, long-time integration of highly subsonic flows compared to compressible approaches.
The recommended range of applicability is for flows where the Mach number does not exceed ~0.1.

In highly subsonic astrophysical phenomena, sound waves carry sufficiently
little energy that they do not significantly affect the convective dynamics of the system.
In many of these flows, modeling long-time convective dynamics are of interest, and numerical
approaches based on compressible hydrodynamics are intractable, even on modern supercomputers.
One approach to this problem is to use low Mach number models. In our
approach, asymptotic model equations are employed that do not contain sound waves.
Our customized low Mach number model retains
compressibilitiy effects due to, e.g., nuclear energy release, large-scale atmospheric stratification,
compositional changes, and thermal diffusion.
When the Mach number (the ratio of the characteristic fluid
velocity over the characteristic sound speed; `Ma = U/c`)
is small, the resulting system can be numerically integrated with much larger time steps than a
compressible model, i.e., at least a factor of `∼ 1/Ma` larger.
Furthermore, MAESTROeX provides a significant advantage over anelastic approaches in that our model
can robustly handle large perturbations in density and temperature.
For a more detailed comparison to other astrophysical approaches to low Mach number flow, see [@maestroex].

``MAESTROeX`` is suitable for modeling full spherical stars as well as planar simulations
of dynamics within localized regions of a star, and can robustly handle several orders of magnitude
of density and pressure stratification.
The code leverages the new AMReX software framework [@AMReX] for block-structured
adaptive mesh refinement (AMR) applications, allowing for scalability
to large fractions of leadership-class machines.
Our approach to parallization uses a hybrid MPI/OpenMP approach, with increasing GPU support.
Microphysics are provided by the Starkiller-Astro libraries [@starkiller].
The code contains documentation through sphinx, a large suite of test problems, and
nightly regression tests on a number of architectures.

Since the series of papers on its predecessor ``MAESTRO`` have been published (see [@maestro]
and references therein), numerous developments have been made to ``MAESTROeX`` to
reduce algorithm complexity and improve parallel performance.
The current numerical algorithm [@maestroex] couples modules for advection (corner transport upwind),
reactions (VODE), thermal diffusion (linear solvers / multigrid),
pressure-projection approaches (linear solvers / multigrid), and spatial mapping routines
used to define and evolve a one-dimensional hydrostatic base state.
``MAESTROeX`` has been used to study a number of problems including convection
in Chandrasekhar mass models for type Ia supernovae, convection in massive stars,
sub-Chandrasekhar white dwarfs, and type I X-ray bursts.

# Acknowledgements

The work at LBNL was supported by the U.S. Department of Energy's Scientific Discovery Through
Advanced Computing (SciDAC) program under contract No. DE-AC02-05CH11231.
The work at Stony Brook was supported by DOE/Office of Nuclear Physics grant DE-FG02-87ER40317
and through the SciDAC program DOE grant DE-SC0017955.
We also acknowledge contributors to the previous, pure F90 implementation of MAESTRO,
John Bell (LBL), Chris Malone (LANL), and Michael Lijewski (LBL).

# References

