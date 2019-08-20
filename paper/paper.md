---
title: 'MAESTROeX: A Massively Parallel Low Mach Number Astrophysical Solver`
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
date: 20 August 2019
bibliography: paper.bib
---

# Summary

``MAESTROeX`` is a massively parallel C++/F90 solver for low Mach number astrophysical flows.
The underlying low Mach number equation set allows for efficient, long-time integration for highly subsonic flows compared to compressible approaches.
``MAESTROeX`` is suitable for modeling full spherical stars as well as well as planar simulations of dynamics within localized regions of a star, and can robustly handle several orders of magnitude of density and pressure stratification.
The code leverages the new AMReX software framework [@AMReX] for block-structured adaptive mesh refinement (AMR) applications, allowing for scalability to large fractions of leadership-class machines.
Microphysics are provided by the Starkiller-Astro libraries [@starkiller].
The code contains unit tests and nightly regression tests.
Parallelization using hybrid MPI/OpenMP is supported with increasing GPU support.

The numerical algorithm is described in a series of papers (see @maestro; @maestroex; and references therein) and has been used to study a number of problems including convection in Chandrasekhar mass models for type Ia supernovae, convection in massive stars, sub-Chandrasekhar white dwarfs, and type I X-ray bursts.

# Acknowledgements

The work at LBNL was supported by the U.S. Department of Energy's Scientific Discovery Through Advanced Computing (SciDAC) program under contract No. DE-AC02-05CH11231.
The work at Stony Brook was supported by DOE/Office of Nuclear Physics grant DE-FG02-87ER40317 and through the SciDAC program DOE grant DE-SC0017955.

# References
