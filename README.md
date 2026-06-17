# RAS_MoM_2D_Point_Matching

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.20737063)

A Fortran 2008 solver for two-dimensional electromagnetic scattering problems,
implementing the **Random Auxiliary Sources (RAS)** method and the
**Method of Moments (MoM) with point matching**, including wideband capability
via **Asymptotic Waveform Evaluation (AWE)**.

> If you use this code in your research, please cite the associated publications
> listed in [CITATION.cff](CITATION.cff) and in the [Citation](#citation) section below.

---

## Table of Contents

- [What the Code Does](#what-the-code-does)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Building](#building)
- [Running an Example](#running-an-example)
- [Input File Reference](#input-file-reference)
- [Output Files](#output-files)
- [Supported Problem Types](#supported-problem-types)
- [Wideband Simulation](#wideband-simulation)
- [Multi-Scatterer Problems](#multi-scatterer-problems)
- [Porting to Other Languages](#porting-to-other-languages)
- [License](#license)
- [Citation](#citation)

---

## What the Code Does

This code solves 2D electromagnetic scattering problems in the frequency domain
for infinitely long cylinders (2D cross-sections) illuminated by a plane wave
or a set of line sources.

Two complementary solvers are implemented:

**Random Auxiliary Sources (RAS)** — an iterative method that places fictitious
line sources (electric or magnetic) inside and/or outside the scatterer to
satisfy electromagnetic boundary conditions. Sources are added group-by-group
until the normalised boundary residual falls below a user-defined threshold.

**Method of Moments (MoM) with point matching** — a classical matrix method
that discretises the boundary integral equations at a set of testing points on
the scatterer surface. The resulting linear system is solved with one of several
LAPACK-based solvers.

Both methods can operate independently or as complementary approaches for
cross-validation. The **AWE** extension accelerates wideband simulations by
expressing the frequency response as a Taylor series around a centre frequency,
then reconstructing the full bandwidth via Padé approximants or Chebyshev
expansions, avoiding a full matrix solve at each frequency point.

### Boundary Conditions Supported

| Code | Boundary type |
|------|---------------|
| 1    | PEC — Perfect Electric Conductor |
| -1   | PMC — Perfect Magnetic Conductor |
| ≥ 2  | Dielectric (arbitrary ε_r, μ_r, including lossy and double-negative) |
| ≤ -2 | Impedance Boundary Condition (IBC) with full 2×2 tensor impedance |

### Geometry

Scatterer cross-sections can be defined as:

- **Superquadrics** — parameterised by centre (cx, cy), semi-axes (a, b), and
  shape exponent g (g=2 gives an ellipse; large g approaches a rectangle).
- **Arbitrary polygons** — provided as a text file of (x, y) boundary points.

---

## Repository Structure

```
RAS_MoM_2D_Point_Matching/
│
├── Source files
│   ├── RAS_MoM_2D_Point_Matching.f08   Main solver (scatterer type, RAS & MoM engines,
│   │                                   field evaluation, far-field & RCS output)
│   ├── Operations.f08                  Mathematical library (vector types, Hankel/Bessel
│   │                                   functions, Simpson quadrature, CJY01 routine)
│   ├── Matrices_Storage_Handling.f08   Matrix storage (linked list), LDL*/SVD/LAPACK solvers
│   ├── sim_par.f08                     Global simulation parameters module
│   ├── constants.f08                   Physical constants (π, η₀, c, ε₀, μ₀)
│   ├── Gauss_Reduction.f03             Gaussian elimination utility
│   └── io_matrix.f03                   Matrix read/write utilities
│
├── Required input files (read by the solver at runtime)
│   ├── Parameters.dat                  Main simulation configuration (fully annotated)
│   ├── materials.dat                   Material properties library (dielectric & IBC)
│   ├── Line_Sources_Model.dat          Line source definitions (used when Source_Model=1)
│   ├── box3.dat                        Example geometry: rectangular PEC cross-section
│   ├── superquad_points.dat            Example geometry: superquadric cross-section
│   └── points_file.dat                 Example geometry: arbitrary boundary point cloud
│
├── Reference example results (one subfolder per validated test case)
│   ├── Case_Ellipse_b_a_0p3_30GHz/    PEC ellipse b/a=0.3 at 30 GHz
│   ├── Case_Ellipse_b_a_0p4_30GHz/    PEC ellipse b/a=0.4 at 30 GHz
│   ├── Case_Ellipse_b_a_0p5_30GHz/    PEC ellipse b/a=0.5 at 30 GHz
│   ├── Case_Ellipse_b_a_0p7_30GHz/    PEC ellipse b/a=0.7 at 30 GHz
│   └── Case_Ellipse_b_a_1_30GHz/      PEC ellipse b/a=1.0 (circle) at 30 GHz
│                                       Each subfolder contains:
│                                         error.dat
│                                         Monostatic_RCS_RAS_AWE.dat
│                                         Monostatic_RCS_MoM_AWE_Chebyshev.dat
│                                         MoM_Formulations_BC_error_Wideband_Frequencies.dat
│                                         solution_spectrum.dat
│                                         solution_spectrum_Pade.dat
│
├── LICENSE                             Apache License 2.0
├── CITATION.cff                        Machine-readable citation metadata
└── README.md                           This file
```

> **Note on output files:** The solver writes all result files (`error.dat`,
> `source_positions.dat`, `far_field_comparison_*.dat`, `Monostatic_RCS_*.dat`,
> `solution_spectrum*.dat`, `currents_plotting_*.dat`, etc.) directly into the
> working directory at runtime. These are not included in the repository —
> use the `Case_Ellipse_*/` subfolders as reference outputs to validate your
> build.

---

## Prerequisites

- A Fortran 2008-compliant compiler:
  - **gfortran** ≥ 7.0 (free, recommended)
  - **Intel Fortran (ifort / ifx)** ≥ 2018
- **LAPACK** and **BLAS** libraries (standard on most HPC systems; on Ubuntu/Debian:
  `sudo apt install liblapack-dev libblas-dev`)
- **OpenMP** support (included with gfortran/ifort by default; used for
  optional parallelism in matrix assembly)

---

## Building

No build system is provided; compilation is straightforward with a single
command. Adjust library flags for your system.

### gfortran

```bash
gfortran -O2 -fopenmp \
    constants.f08 \
    Operations.f08 \
    Gauss_Reduction.f03 \
    io_matrix.f03 \
    sim_par.f08 \
    Matrices_Storage_Handling.f08 \
    RAS_MoM_2D_Point_Matching.f08 \
    -o RAS_MoM_solver \
    -llapack -lblas
```

### Intel Fortran

```bash
ifort -O2 -qopenmp \
    constants.f08 \
    Operations.f08 \
    Gauss_Reduction.f03 \
    io_matrix.f03 \
    sim_par.f08 \
    Matrices_Storage_Handling.f08 \
    RAS_MoM_2D_Point_Matching.f08 \
    -o RAS_MoM_solver \
    -mkl
```

> **Note on compilation order:** The modules must be compiled in dependency
> order as shown above. `constants` → `Operations` → utility files →
> `sim_par` → `Matrices_Storage_Handling` → main program.

---

## Running an Example

The repository includes a ready-to-run example of a **PEC box** illuminated by
a plane wave at **30 GHz**.

1. Compile the code (see [Building](#building)).
2. Ensure `Parameters.dat`, `materials.dat`, and `box3.dat` are in the working
   directory (they are already in the repository root).
3. Run:

```bash
./RAS_MoM_solver
```

The solver reads all configuration from `Parameters.dat` and writes results to
the files listed in [Output Files](#output-files).

Reference results for five ellipse test cases (b/a = 0.3, 0.4, 0.5, 0.7, 1.0)
are provided in the `Case_Ellipse_*` subdirectories for validation.

---

## Input File Reference

### `Parameters.dat`

The main configuration file. Every parameter is documented with inline comments
directly in the file. Key parameters:

| Parameter | Description |
|-----------|-------------|
| `Source_Model` | 0 = plane wave, 1 = line sources |
| `theta_i`, `phi_i`, `alpha_i` | Angles of incidence [degrees] |
| `number_of_scatterers` | Number of objects in the scene |
| `frequency` | Operating frequency [Hz] |
| `simulation_mode` | Single run (0), parameter sweep (2–4), or iterative (5) |
| `N_group` | Sources added per RAS iteration |
| `TOL` | Boundary condition residual convergence threshold |
| `MoM_activation_flag` | Enable/disable MoM and select output type |
| `Matrix_Solution_Method` | 1=Cholesky, 2=LAPACK zhesv, 3=SVD, 4=Preconditioned |
| `Wideband_type` | 0=single frequency, 1=AWE wideband sweep |

See `Parameters.dat` for the complete annotated listing.

### `materials.dat`

Defines all dielectric and IBC materials referenced by ID in `Parameters.dat`.

- **Dielectric** entries: `id  eps_r'  eps_r''  mu_r'  mu_r''`
  where primes denote real part and double-primes the imaginary part.
- **IBC** entries: `id  eta_zz_r  eta_zz_i  eta_zt_r  eta_zt_i  ...`
  (full 2×2 impedance tensor, real and imaginary parts).

Reserved IDs: `0` = free space, `1` = PEC, `-1` = PMC.

### Geometry files (e.g. `box3.dat`, `superquad_points.dat`)

Plain-text files with one boundary point per row: `x  y` [metres].
Points should be ordered consistently (counter-clockwise is conventional).

---

## Output Files

The solver writes results to the working directory. File names are fixed:

| File | Contents |
|------|----------|
| `error.dat` | Normalised boundary condition error vs. iteration number |
| `error_observation.dat` | Residual at field observation points |
| `far_field_comparison_MoM.dat` | Bistatic far-field pattern (MoM) [dBsm vs. angle] |
| `far_field_comparison_RAS.dat` | Bistatic far-field pattern (RAS) |
| `Monostatic_RCS_MoM_AWE_Chebyshev.dat` | Monostatic RCS vs. angle (MoM/AWE Chebyshev) |
| `Monostatic_RCS_RAS_AWE.dat` | Monostatic RCS vs. angle (RAS/AWE) |
| `solution_spectrum.dat` | Wideband frequency response (Taylor AWE) |
| `solution_spectrum_Pade.dat` | Wideband frequency response (Padé AWE) |
| `currents_plotting_RID_1.dat` | Surface current distribution (MoM) |
| `currents_plotting_RID_1_RAS.dat` | Surface current distribution (RAS) |
| `source_positions.dat` | Final RAS source positions |
| `Convergence_rate.dat` | Residual error vs. iteration (convergence history) |
| `MoM_Formulations_BC_error_Wideband_Frequencies.dat` | BC error across frequency (wideband) |

---

## Supported Problem Types

### PEC (material ID = 1)

Enforces the tangential electric field boundary condition (n × E = 0). The
RAS solver places electric line sources inside the scatterer; the MoM solver
uses the EFIE or CFIE formulation.

### PMC (material ID = -1)

Enforces the tangential magnetic field boundary condition (n × H = 0).

### Dielectric (material ID ≥ 2)

Enforces continuity of both tangential E and H across the interface. Sources
are placed both inside (for the exterior field) and outside (for the interior
dielectric field). Set material properties in `materials.dat`.

### Impedance Boundary Condition — IBC (material ID ≤ -2)

Enforces the generalised impedance relation E_t = Z · (n × H). The full
2×2 tensor impedance is supported, enabling anisotropic surfaces (e.g.
corrugated or metasurface-coated objects). Set tensor components in
`materials.dat`. Choose `FormulationType` 1–4 in `Parameters.dat`.

The IBC boundary condition definition and MoM matrix formulation follow
Kishk & Kildal (1995) — see [[0]](#associated-publications) below —
with corrections applied to typographical errors in the original paper.

---

## Wideband Simulation

Set `Wideband_type = 1` in `Parameters.dat` to enable AWE. The solver:

1. Computes Taylor series coefficients of the MoM system matrix and
   excitation vector around the centre frequency.
2. Reconstructs the solution across the band [fl_r, fh_r] × f₀ using either
   Padé approximants (`Pade_L` poles) or Chebyshev expansion.
3. Writes wideband RCS and current spectra to the output files listed above.

The AWE approach avoids re-factorising the system matrix at every frequency,
making it significantly faster than a full sweep for smooth, broadband problems.
See [Hassan & Kishk (2019)](#citation) for algorithmic details.

---

## Multi-Scatterer Problems

Set `number_of_scatterers > 1` and provide a geometry definition line for each
scatterer. Two solution strategies are available via `RAS_solution_method`:

- **Method 1 (Iterative IFB):** Each scatterer is solved independently; fields
  scattered from one object are used as excitation for its neighbours in an
  Iterative Farfield Bouncing loop.
- **Method 2 (Merged system):** All scatterers are assembled into a single
  large linear system and solved at once.

See [Moharram & Kishk (2013, APSURSI)](#citation) for the multi-scatterer
domain decomposition theory.

---

## Porting to Other Languages

This code is published to facilitate reuse and the advancement of research.
If you port this implementation to another language (Python, C++, Julia, etc.),
you are free to do so under the terms of the Apache 2.0 License, subject to
the following attribution requirements:

1. **Retain the copyright notice** from the `LICENSE` file in your derivative
   work's documentation or source code.
2. **Cite the original publications** listed in `CITATION.cff` and below.
   AI-assisted refactoring produces a derivative work of this code; the same
   attribution obligation applies.
3. Include a statement such as:
   > *"This implementation is derived from the original Fortran code by
   > Mohamed A. Moharram Hassan, available at [repository URL]."*

---

## License

Copyright 2024 Mohamed A. Moharram Hassan.

Licensed under the **Apache License, Version 2.0**. See [LICENSE](LICENSE)
for the full text.

---

## Citation

If you use this software, please cite it as:

```bibtex
@software{hassan_ras_mom_2d,
  author  = {Hassan, Mohamed A. Moharram},
  title   = {{RAS\_MoM\_2D\_Point\_Matching}},
  year    = {2024},
  url     = {https://github.com/MAlaaeldin2015/RAS_MoM_2D_Point_Matching},
  license = {Apache-2.0}
}
```

### Associated Publications

Please also cite the relevant paper(s) describing the methods you use.
References are listed in chronological order.

The IBC formulation implemented in this code is based on the following
third-party reference. If you use the IBC functionality, please cite it:

---

**[0] ACES Journal 1995 — IBC formulation basis (third-party reference)**
> A. A. Kishk and P.-S. Kildal, "Electromagnetic scattering from two
> dimensional anisotropic impedance objects under oblique plane wave
> incidence," *Applied Computational Electromagnetics Society Journal*,
> vol. 10, no. 3, pp. 81–92, 1995.

---

**[1] ACES 2013 — Original RAS formulation (PEC scatterers)**
> M. A. Moharram and A. A. Kishk, "Electromagnetic scattering from 2D conducting
> objects using equivalent randomly distributed sources," *The Applied Computational
> Electromagnetic Society (ACES) Conference*, 2013.

---

**[2] EuCAP 2013 — Extension to dielectric scatterers**
> M. A. Moharram and A. A. Kishk, "Electromagnetic scattering from 2D dielectric
> objects using randomly distributed sources," *2013 7th European Conference on
> Antennas and Propagation (EuCAP)*, Gothenburg, Sweden, pp. 1546–1549, 2013.

---

**[3] USNC-URSI 2013 — Efficient frequency-domain RAS technique**
> M. A. Moharram and A. A. Kishk, "Efficient frequency domain technique for
> electromagnetic scattering from arbitrary objects using the Random Auxiliary
> Sources," *2013 USNC-URSI Radio Science Meeting (Joint with AP-S Symposium)*,
> Lake Buena Vista, FL, USA, pp. 143–143, 2013.
> DOI: [10.1109/USNC-URSI.2013.6715449](https://doi.org/10.1109/USNC-URSI.2013.6715449)

---

**[4] IEEE APS/URSI 2013 — Multi-scatterer domain decomposition**
> M. A. Moharram and A. A. Kishk, "Electromagnetic domain decomposition for 2D
> multi-scatterer problem using randomly distributed sources," *2013 IEEE Antennas
> and Propagation Society International Symposium (APSURSI)*, Orlando, FL, USA,
> pp. 1310–1311, 2013.
> DOI: [10.1109/APS.2013.6711315](https://doi.org/10.1109/APS.2013.6711315)

---

**[5] IEEE Antennas and Propagation Magazine 2015 — Primary journal paper (RAS)**
> M. A. Moharram and A. A. Kishk, "Electromagnetic Scattering From Two-Dimensional
> Arbitrary Objects Using Random Auxiliary Sources," *IEEE Antennas and Propagation
> Magazine*, vol. 57, no. 1, pp. 204–216, Feb. 2015.
> DOI: [10.1109/MAP.2015.2397112](https://doi.org/10.1109/MAP.2015.2397112)

---

**[6] IEEE ICCEM 2015 — RAS overview: simplicity, speed, and efficiency**
> M. A. Moharram and A. A. Kishk, "The random auxiliary sources method: A simple,
> fast, and efficient electromagnetic scattering computations approach," *2015 IEEE
> International Conference on Computational Electromagnetics (ICCEM)*, Hong Kong,
> China, pp. 5–7, 2015.
> DOI: [10.1109/COMPEM.2015.7052535](https://doi.org/10.1109/COMPEM.2015.7052535)

---

**[7] URSI GASS 2017 — Wideband RAS-AWE for smooth conducting cylinders**
> M. A. Moharram and A. A. Kishk, "Wideband electromagnetic scattering computations
> for smooth conducting 2D cylinders using the RAS-AWE method," *2017 XXXIInd
> General Assembly and Scientific Symposium of the International Union of Radio
> Science (URSI GASS)*, Montreal, QC, Canada, pp. 1–4, 2017.
> DOI: [10.23919/URSIGASS.2017.8105092](https://doi.org/10.23919/URSIGASS.2017.8105092)

---

**[8] IEEE Transactions on Antennas and Propagation 2019 — Primary journal paper (RAS-AWE)**
> M. A. M. Hassan and A. A. Kishk, "A Combined Asymptotic Waveform Evaluation and
> Random Auxiliary Sources Method for Wideband Solutions of General-Purpose EM
> Problems," *IEEE Transactions on Antennas and Propagation*, vol. 67, no. 6,
> pp. 4010–4021, June 2019.
> DOI: [10.1109/TAP.2019.2902665](https://doi.org/10.1109/TAP.2019.2902665)

