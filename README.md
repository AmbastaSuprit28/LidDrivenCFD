[README.md](https://github.com/user-attachments/files/34645870/README.md)
# LidDrivenCFD

### A from-scratch 2D lid-driven cavity flow solver in Python

A from-scratch Python implementation of the Finite Difference Method, along with the Successive Over-Relaxation (SOR) method, to solve the classic 2D Lid-Driven Cavity problem. The solver employs the vorticity–streamfunction formulation of the incompressible 2D Navier–Stokes equations for low Reynolds numbers. Validated against the Ghia et al. (1982) benchmark at \(Re=100\).

---

## Overview

The lid-driven cavity is a classical benchmark problem in computational fluid dynamics (CFD). A Newtonian fluid occupies a square cavity whose top wall moves at a constant velocity while the remaining walls remain stationary.

This project implements the problem from scratch using finite differences. Pressure is eliminated from the formulation by introducing the streamfunction and vorticity, leaving a Poisson equation for the streamfunction coupled to a transport equation for vorticity.

The implementation is intended as a transparent numerical-methods project: the governing equations, discretisation, iterative Poisson solve, boundary treatment and benchmark comparison are kept explicit.

## Governing equations

For a two-dimensional incompressible flow, the velocity components are defined in terms of the streamfunction \(\psi\):
$$
\[
u = \frac{\partial \psi}{\partial y},
\qquad
v = -\frac{\partial \psi}{\partial x}.
\]
$$
The scalar vorticity is
$$
\[
\omega = \frac{\partial v}{\partial x}
       - \frac{\partial u}{\partial y}
       = -\nabla^2\psi.
\]
$$
Therefore, the streamfunction satisfies the Poisson equation
$$
\[
\nabla^2\psi = -\omega.
\]
$$
The nondimensional vorticity transport equation is
$$
\[
\frac{\partial \omega}{\partial t}
+ u\frac{\partial \omega}{\partial x}
+ v\frac{\partial \omega}{\partial y}
= \frac{1}{Re}\nabla^2\omega,
\]
$$
where
$$
\[
Re = \frac{U_{\mathrm{lid}}L}{\nu}
\]
$$
is the Reynolds number, \(U_{\mathrm{lid}}\) is the lid velocity, \(L\) is the cavity length and \(\nu\) is the kinematic viscosity.

## Numerical methodology

The current implementation uses:

| Component | Method |
|---|---|
| Spatial discretisation | Finite Difference Method |
| Convective terms in the vorticity equation | First-order upwind differences |
| Diffusive terms | Second-order central differences |
| Time integration | Explicit Euler |
| Streamfunction equation | Five-point finite-difference Poisson operator |
| Poisson iteration | Red-black Successive Over-Relaxation (SOR) |
| Wall vorticity | Thom's approximation |

Streamfunction Poisson Equation

For a uniform square grid with spacing $h$, the five-point discretisation of $\nabla^2\psi=-\omega$ gives

$$
\psi_{i,j} =
\frac{
\psi_{i+1,j}+\psi_{i-1,j}
+\psi_{i,j+1}+\psi_{i,j-1}
+h^2\omega_{i,j}
}{4}
$$

SOR relaxes the Gauss–Seidel estimate according to

(1-\alpha)\psi_{i,j}^{,\mathrm{old}}
+\alpha\psi_{i,j}^{,\mathrm{GS}}
$$

where $\alpha$ is the relaxation factor. The implementation updates alternating red and black grid points.

Boundary Conditions

The streamfunction is constant on the cavity walls and is set to zero:

$$
\psi=0 \qquad \text{on all walls}
$$

The top lid moves with the prescribed horizontal velocity $U_{\mathrm{lid}}$; the other three walls are stationary. Wall vorticity is obtained using Thom's approximation.

### Boundary conditions

The streamfunction is constant on the cavity walls and is set to zero:
$$
\[
\psi=0 \qquad \text{on all walls}.
\]
$$
The top lid moves with the prescribed horizontal velocity \(U_{\mathrm{lid}}\); the other three walls are stationary. Wall vorticity is obtained using Thom's approximation.

### Time marching

At each outer time step, the solver applies the boundary conditions, solves the streamfunction Poisson equation, recovers the velocity field, evaluates convection and diffusion, and advances the interior vorticity using explicit Euler time integration.

The outer iteration monitors the maximum change in interior vorticity between successive time steps. The Poisson solver separately monitors its discrete residual.

## Results

The solver has been run for the \(Re=100\) lid-driven cavity on both \(41\times41\) and \(81\times81\) grids.

The final \(81\times81\) run reached the configured outer stopping tolerance of \(10^{-8}\) at iteration 74,705. The reported outer residual was approximately \(1.00\times10^{-8}\).

### Flow visualisations

Add the generated figures to the `results/` directory using these filenames, or update the paths below to match your files.

| Streamfunction contours | Velocity vector field |
|---|---|
| ![Streamfunction contours](results/streamfunction.png) | ![Velocity vector field](results/velocity_vectors.png) |

| Vorticity contours | Velocity magnitude |
|---|---|
| ![Vorticity contours](results/vorticity.png) | ![Velocity magnitude](results/velocity_magnitude.png) |

## Benchmark validation

The \(Re=100\) solution is compared against the centreline velocity data reported by:

**Ghia, U., Ghia, K. N., and Shin, C. T. (1982).** “High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method.” *Journal of Computational Physics*, 48(3), 387–411.

The comparison uses the vertical-centreline horizontal velocity \(u(0.5,y)\) and the horizontal-centreline vertical velocity \(v(x,0.5)\).

![Centreline comparison with Ghia et al. (1982)](validation/results/ghia_81x81_side_by_side.png)

The benchmark comparison is intended to assess whether the computed velocity profiles reproduce the established \(Re=100\) cavity solution. The reference data and comparison script are kept in `validation/ghia_1982/`.

## Repository structure

```text
LidDrivenCFD/
├── README.md
├── requirements.txt
├── src/
│   ├── parameters.py
│   ├── geometry.py
│   ├── initialization.py
│   ├── boundary.py
│   ├── streamfunction.py
│   ├── flow_fields.py
│   ├── solver.py
│   └── visualization.py
├── validation/
│   ├── ghia_1982/
│   │   ├── ghia_Re100.csv
│   │   └── compare.py
│   └── results/
├── results/
├── docs/
└── tests/
```

## Getting started

### Requirements

- Python 3
- NumPy
- Matplotlib

Install the dependencies with:

```bash
pip install -r requirements.txt
```

### Running the solver

The run command will be documented here once the project entry point and final module layout are in place.

## Scope and limitations

This repository focuses on a two-dimensional, incompressible, Newtonian lid-driven cavity at low Reynolds number. It is an educational and portfolio implementation of finite-difference CFD methods; it is not intended to be a general-purpose CFD package.

## Reference

Ghia, U., Ghia, K. N., and Shin, C. T. (1982). “High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method.” *Journal of Computational Physics*, 48(3), 387–411.
