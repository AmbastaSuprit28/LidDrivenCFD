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

## Governing Equations

For a two-dimensional incompressible flow, the velocity components are defined in terms of the streamfunction $\psi$:

$$
u = \frac{\partial \psi}{\partial y},
\qquad
v = -\frac{\partial \psi}{\partial x}
$$

The scalar vorticity is:

$$
\omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
= -\nabla^2 \psi
$$

Therefore, the streamfunction satisfies the Poisson equation:

$$
\nabla^2 \psi = -\omega
$$

The nondimensional vorticity transport equation is:

$$
\frac{\partial \omega}{\partial t}
+ u\frac{\partial \omega}{\partial x}
+ v\frac{\partial \omega}{\partial y}
= \frac{1}{Re}\nabla^2\omega
$$

where:

$$
Re = \frac{U_{\mathrm{lid}}L}{\nu}
$$

is the Reynolds number, $U_{\mathrm{lid}}$ is the lid velocity, $L$ is the cavity length, and $\nu$ is the kinematic viscosity.

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

##Streamfunction Poisson Equation

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
$$
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

The solver has been run for the \(Re=100\) lid-driven cavity on 81 X 81 grid.

The final 81 X 81 run reached the configured outer stopping tolerance of 10^-8 at iteration 74,705. The reported outer residual was approximately 10^-8.

### Flow visualisations
<img width="700" height="600" alt="Figure_1" src="https://github.com/user-attachments/assets/79823c29-9e90-4b25-bf77-3160ba81dc53" />
<img width="700" height="600" alt="2" src="https://github.com/user-attachments/assets/a211016a-777d-4a9b-92e8-bb472a9b3984" />
<img width="700" height="600" alt="3" src="https://github.com/user-attachments/assets/3f05025a-c4cf-4453-b2e0-113911e30e19" />
<img width="700" height="600" alt="4" src="https://github.com/user-attachments/assets/f5f9b75e-1f08-4119-9ba5-2bbc52398181" />
<img width="700" height="500" alt="5" src="https://github.com/user-attachments/assets/0b545a89-b222-43a3-ad2c-e99b5ad6b3f4" />
<img width="700" height="500" alt="6" src="https://github.com/user-attachments/assets/1e7a7d5a-d4ee-4bec-9673-1930711e4aee" />
<img width="1189" height="544" alt="ghia et al 1982 comparison" src="https://github.com/user-attachments/assets/4d1c88cf-dae0-48bc-8859-805d3f3522d9" />

## Benchmark validation

The \(Re=100\) solution is compared against the centreline velocity data reported by:

**Ghia, U., Ghia, K. N., and Shin, C. T. (1982).** “High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method.” *Journal of Computational Physics*, 48(3), 387–411.

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
│   ├── ghia_1982.pdf
│   ├── cavity_Re100.npz
│   └──ghia et al 1982 comparison.png
│   
```
