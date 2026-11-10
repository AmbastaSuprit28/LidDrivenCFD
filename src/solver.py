import numpy as np
import time

from parameters import (
    Re,
    dt,
    max_iterations,
    tolerance
)

from boundary import apply_boundary_conditions
from streamfunction import solve_streamfunction
from flow_fields import calculate_velocity


def solve_cavity(psi, omega, dx, dy):

    start_time = time.perf_counter()

    for iteration in range(max_iterations):

        # Preserve old-time vorticity
        omega_old = omega.copy()

        # Apply boundary conditions
        psi, omega = apply_boundary_conditions(
            psi, omega, dx, dy
        )

        # Solve Poisson equation
        psi = solve_streamfunction(
            psi, omega, dx, dy
        )

        # Calculate velocity
        u, v = calculate_velocity(
            psi, dx, dy
        )

        # Interior vorticity
        w = omega[1:-1, 1:-1]

        # Upwind derivatives in x
        dwdx_backward = (
            omega[1:-1, 1:-1] -
            omega[:-2, 1:-1]
        ) / dx

        dwdx_forward = (
            omega[2:, 1:-1] -
            omega[1:-1, 1:-1]
        ) / dx

        # Upwind derivatives in y
        dwdy_backward = (
            omega[1:-1, 1:-1] -
            omega[1:-1, :-2]
        ) / dy

        dwdy_forward = (
            omega[1:-1, 2:] -
            omega[1:-1, 1:-1]
        ) / dy

        dwdx = np.where(
            u[1:-1, 1:-1] >= 0,
            dwdx_backward,
            dwdx_forward
        )

        dwdy = np.where(
            v[1:-1, 1:-1] >= 0,
            dwdy_backward,
            dwdy_forward
        )

        # Convection
        convection = (
            u[1:-1, 1:-1] * dwdx
            +
            v[1:-1, 1:-1] * dwdy
        )

        # Diffusion
        diffusion = (
            (omega[2:, 1:-1] - 2*w
             + omega[:-2, 1:-1]) / dx**2
            +
            (omega[1:-1, 2:] - 2*w
             + omega[1:-1, :-2]) / dy**2
        )

        # Explicit Euler update
        omega_new = omega.copy()

        omega_new[1:-1, 1:-1] = (
            w + dt * (
                -convection
                + diffusion / Re
            )
        )

        # Replace old-time vorticity
        omega = omega_new

        # Update boundary vorticity
        psi, omega = apply_boundary_conditions(
            psi, omega, dx, dy
        )

        # Outer iteration residual
        error = np.max(
            np.abs(
                omega[1:-1, 1:-1]
                -
                omega_old[1:-1, 1:-1]
            )
        )

        # Progress report
        if iteration % 25 == 0:

            elapsed = time.perf_counter() - start_time

            print(
                f"Iteration: {iteration}, "
                f"Time: {iteration * dt:.3f}, "
                f"Residual: {error:.6e}, "
                f"Elapsed: {elapsed:.2f}s"
            )

              # Steady-state criterion
        if error < tolerance:

            elapsed = time.perf_counter() - start_time

            print(
                f"Steady state reached at iteration {iteration}"
            )
            print(f"Total computation time: {elapsed:.2f}s")

            break

    else:

        elapsed = time.perf_counter() - start_time

        print("Maximum iterations reached.")
        print(f"Total computation time: {elapsed:.2f}s")

    return psi, omega, u, v
