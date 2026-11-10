import numpy as np

from parameters import (
    alpha,
    poisson_tolerance,
    max_poisson_iterations
)


def solve_streamfunction(psi, omega, dx, dy):

    Nx, Ny = psi.shape

    # Checkerboard masks for red-black SOR
    I, J = np.indices((Nx - 2, Ny - 2))

    red = (I + J) % 2 == 0
    black = ~red

    denominator = 2.0 / dx**2 + 2.0 / dy**2

    interior = psi[1:-1, 1:-1]

    for iteration in range(max_poisson_iterations):

        # Calculate Gauss-Seidel estimate for red points
        gs = (
            (psi[2:, 1:-1] + psi[:-2, 1:-1]) / dx**2
            +
            (psi[1:-1, 2:] + psi[1:-1, :-2]) / dy**2
            +
            omega[1:-1, 1:-1]
        ) / denominator

        interior[red] += alpha * (
            gs[red] - interior[red]
        )

        # Recalculate after red-point updates
        gs = (
            (psi[2:, 1:-1] + psi[:-2, 1:-1]) / dx**2
            +
            (psi[1:-1, 2:] + psi[1:-1, :-2]) / dy**2
            +
            omega[1:-1, 1:-1]
        ) / denominator

        interior[black] += alpha * (
            gs[black] - interior[black]
        )

        # Check actual Poisson residual every 10 sweeps
        if iteration % 10 == 0:

            residual = (
                (psi[2:, 1:-1] - 2 * psi[1:-1, 1:-1]
                 + psi[:-2, 1:-1]) / dx**2
                +
                (psi[1:-1, 2:] - 2 * psi[1:-1, 1:-1]
                 + psi[1:-1, :-2]) / dy**2
                +
                omega[1:-1, 1:-1]
            )

            error = np.max(np.abs(residual))

            if error < poisson_tolerance:
                return psi

    print(
        f"Warning: Poisson solver reached maximum iterations. "
        f"Residual = {error:.6e}"
    )

    return psi
