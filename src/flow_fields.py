import numpy as np

from parameters import U_lid


def calculate_velocity(psi, dx, dy):

    Nx, Ny = psi.shape

    u = np.zeros((Nx, Ny))
    v = np.zeros((Nx, Ny))

    # Central differences
    u[1:-1, 1:-1] = (
        psi[1:-1, 2:] - psi[1:-1, :-2]
    ) / (2.0 * dy)

    v[1:-1, 1:-1] = -(
        psi[2:, 1:-1] - psi[:-2, 1:-1]
    ) / (2.0 * dx)

    # Moving lid (excluding corners)
    u[1:-1, -1] = U_lid

    return u, v
