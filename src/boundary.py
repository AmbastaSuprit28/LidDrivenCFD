from parameters import U_lid


def apply_boundary_conditions(psi, omega, dx, dy):

    # Streamfunction boundary conditions
    psi[0, :] = 0.0
    psi[-1, :] = 0.0
    psi[:, 0] = 0.0
    psi[:, -1] = 0.0

    # Thom's wall vorticity approximation
    omega[0, 1:-1] = (
        -2.0 * psi[1, 1:-1] / dx**2
    )

    omega[-1, 1:-1] = (
        -2.0 * psi[-2, 1:-1] / dx**2
    )

    omega[1:-1, 0] = (
        -2.0 * psi[1:-1, 1] / dy**2
    )

    omega[1:-1, -1] = (
        -2.0 * psi[1:-1, -2] / dy**2
        - 2.0 * U_lid / dy
    )

    return psi, omega
