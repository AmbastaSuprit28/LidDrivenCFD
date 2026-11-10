from geometry import X, Y, dx, dy
from initialization import initialize_fields
import numpy as np
from solver import solve_cavity
from visualization import plot_results



# Initialise fields
psi, omega = initialize_fields()

# Solve cavity flow
psi, omega, u, v = solve_cavity(
    psi, omega, dx, dy
)
np.savez(
    "cavity_Re100.npz",
    x=X,
    y=Y,
    u=u,
    v=v,
    psi=psi,
    omega=omega
)


# Visualise results
plot_results(X, Y, psi, omega, u, v)
