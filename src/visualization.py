import numpy as np
import matplotlib.pyplot as plt


def plot_results(X, Y, psi, omega, u, v):

    velocity = np.sqrt(u**2 + v**2)

    # Streamfunction contours
    plt.figure(figsize=(7, 6))
    plt.contour(X, Y, psi, levels=30)
    plt.colorbar(label="Streamfunction")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Streamfunction Contours")
    plt.axis("equal")
    plt.show()

    # Vorticity contours
    plt.figure(figsize=(7, 6))
    plt.contourf(X, Y, omega, levels=40)
    plt.colorbar(label="Vorticity")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Vorticity Contours")
    plt.axis("equal")
    plt.show()

    # Velocity magnitude
    plt.figure(figsize=(7, 6))
    plt.contourf(X, Y, velocity, levels=40)
    plt.colorbar(label="Velocity Magnitude")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Velocity Magnitude")
    plt.axis("equal")
    plt.show()

    # Velocity vectors
    plt.figure(figsize=(7, 6))
    plt.quiver(X[::4, ::4], Y[::4, ::4],
               u[::4, ::4], v[::4, ::4])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Velocity Vector Field")
    plt.axis("equal")
    plt.show()

    # Horizontal centreline velocity
    Nx, Ny = psi.shape
    mid_y = Ny // 2

    plt.figure(figsize=(7, 5))
    plt.plot(X[:, mid_y], u[:, mid_y])
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title("Horizontal Centreline Velocity")
    plt.grid()
    plt.show()

    # Vertical centreline velocity
    mid_x = Nx // 2

    plt.figure(figsize=(7, 5))
    plt.plot(Y[mid_x, :], v[mid_x, :])
    plt.xlabel("y")
    plt.ylabel("v")
    plt.title("Vertical Centreline Velocity")
    plt.grid()
    plt.show()
