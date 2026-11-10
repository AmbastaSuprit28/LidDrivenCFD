import numpy as np
from parameters import L, Nx, Ny

x = np.linspace(0, L, Nx)
y = np.linspace(0, L, Ny)

X, Y = np.meshgrid(x, y, indexing="ij")

dx = x[1] - x[0]
dy = y[1] - y[0]
