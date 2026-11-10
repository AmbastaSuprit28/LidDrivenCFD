import numpy as np
from parameters import Nx, Ny

def initialize_fields():

    psi = np.zeros((Nx, Ny))
    omega = np.zeros((Nx, Ny))

    return psi, omega
