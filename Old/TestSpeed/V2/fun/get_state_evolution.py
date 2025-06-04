import numpy as np
from numba import njit

"""

The purpose of this code is to apply the kinetic scheme
The output is the concentration of each final states

"""

# JIT compile the get_state_evolution function
@njit
def get_state_evolution(S, K, dt):
    return S + dt * (-np.sum(K, axis=1) * S + K.T @ S)
