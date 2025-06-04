import numpy as np

"""

The purpose of this code is to apply the kinetic scheme
The output is the concentration of each final states

"""

def get_state_evolution (S, K, dt):
    
    try:
        
        newS = S + dt * (-np.sum(K, axis=1) * S + K.T @ S) 
        
        return newS
        
    except Exception as e:
        print(f"An error occurred in get_state_evolution function: {e}")
        return None
