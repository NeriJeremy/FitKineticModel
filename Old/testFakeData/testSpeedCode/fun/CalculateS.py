from numba import njit
import numpy as np

@njit
def step_evolution(s_vec, k_mat, dt):
    s_vec = np.ascontiguousarray(s_vec) # contiguous array for njit compatibility
    k_mat = np.ascontiguousarray(k_mat)
    return s_vec + dt * (-k_mat.sum(axis=1) * s_vec + k_mat.T @ s_vec)


@njit
def CalculateS(S, K_L_SEQ_SEQ, K_L_SEQ_KL, K_th, dt,
                         n_cycle, n_phases, N_ON_OFF, N_OFF_ON,
                         phase_order, n_points):
    
    Seq = 0 # Index dt
    half_phases = n_phases // 2 # To allocate phases
    
    for k in range(n_cycle): # Go through all cycles
        for j0 in range(N_ON_OFF): # Go through ON_OFF
            for i0 in range(half_phases): # Go through ON_OFF phases only
                i = phase_order[i0] # Set correct phase
                N = n_points[i] # Put in N the number of point of the selected phase
                if N > 0: # Continue till there is no more point in the selected phase
                    for j in range(N):
                        if K_L_SEQ_SEQ[i0] == 1: # If there is a constant rate matrix during the phase (see K_L_SEQ matrix)
                            K = K_th + K_L_SEQ_KL[i0, :, :, 0] # Add thermal relaxation rate to K_L_SEQ
                        else:
                            K = K_th + K_L_SEQ_KL[i0, :, :, j] # Add K_th to each constant rate matrix of the phase
                        S[:, Seq + 1] = step_evolution(S[:, Seq], K, dt) # Apply kinetic scheme
                        Seq += 1

        for j0 in range(N_OFF_ON):
            for i0 in range(half_phases, n_phases):
                i = phase_order[i0]
                N = n_points[i]
                if N > 0:
                    for j in range(N):
                        if K_L_SEQ_SEQ[i0] == 1:
                            K = K_th + K_L_SEQ_KL[i0, :, :, 0]
                        else:
                            K = K_th + K_L_SEQ_KL[i0, :, :, j]
                        S[:, Seq + 1] = step_evolution(S[:, Seq], K, dt)
                        Seq += 1
    return S