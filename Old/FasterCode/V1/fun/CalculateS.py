from numba import njit
import numpy as np
from fun.get_state_evolution import get_state_evolution

#@njit
def CalculateS(S, K_L_SEQ_SEQ, K_L_SEQ_KL, K_th, dt,
                         n_cycle, n_phases, N_ON_OFF, N_OFF_ON,
                         phase_order, n_points):
    
    Seq = 0
    half_phases = n_phases // 2

    for k in range(n_cycle):
        for j0 in range(N_ON_OFF):
            for i0 in range(half_phases):
                i = phase_order[i0]
                N = n_points[i]
                if N > 0:
                    for j in range(N):
                        if K_L_SEQ_SEQ[i0] == 1:
                            K = K_th + K_L_SEQ_KL[i0, :, :, 0]
                        else:
                            K = K_th + K_L_SEQ_KL[i0, :, :, j]
                        #print(f"[ON-OFF] Cycle={k}, PhaseIdx={i0}, j={j}, Seq={Seq}, sum(K)={np.sum(K):.2e}, sum(S)={np.sum(S[:, Seq]):.2e}", flush=True)
                        S[:, Seq + 1] = S[:, Seq] + dt * (-K.sum(axis=1) * S[:, Seq] + K.T @ S[:, Seq])
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
                        #print(f"[OFF-ON] Cycle={k}, PhaseIdx={i0}, j={j}, Seq={Seq}, sum(K)={np.sum(K):.2e}, sum(S)={np.sum(S[:, Seq]):.2e}", flush=True)
                        S[:, Seq + 1] = S[:, Seq] + dt * (-K.sum(axis=1) * S[:, Seq] + K.T @ S[:, Seq])
                        Seq += 1
    return S