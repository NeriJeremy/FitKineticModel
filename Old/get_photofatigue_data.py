import numpy as np
from get_photofatigue_model import get_photofatigue_model

"""
This part of the code is calculating the different rates for the following model:
    DARK => ON <=> OFF
"""

def get_photofatigue_data(start, Switch_params):
    
    try:
        
        #Calculate lasers_i for each laser
        lasers_i = []
        
        # Iterate over each laser dictionary in lasers list by index
        for idx, laser in enumerate(Switch_params['lasers']):
            if isinstance(laser['l'], (int, float)) and isinstance(laser['pd'], (int, float)):
                lasers_i.append((laser['pd'] * laser['l']) / 1.9846e-16)
        lasers_i = np.array(lasers_i)
        
        #Create an array that takes into account the number of states and the number of lasers
        eps = np.zeros((Switch_params['n_states'], Switch_params['n_lasers']))
        eps[0,:] = [Switch_params['EpsOn_actOn'], Switch_params['EpsOn_actOff'], Switch_params['EpsOn_RO']] #Attribute Eps to On state
        eps[1,:] = [Switch_params['EpsOff_actOn'], Switch_params['EpsOff_actOff'], Switch_params['EpsOff_RO']] #Attribute Eps to Off state 
        eps[2,:] = [start['EpsDark_actOn'], start['EpsDark_actOff'], start['EpsDark_RO']] #Attribute Eps to Dark state 
        
        #Calculate the excitation rate for each light sensisive state
        light_sensitive_states = [0, 1, 2]  # Choose the light sensistive state in eps array
        
        k_ex = np.zeros((Switch_params['n_states'], Switch_params['n_lasers']))
        
        # Loop through the light-sensitive states and compute excitation rates
        for state in light_sensitive_states:
            k_ex[state, :] = lasers_i * 3.82e-21 * eps[state, :]  
        
        """
        Considering adding the quadratic excitation rates
        """
        #Create an array with all the Thermal relaxation rates
        K_th = np.zeros((Switch_params['n_states'], Switch_params['n_states']))
        #Asign ThermalK 
        K_th[0,1] = start['Th_OnOff']
        K_th[2,1] = start['Th_DarkOff']
        
        #Create an array with all the Qy
        QuantumY = np.zeros((Switch_params['n_states'], Switch_params['n_states']))        
        #Asign Qy 
        QuantumY[0,1] = start['qOnOff']
        QuantumY[1,0] = start['qOffOn']
        QuantumY[2,0] = start['qDarkOn']
        
        #Compute the light induced rates
        K_l = np.zeros((Switch_params['n_states'], Switch_params['n_states'], Switch_params['n_lasers'])) #(Nb arrays, rows, columns)
        # Iterate over states and lasers to compute k_l
        #For each state, create a 2D array, multiply the switching Qy by the corresponding k_ex (here row 0 = k_ex for Onstate laser, row 1 = k_ex for Offstate laser...)
        for i in range(Switch_params['n_states']):  # starting state
        
            for j in range(Switch_params['n_states']):  # ending state
            
                for k in range(Switch_params['n_lasers']):  # iterate over lasers
                
                    K_l[i, j, k] = k_ex[i, k] * QuantumY[i, j]
        
        #Compute fluo rate
        K_fluo = np.zeros(Switch_params['n_lasers'])
        
        # Iterate over all lasers (phases) to compute k_fluo
        for k in range(Switch_params['n_lasers']):  # go over all phases
        
            K_fluo[k] = k_ex[0, k] * Switch_params['Q_Fluo']
        
        #Store the calculated rate in a dictionnary
        rate_array = {
            'K_th' : K_th,
            'K_l' : K_l,
            'K_fluo' : K_fluo
            }
        
        FL, S_RECORDED, S, T_DET, used_dt = get_photofatigue_model(rate_array, Switch_params, start)
        
        # Find indices where T_DET is within the fitting range
        T_select = np.where((T_DET >= Switch_params['fit_range'][0]) & (T_DET <= Switch_params['fit_range'][1]))[0]
        
        # Select only the values within the fit range
        T_DET_IN_FIT_RANGE = T_DET[T_select]
        FL_IN_FIT_RANGE = FL[T_select]
        S_RECORDED_IN_FIT_RANGE = S_RECORDED[:, T_select]
        
        return FL, S_RECORDED, S, T_DET, FL_IN_FIT_RANGE, S_RECORDED_IN_FIT_RANGE, T_DET_IN_FIT_RANGE, used_dt
    
    except Exception as e:
        print(f"An error occurred in get_photofatigue_data function: {e}")
        return None
    