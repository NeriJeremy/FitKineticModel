import numpy as np
import math
from tqdm import tqdm
from fun.check_lasers_on import check_lasers_on
from fun.CalculateS import CalculateS


"""
This part of the code is generating the kinetic scheme for the following model:
    DARK => ON <=> OFF
"""

def get_photofatigue_model(rate_array, Switch_params, start):
    
    try:
        K_l = rate_array['K_l']
        K_th = rate_array['K_th']
        
        Nb_cycles_Offswitch = Switch_params['T_ON-OFF']/(Switch_params['Frametime']+Switch_params['AddTime'])
        Nb_cycles_Onswitch = Switch_params['T_OFF-ON']/(Switch_params['Frametime']+Switch_params['AddTime'])
        print(f'Number of full switching cycles: {Switch_params['n_cycle']}')
        print(f'Number of full on-off cycles: {Nb_cycles_Offswitch}')
        print(f'Number of full off-on cycles: {Nb_cycles_Onswitch}')

    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during Nb of cycles calculation: {e}")

#######################################################################################################################################
############################ Calculate number of dt and time steps ################################################################################# 
        
    try:
        
        # Create an array where rows are lasers, and columns are phases
        lasers_on = np.array([0 < np.array(laser['dc']) for laser in Switch_params['lasers']])
        
        # Calculate maximum excitation rate over all phases
        max_rate = 0  # Initialize max_rate
        # Check if max_rates has to be an array (one value per phase)
        # Iterate over all phases to calculate max excitation rate only when the laser is on
        for i in range(1, Switch_params['n_phases'] + 1):
            K = K_th.copy()
            # Iterate over lasers
            for j in range(1, Switch_params['n_lasers'] + 1):
                # lasers_on[j] gives the status of the j-th laser for phase i
                laser_status = lasers_on[j -1, i -1]
                
                # Multiply K_l by the laser status for the current phase
                K += K_l[:, :, j-1] * laser_status
    
            # Update max_rate with the maximum value of the current K matrix
            max_rate = max(np.max(K), max_rate)  # Compare max(K) with max_rate and update
        
        #Compute the recommended dt, correspond to dt in [s] make 10 steps for process the highest rate
        recommended_dt = 1/(max_rate*10)
        
        #If the recommended dt is giving less than 'n_min_samples' per frame, increase the recommended dt
        if Switch_params['AddTime'] == 0:
            recommended_dt = min(recommended_dt, Switch_params['Frametime'] /
                                 Switch_params['n_min_samples'])  # get at minimum 'n_min_samples' per frame
            
            # we need that dt be an integer # of the frametime
            recommended_dt = Switch_params['Frametime'] / math.ceil(1e+3 * Switch_params['Frametime'] /
                                                                    (1e+3 * recommended_dt))
        else:
            recommended_dt = min(recommended_dt, Switch_params['Frametime'] / Switch_params['n_min_samples'], Switch_params['AddTime'] / Switch_params['n_min_samples'])  # get at minimum 20 points per frame
            # we need that dt be an integer # of the frametime & addtime
            recommended_dt = 1e-3 * math.gcd(int(1e+3 * Switch_params['Frametime']), int(1e+3 * Switch_params['AddTime'])) / math.ceil(
                math.gcd(int(1e+3 * Switch_params['Frametime']), int(1e+3 * Switch_params['AddTime'])) / (1e+3 * recommended_dt))
            
        print(f'Minimum recommended dt [µs]: {recommended_dt*1e+06}')
        
        if Switch_params['fixed_dt'] == -1:
            dt = recommended_dt
        else:
            if Switch_params['AddTime'] == 0:
                #math.floor round the number to the nearest higher integer
                dt = Switch_params['Frametime'] / math.ceil(1e+3 * Switch_params['Frametime'] / (1e+3 * Switch_params['fixed_dt']))
            else:
                dt = 1e-3 * math.gcd(int(1e+3 * Switch_params['Frametime']), int(1e+3 * Switch_params['AddTime'])) / math.ceil(math.gcd(int(1e+3 * Switch_params['Frametime']), int(1e+3 * Switch_params['AddTime'])) / (1e+3 * Switch_params['fixed_dt']))
    
            if dt > recommended_dt:
                print('****ERROR***** Time step for calculation looks TOO BIG!')

        print(f'Time step for calculation (dt) [µs]: {dt * 1e+06}')
        
        #math.floor round the number to the nearest lower integer
        n_points_per_frametime = math.floor(Switch_params['Frametime'] / dt)
        n_points_per_addtime = math.floor(Switch_params['AddTime'] / dt)
        
        if Switch_params['f_frametime'] == 1:
            print(f'Number of points per EMCCD frame: {n_points_per_frametime}')
        else:
            # We need that f_frametime be an integer multiple of the frametime
            Switch_params['f_frametime'] = math.ceil(n_points_per_frametime * Switch_params['f_frametime']) / n_points_per_frametime
            print(f'Number of points per collected frame: {n_points_per_frametime * Switch_params["f_frametime"]}')
            print(f'Effective frametime set to [s]: {Switch_params["Frametime"] * Switch_params["f_frametime"]}')
            
        # resolution, we might not have exactly correct t_ON_OFF and t_OFF_ON times
        t_ON_OFF_real = dt * math.floor(Switch_params['T_ON-OFF'] / dt)
        t_OFF_ON_real = dt * math.floor(Switch_params['T_OFF-ON'] / dt)
        print(f'Used t_ON_OFF due to integer # of elementary dt  : {t_ON_OFF_real}')
        print(f'Used t_OFF_ON due to integer # of elementary dt  : {t_OFF_ON_real}')
        
        # Calculate the number of frames
        n_frames = math.floor(Switch_params['n_cycle'] * (Switch_params['T_ON-OFF'] + Switch_params['T_OFF-ON']) / (Switch_params['Frametime'] + Switch_params['AddTime']))
        print(f'Number of treated EMCCD frames during experiment: {n_frames}')
        
        # Number of points per phase
        n_points = [0] * Switch_params['n_phases']
        for i in range(0, Switch_params['n_phases'], 2):  # Odd indexed phases get n_points_per_frametime
            n_points[i] = n_points_per_frametime
        for i in range(1, Switch_params['n_phases'], 2):  # Even indexed phases get n_points_per_addtime
            n_points[i] = n_points_per_addtime
        
        # Number of points during frametime + addtime
        N_AF = n_points_per_frametime + n_points_per_addtime
        
        # Number of sub-cycles during OFF_ON
        N_ON_OFF = round(Switch_params['T_ON-OFF'] / (dt * N_AF))
        N_OFF_ON = round(Switch_params['T_OFF-ON'] / (dt * N_AF))
        
        # Total number of points per cycle
        tot_n_points_cycle = (N_ON_OFF + N_OFF_ON) * N_AF

    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during number of dt and time steps calculations: {e}")
        
######################################################################################################################################
#################################### Define concentrations for each states ############################################################
    
    try:
        
        #Define states concentrations over each cycle, rows are for each state, columns are for each point
        S = np.zeros((Switch_params['n_states'], (Switch_params['n_cycle'] * tot_n_points_cycle)+1))
        
        S[0,0] = start['FP_init'] * start['On_init']
        S[1,0] = start['FP_init'] * start['Off_init']
        S[2,0] = start['FP_init'] * start['Dark_init']
        
    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during starting concentrations assignement : {e}")

######################################################################################################################################
########################################################### Define Phases ############################################################

    try:

        #Set phase order
        phase_order = np.arange(1, Switch_params['n_phases'] + 1)
        
        if Switch_params['AddtimeBeforeFrametime']:  # If True, invert order: addtime before frametime
            tmp = phase_order.copy()
            phase_order[0::2] = tmp[1::2]
            phase_order[1::2] = tmp[0::2]
            
    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during phase order adjustment : {e}")

######################################################################################################################################
########################### Calculate the illumination sequence for light induced processes ############################################################    
        
    try:
        
        #Create an array with one row for each laser and one column for each phase with the corresponding phase delay
        lasers_phi = np.array([np.array(laser['phi']) for laser in Switch_params['lasers']])
        #Create an array with one row for each laser and one column for each phase with the corresponding duty cycle
        lasers_dc = np.array([np.array(laser['dc']) for laser in Switch_params['lasers']])
        
        n_cycle = Switch_params['n_cycle']
        n_phases = Switch_params['n_phases']
        n_states = Switch_params['n_states']
        
        # Pre-allocate output arrays
        max_points = max(n_points)
        K_L_SEQ_SEQ = np.zeros(n_phases, dtype=np.int32)
        K_L_SEQ_KL = np.zeros((n_phases, n_states, n_states, max_points))

        # Fill the K_L_SEQ arrays
        for k in range(n_phases):
           i = phase_order[k]
           dc = lasers_dc[:, k]
           N = n_points[i - 1]
           phase_duration = dt * N
       
           w_frametime = next((idx for idx, val in enumerate(dc) if val != 0 and val != 1), None)
       
           if w_frametime is None:
               L_ON = check_lasers_on(lasers_phi, lasers_dc, i, 0, phase_duration)
               K_sum = np.sum(K_l[:, :, L_ON > 0], axis=2)
               for t in range(N):
                   K_L_SEQ_KL[k, :, :, t] = K_sum
               K_L_SEQ_SEQ[k] = 1
       
           else:
               K_L_SEQ_SEQ[k] = 0
               for j in range(N):
                   t = j * dt
                   L_ON = check_lasers_on(lasers_phi, lasers_dc, i, t, phase_duration)
                   K_L_SEQ_KL[k, :, :, j] = np.sum(K_l[:, :, L_ON > 0], axis=2)   

    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during calculation of light induced processes sequences : {e}")    

######################################################################################################################################
########################### Perform the calculation through all cycles #############################################################
    
    try:

        S = CalculateS(S, K_L_SEQ_SEQ, K_L_SEQ_KL, K_th, dt,
                                 n_cycle, n_phases, N_ON_OFF, N_OFF_ON,
                                 phase_order, n_points)
    
    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during calculation through all cycles : {e}")

######################################################################################################################################
########################################## Compute the fluorescent signal #############################################################
    
    try:
        
        k_fluo = rate_array['K_fluo']
        
        # Select the phases where fluorescence is measured (with readout laser)
        read_frames = phase_order[::2] 
       
        # Get fraction of time lasers are on during read_frames
        DC = [laser['dc'] for laser in Switch_params['lasers']]
        DC = np.array(DC).T  # Convert to NumPy array and transpose if needed to match MATLAB's dimension
        DC_read = DC[read_frames-1, :]  # Select rows corresponding to read_frames
        
        F_I = np.zeros(len(read_frames)) #Define new array to store fluo signal
        
        for k in range(len(read_frames)): #Define starting fluo intensity for each readout phase
            F_I[k] = Switch_params['Frametime'] * np.sum(DC_read[k, :] * k_fluo)
        
        # To compute the recorded signal, remove the signal during addtime
        
        #Put in w_frametime the points that are only corresponding to the frametime
        if Switch_params['AddTime'] > 0:
            IND1 = np.arange(1, Switch_params['n_cycle'] * tot_n_points_cycle + 1)  # Create an array the size of the number of cycles
        
            if not Switch_params['AddtimeBeforeFrametime'] :
                # Find indices where mod(IND1, N_AF) is in (0, n_points_per_frametime]
                w_frametime = np.where((IND1 % N_AF <= n_points_per_frametime) & (IND1 % N_AF > 0))[0]
            else:
                # Find indices where mod(IND1, N_AF) > n_points_per_addtime or mod(IND1, N_AF) == 0
                w_frametime = np.where((IND1 % N_AF > n_points_per_addtime) | (IND1 % N_AF == 0))[0]
        
        # Look which points correspond to ON-OFF and OFF-ON cycles
        IND2 = np.arange(1, n_frames + 1)  # Create an array with total number of frames
        
        # Find ON-OFF frame indices: where mod(IND2, N_ON_OFF + N_OFF_ON) is in (0, N_ON_OFF]
        w_on_off = np.where((IND2 % (N_ON_OFF + N_OFF_ON) <= N_ON_OFF) & (IND2 % (N_ON_OFF + N_OFF_ON) > 0))[0]
        
        # Find OFF-ON frame indices: where mod(IND2, N_ON_OFF + N_OFF_ON) > N_ON_OFF or == 0
        w_off_on = np.where((IND2 % (N_ON_OFF + N_OFF_ON) > N_ON_OFF) | (IND2 % (N_ON_OFF + N_OFF_ON) == 0))[0]
        
        # Get emccd_shift as an integer number of dt steps
        used_emccd_shift = round(Switch_params['emccd_shift'] / dt)
        
        # Remove last point from S 
        S = S[:, :-1]
        
        # Use only frametime points if addtime > 0
        if Switch_params['AddTime'] > 0:
            S_FRAMETIME = S[:, w_frametime]
        else:
            S_FRAMETIME = S
        
        # Shift S_FRAMETIME if used_emccd_shift is not zero
        if used_emccd_shift != 0:
            S_FRAMETIME = np.roll(S_FRAMETIME, shift=used_emccd_shift, axis=1)
        
        # Initialize array to store averaged states over each frametime
        S_RECORDED = np.zeros((Switch_params['n_states'], n_frames))
        
        # Take the mean over each frametime
        for k in range(n_frames):
            start_idx = k * n_points_per_frametime
            end_idx = (k + 1) * n_points_per_frametime
            S_RECORDED[:, k] = np.mean(S_FRAMETIME[:, start_idx:end_idx], axis=1)
        
        # This is the fluorescence signal
        FL = S_RECORDED[0, :].copy()  # Copy the first state's time trace
        
        # Scale fluorescence signal for ON-OFF and OFF-ON frames
        FL[w_on_off] *= F_I[0]
        FL[w_off_on] *= F_I[1]
        
        # We remove 1% of the frames and at least one frame because the last frame may not be entirely sampled
        total_duration = Switch_params['n_cycle'] * (Switch_params['T_ON-OFF'] + Switch_params['T_OFF-ON'])

        frame_unit = Switch_params['Frametime'] + Switch_params['AddTime']

        raw_n_frames = total_duration / frame_unit

        n_frames_reduced = int(np.floor(raw_n_frames - max(1, 0.01 * raw_n_frames)))
        
        # Take only the reduced number of frames and apply the frametime fraction scaling
        FL = Switch_params['f_frametime'] * FL[:n_frames_reduced]
        
        FL = start['FP_init'] * start['On_init'] * FL / FL[0]  # Rescale
        
        # Reduce S_RECORDED accordingly
        S_RECORDED = S_RECORDED[:, :n_frames_reduced]
        
        # Get the detection times in seconds at end of EMCCD integration time
        frame_indices = np.arange(1, len(FL) + 1)  # MATLAB 1:numel(FL)
        frame_length = dt * (n_points_per_frametime + n_points_per_addtime)
        
        if not Switch_params['AddtimeBeforeFrametime'] :
            # Invert order: addtime before frametime
            T_DET = frame_length * frame_indices - dt * n_points_per_addtime + used_emccd_shift
        else:
            T_DET = frame_length * frame_indices + used_emccd_shift
        
        used_dt = dt  # Return this value if needed later
        
        return FL, S_RECORDED, S, T_DET, used_dt
        
    except Exception as e:
        print(f"An error occurred in get_photofatigue_model def during computation of fluo signal : {e}")