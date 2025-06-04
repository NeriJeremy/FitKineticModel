import numpy as np


def AssignParams(Switch_params):
    
    try:
    
        lasers = []
        
        # Create n_lasers laser objects with default values
        for _ in range(Switch_params['n_lasers']):
            lasers.append({
                'l': [],          # Wavelength
                'pd': [],         # Power density
                'phi': np.zeros(Switch_params['n_phases']),  # Phase delay 
                'dc': np.zeros(Switch_params['n_phases'])    # Duty cycle
            })       
    
        # Assigning wavelengths to each laser
        lasers[0]['l'] = Switch_params['LambdaOnActinic']
        lasers[1]['l'] = Switch_params['LambdaOffActinic']
        lasers[2]['l'] = Switch_params['LambdaReadout']
        
        # Assigning phase delay to each laser
        # Phase delay for Off switching
        lasers[0]['phi'][0] = Switch_params['t2_actinic_during_ON_OFF_frametime'][0]
        lasers[1]['phi'][0] = Switch_params['t1_actinic_during_ON_OFF_frametime'][0]
        lasers[2]['phi'][0] = Switch_params['t_readout_during_ON_OFF_frametime'][0]
        
        lasers[0]['phi'][1] = Switch_params['t2_actinic_during_ON_OFF_addtime'][0]
        lasers[1]['phi'][1] = Switch_params['t1_actinic_during_ON_OFF_addtime'][0]
        lasers[2]['phi'][1] = Switch_params['t_readout_during_ON_OFF_addtime'][0]
        
        # Phase delay for On switching
        lasers[0]['phi'][2] = Switch_params['t2_actinic_during_OFF_ON_frametime'][0]
        lasers[1]['phi'][2] = Switch_params['t1_actinic_during_OFF_ON_frametime'][0]
        lasers[2]['phi'][2] = Switch_params['t_readout_during_OFF_ON_frametime'][0]
        
        lasers[0]['phi'][3] = Switch_params['t2_actinic_during_OFF_ON_addtime'][0]
        lasers[1]['phi'][3] = Switch_params['t1_actinic_during_OFF_ON_addtime'][0]
        lasers[2]['phi'][3] = Switch_params['t_readout_during_OFF_ON_addtime'][0]
    
        # Assigning duty cycle values to each laser
        # Duty cycles for Off switching
        lasers[0]['dc'][0] = Switch_params['t2_actinic_during_ON_OFF_frametime'][1]
        lasers[1]['dc'][0] = Switch_params['t1_actinic_during_ON_OFF_frametime'][1]
        lasers[2]['dc'][0] = Switch_params['t_readout_during_ON_OFF_frametime'][1]
        
        lasers[0]['dc'][1] = Switch_params['t2_actinic_during_ON_OFF_addtime'][1]
        lasers[1]['dc'][1] = Switch_params['t1_actinic_during_ON_OFF_addtime'][1]
        lasers[2]['dc'][1] = Switch_params['t_readout_during_ON_OFF_addtime'][1]
        
        # Duty cycles for On switching
        lasers[0]['dc'][2] = Switch_params['t2_actinic_during_OFF_ON_frametime'][1]
        lasers[1]['dc'][2] = Switch_params['t1_actinic_during_OFF_ON_frametime'][1]
        lasers[2]['dc'][2] = Switch_params['t_readout_during_OFF_ON_frametime'][1]
        
        lasers[0]['dc'][3] = Switch_params['t2_actinic_during_OFF_ON_addtime'][1]
        lasers[1]['dc'][3] = Switch_params['t1_actinic_during_OFF_ON_addtime'][1]
        lasers[2]['dc'][3] = Switch_params['t_readout_during_OFF_ON_addtime'][1]
    
        # Assigning power densities to each laser
        lasers[0]['pd'] = Switch_params['POnActinic']
        lasers[1]['pd'] = Switch_params['POffActinic']
        lasers[2]['pd'] = Switch_params['PReadout']
        
        #Add the array in the Switch_params dico
        Switch_params['lasers'] = lasers
        
        return Switch_params

    except Exception as e:
        print(f"An error occurred in Assign params function: {e}")
        return None