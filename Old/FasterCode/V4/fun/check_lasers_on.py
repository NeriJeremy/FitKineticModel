

"""

This part of the code is to check if the laser is on during a certain phase at time t

"""

def check_lasers_on(phi, dc, phase, t, T):
    
    try:

        laser_start_times = T * phi[:, phase -1]
        
        laser_end_times = laser_start_times + T * dc[:, phase -1]
        
        L_ON = (laser_start_times <= t) & (laser_end_times > t)
        
        return L_ON
    
    except Exception as e:
        print(f"An error occurred in check_lasers function: {e}")
        return None