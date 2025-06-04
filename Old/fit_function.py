from get_photofatigue_data import get_photofatigue_data

def fit_function(params, x, start, w_variables, fitting_cycle):
    
    try:
    
        print(f"Fitting cycle: {fitting_cycle}")
        print(f"Fitting parameters: {start}")
    
        # Copy the 'start' dictionary
        a2 = start.copy()
        
        # Get list of field names
        fields_start = list(start.keys())
        
        # Assign new values to the specified fields using w_variables
        for i, idx in enumerate(w_variables):
            field_name = fields_start[idx]
            a2[field_name] = start[i]
    
        # Call your data function
        # Assuming get_photofatigue_data_vsn14 returns a tuple where the fifth element is 'y'
        _, _, _, _, y, _, _, _ = get_photofatigue_data(a2, Switch_params)
        
        # Increment the fitting cycle
        fitting_cycle += 1
        
        return y
    
    except Exception as e:
        print(f"An error occurred in fit_function : {e}")
        return None