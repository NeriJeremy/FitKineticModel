from DARKONOFFswitch_model import DARKONOFFswitch_model

def DARKONOFFswitch_fit_function(params, x, start, w_variables, ONphasefirst, ONphaselast, OFFphasefirst, OFFphaselast):
                                                                 
    # Update 'start' values with the new parameters 'a'
    a2 = start.copy()# Make a copy of the start array
    a2[w_variables] = params  # Update specified indices with the new parameters
    
    #retain only the second parameters (y)
    y,_,_ = DARKONOFFswitch_model(a2, x, ONphasefirst, ONphaselast, OFFphasefirst, OFFphaselast)
    
    return y