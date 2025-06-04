import numpy as np

"""

d[ON]/d(t)= -(kAB+kOnDark)*[ON]+kBA[OFF]+kDarkOn[DARK]
d[OFF]/d(t)= -kAB[OFF]+kBA[ON]
d[DARK]/d(t)= -kDarkOn[DARK]+kOnDark[ON]

k=(a.E.Qs.P.lambda)

"""

def DARKONOFFswitch_model(a2, x, ONphasefirst, ONphaselast, OFFphasefirst, OFFphaselast):
    """
    For fitting, returns A(x) in a kinetic model C <=> A <=> B
    
    Parameters:
    - a: list or array of parameters, where:
        - a[0] is the total initial concentration of A
        - a[1] is the rate constant kABon
        - a[2] is the rate constant kBAon
        - a[3] is the rate constant kABoff
        - a[4] is the rate constant kBAoff
        - a[5] is the rate constant kOnDarkOn
        - a[6] is the rate constant kOnDarkOff
        - a[7] is the rate constant kDarkOnOn
        - a[8] is the rate constant kDarkOnOff
    
    """
    xsize = len(x)

    #create arrays for each concentration of sa and sb (the size of xlength)
    sa = np.zeros(xsize)
    sb = np.zeros(xsize)
    sc = np.zeros(xsize)
    
    
    sa[0] = a2[0]  # Initial concentration of A
    sb[0] = a2[0]  # Initial concentration of B
    sc[0] = a2[0]  # Initial concentration of C
    
    # Define the time increments (dx)
    dx = x - np.roll(x, 1)  # Shift x values by 1 for dx, mimicking circshift
    dx[0] = 0  # Set dx[0] to 0 because it was shifted
    
    kABon = a2[1]
    kBAon = a2[2]
    kABoff = a2[3]
    kBAoff = a2[4]
    kOnDarkOn = a2[5]
    kOnDarkOff = a2[6]
    kDarkOnOn = a2[7]
    kDarkOnOff = a2[8]

    # Loop to compute values for sa, sb and sc
    for i in range(1, xsize):
        
        if ONphasefirst <= i <= ONphaselast:
            sa[i] = sa[i-1] + (-(kABon + kOnDarkOn) * sa[i-1] + kBAon * sb[i-1] + kDarkOnOn * sc[i-1]) * dx[i]
            sb[i] = sb[i-1] + (- kBAon * sb[i-1] + kABon * sa[i-1] ) * dx[i]
            sc[i] = sc[i-1] + (-kDarkOnOn * sc[i-1] + kOnDarkOn * sa[i-1]) * dx[i]
        
        elif OFFphasefirst <= i <= OFFphaselast:
            sa[i] = sa[i-1] + (-(kABoff + kOnDarkOff) * sa[i-1] + kBAoff * sb[i-1] + kDarkOnOff * sc[i-1]) * dx[i]
            sb[i] = sb[i-1] + (- kBAoff * sb[i-1] + kABoff * sa[i-1] ) * dx[i]
            sc[i] = sc[i-1] + (-kDarkOnOff * sc[i-1] + kOnDarkOff * sa[i-1]) * dx[i]      

        
    return sa, sb, sc
