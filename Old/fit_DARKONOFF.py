import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from get_photofatigue_data import get_photofatigue_data
from fit_function import fit_function

def fit_DARKONOFF(xtofit, ytofit, Use_Start_Values, Fit_params, Save_csv, directory, Expname, Save_dir, Fit, check_fit, Switch_params, start):
    try:
        
        if Fit:
            
            #init fitting cycle variable
            fitting_cycle=1
            
            fixed = start[:, 1].astype(bool)    
            #indices the non fixed variables
            w_variables = np.where(~fixed)[0]
            
            InitValues = start[:,0]
            #Get the starting values for the fitting variables
            startvalues = InitValues[w_variables]
                   
        else:
            
            print("""
            Getting results from starting values only
            """)
            
        # Find the indices where xtofit is within the specified range
        w_data = np.where((xtofit > Switch_params['fit_range'][0]) & (xtofit < Switch_params['fit_range'][1]))[0]                    
        #Define data to fit within the specified range
        xdata = xtofit[w_data]
        ydata = ytofit[w_data]
        
        # Get the fitted result with starting values only
        yfit_start_full_range, _, _, xfit_start_full_range, yfit_start, _, xfit_start, _ = get_photofatigue_data(Fit_params, Switch_params)
        
        if check_fit:
            
            # Plot experimental data and fitted results (with starting values only)
            plt.figure()
            plt.plot(xdata, ydata, label='Experimental Data')  
            plt.plot(xfit_start, yfit_start, label='Fitted Data')  
            
            # Compute the correlation between experimental and modeled data
            min_len = min(len(yfit_start), len(ydata)) # get min length of the smallest array
            c_exp_m = np.corrcoef(yfit_start[:min_len], ydata[:min_len])[0, 1]
            
            print(f"""
Final correlation between experimental and modeled data: {c_exp_m}
                  """)
            
            # Titles and labels
            xtitle = '[Frame #]'
            ytitle = 'Fluorescence [AU]'
            plt.xlabel(xtitle, fontsize=10)
            plt.ylabel(ytitle, fontsize=10)
            plt.title('Fit between exp data and non fitted fitting parameters', fontsize=16)

            plt.legend()
            plt.show()
        
            # Input from user
            ok = input("""
Is matching between experimental data and starting model data OK ? (y/n): 
                        """)
                           
            if ok.lower() != 'y':
                print("Exiting the program due to user input.")
                sys.exit()  # Stop the code if user input != y
        
        if Fit:
            print("""
Fitting of the model
                  """)
            
            # Plotting initial data and start fit line
            min_len = min(len(ydata), len(yfit_start))

            plt.plot(ydata[:min_len], linewidth=1.5, label='Original Data')
            plt.plot(range(1, min_len + 1), yfit_start[:min_len], color='red', linewidth=1.5, label='Fit Line')
            plt.legend()
            plt.show()
            
            # Fit the model
            popt, pcov = curve_fit(lambda x, *params: fit_function(params, x, start, w_variables), np.zeros_like(ydata), ydata, p0=startvalues)
            
            # Compute standard error (95% CI approximation assuming normality)
            perr = np.sqrt(np.diag(pcov))
            ci_err = 1.96 * perr  # 95% confidence interval error (assuming normal distributio
            print(ci_err)
            
            results = start
        
        
        
        
        return _,_,_

    except Exception as e:
        print(f"An error occurred in fit_DARKONOFF function: {e}")
   
        