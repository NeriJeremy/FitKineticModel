import os
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from fun.get_photofatigue_model import get_photofatigue_model

class FitFun:
    
    def __init__(self, Switch_params, directory, x_scale_factor, x_range, Use_Start_Values,
                 Fit_params, Fit, check_fit, start, ShowFullModel):
        
        self.Switch_params = Switch_params
        self.AssignParams() #Add lasers to Switch_params
        self.directory = directory
        self.x_scale_factor = x_scale_factor
        self.x_range = x_range
        self.Use_Start_Values = Use_Start_Values
        self.Fit_params = Fit_params
        self.Fit = Fit
        self.check_fit = check_fit
        self.start = start
        self.fitting_cycle = 0
        self.w_variables = 0
        self.ShowFullModel = ShowFullModel
        
    def Open_csv_fit(self):

        files = os.listdir(self.directory)
        
        # Find the .csv files
        csv_files = [file for file in files if file.endswith('.csv')]
        
        if csv_files:

            if len(csv_files) == 1:   
                filename = csv_files[0]  # First .csv file
                file_path = os.path.join(self.directory, filename)
                
                # Read the .csv file into a DataFrame
                df = pd.read_csv(file_path, sep=',', engine='python')
                
                #Convert the X data in sec and put X and Y data in a variable
                xtofit = self.x_scale_factor * df['Slice']
                ytofit = df['Mean Intensity']
                
                #Convert dfs into numpy arrays
                xtofit = np.array(xtofit)
                ytofit = np.array(ytofit)     
                
                # Restrict to specified x_range
                mask = (xtofit >= self.x_range[0]) & (xtofit <= self.x_range[1])
                ytofit = ytofit[mask]
                xtofit = xtofit[mask]
                
                # Shift xtofit data so that first frame corresponds to t=0
                indices = np.arange(1, len(xtofit) + 1)  
                xtofit = (xtofit - xtofit[0] + self.Switch_params['Frametime'] + self.Switch_params['AddTime'] *
                          (indices - int(self.Switch_params['AddtimeBeforeFrametime'])))
                
                # Shift ytofit to correct for experimental offset
                ytofit = ytofit - self.Switch_params['Fluo_Offset']
                
                # Expected number of frames during specified number of cycles
                n_frames = int(np.floor(self.Switch_params['n_cycle'] * (self.Switch_params['T_ON-OFF'] + self.Switch_params['T_OFF-ON']) /
                                        (self.Switch_params['Frametime'] + self.Switch_params['AddTime'])))

                # Check if sufficient data are there
                if n_frames <= len(xtofit):  # OK in this case
                    xtofit = xtofit[:n_frames]
                    ytofit = ytofit[:n_frames]
                else:
                    # If not enough data, add some zeros at the end
                    zeropad2 = n_frames - len(xtofit)
                    
                    if zeropad2 > 0:
                        # Extend xtofit with extra time points
                        new_xtofit = xtofit[-1] + (np.arange(1, zeropad2 + 1) * (self.Switch_params['Frametime'] + self.Switch_params['AddTime']))
                        xtofit = np.concatenate([xtofit, new_xtofit])
                        
                        # Extend ytofit with zeros
                        ytofit = np.concatenate([ytofit, np.zeros(zeropad2)])
                
                return xtofit, ytofit   
            else:            
                raise ValueError("Too many .csv found in the directory")
        else:
            raise FileNotFoundError("No CSV file found in the directory")
    
    def AssignParams(self):
        
        try:
        
            lasers = []
            
            # Create n_lasers laser objects with default values
            for _ in range(self.Switch_params['n_lasers']):
                lasers.append({
                    'l': [],          # Wavelength
                    'pd': [],         # Power density
                    'phi': np.zeros(self.Switch_params['n_phases']),  # Phase delay 
                    'dc': np.zeros(self.Switch_params['n_phases'])    # Duty cycle
                })       
        
            # Assigning wavelengths to each laser
            lasers[0]['l'] = self.Switch_params['LambdaOnActinic']
            lasers[1]['l'] = self.Switch_params['LambdaOffActinic']
            lasers[2]['l'] = self.Switch_params['LambdaReadout']
            
            # Assigning phase delay to each laser
            # Phase delay for Off switching
            lasers[0]['phi'][0] = self.Switch_params['t2_actinic_during_ON_OFF_frametime'][0]
            lasers[1]['phi'][0] = self.Switch_params['t1_actinic_during_ON_OFF_frametime'][0]
            lasers[2]['phi'][0] = self.Switch_params['t_readout_during_ON_OFF_frametime'][0]
            
            lasers[0]['phi'][1] = self.Switch_params['t2_actinic_during_ON_OFF_addtime'][0]
            lasers[1]['phi'][1] = self.Switch_params['t1_actinic_during_ON_OFF_addtime'][0]
            lasers[2]['phi'][1] = self.Switch_params['t_readout_during_ON_OFF_addtime'][0]
            
            # Phase delay for On switching
            lasers[0]['phi'][2] = self.Switch_params['t2_actinic_during_OFF_ON_frametime'][0]
            lasers[1]['phi'][2] = self.Switch_params['t1_actinic_during_OFF_ON_frametime'][0]
            lasers[2]['phi'][2] = self.Switch_params['t_readout_during_OFF_ON_frametime'][0]
            
            lasers[0]['phi'][3] = self.Switch_params['t2_actinic_during_OFF_ON_addtime'][0]
            lasers[1]['phi'][3] = self.Switch_params['t1_actinic_during_OFF_ON_addtime'][0]
            lasers[2]['phi'][3] = self.Switch_params['t_readout_during_OFF_ON_addtime'][0]
        
            # Assigning duty cycle values to each laser
            # Duty cycles for Off switching
            lasers[0]['dc'][0] = self.Switch_params['t2_actinic_during_ON_OFF_frametime'][1]
            lasers[1]['dc'][0] = self.Switch_params['t1_actinic_during_ON_OFF_frametime'][1]
            lasers[2]['dc'][0] = self.Switch_params['t_readout_during_ON_OFF_frametime'][1]
            
            lasers[0]['dc'][1] = self.Switch_params['t2_actinic_during_ON_OFF_addtime'][1]
            lasers[1]['dc'][1] = self.Switch_params['t1_actinic_during_ON_OFF_addtime'][1]
            lasers[2]['dc'][1] = self.Switch_params['t_readout_during_ON_OFF_addtime'][1]
            
            # Duty cycles for On switching
            lasers[0]['dc'][2] = self.Switch_params['t2_actinic_during_OFF_ON_frametime'][1]
            lasers[1]['dc'][2] = self.Switch_params['t1_actinic_during_OFF_ON_frametime'][1]
            lasers[2]['dc'][2] = self.Switch_params['t_readout_during_OFF_ON_frametime'][1]
            
            lasers[0]['dc'][3] = self.Switch_params['t2_actinic_during_OFF_ON_addtime'][1]
            lasers[1]['dc'][3] = self.Switch_params['t1_actinic_during_OFF_ON_addtime'][1]
            lasers[2]['dc'][3] = self.Switch_params['t_readout_during_OFF_ON_addtime'][1]
        
            # Assigning power densities to each laser
            lasers[0]['pd'] = self.Switch_params['POnActinic']
            lasers[1]['pd'] = self.Switch_params['POffActinic']
            lasers[2]['pd'] = self.Switch_params['PReadout']
            
            #Add the array in the Switch_params dico
            self.Switch_params['lasers'] = lasers

        except Exception as e:
            print(f"An error occurred in Assign params function: {e}")
            return None

    def get_photofatigue_data(self, params):
        
        """
        This part of the code is calculating the different rates for the following model:
            DARK => ON <=> OFF
        """
        
        try:
            
            #Calculate lasers_i for each laser
            lasers_i = []
            
            # Iterate over each laser dictionary in lasers list by index
            for idx, laser in enumerate(self.Switch_params['lasers']):
                if isinstance(laser['l'], (int, float)) and isinstance(laser['pd'], (int, float)):
                    lasers_i.append((laser['pd'] * laser['l']) / 1.9846e-16)
            lasers_i = np.array(lasers_i)
            
            #Create an array that takes into account the number of states and the number of lasers
            eps = np.zeros((self.Switch_params['n_states'], self.Switch_params['n_lasers']))
            eps[0,:] = [self.Switch_params['EpsOn_actOn'], self.Switch_params['EpsOn_actOff'], self.Switch_params['EpsOn_RO']] #Attribute Eps to On state
            eps[1,:] = [self.Switch_params['EpsOff_actOn'], self.Switch_params['EpsOff_actOff'], self.Switch_params['EpsOff_RO']] #Attribute Eps to Off state 
            eps[2,:] = [params['EpsDark_actOn'], params['EpsDark_actOff'], params['EpsDark_RO']] #Attribute Eps to Dark state 
            
            #Calculate the excitation rate for each light sensisive state
            light_sensitive_states = [0, 1, 2]  # Choose the light sensistive state in eps array
            
            k_ex = np.zeros((self.Switch_params['n_states'], self.Switch_params['n_lasers']))
            
            # Loop through the light-sensitive states and compute excitation rates
            for state in light_sensitive_states:
                k_ex[state, :] = lasers_i * 3.82e-21 * eps[state, :]  
            
            """
            Considering adding the quadratic excitation rates
            """
            #Create an array with all the Thermal relaxation rates
            K_th = np.zeros((self.Switch_params['n_states'], self.Switch_params['n_states']))
            #Asign ThermalK 
            K_th[0,1] = params['Th_OnOff']
            K_th[2,1] = params['Th_DarkOff']
            
            #Create an array with all the Qy
            QuantumY = np.zeros((self.Switch_params['n_states'], self.Switch_params['n_states']))        
            #Asign Qy 
            QuantumY[0,1] = params['qOnOff']
            QuantumY[1,0] = params['qOffOn']
            QuantumY[2,0] = params['qDarkOn']
            
            #Compute the light induced rates
            K_l = np.zeros((self.Switch_params['n_states'], self.Switch_params['n_states'], self.Switch_params['n_lasers'])) #(Nb arrays, rows, columns)
            # Iterate over states and lasers to compute k_l
            #For each state, create a 2D array, multiply the switching Qy by the corresponding k_ex (here row 0 = k_ex for Onstate laser, row 1 = k_ex for Offstate laser...)
            for i in range(self.Switch_params['n_states']):  # starting state
            
                for j in range(self.Switch_params['n_states']):  # ending state
                
                    for k in range(self.Switch_params['n_lasers']):  # iterate over lasers
                    
                        K_l[i, j, k] = k_ex[i, k] * QuantumY[i, j]
            
            #Compute fluo rate
            K_fluo = np.zeros(self.Switch_params['n_lasers'])
            
            # Iterate over all lasers (phases) to compute k_fluo
            for k in range(self.Switch_params['n_lasers']):  # go over all phases
            
                K_fluo[k] = k_ex[0, k] * self.Switch_params['Q_Fluo']
            
            #Store the calculated rate in a dictionnary
            rate_array = {
                'K_th' : K_th,
                'K_l' : K_l,
                'K_fluo' : K_fluo
                }
            
            FL, S_RECORDED, S, T_DET, used_dt = get_photofatigue_model(rate_array, self.Switch_params, params)
            
            # Find indices where T_DET is within the fitting range
            T_select = np.where((T_DET >= self.Switch_params['fit_range'][0]) & (T_DET <= self.Switch_params['fit_range'][1]))[0]
            
            # Select only the values within the fit range
            T_DET_IN_FIT_RANGE = T_DET[T_select]
            FL_IN_FIT_RANGE = FL[T_select]
            S_RECORDED_IN_FIT_RANGE = S_RECORDED[:, T_select]
            
            return FL, S_RECORDED, S, T_DET, FL_IN_FIT_RANGE, S_RECORDED_IN_FIT_RANGE, T_DET_IN_FIT_RANGE, used_dt
        
        except Exception as e:
            print(f"An error occurred in get_photofatigue_data function: {e}")
            return None
        
    def fit_function(self, params, x):
        
        try:
        
            print(f"Fitting cycle: {self.fitting_cycle}")
            print(f"Fitting parameters: {params}")

            a2 = self.Fit_params.copy()

            fit_param_keys = list(self.Fit_params.keys())

            for i, idx in enumerate(self.w_variables):
                if idx < len(fit_param_keys):
                    field_name = fit_param_keys[idx]
                    a2[field_name] = params[i]

            # Call your data function
            # Assuming get_photofatigue_data_vsn14 returns a tuple where the fifth element is 'y'
            _, _, _, _, y, _, _, _ = self.get_photofatigue_data(a2)
            
            # Increment the fitting cycle
            self.fitting_cycle += 1
            
            # Ensure output matches the expected size
            y = y[:len(self.ydata)]  # truncate if too long
            if len(y) < len(self.ydata):
                # pad with zeros if too short
                y = np.pad(y, (0, len(self.ydata) - len(y)), mode='constant')
            
            return y
        
        except Exception as e:
            print(f"An error occurred in fit_function : {e}")
            return None   
        
    def fit_DARKONOFF(self, xtofit, ytofit):
        try:
            
            if self.Fit:
                
                #init fitting cycle variable
                self.fitting_cycle=1
                
                #indices the non fixed variables
                self.w_variables = [i for i, param in enumerate(self.start) if param[1] == True]

                InitValues = self.start[:,0]
                #Get the starting values for the fitting variables
                startvalues = InitValues[self.w_variables]
                       
            else:
                
                print("""
                Getting results from starting values only
                """)
                
            # Find the indices where xtofit is within the specified range
            w_data = np.where((xtofit > self.Switch_params['fit_range'][0]) & (xtofit < self.Switch_params['fit_range'][1]))[0]                    
            #Define data to fit within the specified range
            xdata = xtofit[w_data]
            ydata = ytofit[w_data]
            
            parameters = self.Fit_params
            
            # Get the fitted result with starting values only
            yfit_start_full_range, _, _, xfit_start_full_range, yfit_start, _, xfit_start, _ = self.get_photofatigue_data(self.Fit_params)
            
            if self.check_fit:
                
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
            
            if self.Fit:
                print("""
    Fitting of the model
                      """)
                
                # Plotting initial data and start fit line
                min_len = min(len(ydata), len(yfit_start))

                plt.plot(ydata[:min_len], linewidth=1.5, label='Original Data')
                plt.plot(range(1, min_len + 1), yfit_start[:min_len], color='red', linewidth=1.5, label='Fit Line')
                plt.legend()
                plt.show()
                
                # Save ydata for use in fit_function
                self.ydata = ydata
                
                # Fit the model
                popt, pcov = curve_fit(lambda x, *params: self.fit_function(params, x), np.zeros_like(ydata), ydata, p0=startvalues)
                
                # Compute standard error (95% CI approximation assuming normality)
                perr = np.sqrt(np.diag(pcov))
                ci_err = 1.96 * perr  # 95% confidence interval error (assuming normal distributio
                
                print(f"Fitted parameters (popt): {popt}")

                # Initializing fitted_parameters with a copy of self.start
                fitted_parameters = np.copy(self.start)
                print(f"Initial fitted_parameters: {fitted_parameters}")
                
                # Ensure self.w_variables contains valid indices
                print(f"Indices for updating (self.w_variables): {self.w_variables}")
                
                # Check if the length of popt matches the expected size of indices in w_variables
                if len(self.w_variables) != len(popt):
                    print(f"Warning: Length of self.w_variables ({len(self.w_variables)}) does not match length of popt ({len(popt)})")
                
                # Update fitted_parameters with values from popt at the specified indices
                fitted_parameters[self.w_variables] = popt
                print(f"Updated fitted_parameters: {fitted_parameters}")
            else:
                
                fitted_parameters = np.copy(self.start)
        
            # Get the fitted values with starting parameters or fitted parameters
            yfit_all, _, _, xfit_all, yfit, S_RECORDED_IN_FIT_RANGE, xfit, used_dt = self.get_photofatigue_data(fitted_parameters)

            # Plot fitted values and experimental values
            xtitle = 'Time [s]'
            ytitle = 'Concentration [AU]'
            plt.xlabel(xtitle)
            plt.ylabel(ytitle)
            plt.title('PhotoFatigue Model')
            plt.plot(xfit, yfit, linewidth=1.5, label = 'Fitted values')
            plt.plot(xdata, ydata, linewidth=1.5, label = 'Experimental values')
            # Show plot
            plt.show()
            
            if self.ShowFullModel:
        
                # Plot full experimental data
                plt.plot(xtofit, ytofit, color='blue', linewidth=1.5)
                plt.plot(xdata, ydata, color='blue', linewidth=2.5)
        
                # Optional: plot specific data points (commented out in MATLAB)
                # plt.plot(xtofit, ytofit, '*b')
        
                # Labeling
                fig_fontsize = 16
                xtitle = 'Time [s]'
                ytitle = 'Concentration [AU]'
                plt.xlabel(xtitle, fontsize=fig_fontsize, fontweight='bold')
                plt.ylabel(ytitle, fontsize=fig_fontsize, fontweight='bold')
                plt.title('PhotoFatigue Model', fontsize=fig_fontsize, fontweight='bold')
        
                # Plot full and partial fitted curves
                plt.plot(xfit_all, yfit_all, linewidth=1.5)
                plt.plot(xfit, yfit, linewidth=2.5)
        
                # Add legend
                plt.legend([
                    'EXPERIMENTAL DATA FULL',
                    'EXPERIMENTAL DATA FITTED PART',
                    'FITTED DATA FULL',
                    'FITTED DATA WHERE FITTED'
                ])
        
                plt.show()    
            
            return fitted_parameters

        except Exception as e:
            print(f"An error occurred in fit_DARKONOFF function: {e}")
        
        
        
        
