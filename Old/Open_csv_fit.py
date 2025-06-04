import os
import pandas as pd
import numpy as np

def Open_csv_fit(directory, x_scale_factor, x_range, Switch_params):

    files = os.listdir(directory)
    
    # Find the .csv files
    csv_files = [file for file in files if file.endswith('.csv')]
    
    if csv_files:

        if len(csv_files) == 1:   
            filename = csv_files[0]  # First .csv file
            file_path = os.path.join(directory, filename)
            
            # Read the .csv file into a DataFrame
            df = pd.read_csv(file_path, sep=',', engine='python')
            
            #Convert the X data in sec and put X and Y data in a variable
            xtofit = x_scale_factor * df['Slice']
            ytofit = df['Mean Intensity']
            
            #Convert dfs into numpy arrays
            xtofit = np.array(xtofit)
            ytofit = np.array(ytofit)     
            
            # Restrict to specified x_range
            mask = (xtofit >= x_range[0]) & (xtofit <= x_range[1])
            ytofit = ytofit[mask]
            xtofit = xtofit[mask]
            
            # Shift xtofit data so that first frame corresponds to t=0
            indices = np.arange(1, len(xtofit) + 1)  
            xtofit = (xtofit - xtofit[0] + Switch_params['Frametime'] + Switch_params['AddTime'] *
                      (indices - int(Switch_params['AddtimeBeforeFrametime'])))
            
            # Shift ytofit to correct for experimental offset
            ytofit = ytofit - Switch_params['Fluo_Offset']
            
            # Expected number of frames during specified number of cycles
            n_frames = int(np.floor(Switch_params['n_cycle'] * (Switch_params['T_ON-OFF'] + Switch_params['T_OFF-ON']) / (Switch_params['Frametime'] + Switch_params['AddTime'])))

            # Check if sufficient data are there
            if n_frames <= len(xtofit):  # OK in this case
                xtofit = xtofit[:n_frames]
                ytofit = ytofit[:n_frames]
            else:
                # If not enough data, add some zeros at the end
                zeropad2 = n_frames - len(xtofit)
                
                if zeropad2 > 0:
                    # Extend xtofit with extra time points
                    new_xtofit = xtofit[-1] + (np.arange(1, zeropad2 + 1) * (Switch_params['Frametime'] + Switch_params['AddTime']))
                    xtofit = np.concatenate([xtofit, new_xtofit])
                    
                    # Extend ytofit with zeros
                    ytofit = np.concatenate([ytofit, np.zeros(zeropad2)])
            
            return xtofit, ytofit   
        else:            
            raise ValueError("Too many .csv found in the directory")
    else:
        raise FileNotFoundError("No CSV file found in the directory")       