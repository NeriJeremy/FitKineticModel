import numpy as np
from AssignParams import AssignParams
from Open_csv_fit import Open_csv_fit
from fit_DARKONOFF import fit_DARKONOFF
from Plot_fit import Plot_fit
from DARKONOFFCalcQs import CalcQs

"""

Requirements :Put in the directory one .csv file. The headers should be the following: 'Slice' for the different timepoints and 'Mean Intensity' for fluo data
(You can obtain such a file with PALM data preprocessed with 'show_data' python code)
    
fitting model = DARK=>ON<->OFF

"""
global start, Switch_params
###############################################################################################################################################
#################################################### INPUT PARAMETERS #########################################################################

# Choose the directory containing the .txt files
directory='C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/'

#To check if the starting values are matching with the experimental data
check_fit = False

#If True, fit the starting values or get results from starting values
Fit = True

# Only use the starting values for the fit
Use_Start_Values = True

# If True, Switching quantum yields will be calculated
OptainQy=True

# Save data as csv file and/or as png
Save_csv=False
Save_png=False
Save_dir='C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/'
Expname='Fit_with_DARK_cycles'

# Define range to fit 
x_scale_factor=0.02 # Set to frametime [s]
x_range = [0, np.inf] # Set the limit of data to fit [sec,sec]

#################################################################################################################################################
# Define switching cycles parameters

# Put all the Experiment parameters in a dico
Switch_params = {
    'T_OFF-ON' : 40, #On switching cycle time [s]
    'T_ON-OFF' : 40, #Off switching cycle time [s]
    'Frametime' : 0.02, #Frame duration [s] EMCCD acquisition
    'AddTime' : 0.08, #Add time in between frames [s] For other lasers (actinic)
    'AddtimeBeforeFrametime' : False, #Set to True if addtime is before frametime
    'f_frametime' : 1, #Fraction of frametime effectively detected: set to 1 for continuous EMCCD recording [s] if laser ON shorter than frametime
    'n_cycle' : 10, #Full switching cycles
    'emccd_shift' : 0, #Shift of EMCCD time-integration window relative to t=0; should be limited to a few units of dt [s]
    'offset_to_exp_data' : 200, #Offset to be subtracted to experimental data due to sources of background [Counts]
    'exp_frame_shift' : 0, #Offset to apply to experimental data. <0 value will shift exp data to the left of the plots [# of frames]
    'n_min_samples' : 100, #Minimum # of time points during frametime and addtime
    'fixed_dt' : -1, # [s] set to -1 if automatic determination of dt at every step; otherwise set to a fixed value
    'n_lasers' : 3, #Nb lasers
    'n_phases' : 4, #Nb of phases per cycle (frametime ON-OFF, addtime ON-OFF, frametime OFF-ON, addtime OFF-ON)
    'n_states' : 3, #Nb of states the protein has
    'Q_Fluo' : 0.2,
    'Fluo_Offset' : 100, # Remove offset from experimental fluo data
    'fit_range' : [0,np.inf], # [s]
    
    # Define lasers parameters
    'LambdaOnActinic' : 405, #(nm) On switching laser
    'LambdaOffActinic' : 561, #(nm) Off switching laser
    'LambdaReadout' : 561, #(nm)
    'POnActinic' : 1, #(W.cm-2)
    'POffActinic' : 30, #(W.cm-2)
    'PReadout' : 9, #(W.cm-2)
    
    # Define protein Epsilons 
    'EpsOn_actOn' : 50000, #On state Epsilon during On switch phases for actinic laser (405nm)
    'EpsOn_actOff' : 3000, #On state Epsilon during Off switch phases for actinic laser (561nm)
    'EpsOn_RO' : 3000, #On state Epsilon for readout laser (561nm)
    'EpsOff_actOn' : 2000,
    'EpsOff_actOff' : 30000,
    'EpsOff_RO' : 3000,
    
    # Readout laser informations
    't_readout_during_ON_OFF_frametime' : [0, 1], #[AU] Phase delay and duty cycle of the readout laser during ON_OFF frametime
    't_readout_during_ON_OFF_addtime' : [0, 0], #[AU] Phase delay and duty cycle of the readout laser during ON_OFF addtime
    't_readout_during_OFF_ON_frametime' : [0, 1], #[AU] Phase delay and duty cycle of the readout laser during OFF_ON frametime    
    't_readout_during_OFF_ON_addtime' : [0, 0], #[AU] Phase delay and duty cycle of the readout laser during OFF_ON addtime
    
    # Actinic off switching laser informations
    't1_actinic_during_ON_OFF_frametime' : [0, 0], #[AU] Phase delay and duty cycle of the actinic off switching laser during ON_OFF frametime
    't1_actinic_during_ON_OFF_addtime' : [0.0625, 0.875], #[AU] Phase delay and duty cycle of off switching laser during ON_OFF addtime
    't1_actinic_during_OFF_ON_frametime' : [0, 0], #[AU] Phase delay and duty cycle of off switching laser during OFF_ON frametime    
    't1_actinic_during_OFF_ON_addtime' : [0, 0], #[AU] Phase delay and duty cycle of off switching laser during OFF_ON addtime
    
    # Actinic on switching laser informations
    't2_actinic_during_ON_OFF_frametime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during ON_OFF frametime
    't2_actinic_during_ON_OFF_addtime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during ON_OFF addtime
    't2_actinic_during_OFF_ON_frametime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during OFF_ON frametime
    't2_actinic_during_OFF_ON_addtime' : [0.0625, 0.875], #[AU] Phase delay and duty cycle of on switching laser during OFF_ON addtime

}

###################################################################################################################################################

#################################################################################################################################################

# Put all the Experiment parameters in a dico
Fit_params = {
    
    # Define protein Epsilons to fit
    'EpsDark_actOn' : 1000,
    'EpsDark_actOff' : 1000,
    'EpsDark_RO' : 1000,

    # Define initial concentration of each state, no concentration should be equal to zero
    'FP_init' : 200, #Initial level of fluorescence
    'On_init' : 0.21, #Fraction of On state proteins
    'Off_init' : 0.11, #Fraction of Off state proteins
    'Dark_init' : 0.68, #Fracion of Dark state proteins
    
    # Define the switching quantum yields
    'qOnOff' : 0.004,
    'qOffOn' : 0.004,
    'qDarkOn' : 0.00005,
    
    # Define the thermal recovery rates
    'Th_OnOff' : 0,
    'Th_DarkOff' : 0,
    
}

# Define the fixed parameters
start = np.array([[Fit_params['EpsDark_actOn'], False],
                  [Fit_params['EpsDark_actOff'], True],
                  [Fit_params['EpsDark_RO'], True],
                  [Fit_params['FP_init'], True],
                  [Fit_params['On_init'], True],
                  [Fit_params['Off_init'], True],
                  [Fit_params['Dark_init'], True],
                  [Fit_params['qOnOff'], True],
                  [Fit_params['qOffOn'], True],
                  [Fit_params['qDarkOn'], True],
                  [Fit_params['Th_OnOff'], True],
                  [Fit_params['Th_DarkOff'], True],
                  ])

###################################################################################################################################################

##################################################END OF INPUT PARAMETERS##########################################################################
###################################################################################################################################################

Switch_params = AssignParams(Switch_params)

xtofit, ytofit = Open_csv_fit(directory, x_scale_factor, x_range, Switch_params)

fitted_df, fitted_parameters, Constants = fit_DARKONOFF(xtofit, ytofit, Use_Start_Values, Fit_params, Save_csv, directory, Expname, Save_dir, Fit, check_fit, Switch_params, start)

plt = Plot_fit(fitted_df, Save_png, Save_dir, Expname)

QySwitching = CalcQs(Exp_params, Constants, Save_csv, Save_dir, Expname, fitted_parameters, OptainQy)

print(fitted_parameters)
print(QySwitching)