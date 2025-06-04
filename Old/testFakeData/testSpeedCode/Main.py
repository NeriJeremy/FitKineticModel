import numpy as np
from fun.FitModule import FitFun

"""

Requirements :Put in the directory one .csv file. The headers should be the following: 'Slice' for the different timepoints and 'Mean Intensity' for fluo data
(You can obtain such a file with PALM data preprocessed with 'show_data' python code)
    
fitting model = DARK=>ON<->OFF

"""

###############################################################################################################################################
#################################################### INPUT PARAMETERS #########################################################################
Global_Params = {
    # Choose the directory containing the .txt files
    'directory' : 'C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/testFakeData',

    #To check if the starting values are matching with the experimental data
    'check_fit' : False,

    #If True, fit the starting values or get results from starting values
    'Fit' : True,

    # Only use the starting values for the fit
    'Use_Start_Values' : False,

    # Plot fitted values
    'PlotFittedValues' : True,

    # Plot percentage of each population
    'ShowFullModel' : True,

    # Save data as csv file and/or as png
    'SaveFitParam' : False, # The output is .txt file
    'SavePlots' : False, # The outptut is .png
    'Save_dir' : 'C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/testFakeData/testSpeedCode',
    'Expname' : 'Fit_with_DARK_cycles',

    # Define range to fit 
    'x_scale_factor' : 1, # Set to frametime [s]
    'x_range' : [0, np.inf], # Set the limit of data to fit [sec,sec]
    
    }
#################################################################################################################################################
# Define switching cycles parameters

# Put all the Experiment parameters in a dico
Switch_params = {
    'T_OFF-ON' : 600, #On switching cycle time [s]
    'T_ON-OFF' : 1742, #Off switching cycle time [s]
    'Frametime' : -1, #Frame duration [s] EMCCD acquisition -1 for automated frametime determination
    'AddTime' : 0.025, #Add time in between frames [s] For other lasers (actinic)
    'AddtimeBeforeFrametime' : False, #Set to True if addtime is before frametime
    'f_frametime' : 1, #Fraction of frametime effectively detected: set to 1 for continuous EMCCD recording [s] if laser ON shorter than frametime
    'n_cycle' : 1, #Full switching cycles
    'emccd_shift' : 0, #Shift of EMCCD time-integration window relative to t=0; should be limited to a few units of dt [s]
    'offset_to_exp_data' : 10, #Offset to be subtracted to experimental data due to sources of background [Counts]
    'exp_frame_shift' : 0, #Offset to apply to experimental data. <0 value will shift exp data to the left of the plots [# of frames]
    'n_min_samples' : 10, #Minimum # of time points during frametime and addtime
    'fixed_dt' : -1, # [s] set to -1 if automatic determination of dt at every step; otherwise set to a fixed value
    'n_lasers' : 3, #Nb lasers
    'n_phases' : 4, #Nb of phases per cycle (frametime ON-OFF, addtime ON-OFF, frametime OFF-ON, addtime OFF-ON)
    'n_states' : 3, #Nb of states the protein has
    'Q_Fluo' : 0.7,
    'Fluo_Offset' : 10, # Remove offset from experimental fluo data
    'fit_range' : [0,np.inf], # [s]
    
    # Define lasers parameters
    'LambdaOnActinic' : 405, #(nm) On switching laser
    'LambdaOffActinic' : 488, #(nm) Off switching laser
    'LambdaReadout' : 488, #(nm)
    'POnActinic' : 1, #(W.cm-2)
    'POffActinic' : 0, #(W.cm-2)
    'PReadout' : 1, #(W.cm-2)
    
    # Define protein Epsilons 
    'EpsOn_actOn' : 0, #On state Epsilon during On switch phases for actinic laser (405nm)
    'EpsOn_actOff' : 0, #On state Epsilon during Off switch phases for actinic laser (561nm)
    'EpsOn_RO' : 30000, #On state Epsilon for readout laser (561nm)
    'EpsOff_actOn' : 30000,
    'EpsOff_actOff' : 0,
    'EpsOff_RO' : 0, 
    
    # Readout laser informations
    't_readout_during_ON_OFF_frametime' : [0, 1], #[AU] Phase delay and duty cycle of the readout laser during ON_OFF frametime
    't_readout_during_ON_OFF_addtime' : [0, 0], #[AU] Phase delay and duty cycle of the readout laser during ON_OFF addtime
    't_readout_during_OFF_ON_frametime' : [0, 1], #[AU] Phase delay and duty cycle of the readout laser during OFF_ON frametime    
    't_readout_during_OFF_ON_addtime' : [0, 0], #[AU] Phase delay and duty cycle of the readout laser during OFF_ON addtime
    
    # Actinic off switching laser informations
    't1_actinic_during_ON_OFF_frametime' : [0, 0], #[AU] Phase delay and duty cycle of the actinic off switching laser during ON_OFF frametime
    't1_actinic_during_ON_OFF_addtime' : [0, 0], #[AU] Phase delay and duty cycle of off switching laser during ON_OFF addtime
    't1_actinic_during_OFF_ON_frametime' : [0, 0], #[AU] Phase delay and duty cycle of off switching laser during OFF_ON frametime    
    't1_actinic_during_OFF_ON_addtime' : [0, 0], #[AU] Phase delay and duty cycle of off switching laser during OFF_ON addtime
    
    # Actinic on switching laser informations
    't2_actinic_during_ON_OFF_frametime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during ON_OFF frametime
    't2_actinic_during_ON_OFF_addtime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during ON_OFF addtime
    't2_actinic_during_OFF_ON_frametime' : [0, 0], #[AU] Phase delay and duty cycle of on switching laser during OFF_ON frametime
    't2_actinic_during_OFF_ON_addtime' : [0.03, 0.96], #[AU] Phase delay and duty cycle of on switching laser during OFF_ON addtime

}

###################################################################################################################################################

#################################################################################################################################################

# Define protein Epsilons to fit
EpsDark_actOn = 0
EpsDark_actOff = 0
EpsDark_RO = 0

# Define initial concentration of each state, no concentration should be equal to zero
FP_init = 2.92943106e+03 #Initial level of fluorescence
On_init = 1#8.86542249e-04 #Fraction of On state proteins
Off_init = 0#8.28699531e-01 #Fraction of Off state proteins
Dark_init = 0 #Fracion of Dark state proteins
    
# Define the switching quantum yields
qOnOff = 1.22405067e-04
qOffOn = 5.14903077e-03
qDarkOn = 0
    
# Define the thermal recovery rates
Th_OnOff = 2.48729169e-07
Th_DarkOff = 0

# Define the fixed parameters
start = np.array([[EpsDark_actOn, False, 0, np.inf], # [Value to Fit, Fit ?, lower value to fit, upper value to fit]
                  [EpsDark_actOff, False, 0, np.inf],
                  [EpsDark_RO, False, 0, np.inf],
                  [FP_init, True, 0, np.inf],
                  [On_init, False, 0, 1],
                  [Off_init, False, 0, 1],
                  [Dark_init, False, 0, 1],
                  [qOnOff, True, 0, 1],
                  [qOffOn, True, 0, 1],
                  [qDarkOn, False, 0, 1],
                  [Th_OnOff, True, 0, 1],
                  [Th_DarkOff, False, 0, 1],
                  ])

###################################################################################################################################################

##################################################END OF INPUT PARAMETERS##########################################################################
###################################################################################################################################################

data = FitFun(Switch_params, Global_Params, start)

xtofit, ytofit = data.Open_csv_fit()

coef = data.fit_DARKONOFF(xtofit, ytofit)

