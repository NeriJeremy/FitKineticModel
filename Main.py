import numpy as np
from fun.FitModule import FitFun

"""

Requirements :Put in the directory one .csv file. The headers should be the following: 'Slice' for the different timepoints and 'Mean Intensity' for fluo data
(You can obtain such a file with PALM data preprocessed with 'show_data' python code)

V1: _ This version takes calculate the kinetic model for Pink taking into account
the frametime (Cmos recording time) and the the addtime (time without fluorescence recording)
_ Considering a revision for the thermal recovery calculation
_ This version fit the following model : DARK=>ON<->OFF
_ Versionning on GitHub 


"""

###############################################################################################################################################
#################################################### INPUT PARAMETERS #########################################################################
Global_Params = {
    # Choose the directory containing the .txt files
    'directory' : 'C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/',

    # To check if the starting values are matching with the experimental data
    'check_fit' : False,

    # If True, fit the starting values or get results from starting values
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
    'Save_dir' : 'C:/Users/jneri/Desktop/Local/Code/fit_OnOffDarkCycles/TestSpeed/V6/output',
    'Expname' : 'Test_NewEps',

    # Define range to fit 
    'x_scale_factor' : 0.02, # Set to frametime [s]
    'x_range' : [0, 799.5], # Set the limit of data to fit [sec,sec]
    
    }
#################################################################################################################################################
# Define switching cycles parameters

# Put all the Experiment parameters in a dico
Switch_params = {
    'T_OFF-ON' : 40, #On switching cycle time [s]
    'T_ON-OFF' : 40, #Off switching cycle time [s]
    'Frametime' : -1, #Frame duration [s] EMCCD acquisition -1 for automated frametime determination
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
    'Q_Fluo' : 0.1,
    'Fluo_Offset' : 300, # Remove offset from experimental fluo data
    'fit_range' : [0,799.5], # [s]
    
    # Define lasers parameters
    'LambdaOnActinic' : 405, #(nm) On switching laser
    'LambdaOffActinic' : 561, #(nm) Off switching laser
    'LambdaReadout' : 561, #(nm)
    'POnActinic' : 1, #(W.cm-2)
    'POffActinic' : 30, #(W.cm-2)
    'PReadout' : 9, #(W.cm-2)
    
    # Define protein Epsilons 
    'EpsOn_actOn' : 0, #On state Epsilon during On switch phases for actinic laser (405nm)
    'EpsOn_actOff' : 25000, #On state Epsilon during Off switch phases for actinic laser (561nm)
    'EpsOn_RO' : 25000, #On state Epsilon for readout laser (561nm)
    'EpsOff_actOn' : 15450,
    'EpsOff_actOff' : 0,
    'EpsOff_RO' : 0,
    
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

# Define protein Epsilons to fit
EpsDark_actOn = 10950
EpsDark_actOff = 0
EpsDark_RO = 0

# Define initial concentration of each state, no concentration should be equal to zero
FP_init = 879 # Initial level of fluorescence
On_init = 0.2 # Fraction of On state proteins
Off_init = 0.1 # Fraction of Off state proteins
Dark_init = 0.7 # Fracion of Dark state proteins
    
# Define the switching quantum yields
qOnOff = 0.914417594273859e-04
qOffOn = 2.433273959647497e-02
qDarkOn = 0.00005012302919827813
    
# Define the thermal recovery rates
Th_OnOff = 8.879256196310384e-28
Th_DarkOff = 0

# Define the fixed parameters
start = np.array([[EpsDark_actOn, True, 0, np.inf], # [Value to Fit, Fit ?, lower value to fit, upper value to fit]
                  [EpsDark_actOff, False, 0, np.inf],
                  [EpsDark_RO, False, 0, np.inf],
                  [FP_init, True, 0, np.inf],
                  [On_init, True, 0, 1],
                  [Off_init, True, 0, 1],
                  [Dark_init, False, 0, 1],
                  [qOnOff, True, 0, 1],
                  [qOffOn, True, 0, 1],
                  [qDarkOn, True, 0, 1],
                  [Th_OnOff, False, 0, 1],
                  [Th_DarkOff, False, 0, 1],
                  ])

###################################################################################################################################################

##################################################END OF INPUT PARAMETERS##########################################################################
###################################################################################################################################################

data = FitFun(Switch_params, Global_Params, start)

xtofit, ytofit = data.Open_csv_fit()

coef = data.fit_DARKONOFF(xtofit, ytofit)

