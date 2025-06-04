import pandas as pd

"""

kOn switching phase:
    kABon = a.E(fluo561).QsOff.PfluoReadout.lambda + a.E(fluo405).QsOff.Pswitchon.lambda
    kBAon = a.E(off561).QsOn.PfluoReadout.lambda + a.E(off405).QsOn.Pswitchon.lambda
    kOnDarkOn = a.(Efluo561).QsOnDark.PfluoReadout.lambda + a.(Efluo405).QsOnDark.Pswitchon.lambda
    kDarkOnOn = a.(Edark561).QsDarkOn.PfluoReadout.lambda + a.(Edark405).QsDarkOn.Pswitchon.lambda

kOFF switching phase:
    kABoff = a.E(fluo561).Qsoff.PfluoReadout.lambda + a.E(fluo561).Qsoff.Pswitchoff.lambda
    KBAoff = a.E(off561).Qson.PfluoReadout.lambda + a.E(off561).Qson.Pswitchoff.lambda
    kOnDarkOff = a.(Efluo561).QsOnDark.PfluoReadout.lambda + a.(Efluo561).QsOnDark.Pswitchoff.lambda
    kDarkOnOff= a.(Edark561).QsDarkOn.PfluoReadout.lambda + a.(Edark561).QsDarkOn.Pswitchoff.lambda

"""

def CalcQs(Exp_params, Constants, Save_csv, Save_dir, Expname, fitted_parameters, OptainQy):
    
    if OptainQy:

        #Calculate QySwitching for OnSwitching phase
        QsOnOn = Constants['kBAon'] / (Exp_params['aMulti']
                                       * (Exp_params['Eoff561'] * Exp_params['P561onReadout'] * Exp_params['LambdaOnReadout']
                                          + Exp_params['Eoff405'] * Exp_params['P405onActinic'] * Exp_params['LambdaOnActinic']))
        
        QsOffOn = Constants['kABon'] / (Exp_params['aMulti']
                                        * (Exp_params['Efluo561'] * Exp_params['P561onReadout'] * Exp_params['LambdaOnReadout']
                                         + Exp_params['Efluo405'] * Exp_params['P405onActinic'] * Exp_params['LambdaOnActinic']))
        
        QsOnDarkOn = Constants['kOnDarkOn'] / (Exp_params['aMulti']
                                               * (Exp_params['Efluo561'] * Exp_params['P561onReadout']* Exp_params['LambdaOnReadout']
                                                  + Exp_params['Efluo405'] * Exp_params['P405onActinic'] * Exp_params['LambdaOnActinic']))
        
        QsDarkOnOn = Constants['kDarkOnOn'] / (Exp_params['aMulti']
                                               * (Exp_params['Edark561'] * Exp_params['P561onReadout'] * Exp_params['LambdaOnReadout']
                                                  + Exp_params['Edark405'] * Exp_params['P405onActinic'] * Exp_params['LambdaOnActinic']))
        
        #Calculate QySwitching for OffSwitching phase
        QsOnOff = Constants['kBAoff'] / (Exp_params['aMulti']
                                         * Exp_params['Eoff561'] * (Exp_params['P561offReadout'] * Exp_params['LambdaOffReadout']
                                                                    + Exp_params['P561offActinic'] * Exp_params['LambdaOffActinic']))
        
        QsOffOff = Constants['kABoff'] / (Exp_params['aMulti']
                                          * Exp_params['Efluo561'] * (Exp_params['P561offReadout'] * Exp_params['LambdaOffReadout']
                                                                      + Exp_params['P561offActinic'] * Exp_params['LambdaOffActinic']))
        
        QsOnDarkOff = Constants['kOnDarkOff'] / (Exp_params['aMulti']
                                                 * Exp_params['Efluo561'] * (Exp_params['P561offReadout']* Exp_params['LambdaOffReadout']
                                                                             + Exp_params['P561offActinic'] * Exp_params['LambdaOffActinic']))
        
        QsDarkOnOff = Constants['kDarkOnOff'] / (Exp_params['aMulti']
                                                 * Exp_params['Edark561'] * (Exp_params['P561offReadout'] * Exp_params['LambdaOffReadout']
                                                                             + Exp_params['P561offActinic'] * Exp_params['LambdaOffActinic']))
        
        #Create a df with the calculated QySwitch
        QySwitching = pd.DataFrame(
            {'QyOnSwitching' : [QsOnOn, QsOffOn, QsOnDarkOn, QsDarkOnOn],
             'QyOffSwitching' : [QsOnOff, QsOffOff, QsOnDarkOff, QsDarkOnOff]},
            index = ['QsOn', 'QsOff', 'QsOnDark', 'QsDarkOn'])
        
        if Save_csv:
            QySwitching.to_csv(Save_dir + Expname + 'QySwitching' + '.csv')
    
        return QySwitching
    
    if not OptainQy:
        return None