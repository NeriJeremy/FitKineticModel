import matplotlib.pyplot as plt

def Plot_fit(fitted_df, Save_png, Save_dir, Expname):

    plt.plot(fitted_df['Slice'], fitted_df['Mean Intensity'], label='Experimental Data')
    plt.plot(fitted_df['Slice'], fitted_df['Fitted Intensity'], label='Fitted Curve', linestyle='--')
    plt.xlabel('Slices')
    plt.ylabel('Intensity')
    plt.legend()
        
    
    if Save_png:    
        #Save figure
        plt.savefig(Save_dir + Expname + '.png')
    
    plt.show()   
    
    return plt