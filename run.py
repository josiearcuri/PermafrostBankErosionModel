import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import Image

from src import permabankerosionmodel, boundaryconditionstools, watertemperaturemodel

def get_params():
    """
    define inputs in a dictionary
    
    Return
    ------
    options: dict
        contains all parameters for defining boundary conditions, running the water temperature model, and running the bank erosion model
        
    """
    options = {}

    # Define input directories, this date is used in the accompanying manuscript
    options["dischargedata"] = 'data/primarydata.csv'
    options["metdata"] = 'data/CanningRiverAll_JA.csv'
    options["metdata2"] = 'data/Mfile2.csv'

    # Simulation name
    options['identifier'] = "default"

    # Ablation model "gravelroughness" (Costard et al. 2003) or "bankroughness" (Douglas and Lamb, 2024)
    options['ablationmodel'] = "gravelroughness"  

    # Spacing and timing
    options['dt'] = 60*60*3
    options['dz'] = 0.01
    options['runtime'] =60*60*24*365*2

    ##### River Ice Conditions #####

    # river ice decay starting date,  in days from January 1st
    options['todoy']= 151
    # river ice decay timescale, in days
    options['tau'] = 14
    #day ice starts to accumulate
    options['iceondoy'] = 274 #october 1st 


    ##### River bank material properties #####
    #permafrost temperature
    options["T_permafrost_C"] = -8

    # grain size (84th percentile) of the river bed gravel
    options['D84bed'] = .064
    # grain size (median) of the river bank sediment
    options['D50bank'] = 50e-6

    #height of river bank exposed above gravel bank toe
    options['bh']  = 1.7

    #gravel thickness on river bank toe above the river bed
    options['gravelz'] = .5

    #bankfull river bank height - smallest along the cross section runnign perpendicular to the bank
    options['bfbh']  =   1.6

    #bankfull channel wetted perimeter, width at high stage
    options['width']  = 790

    #ice wedge spacing / distance between troughs on surface, thiis sets the overhanging block dimensions, and how often collapses occur
    options['Xiw'] = 35

    #river bank ice content, in percent of volume, ranging between 40 and 100 %
    options['Fi'] = 75

    #starting river water depth. values below 'starth' are set equal to 'starth', because the river water temperature model is innacurate when flow is shallow. 
    options['starth'] = options['gravelz'] 

    #channel bed slope in m/m
    options['bedslope'] = .003

    #for scenarios without river ice
    options['nullflag'] = False
    return options

def run_model_scenario(params):
    """
    
    run one model scenario based on input dictionary and gather results in a pandas dataframe
    
        
    Parameters
    ----------
    params: dict
        all inputs for one simulation, uses the default values if empty
    
    Return 
    ----------
    results: Pandas DataFrame
        all inputs for one simulation and the resulting mean annual erosion rate and nondimensional number
    """
        
    # define wiver water stage and meteorological boundary conditions in a dictionary
    BCT= boundaryconditionstools.BoundaryConditionsTools(params)
    bcdict = BCT.get_bcs()

    bcdict= watertemperaturemodel.River_temp_model(bcdict)
    
    # run bank erosion model
    BEM = permabankerosionmodel.BankErosionModel(params, bcdict)
    BEM.run()
    erosionrate = BEM.get_erosion()
    resultsdict = {"meanannualerosionrate_mperyr":erosionrate}
    
    #add parameters and results to pandas dataframe
    df = pd.DataFrame(params|resultsdict, index = [0])
    
    #print and plot results
    fig = BEM.plot_bankprofile(title = str("E = " + str(np.round(erosionrate,2))+ " meters/year"))
    plt.savefig("outputs/"+str(params["identifier"])+"_finalprofile.jpg", dpi = 300)
    plt.close()
    return df
#define parameters
params = get_params()

#run one scenario
results = run_model_scenario(params)

#export results
results.to_csv("outputs/results.csv") 