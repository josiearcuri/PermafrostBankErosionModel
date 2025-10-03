import pandas as pd
import numpy as np

class BoundaryConditionsTools:
    """
    functions for generating boundary conditions as inputs for the permafrost bank erosion model
    
    JA 2025
    """
    def __init__(self, params):
        self.params = params
        self.time_s = np.arange(0, int(365*24*60*60), params["dt"])

        return
        
    def get_bcs(self):
        """
        
        Return
        ------
        bcdict: dict
            concatenated hyrological and thermal boundary conditions used as inputs for permafrost bank erosion model
        """

        #met inputs
        metdict = self.prep_metdata()
        
        
        #hydro inputs
        hydrodf, nondims = self.prep_hydrodata()

        bcdict = {"t_s": self.time_s,'H_m': hydrodf.H , 'Qw_m3s': hydrodf.Q, 'Ai_m2': hydrodf.Ai}
        
        
        bcdict["todoy"] = self.params["todoy"]
        bcdict["iceon"] = self.params["iceondoy"]
        bcdict["tau"] = self.params['tau']
        bcdict["starth"] = self.params["starth"]
        bcdict["gravelz"] = self.params["gravelz"]

        return bcdict|metdict|nondims
    def prep_metdata(self):
        """
        read in meteorological data and interpolate arrays
        
        """
        options = self.params
        dt = options['dt']
        nits = len(self.time_s)
        #USGS
        #met station
        Mdata = pd.read_csv(options['metdata'], header=0)
        Mdata['TIMESTAMP']  = pd.to_datetime(Mdata['TIMESTAMP'])
        Mdata['doy'] = Mdata['TIMESTAMP'].dt.dayofyear #
        Mdata['year']= Mdata['TIMESTAMP'].dt.year.values
        Mdata['hr'] = Mdata['TIMESTAMP'].dt.hour.values
        Mdata['ddoy'] = Mdata['doy'].values + Mdata['hr'].values/24 
        Mdata['relsecs']= Mdata['ddoy']*24*60*60
        
        Mdata= Mdata[Mdata['year'] == np.max(Mdata['year']) -1]
        
        #UAF 
        Mdata2 = pd.read_csv(options['metdata2'], header=0)#.fillna(0)
        Mdata2['ddoy'] =  Mdata2['DOY']+.5
        Mdata2['hr'] = 12
        Mdata2['relsecs']= Mdata2['ddoy']*24*60*60

        # His = np.zeros(nits)
        # P = np.zeros(nits)
        # Uw = np.zeros(nits)
        # rh = np.zeros(nits)
        # Ta = np.zeros(nits)

#         for i in range(nits):
#             tstart = dt*i
#             tend = tstart+dt

        His = np.interp(self.time_s, Mdata['relsecs'], Mdata['SlrkW_m2_Avg_down'])*1000
        P = np.interp(self.time_s, Mdata['relsecs'], Mdata["P_kpa"])*1000
        Uw = np.interp(self.time_s, Mdata['relsecs'], Mdata["WS_ms_Avg"])
        rh = np.interp(self.time_s, Mdata['relsecs'],Mdata["RelHum_per"])/100
        Ta = np.interp(self.time_s, Mdata['relsecs'],Mdata["AirTC_Avg"])
        Hdl= np.interp(self.time_s, Mdata2['relsecs'], Mdata2['Hdl'])


        inputsdf = {'Hdl': Hdl, 'His':His, 'Ta':Ta, 'Uw':Uw, 'rh': rh, 'P': P, "Ta_C":Ta}

        return inputsdf
    
    
    def prep_hydrodata(self):
        """
        read in river water discharge time series and interpolate arrays for stage
        Return
        ------
        hydrodf: Pandas DataFrame
            hydrological and thermal inputs for permafrost bank erosion model
        cumulativeexposure: float
            total river water heat flux over 1 year, kW^2
        cumulativeiceblockage: float
            annual cumulative river channel area blockage by river ice
        cumulativedischarge: float
            annual cumulative river water dicharge 
            
        """
        
        #Read usgs daily mean discharge
        options = self.params
        data = pd.read_csv(options['dischargedata'], skiprows= range(33), usecols = [2,4,5], names = ["datetime","discharge_cfs","code"], parse_dates = ["datetime"])

        #number of iterations
        nits = len(self.time_s)

        #seconds to decimal days
        ddoy = self.time_s/(60*60*24)#)%(365) 
        # aconvert to metric 
        data["discharge_m3s"] = data["discharge_cfs"] * 0.028316847
        
        #add some new columns to interpolating discharge
        data['doy'] = data["datetime"].dt.dayofyear
        data["ddoy"] = data["datetime"].dt.dayofyear + (data["datetime"].dt.hour +data["datetime"].dt.minute/60)/24
        data["month"] = data["datetime"].dt.month
        data["year"] = data["datetime"].dt.year

        #Bankfull Channel Area
        Abf = options["width"]*options["bfbh"] 

        #Proportion of channel area filled by river ice
        #ice area based on exponential decay
        Ai= 0*ddoy
        Ai[ddoy>=options["todoy"]] = Abf*(np.exp(-1*((ddoy[ddoy >=options["todoy"]] -options["todoy"]))/(options["tau"])))

        #River Water Discharge
        #daily mean from 2009 
        #apex initial discharge
        Qtot = np.interp(ddoy, data.ddoy[data.year == 2009], data.discharge_m3s[data.year == 2009])

        Q= Qtot*(1-(Ai/ Abf))
        Q[(ddoy<options["todoy"])&(ddoy<151)] = 0
        
        #Beginning of record is not avaiable
        Q[ddoy>options["iceondoy"]] = 0

        Q[np.isnan(Q)] = 0

        n = 0.05
        H = np.zeros_like(Q)
        
        #Prowse model for mean water height during breakup
        H=1.32*((Q*n)/(1.49*options["width"]*options["bedslope"]**0.5))**(3/5) +917*(Ai/options["width"])/1000 

       #no water before breakup starts, or after freezing starts
        H[ddoy<options["todoy"]] =H[options["todoy"]]

        H[ddoy>=options["iceondoy"]] = 0

        hydrodf = pd.DataFrame({"Q":Q,"H":H,"Ai": Ai})
        
        cumulativeexposure = np.nansum(H/options["bfbh"])*options["dt"]
        cumulativeiceblockage =np.nansum(Ai/Abf)*options["dt"]
        cumulativedischarge = np.nansum(Q/Qtot)*options["dt"]
    
        return hydrodf, {"Hstar": cumulativeexposure, "Astar":cumulativeiceblockage, "Qstar":cumulativedischarge}
