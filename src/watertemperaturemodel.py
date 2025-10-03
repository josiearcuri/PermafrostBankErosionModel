import pandas as pd
import numpy as np

def River_temp_model(params):
    """
        
    river water temperature model from Zheng, Overeem and Clow (2011)
    
    Parameters
        ----------
        params: dict
            meteorological inputs as time series
        
    Return
        ----------
        params: dict
            input dictionary with water temperature and thawing degree days appended
    """
    
    His = params["His"]
    Hdl = params["Hdl"]
    Ta = params["Ta"]
    Uw =  params["Uw"]
    rh = params["rh"]
    P = params["P"]
    hin = params["H_m"]
    iceoffstart = params["todoy"]
    iceon = params["iceon"]
    tau = params['tau']
    starth = params["gravelz"]
    
    
    dt = 3*60*60
    # --------------------Declare constants----------------------
    c      = 4182  # Heat capacity of fresh water (J kg-1 K-1)
    rho    = 1000  # Density of fresh water (kg m-3)
    R      = 0.2 # White water albedo-McMahon et al., 2017
    Ri = 0.8
    sigma  = 5.67e-8  # Stefan-Boltzmann constant (W m-2 K-4)
    #               ----Riverbed module----
    dz     = 1  # Each depth step is 1 meter
    Nz     = 30  # Choose the number of depth steps
    Nt     = len(params['t_s'])  # Choose the number of time steps
  #  dt     = int((365*24*60*60)/Nt)  # Length of each time step in seconds
    K      = 1.3  # Thermal conductivity of sediment ( W m-1 K-1 )
    density= 2650  # Density of sediment(kg m-3)
    C      = 890 # Heat capacity of sediment
    dsoil  = 30  # Depth without soil temperature variations (m)
    Tc     = -8  # Soil temperature keeps -8°C at dsoil (30 m) same as MAAT
    T      = Tc*np.ones([Nz, Nt-1])  # Create ground temperature matrix with Nz+1 rows, and Nt+1 columns
    r      = 0.6  # fraction of shortwave absorbed in water surface layer
    d      = 0.2  # fraction of shortwave reflected by riverbed-Web and Zhang, 1997
    f      = 0.05  # Attenuation coefficient (m-1)-Web and Zhang, 1997

    # -------------------River inundation timing------------------
    #number of iterations
    iceoffend = iceoffstart+tau
    #seconds 
    secs = params['t_s']
    ddoy= (secs/(60*60*24))%(365)    
 

    DOY = np.floor(ddoy)
    riverice = np.zeros(Nt)
    riverice_index = np.where((ddoy>=iceoffend)&(ddoy<=iceon-1))[0]
    Rs = np.ones(Nt)*Ri
    Rs[riverice_index] = R   
    riverice[np.asarray(riverice_index)] = 1
    h = hin.copy()
    
    startidx = np.where((DOY> iceoffend))[0][0]+1# int(Nt *iceoff/365)/ (60*60*24)#int(iceoff*(24*60*60)/dt)-1
    # water level must allways be => minimum
    h[startidx]= starth
    h[h<starth] = starth
    # -------------------DAILY WATER TEMPERATURE-----------------
    Tw     = np.zeros(Nt)
    Hsr    = np.zeros(Nt)
    Hlr    = np.zeros(Nt)
    Hlh    = np.zeros(Nt)
    Hsc    = np.zeros(Nt)
    Hht    = np.zeros(Nt)
    Hba    = np.zeros(Nt)
    Hbc    = np.zeros(Nt)
    DeltaH = np.zeros(Nt)  # Heat balance
    wt = 0
    endidx = np.where((DOY< iceon))[0][-1]#int(Nt *iceon/365)*dt / (60*60*24)

    for i in range(startidx,endidx):
        
        # -------------------Solar radiation heat gain-------------------
        Hsr[i-1] = (1-R)*His[i-1]

        # -------Longwave radiation heat (Gao and Merrick, 1996)---------
        Hlr[i-1] = Hdl[i-1]-((Tw[i-1]+273.15)**4)*0.97*sigma

        # ----------Evaporation heat gain (Hebert et al., 2011)----------
        # Saturated vapor pressure at the water temperature (mm Hg)
        Es = 25.374*np.exp(17.62-5271/(Tw[i-1]+273.15))
        # Atmospheric water vapour pressure (mm Hg)
        Ea = rh[i-1]*25.374*np.exp(17.62-5271/(Ta[i-1]+273.15))
        Hlh[i-1] = (3*Uw[i-1]+6)*(Ea-Es)

        # ----------Convective heat gain (Hebert et al., 2011)-----------
        Hsc[i-1] = (3.66+1.83*Uw[i-1])*(P[i-1]*0.0075/1000)*(Ta[i-1]-Tw[i-1])

        # ---------riverbed heat conduction (Web and Zhang, 1997)--------
        T[0,i-1] = Tw[i-1]  # First layer temperature (same as water temperature)
        # ------------Finite difference approximations--------------------
        depth_2D = (T[0:-3,i-2]-2*T[1:-2,i-2]+T[2:-1,i-2])/dz**2
        time_1D  = (K/density/C)*(depth_2D)
        T[1:-2,i-1] = time_1D*dt+T[1:-2,i-2]
        # ----------------------------------------------------------------
        T[int(dsoil/dz)-1,i]= Tc  # Enforce bottom BC
        # --------Riverbed heat transfer-------
        Hht[i-1] = -K*(T[0,i-1]-T[1,i-1])/ dz
        # -----Riverbed absorbed radiation-----
        Hba[i-1] = -(1-Rs[i-1])*(1-r)*(1-d)*Hsr[i-1]*np.exp(-f* h[i-1])
        # --------Riverbed heat balance---------
        Hbc[i-1] = Hht[i-1]+Hba[i-1]

        # -----------------------Total heat gain-------------------------
        
        DeltaH[i-1] = Hsr[i-1]+Hlr[i-1]+Hlh[i-1]+Hsc[i-1]+Hbc[i-1]

        # -----------------------Water temperature-----------------------
        dwt = (1/((h[i-1]*rho*c)))*DeltaH[i-1]
        wt = wt+(dwt*dt)
        
        if riverice[i-1]==0:
            wt  = 0
        if h[i-1]==0: 
            wt = wt[i-1]
        wt = max(0,wt)

        Tw[i] = np.max((wt,0))

    DeltaH[(DOY >= iceon) ] = np.nan
    Tw[(    DOY >= iceon)] = np.nan
    Hlh[(   DOY >= iceon) ] = np.nan
    Hsc[(   DOY >= iceon) ] = np.nan
    Hsr[(   DOY >= iceon) ] = np.nan
    Hlr[(   DOY >= iceon) ] = np.nan
    Hbc[(   DOY >= iceon) ] = np.nan
    Tw[(DOY < np.max((iceoffstart,151)))] = np.nan
    #Tw = np.interp(params["t_s"], params["t_s"]//(60*60*24),Tw)
    TDD = np.nansum(Tw*dt)/(60*60*24)


    #ddf_dailymean = np.asarray([np.nanmean(Tw[DOY==i]) for i in range(1, 366)])
    params["Tw_C"] = Tw
  #  params["time_d"] = ddoy_full

    params["thawingdegreedays"] = TDD
    return params