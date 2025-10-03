"""
permafrostbankerosionrmodel.py

JA 2025

"""
import numpy as np
import matplotlib.pyplot as plt
import os
#from scipy import interpolate
import pandas as pd
import warnings
import datetime

warnings.filterwarnings("ignore")

class BankProfile:
    """class for bank profile objects"""

    def __init__(self, x, time, Tw, Ta, stage, D50block):
        """
        initialize profile object

        Parameters
        ----------
        x: array
            x-coordinates of nodes in m
        time: float
            representative time - or age - of profile in seconds
        Tw: float
            water temperature in degrees K
        Ta: float
            air temperature in degrees K
        stage: float
            river water stage in m
        D50block: float
           diameter of collapsed block in m
        """
        self.x = x
        self.time = time
        self.Tw = Tw
        self.Ta = Ta
        self.stage = stage
        self.D50block = D50block

class BankErosionModel:
    """permafrost bank erosion model"""
    def __init__(self, params, bcdf):
        """
        initialize model

        Parameters
        ----------
        params: dict
            all model parameters
        bcdict: dict
            all boundary conditions
        """
        self.D84bed = params["D84bed"] #fluvial gravel 84th percentile diameter
        self.D50bank = params["D50bank"]#bank sediment 50th percentile diameter
        
        self.run_time = params["runtime"]
        self.time = 0
        
        self.dz = params["dz"]
        self.dt = params["dt"]
    
        ### river temperature and seasonality
        self.river_temp = 273.15
        self.To = 273.15 
        self.tensileice = 1e5 #Pa         

        self.tensilepf = 2e5 #Pa         
        self.D50block =  self.D50bank
        self.width = params["width"]
        ### material properties
        self.Fi = params["Fi"] / 100  # .5
        self.L_ice = 334000  # J/kg
        self.g = 9.81
        #heat capacity
        Cs = 1000
        Ci = 2100
        self.Cwater = 4184 #J/kg/degreesC
        self.C_perma =Cs*(1-self.Fi) + Ci*(self.Fi) 

        self.rho_ice = 917  # ice density, kg/m3
        self.k_water =  0.598  # water thermal conductivity, J/smK
        self.k_s = 1 # substrate thermal conductivity J/smK
        self.k_ice=  0.54  # ice therml conductivity J/smK
        self.rho_w = 1000  # water density kg/m3
        self.rho_soil = (2650*(1-self.Fi))+(917*(self.Fi))# substrate density kg/m3
        self.rho_s = 2650 # gravel density kg/m3
        self.rho_block =self.rho_soil# gravel density kg/m3
        self.gravelz = params["gravelz"]
        self.bedslope = params["bedslope"] #longitudnal river bed slope
        self.bankheight = params["bh"] # bank height above fuvial gravel
        self.bankslope_init = 90#params["initial_bankslope"] # 
        ### load external forcing timeseris
        self.T_perma = 273.15 + params["T_permafrost_C"]  #bank substrate temperature, C
        self.tfreeze = params["iceondoy"] #day of year water begins to freeze
        ##ice wedge geometry
        self.iw_spacing = params["Xiw"] # ice wedge spacing, m
        self.iw_depth = self.bankheight
        self.load_bcs(bcdf)
        
        #choose ablation model
        self.ablationmodel = params["ablationmodel"] #or "bankroughness"
        
        # set initial river bank profile
        self.profile = self.create_initial_profile() 

        self.cumulativeheatflux = 0

        return

    def run(self):
        """Run simulation from current profile to specified end time"""
        self.final_time= self.time +self.run_time
                                       
        if str(self.ablationmodel) =='gravelroughness':
            while self.time <= self.final_time:
                self.time += self.dt      
                doy = (self.time / (60 * 60 * 24)) % 365
                    #minimize useless profiles
                lateralerosion = np.zeros_like(self.profile.x)
                collapsedeltax = 0

                self.update_bcs()
        
                if (doy >= 151) & (doy <self.tfreeze):
                    Va, heatflux  = self.get_ablation()

#                     blockVa = Va*np.min(((self.h_stage)/(self.D50block), 1))
#                     self.D50block =np.nanmax((self.D50block - blockVa*self.dt, self.D50bank))

#                     # entrainment = self.entrain()    
                    #                     #block collapse
                    #                     # if entrainment== True:
                    lateralerosion[self.z<= self.h_stage - self.gravelz] = Va*self.dt

                    collapsedeltax, collapsearea, collapsecount = self.solve_collapse() 

                    lateralerosion = lateralerosion +collapsedeltax
                    xnew = self.profile.x+ lateralerosion
                    newprofile = BankProfile(xnew,self.time,self.Tw, self.Ta, self.h_stage- self.gravelz, self.D50block)
                    self.cumulativeheatflux = np.nansum([heatflux*self.dt, self.cumulativeheatflux])
                    self.profile = newprofile
                        #  print("entrained")
        if str(self.ablationmodel) =='bankroughness':
            while self.time <= self.final_time:
                self.time += self.dt      
                doy = (self.time / (60 * 60 * 24)) % 365
                    #minimize useless profiles


                self.update_bcs()
        
                if (self.Tw>self.To) & (doy <self.tfreeze):
                    lateralerosion = np.zeros_like(self.profile.x)
                    collapsedeltax = 0
                    Va, heatflux  = self.get_ablation_bankgrainroughness()

#                     blockVa = Va*np.min(((self.h_stage)/(self.D50block), 1))
#                     self.D50block =np.nanmax((self.D50block - blockVa*self.dt, self.D50bank))

#                     # entrainment = self.entrain()    
                    #                     #block collapse
                    #                     # if entrainment== True:
                    lateralerosion[self.z<= self.h_stage] = Va*self.dt

                    collapsedeltax, collapsearea, collapsecount = self.solve_collapse() 

                    lateralerosion = lateralerosion +collapsedeltax                      #  print("entrained")
                    xnew = self.profile.x+ lateralerosion
                    newprofile = BankProfile(xnew,self.time,self.Tw, self.Ta, self.h_stage, self.D50block)
                    self.cumulativeheatflux = np.nansum([heatflux*self.dt, self.cumulativeheatflux])
                    self.profile = newprofile
  
        return
    def entrain(self):
        """check if block is entrained and reset size if entrained
        Return
        ----------
        excess: byte
            flag for if block/gravel is entrained or not
        """
        excess = True
        phi = 0.05#2.12*self.bedslope + 0.02
        g = 9.81

        U = self.get_Ustar()
        Uc =np.sqrt(phi*g*(self.rho_block - self.rho_w)*self.D50block/self.rho_w)
        excessu = U - Uc

        if (excessu>0)== True:
            excess = True
            self.rho_block = self.rho_s
            self.D50block = self.D50bank
        return excess
# shear velocity
    def get_Ustar(self):
        """get river water shear velocity
        Return
        ------
        ustar: float
            river water shear velocity
        """
   
        ustar = np.sqrt(self.g * self.h_stage* self.bedslope)
        return ustar
    
    def solve_collapse(self):
        """calculate erosion from upper bank collapse, based on Barnhart et al. (2003)
            Return
            -------
            C: array
                erosion for each node in the bank profile
            Carea: float
                area of collapsed block
            Ccount: int
                number of collapses in timestep
            
        """
        x = self.profile.x
        Carea = 0    
        C= np.zeros_like(x)
        Ccount = 0
       # if (np.max(x)< x[-1]):
            #pivot point
        idx =  np.where(x==np.nanmax(x))[0][0]
        xpiv =x[idx]
        zpiv = self.z[idx]
        blockidx = np.where(self.z >zpiv)[0][:]

            #distance to ice wedge
        xiw = self.xiw
            
            #block height
        bheight = np.max(self.z) - zpiv
            
            #bulk density
        rhobulk = self.rho_soil
        lend = (xpiv - x[blockidx])
        lenr = xiw - xpiv 
            
            #driving torque
        blockarea = np.sum((xpiv - x[blockidx])*self.dz)
        Td = np.sum(self.dz*rhobulk*self.g*(xpiv - x[blockidx])) *np.mean((xpiv - x[blockidx]))
             #resisting force * length
        Tiw = np.sum(self.dz*(self.tensileice)*(xpiv - xiw)*(x[blockidx] - xpiv))/bheight
            
            #resisting force * length
        Tpf = ((self.tensilepf)*bheight*(xpiv - x[-1]))          
            #resisting torque
        Tr = bheight*rhobulk*self.g*(xiw - xpiv) + Tpf+Tiw           

            
        if Td>(Tr):
                # view bank profile just before collapse
            C[self.z>=zpiv]=  (xpiv - x[self.z>=zpiv])
            Carea += blockarea
            Ccount += 1
            self.rho_block = self.rho_soil

            self.D50block = blockarea/2#/np.sin(np.pi*35/180)#(blockarea/bheight +  bheight)/2

            while xiw >= self.xiw:          
                self.xiw += self.iw_spacing

        return C, Carea, Ccount
    
    def get_ablation(self):
        """Get river water shear velocity using the model of Costard et al. (2003). This model calculates hydraulic roughness based on bed gravel roughness.
        
        Return
        ------
        Va: float
            instantaneous ablation velocity, m/s
        heatflux: float
            instantaneous heatflux, W/m2       
        """
     
        #Costard et al. (2003)
        H = self.h_stage 
      
        ksgravel =  self.D84bed*2 #gravel roughness length scale
        
        k_w = 0.598
        visc = 1.56e-6
        g = 9.81
        
        Cfbank = ((1/8.1)*(H/(ksgravel))*(-1/6))**2   
        U = self.get_Ustar()/np.sqrt(Cfbank)
        

        #flow reynolds number
        Re = U*H / visc
        
        #flow Prandtl number
        Pr = visc*self.rho_w*self.Cwater / k_w

        #heat transfer coefficient
        Ch = (0.00249*k_w*(Pr**(1/3)) * (Re**(1.0552))) / (H)
     
        #Ablation Velocity       
        denom = ((self.rho_ice*self.L_ice*self.Fi+ self.C_perma*self.rho_soil*(self.To - self.T_perma)))
        Va = Ch*(self.Tw - self.To)/denom
        heatflux= Va*denom
   
        return np.max((0, Va)), heatflux
    def get_ablation_bankgrainroughness(self):
        
        """Get river water shear velocity using the model of Douglas and Lamb (2024). This model calculates hydraulic roughness based on bank sediment roughness.
        
        Return
        ------
        Va: float
            instantaneous ablation velocity, m/s
        heatflux: float
            instantaneous heatflux, W/m2
               
        """
        #modified approach from Yaglom and Kader (1974), Douglas and Lamb (2024)
        H = self.h_stage 
        #bank friction coefficient
        ksbank=  self.D50bank*2.2*2.5
        k_w = 0.598
        visc = 1.56e-6
        g = 9.81

        

        #bank grain roughness   
        Reb = self.get_Ustar()*ksbank/visc
        #fluvial heat transefr coefficient
        Pr = visc*self.rho_w*self.Cwater / k_w

        Betar = np.sqrt(Reb)*(0.55*Pr**(2/3) - 1/11) +0.95
        Betas = 12.5*(Pr**(2/3)) - 6
        Betatrans = (Reb/100)*Betar + (1 - Reb/100)*Betas
        Bt = Betatrans
        
        
        if Reb>100:
             Bt = Betar
    
        Beta1 = 0.5
        alpha2= 2.12
        # heat transfer coefficient
        Ch = self.get_Ustar()/(-alpha2*np.log(ksbank/H) +Bt +Beta1)      
        #Ablation Velocity
        denom = ((self.rho_ice*self.L_ice*self.Fi+ self.C_perma*self.rho_soil*(self.To - self.T_perma)))
        Va = Ch*self.Cwater*self.rho_w*(self.Tw - self.To)/denom
        heatflux= Va*denom
        return np.max((0, Va)), heatflux

    ### Boundary Conditions  ########
    def load_bcs(self, bcdf):
        """read river water conditions from data records"""
        self.bc_record_H =bcdf['H_m']
        self.bc_record_Qw = bcdf['Qw_m3s']
        self.bc_record_Tw = bcdf['Tw_C']+273.15
        self.bc_record_t = bcdf['t_s']
        self.bc_record_Ta = bcdf['Ta_C']+273.15
        
        self.h_stage = 0
        self.discharge = 0
        self.Ta = 273.15
        self.Tw = 273.15

        return 
    
    def update_bcs(self):
        """interpolate river water conditions from data records"""
        time = self.time  % (60*60*24*365)
        self.h_stage = np.interp(time, self.bc_record_t, self.bc_record_H) #self.bc_record_Qw[self.bc_record_t>=self.time][0]
        self.discharge = np.interp(time, self.bc_record_t, self.bc_record_Qw)
        self.Tw = np.interp(time, self.bc_record_t, self.bc_record_Tw)
        self.Ta = np.interp(time, self.bc_record_t, self.bc_record_Ta)

        return 
    #Inital conditions    
                               
    def create_initial_profile(self):
        """
        Return
        ------
        profile: Object
            initial river bank vertical profile
        
        """
        initslope = self.bankslope_init
        self.xiw = self.iw_spacing
        topz = self.bankheight
        self.z = np.arange(0, topz, self.dz)
        x = self.z *(-1/np.tan(initslope/(2*np.pi)))
        profile= BankProfile(x, 0, self.Tw, self.Ta, self.h_stage, self.D50block)
        return profile
    
    def get_erosion(self):
        """
        Return
        ------
        meanerosionrate_mperyear: float
            mean annual erosion rate per year, averaged over the entire simulation
        
        """
        dz = self.dz
        dt = self.dt
        nyears = np.ceil(self.final_time/(60*60*24*365))
        z_interp = self.z[1:-1]+self.dz/2
        meanerosionrate_mperyear=self.dz*np.nanmean(np.interp(z_interp, self.z, self.profile.x))/(nyears*self.bankheight)
        return meanerosionrate_mperyear
    
    def plot_bankprofile(self, title):
        """
        Parameters
        ------
        title: str
            plot title
            
        Return
        ------
        fig: matplotlib figure
            visual represntation of river bank vertical profile
        """
        x = np.concatenate((np.asarray([np.round(np.min(self.profile.x)-10)]), self.profile.x, np.asarray([np.max(self.profile.x)+10])))
        z = np.concatenate((-1.25*np.asarray([self.gravelz]),np.asarray(self.z), np.asarray([np.max(self.z)])))

        fig, ax = plt.subplots(1,1, figsize = [8, 3])

        plt.fill_between(x, z, np.nanpercentile(self.bc_record_H, 100), color = "lightblue", edgecolor = 'blue')

        plt.fill_between(x, z, self.gravelz*-1, color = "grey", edgecolor = 'k')
        plt.scatter(x, z, s = 2, c= 'k')
        plt.title(title)
        ax.set_ylabel("z (m)")
        ax.set_ylabel("x (m)")
        ax.set_xticks(np.arange(0, np.max(self.profile.x)+20, 2))
        ax.set_xticks(np.arange(0, np.max(self.profile.x)+20, .5), minor = True)

        ax.set_yticks(np.arange(-self.gravelz, np.max(self.z)+5, 1))
        ax.set_yticks(np.arange(-self.gravelz, np.max(self.z)+5, .5), minor = True)
        ax.set_xlim((np.round(np.asarray([np.min(self.profile.x)-2])), np.asarray([np.max(self.profile.x)+2])))
        ax.set_ylim((-1*np.asarray([self.gravelz]), 1.25*np.asarray([self.bankheight])))
        return fig
  