# PermafrostBankErosionModel

**Author:** Josie Arcuri  
**Date:** November 2025  

PermafrostBankErosionModel is a Python-based simulation tool for modeling riverbank erosion in permafrost-dominated environments, where thermal and mechanical processes interact seasonally. The model simulates river bank ablation and collapse under fluvial and thermal forcing.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Input Files](#input-files)
- [Running the Model](#running-the-model)
- [Output](#output)
- [Dependencies](#dependencies)
- [License](#license)
- [Acknowledgments](#acknowledgments)
---

## Features

- Simulates thermal ablation and mechanical collapse of permafrost banks
- Uses realistic forcing from boundary condition time series
- Supports customizable input parameters for regional prmafrost river characteristics
- Outputs daily and annual erosion rates
---


## Installation

Clone the repository and install the required dependencies:

git clone https://github.com/yourusername/PermafrostBankErosionModel.git

cd PermafrostBankErosionModel

### Option 1: Using pip

Install Python dependencies directly:

pip install -r requirements.txt

### Option 2: Using conda

conda env create -f environment.yml
conda activate permafrostbankerosionmodel

---


## Input Files
The model requires two inputs:

Parameter dictionary (params) with keys:

"run_duration" (in seconds)

"dt" (timestep in seconds)

"dz" (vertical resolution in meters)

"icepercent" (volumetric ice content in percent)

"gravelz" ( in meters)

"surface_d84" (84th percentile grain diameter of gravel on the bank toe in meters)

"bank_d50" (median grain size of inorganic particles in bank material in meters)

"bedslope" (river bed slope in m/m)

"bankheight" (bank height in meters)

"initial_bankslope" ()

"iceon" (day of year freezing starts)

"IWspacing" (ice wedge spacing in m)

"ablationmodel" (currently supports 'bank')

Boundary condition dataframe (bcdf) in pandas.DataFrame format with the following columns:

't_s': time in seconds

'Qw_m3s': water discharge (m³/s)

'H_m': stage height (m)

'Tw_C': water temperature (°C)

'Ta_C': air temperature (°C)

Other datasets that are used as input can be found in the data folder and are:

    Mfile1.csv
* 2023 Meteorological data from the Canning River, AK collected with instruments maintained by Frank Urban.

    Mfile2.csv
* 2023 Meteorological data from Deadhorse, AK made available by the Permafrost Laboratory Geophysical Institute (Deadhorse, 70.1613,-148.4653), available from www.permafrostwatch.org

   dmddischarge.csv
* 2009 daily mean discharge from the Canning River above Staines River near Deadhorse station (15955000, 69.8817, -146.3887) (US Geological Survey, 2023).


---

## Running the Model

from src import permafrostbankerosionmodel as pbem
import pandas as pd

### Define parameters
params = {
    "run_duration": 60 * 60 *24 * 365 * 3,  # 3 years in seconds
    "dt": 86400,                      # daily timestep (seconds)
    "dz": 0.05,                       # vertical resolution (m)
    "icepercent": 50,                 # volumetric ice content of bank material in %
    "gravelz": 0.3,                   # elevation of the uppermost surface of the bank toe above the river bed (m)
    "surface_d84": 0.064,             # 84th percentile grain diameter of the clasts on bank toes (m)
    "bank_d50": 0.01,                 # 50th percentile grain diameter of the inorganic sediment in bank soil (m)
    "bedslope": 0.003,                # river bed slope (m/m)
    "bankheight": 2.0,                # river bank surface height above the river bed (m)
    "iceon": 300,                     # end of the open water season (julian day)
    "IWspacing": 10,                  # distance between ice wedge center, or troughs in the floodplain (m)
}

### Load boundary condition data
bcdataframe = pd.read_csv("boundary_conditions.csv")

### Initialize and run model
model = pbem.BankErosionModel(params, bcdataframe)
model.run()

### Get erosion results
mean_annual_erosion_rate = model.get_erosion()
print("Mean annual erosion rate (m/year):", mean_annual_erosion_rate)

---


## Output
default_finalprofile: Bank profile at the end of the simulation. this is a 1-d array of horizantal distance from x = 0.
results: Model inputs and outputs in a data table


---


## Dependencies
numpy
matplotlib
pandas
sphinx
sphinx-autodoc-typehints
furo

---


## License
GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 (GPLv3)

---


## Acknowledgments
Developed by Josie Arcuri (2025). For questions or collaboration inquiries, please contact [Josephine.Arcuri@colorado.edu].

---
