# PermafrostBankErosionModel

**Author:** Josie Arcuri  
**Date:** October 2025  

PermafrostBankErosionModel is a Python-based simulation tool for modeling riverbank erosion in permafrost-dominated environments, where thermal and mechanical processes interact seasonally. The model simulates river bank ablation and collapse under fluvial and thermal forcing.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Input Files](#input-files)
- [Running the Model](#running-the-model)
- [Output](#output)
- [Dependencies](#dependencies)
- [License](#license)
- [Acknowledgments] (#acknowledgments)
---

## Features

- Simulates thermal ablation and mechanical collapse of permafrost banks
- Uses realistic fluvial forcing from boundary condition time series
- Supports customizable input parameters for river and permafrost characteristics
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
conda activate permafrostbankerosionrmodel

---


## Input Files
The model requires two inputs:

Parameter dictionary (params) with keys such as:

"run_duration" (in seconds)

"dt" (timestep in seconds)

"dz" (vertical resolution in meters)

"icepercent", "gravelz", "surface_d84", "bank_d50"

"bedslope", "bankheight", "initial_bankslope"

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

    CanningRiverAll_JA.csv
* Met data from the USGS Canning River Meterological Station

    Mfile1.csv
* Met data from XXXXXX 

    Mfile2.csv
* Met data from XXXXXXX 

    primarydata.csv
* USGS daily mean discharge


---


## Running the Model

from permafrostbankerosionrmodel import permafrostbankerosionrmodel
import pandas as pd

### Define parameters
params = {
    "run_duration": 86400 * 365 * 3,  # 3 years in seconds
    "dt": 86400,                      # daily timestep
    "dz": 0.05,                       # vertical resolution (m)
    "icepercent": 50,
    "gravelz": 0.3,
    "surface_d84": 0.05,
    "bank_d50": 0.01,
    "bedslope": 0.001,
    "bankheight": 2.0,
    "initial_bankslope": 0.5,
    "iceon": 300,
    "IWspacing": 1.5,
    "ablationmodel": "bank"
}

### Load boundary condition data
bcdf = pd.read_csv("boundary_conditions.csv")

### Initialize and run model
model = IcyRiver(params, bcdf)
model.run()

### Get erosion results
daily_rates, yearly_mean, erosion_series = model.get_erosion()
print("Mean yearly erosion rate (m/year):", yearly_mean)

---


## Output
default_finalprofile: Bank profile at the end of the simulation
results: Model input and output data


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
