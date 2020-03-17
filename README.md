# hamster

**HAMSTER - a Heat And MoiSturE Tracking fRamework**

**HAMSTER** is an open source software framework to trace heat and moisture through the atmosphere and establish (bias-corrected) source-sink relationships, using output from a Lagrangian model. 

- - - -
## Getting started. 

This sections describes the prerequisites required to run HAMSTER, and describes the steps to install it. 

### Prerequisites
To run HAMSTER, you need 
* python 3
* [anaconda](https://www.anaconda.com/) (or similar to manage python packages)
* [git](https://git-scm.com/)

### Installation
1. Clone the repository
```
git clone https://github.ugent.be/jkeune/hamster.git
cd hamster
```
2. Make an anaconda environment with the necessary python packages
```
conda create -n _newenvironment_ --file requirements.txt
```
or install the packages listed in requirements.txt in your local environment. 

### Contact and support
Dominik Schumacher and Jessica Keune

### License
Copyright 2019 Dominik Schumacher and Jessica Keune

This software is published under the GPLv3 license. This means: 
1. Anyone can copy, modify and distribute this software. 
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. If you dare build your business solely from this code, you risk open-sourcing the whole code base.
6. If you modify it, you have to indicate changes made to the code.
7. Any modifications of this code base MUST be distributed with the same license, GPLv3.
8. This software is provided without warranty.
9. The software author or license can not be held liable for any damages inflicted by the software.

- - - - 
## HAMSTER: modules.
**HAMSTER** consists of 3 modules, 
1. Diagnosis
2. Attribution
3. Bias-correction

which build up on each other. It is suggested to run them sequentially to obtain the most efficient and informative workflow. 

### 1. Diagnosis
The diagnosis part of **HAMSTER** identifies atmospheric fluxes of humidity (precipitation and evaporation) or heat (sensible heat flux) using output from a Lagrangian model. There are several thresholds and criteria that can be set (see docs) to reduce the bias, increase the probability of detection and reduce the probability of false detection. The output from this part can be used to bias correct source-sink relationships. 

### 2. Attribution
The attribution part of **HAMSTER** constructs mass- and energy-conserving trajectories of heat and moisture (e.g. using a linear discounting of changes en route), and establishes a first (biased) source-sink relationship. Multiple options to ensure mass- and energy conservation along trajectories are available (see docs). Various time and space-scales for attribution are possible (see docs). 

### 3. Bias-correction
The last module of **HAMSTER** uses information from the former two modules to bias-correct source-receptor relationships. Multiple options for bias-correction are available (see docs). 

## Running HAMSTER.
### Prerequisites
To execute the full chain (all 3 modules) of **HAMSTER**, the only prerequisites are: 
* Output from a Lagrangian model that traces air parcels and their properties (driven with a reanalysis or output from a GCM/RCM)
* Benchmarking data; e.g., the reanalysis used to run FLEXPART and track parcels
* A file paths.txt which lists the paths where the above data is found and where output will be stored.

The file paths.txt is not part of **HAMSTER**. Users have to create the file theirselves. The order in this file is arbitrary, but it has to contain paths for diagnosis, attribution and biascorrection and reference (benchmark) data: 
```
# This file contains all required paths to run hamster
# INPUT paths
ipath_ATR = "./data/FLEXPART/era_global/"
ipath_DGN = "./data/FLEXPART/era_global/"

# INPUT file name base
ibase_ATR = ["pom_ecoreg1_traj10d_AUXTRAJ_"]
ibase_DGN = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"]

# OUTPUT paths
opath_DGN = "./flexpart_data/hamster/01_diagnosis"
opath_ATR = "./flexpart_data/hamster/02_attribution"
opath_BIA = "./flexpart_data/hamster/03_biascorrection"
```

### Run and settings.
To run **HAMSTER**, do
```
python main.py
```

Note that - without any flags - main.py is run with default values. Use 
```
python main.py -h
```
for more details on setting dates, thresholds and other options. All user-specific paths are set in paths.txt. 
