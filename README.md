# HAMSTER – a Heat And MoiSture Tracking framEwoRk

**HAMSTER** is an open source software framework to trace heat and moisture through the atmosphere and establish (bias-corrected) source–receptor relationships, using output from a Lagrangian model. It has been developed within the DRY-2-DRY project at the Hydro-Climatic Extremes Laboratory (H-CEL) at Ghent University. 

- - - -
## Getting started. 

This section describes the prerequisites required to run HAMSTER, as well as the steps to install it. 

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

- - - - 
## HAMSTER: modules.
**HAMSTER** consists of 4 modules, 

0. flex2traj
1. Diagnosis
2. Attribution
3. Bias-correction

which build up on each other. It is suggested to run them sequentially to obtain the most efficient and informative workflow. 

### 0. flex2traj
This module of **HAMSTER** reads in the instantaneous binary FLEXPART files, filters for a specific region (using a netcdf mask), constructs trajectories and writes them to a file.

### 1. Diagnosis
The diagnosis part of **HAMSTER** identifies atmospheric fluxes of humidity (precipitation and evaporation) or heat (sensible heat flux) using trajectories constructed from FLEXPART binary data. There are several thresholds and criteria that can be set (see docs) to reduce the bias, increase the probability of detection and reduce the probability of false detection. The output from this part can be used to bias correct source–receptor relationships. 

### 2. Attribution
The attribution part of **HAMSTER** constructs mass- and energy-conserving trajectories of heat and moisture (e.g. using a linear discounting of changes en route, or applying a random attribution for moisture), and establishes a first (biased) source–receptor relationship. Multiple options to ensure mass- and energy conservation along trajectories are available (see docs). Various time and space-scales for attribution are possible (see docs). 

### 3. Bias-correction
The last module of **HAMSTER** uses information from the former two modules to bias-correct source–receptor relationships. Multiple options for bias-correction are available (see docs). 

## Running HAMSTER.
### Prerequisites
To execute the full chain (all 3 modules) of **HAMSTER**, the only prerequisites are: 
* Output from a Lagrangian model that traces air parcels and their properties (driven with a reanalysis or output from a GCM/RCM)
* Benchmarking data; e.g., the reanalysis used to run FLEXPART and track parcels
* A file paths.txt which lists the paths where the above data is found and where output will be stored.

The file paths.txt is not part of **HAMSTER**. Users have to create the file themselves. The order in this file is arbitrary, but it has to contain paths for diagnosis, attribution and biascorrection and reference (benchmark) data: 
```
# This file contains all required paths and file names to run hamster; the order doesn't matter and paths can also be empty (if, e.g., not used)

# MASK
maskfile  = "./flexpart_data/masks/mask.nc"

# INPUT paths
ipath_f2t = "/user/data/gent/gvo000/gvo00090/EXT/data/FLEXPART/era_global/orig_untar" 
ipath_DGN = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t2"
ipath_ATR = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t62/MON"
ipath_REF = "/data/gent/vo/000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/1x1"

# INPUT file name base
ibase_DGN = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"]
ibase_ATR = ["pom_ecoreg1_traj10d_AUXTRAJ_"]

# OUTPUT paths
opath_f2t = "./flexpart_data/hamster/00_eraglobal"
opath_DGN = "./flexpart_data/hamster/01_diagnosis"
opath_ATR = "./flexpart_data/hamster/02_attribution"
opath_BIA = "./flexpart_data/hamster/03_biascorrection"
```

The sample paths provided here are (mostly) accessible for members of the virtual organization (VO00090) from H-CEL at the HPC @ Gent. Note, however, that the binary FLEXPART data needs to be untarred. 

### Run and settings.
To run **HAMSTER**, run
```
python main.py
```

Note that — without any flags — main.py is run with default values. Use 
```
python main.py -h
```
for more details on setting dates, thresholds and other options. All user-specific paths are set in paths.txt. 


#### The most important settings are: 

- `--steps` to select the part of hamster that is being executed (e.g., `--steps 1` runs the diagnosis, `--steps 2` performs the attribution, ...)
- `--ayyyy` and `--am` to select the analysis year and month (e.g., `--ayyyy 2002 --am 1`)
- `--tdiagnosis` to select a diagnosis methodology with (default) criteria (e.g., `--tdiagnosis SOD` for constant thresholds of moisture for evaporation and potential temperature for sensible heat, restricted to the planetary boundary layer (PBL), following Sodemann et al. (2008) and Schumacher et al. (2019); `--tdiagnosis SOD2` to use all moisture uptakes above a threshold to detect and quantify evaporation, following Sodemann (2020); `--tdiagnosis KAS` to select a Clausius-Clapeyron dependent threshold for evaporation and sensible heat fluxes, i.e. define moisture thresholds based on the temperature and temperature-thresholds based on the moisture)
- `--expid` to name a setting (e.g., `--expid Ghent_SOD`)
- `--ctraj_len` to determine the maximum length of a trajectory for evaluation (e.g., `--ctraj_len 15` to select 15 days)
- `--mattribution` to determine the attribution method for precipitation (e.g., `--mattribution random` uses the random attribution to attribute moisture for precipitation along trajectories – it keeps linear discounting for heat though)
- `--maskval` to filter for a value other than 1 using the maskfile from `paths.txt`(e.g., `--maskval 5001`)


#### Fine-tuning of detection criteria can be done using, e.g., 
- `--cpbl_strict` to determine the 'strictness' of the PBL criteria (`--cpbl_strict 1` requires both instances to be within the maximum PBL, `--cpbl_strict 2` requires only one instance to be within the maximum PBL; `--cpbl_strict 0` does not filter for the PBL at all)
- `-–cprec_dqv` and `–-cprec_rh` to adjust the detection of preciptation
- `--cevap_c` and `--cheat_cc` to adjust the Clausius-Clapeyron criteria for evaporation and sensible heat, respectively
- `--cevap_hgt`, etc., to filter for specific heights
- `--fjumps` and `--fjumpsfull` to filter for jumps larger than `--cjumps` at the beginning of the trajectory or along the full trajectory
- ... among others. 


#### A few more notes on flags...
- Short flags available! See `python main.py -h` for details (e.g., `-–ayyyy`can be replaced with `-ay` and `--tdiagnosis` can be replaced with `-dgn`). 
- Analysis is performed on a monthly basis: for an independent analysis of months, the flag `--memento` is incorporated (default: True) and requires additional data for the previous month in 02_attribution. 
- The `expid` has to be used consistently for the settings between steps 1-2-3. Otherwise, source-sink relationships may be bias-corrected with other criteria (DANGER!). There is no proper check for this – the user has to make sure they are using everything correctly. Various regions or attribution methods can be run using separate directories. 

## Miscellaneous notes
- Everything is more or less hard-coded for (global) FLEXPART–ERA-Interim simulations with a 6-hourly time step and a maximum of ~2 million parcels. Any changes in resolution or input data require code adjustments!
- In this context, the bias correction is currently implemented for the driving ERA-Interim data only (again, using a hard-coded structure of that data). This data can, however, be easily substituted with other data sets, but it requires changes in the code. 
- 'flex2traj' is the python replacement for *particle-o-matic*. Note that 'flex2traj' is currently under development and not fully integrated yet. All other modules (01_diagnosis and 02_attribution) currently only read data from *particle-o-matic* (dat-files). Once fully integrated, `ipath_ATR` should be set identical to `opath_f2t`.
- Directories are currently assumed to have an annual structure (e.g., ipath_ATR + "/2002")
- The 'minimum' time scale for steps 1-2-3 is daily, which we assumed to be a reasonable limit for the FLEXPART–ERA-Interim simulations with 6-hourly time steps. This could be adjusted and tested though...  


## Epilogue
Keep in mind that... 
- **This code is not bug-free.** Please report any bugs through 'Issues' on https://github.ugent.be/jkeune/hamster/issues. 
- **This code is not intended to cover specific research-related tasks.** This code is intended to serve as a common base for the analysis of (FLEXPART) trajectories. Every user may create their own branch and adjust the code accordingly. Features of broad interest may be merged in future releases. 

### Contact and support
Dominik Schumacher and Jessica Keune

### License
Copyright 2019 Dominik Schumacher, Jessica Keune, Diego Miralles. 

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

