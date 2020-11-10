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
To execute the full chain (all 4 modules) of **HAMSTER**, the only prerequisites are: 
* Output from a Lagrangian model that traces air parcels and their properties (driven with a reanalysis or output from a GCM/RCM)
* Benchmarking data; e.g., the reanalysis used to run FLEXPART and track parcels
* A file paths.txt which lists the paths where the above data is found and where output will be stored.

The file paths.txt is not part of **HAMSTER**. Users have to create the file themselves. The order in this file is arbitrary, but it has to contain paths for diagnosis, attribution and biascorrection and reference (benchmark) data: 
```
# This file contains all required paths and file names to run hamster; the order doesn't matter and paths can also be empty (if, e.g., not used)

# MASK
maskfile  = "./flexpart_data/masks/mask.nc"

# INPUT paths
ipath_f2t = "./flexpart_data/orig"
ipath_DGN = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/flex2traj_t2"
ipath_ATR = ".flexpart_data/00_flex2traj/myregion"
ipath_REF = "/data/gent/vo/000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/1x1"

# INPUT file name base
ibase_f2t = "bahamas"
ibase_DGN = "global"
ibase_ATR = "bahamas"

# OUTPUT paths
opath_f2t = "./flexpart_data/hamster/00_eraglobal"
opath_DGN = "./flexpart_data/hamster/01_diagnosis"
opath_ATR = "./flexpart_data/hamster/02_attribution"
opath_BIA = "./flexpart_data/hamster/03_biascorrection"
```

The sample paths provided here are (mostly) accessible for members of the virtual organization (VO00090) from H-CEL at the HPC @ Gent. Note, however, that the binary FLEXPART data needs to be untarred from the archive.

### Run and settings.
To run **HAMSTER**, run
```python
python main.py
```

Note that — without any flags — main.py is run with default values. Use 
```python
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
- ... among a lot of other options. 


#### A few more notes on flags...
- Short flags available! See `python main.py -h` for details (e.g., `-–ayyyy`can be replaced with `-ay` and `--tdiagnosis` can be replaced with `-dgn`). 
- Analysis is performed on a monthly basis: for an independent analysis of months, the flag `--memento` is incorporated (default: True) and requires additional data for the previous month in 02_attribution. 
- The `expid` has to be used consistently for the settings between steps 1-2-3. Otherwise, source-sink relationships may be bias-corrected with other criteria (DANGER!). There is no proper check for this – the user has to make sure they are using everything correctly. Various regions or attribution methods can be run using separate directories. 
- There are quite a few flags for 02_attribution (e.g., refering to settings concerning the random attribution) and 03_biascorrection (e.g., refering to the applied time scale and the aggregation of the output) available. Please use the help option for details for now. 
- While the output of flex2traj could be adjusted through modifications in 00_flex2traj.py, currently, all other steps require the following 9 variables (and in that specific order): `parcel id`, `lon`, `lat`, `ztra1`, `topo`, `qvi`, `rhoi`, `hmixi`, `tti`.
- If `--writestats True` is set for `--steps 2`, then the attribution statistics are written to a file `*_stats.csv` (absolute fraction of attributed precipitation, etc.). If `--writestats True` is set for `--steps 3`, then the validation statistics are written to a file `*_stats.csv` (bias in the sink region, the probability of detection etc.).  
- Use `--maskval -999` (or set maskfile="" in paths.txt) in combination with `--ctraj_len 0` to extract global 2-step trajectories for a global 'diagnosis' with flex2traj (data already available on the HPC).

#### One important note. 
- The mask for 00_flextraj has to be bigger than the mask for 02_attribution. This is especially important for precipitation, as we use the midpoint to detect precipitation over a certain region. If the same mask is used, your results of 02_attribution and 03_biascorrection won't be representative. 

#### A very basic example. 
1. Create a (global) netcdf file with a mask (value=1) for a specific region of interest, e.g., the Bahamas. Make two masks: one that represents the area you're really interested in (for the attribution), and one that is +4-5° bigger on each side (used for flex2traj). 
2. Adjust the maskfile in `paths.txt` (first, for flex2traj, i.e. the extended mask). 
3. Construct trajectories arriving at the Bahamas (don't forget to untar the binary FLEXPART simulations for this and the previous month):
  ```python
  python main.py --steps 0 --ayyyy 2000 --am 6 --ctraj_len 16 --maskval 1 
  ```
4. Perform a global analysis of fluxes (and the previous month), and evaluate the bias and the reliability of detection for your region of interest and its (potential) source region, possibly selecting various diagnosis methods and fine tuning detection criteria (using the already available global data set on the VO), e.g.,  
  ```python
  python main.py --steps 1 --ayyyy 2000 --am 6 --tdiagnosis SOD --cprec_rh 70 --expid "SOD_prh-70"
  ...
  python main.py --steps 1 --ayyyy 2000 --am 6 --tdiagnosis KAS --cprec_rh 70 --cpbl_strict 2 --cevap_cc 0.9 --expid "KAS_prh70_cpbl2_cevapcc0.9"
  ```
5. Adjust the mask in `paths.txt` to your region of interest (the smaller one) to continue with the attribution.
6. Once you have fine-tuned your detection criteria, perform a first backward analysis considering a trajectory length of 15 days, e.g.
  ```python
  python main.py --steps 2 --ayyyy 2000 --am 6 --tdiagnosis KAS --cprec_rh 70 --cpbl_strict 2 --cevap_cc 0.9 --ctraj_len 15 --expid "KAS_prh70_cpbl2_cevapcc0.9"
  ```
7. Bias-correct the established source and aggregate the results over the backward time dimension
  ```python
  python main.py --steps 3 --ayyyy 2000 --am 6 --expid "KAS_prh70_cpbl2_cevapcc0.9" --bc_aggbwtime True
  ```
  The final netcdf file, `KAS_prh70_cpbl2_cevapcc0.9_biascor-attr_r02_2002-06.nc` then contains all the source regions of heat and precipitation, both the raw and bias-corrected version (i.e., Had and Had_Hs, and E2P, E2P_Es, E2P_Ps, and E2P_EPs).  


## Miscellaneous notes
- Everything is coded for a **backward** analysis (Where does the heat come from? What is the source region of precipitation?). Adjustments for a forward analysis can be easily made, but require code changes.
- Note that, however, flex2traj writes out data in a forward format (startdate --> enddate; but still filtering for the last step, i.e. in a backward manner), but that the time axis is swapped when reading this data in (enddate <-- startdate, see function `readtraj`) for all the remaining analysis steps.
- Everything is more or less hard-coded for (global) FLEXPART–ERA-Interim simulations with a 6-hourly time step and a maximum of ~2 million parcels. Any changes in resolution or input data require code adjustments!
- The bias correction is currently implemented for the driving ERA-Interim data only (again, using a hard-coded structure of that data). This data can, however, be easily substituted with other data sets, but it requires changes in the code. 
- Note that regardless of the sink region size, 'flex2traj' reads in and temporarily stores data from all parcels during the backward analysis time period; in case of 15-day trajectories and 9 variables of interest, this translates to a numpy array with a size of ~ 7.2 GB (62 x 2e6 x 9 x 64 bit). For a small sink region with ~13'000 parcels (trajectory array: 62 x 13'000 x 9 x 64 bit ~ 0.5 GB), a total of 10 GB RAM is recommended to safely run flex2traj with a trajectory length of 15 days.
- In paths.txt, use a single string to describe the filename base for `ibase_DGN`, `ìbase_ATR`, `ìbase_f2t` (before, `ibase_DGN` and `ibase_ATR` used to be lists of string and had to have underscores at the end --> this is not needed anymore!) -- so you can only read in **one** file from now on. 
- flex2traj-related directories are currently assumed to have an annual structure (e.g., ipath_ATR + "/2002") - these are created automatically.
- The 'minimum' time scale for steps 1-2-3 is daily, which we assumed to be a reasonable limit for the FLEXPART–ERA-Interim simulations with 6-hourly time steps. This could be adjusted and tested though...  
- An additional file `*_warning.txt` is written, if a monthly bias-correction was required and daily data cannot be trusted (this is the case if, e.g., the reference data set contains precipitation for a specific day, but precipitation was not detected using FLEXPART and the selected detection criteria; and hence no trajectories were evaluated and no attribution for that specific day was performed, but the contribution of other precipitation days was upscaled to match the monthly precipitation amount). 

## Epilogue
Keep in mind that... 
- **This code is not bug-free.** Please report any bugs through 'Issues' on https://github.ugent.be/jkeune/hamster/issues. 
- **This code is not intended to cover specific research-related tasks.** This code is intended to serve as a common base for the analysis of (FLEXPART) trajectories. Every user may create their own branch and adjust the code accordingly. Features of broad interest may be merged in future releases. 

### Contact and support
Dominik Schumacher and Jessica Keune

### License
Copyright 2019 Dominik Schumacher, Jessica Keune, Diego G. Miralles. 

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

