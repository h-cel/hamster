# HAMSTER – a Heat And MoiSture Tracking framEwoRk

**HAMSTER** is an open source software framework to trace heat and moisture through the atmosphere and establish (bias-corrected) source–receptor relationships, using output from a Lagrangian model. It has been developed within the DRY-2-DRY project at the Hydro-Climatic Extremes Laboratory (H-CEL) at Ghent University. 

- - - -
## Getting started. 

This section describes the prerequisites required to run HAMSTER, as well as the steps to install it. 

### Prerequisites
To run HAMSTER, you need 
* python 3
* [git](https://git-scm.com/)

and
* [anaconda](https://www.anaconda.com/) (or similar to manage python packages)

or
* python 3 and the required modules on a cluster

### Installation
1. Clone the repository
    ```
    git clone https://github.ugent.be/jkeune/hamster.git
    cd hamster
    ```
2. Load the following python 3.7 environment: 
    ```
    module load h5py/2.10.0-intel-2019b-Python-3.7.4
    module load netcdf4-python/1.5.3-intel-2019b-Python-3.7.4
    ```
Alternatively, make an anaconda environment with the necessary python packages
    ```
    conda create -n _newenvironment_ --file requirements.txt
    ```
or install the packages listed in requirements.txt in your local environment. Note, however, that due to different versions, some errors might occur. It is thus recommended to load preinstalled environments, such as the one above. 

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

# location of original flexpart files (untarred)
## (untar the original FLEXPART partposit_* files to this directory)
path_orig = "./flexpart_data/orig"

# location of the reference data used for bias correction (e.g., ERA-Interim)
## (available on the VSC Tier-2 HPC clusters)
path_ref  = "/data/gent/vo/000/gvo00090/EXT/data/ERA-INTERIM/by_var_nc/1x1"

# path and base name for global parcel diag data (t2)
## (available on the VSC Tier-2 HPC clusters)
base_f2t_diag = "global"
path_f2t_diag = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/flex2traj_t2"

# path and base name for parcel trajectory data
base_f2t_traj = "bahamas_15d"
path_f2t_traj = "./flexpart_data/bahamas/f2t_traj"

# paths for processed data
path_diag = "./flexpart_data/bahamas/diag"
path_attr = "./flexpart_data/bahamas/attr"
path_bias = "./flexpart_data/bahamas/bias"
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

- `--steps` to select the part of hamster that is being executed (e.g., `--steps 0` runs flex2traj, `--steps 1` runs the diagnosis, `--steps 2` performs the attribution, ...)
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
- Analysis is performed on a monthly basis: for an independent analysis of months, the flag `--memento` is incorporated (default: True). This flag requires additional data of the previous month, that is extracted from the trajectories or, if not available from the trajectories, read in from the binary FLEXPART files in a preloop. Hence, the smoothest workflow is obtained if flex2traj dumps trajectories that are one day longer than what is needed in 02_attribution (e.g., run `python main.py --steps 0 --ctraj_len 16` but `python main.py --steps 0 --ctraj_len 15` -- the preloop is not needed in this case). 
- The `expid` has to be used consistently for the settings between steps 1-2-3. Otherwise, source-sink relationships may be bias-corrected with other criteria (DANGER!). There is no proper check for this – the user has to make sure they are using everything correctly. Various regions or attribution methods can be run using separate directories. 
- There are quite a few flags for 02_attribution (e.g., refering to settings concerning the random attribution) and 03_biascorrection (e.g., refering to the applied time scale and the aggregation of the output) available. Please use the help option for details for now. 
- While the output of flex2traj could be adjusted through modifications in 00_flex2traj.py, currently, all other steps require the following 9 variables (and in that specific order): `parcel id`, `lon`, `lat`, `ztra1`, `topo`, `qvi`, `rhoi`, `hmixi`, `tti`.
- If `--writestats True` is set for `--steps 2`, then the attribution statistics are written to a file `*_stats.csv` (absolute fraction of attributed precipitation, etc.). If `--writestats True` is set for `--steps 3`, then the validation statistics are written to a file `*_stats.csv` (bias in the sink region, the probability of detection etc.).  
- Use `--maskval -999` (or set maskfile=None in paths.txt) in combination with `--ctraj_len 0` to extract global 2-step trajectories for a global 'diagnosis' with flex2traj (data already available on the HPC).
- To filter for all PBL processes without any thresholds, set `--fallingdry False`and a high value for the CC criterion, e.g. `--cevap_cc 1000`. If, in addition `--cpbl_strict 0`, then also above PBL parcels are evaluated.

#### A very basic example. 
1. Create a (global) netcdf file with a mask (value=1) for a specific region of interest, e.g., the Bahamas. 
2. Adjust the maskfile in `paths.txt`.
3. Construct trajectories for parcels whose arrival+midpoints are over the Bahamas (don't forget to untar the binary FLEXPART simulations for this and the previous month and the first day of the following month):
  ```python
  python main.py --steps 0 --ayyyy 2000 --am 6 --ctraj_len 16 --maskval 1 
  ```
4. Perform a global analysis of fluxes (and the previous month), and evaluate the bias and the reliability of detection for your region of interest and its (potential) source region, possibly selecting various diagnosis methods and fine tuning detection criteria (using the already available global data set on the VO), e.g.,  
  ```python
  python main.py --steps 1 --ayyyy 2000 --am 6 --tdiagnosis SOD --cprec_rh 70 --expid "SOD_prh-70"
  ...
  python main.py --steps 1 --ayyyy 2000 --am 6 --tdiagnosis KAS --cprec_rh 70 --cpbl_strict 2 --cevap_cc 0.9 --expid "KAS_prh70_cpbl2_cevapcc0.9"
  ```
5. Once you have fine-tuned your detection criteria, perform a first backward analysis considering a trajectory length of 15 days, e.g.
  ```python
  python main.py --steps 2 --ayyyy 2000 --am 6 --tdiagnosis KAS --cprec_rh 70 --cpbl_strict 2 --cevap_cc 0.9 --ctraj_len 15 --expid "KAS_prh70_cpbl2_cevapcc0.9"
  ```
6. Bias-correct the established source and aggregate the results over the backward time dimension
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
- flex2traj-related directories are currently assumed to have an annual structure (e.g., path_f2t_diag + "/2002") - these are created automatically.
- The 'minimum' time scale for steps 1-2-3 is daily, which we assumed to be a reasonable limit for the FLEXPART–ERA-Interim simulations with 6-hourly time steps. This could be adjusted and tested though...  
- An additional file `*_warning.txt` is written, if a monthly bias-correction was required and daily data cannot be trusted (this is the case if, e.g., the reference data set contains precipitation for a specific day, but precipitation was not detected using FLEXPART and the selected detection criteria; and hence no trajectories were evaluated and no attribution for that specific day was performed, but the contribution of other precipitation days was upscaled to match the monthly precipitation amount). 

## Known errors. 
- **Remote I/O error** when reading/writing a netcdf or h5 file: this only seems to happen with specific h5py / netCDF4 module versions (used in anaconda) and only on some clusters (victini). If you're on the UGent HPC cluster, please try to use the HPC modules listed above - then the error should disappear... 
- **Resolutions other than 1°**: the current version only supports grids of 1 degree - for ALL data (mask, reference data, hamster output) - the code should abort if that is not the case, stating that the grids are not identical, but it may not necessarily do so... so check your outputs carefully, if its runs through anyhow (at least step 3 - bias correction - should be erroneous!)

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

