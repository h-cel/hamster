# HAMSTER – a Heat And MoiSture Tracking framEwoRk

**HAMSTER** is an open source software framework to trace heat and moisture through the atmosphere and establish (bias-corrected) source–receptor relationships, using output from a Lagrangian model. It has been developed within the DRY-2-DRY project at the Hydro-Climatic Extremes Laboratory (H-CEL) at Ghent University. 

[![DOI](https://zenodo.org/badge/372759352.svg)](https://zenodo.org/badge/latestdoi/372759352)

- - - -
## What is HAMSTER?
**HAMSTER** is a heat- and moisture tracking framwork to evaluate air parcel trajectories from a Lagrangian model, such as FLEXPART (Stohl et al., 2005) and to establish source–receptor relationships, such as the source regions of precipitation or heat. The current version of **HAMSTER** has been built using simulations from FLEXPART driven with ERA-Interim reanalysis data, but other simulations may be supported as well. 

**HAMSTER** consists of 4 modules, 

0. flex2traj
1. Diagnosis
2. Attribution
3. Bias-correction

which build up on each other. It is suggested to run them sequentially to obtain the most efficient and informative workflow. In the following, a short description of each module (step) is given. 

### 0. flex2traj
This module of **HAMSTER** reads in the instantaneous binary FLEXPART files, filters for a specific region (using a netcdf mask), constructs trajectories and writes them to a file.

### 1. Diagnosis
The diagnosis part of **HAMSTER** identifies atmospheric fluxes of humidity (precipitation and evaporation) or heat (sensible heat flux) using trajectories constructed from FLEXPART binary data. There are several thresholds and criteria that can be set to reduce the bias and to increase the probability of detection for each flux. The output from this diagnosis step is used to bias correct source–receptor relationships. 

### 2. Attribution
The attribution part of **HAMSTER** constructs mass- and energy-conserving trajectories of heat and moisture (e.g. using a linear discounting of changes en route, or applying a random attribution for moisture), and establishes a first (biased) source–receptor relationship. The output of this step are, e.g., 4D netcdf files that illustrate the spatio-temporal origins of precipitation or *heat* arriving in the receptor region. Multiple options to construct these relationships are available.

### 3. Bias-correction
The last module of **HAMSTER** uses information from the former two steps to bias-correct source–receptor relationships. Multiple options for bias-correction are available. 

- - - -
## What do I need to get and run HAMSTER?

### Prerequisites and Installation
This section describes the prerequisites required to run HAMSTER, as well as the steps to install it. 

#### Prerequisites
To run HAMSTER, you need 
* python 3
* [git](https://git-scm.com/)

and
* [anaconda](https://www.anaconda.com/) (or similar to manage python packages)

or
* python 3 and the required modules on a cluster

The main packages required to run **HAMSTER** are: 
```bash
import argparse
import calendar
import csv
import fnmatch
import gzip
import imp
import math
import os
import random
import re
import struct
import sys
import time
import timeit
import warnings
import datetime
import functools

import h5py
import netCDF4
import numpy
import pandas
import dateutil
```

In addition, to execute the full chain (all 4 modules) of **HAMSTER**, the following data sets are needed: 
* **Output from a Lagrangian model** that traces air parcels and their properties (driven with a reanalysis or output from a GCM/RCM)
* **Reference data**; e.g., the reanalysis used to run FLEXPART and track parcels, or any other reference data set used for bias-correction

To run **HAMSTER**, you will also need to create one file:
* `paths.txt` 

which lists the paths where the above data is found and where output will be stored.

The file `paths.txt` is not part of **HAMSTER**. The order in this file is arbitrary, but it has to contain paths for diagnosis, attribution and biascorrection and reference data: 
```
# This file contains all required paths and file names to run hamster; the order doesn't matter and paths can also be empty (if, e.g., not used)

# MASK
maskfile  = "./flexpart_data/masks/mask.nc"

# location of original flexpart files (untarred)
## (untar the original FLEXPART partposit_* files to this directory)
path_orig = "./flexpart_data/orig"

# location of the reference data used for bias correction (e.g., ERA-Interim)
## for each variable (P, E, H)
path_ref_p  = "./ERA-INTERIM/1x1/tp_12hourly"
path_ref_e  = "./ERA-INTERIM/1x1/evap_12hourly"
path_ref_h  = "./ERA-INTERIM/1x1/sshf_12hourly"

# path and base name for global parcel diag data (t2)
base_f2t_diag = "global"
path_f2t_diag = "./flexpart_data/global/f2t_diag"

# path and base name for parcel trajectory data
base_f2t_traj = "bahamas_10d"
path_f2t_traj = "./flexpart_data/bahamas/f2t_traj"

# paths for processed data
path_diag = "./flexpart_data/global/diag"
path_attr = "./flexpart_data/bahamas/attr"
path_bias = "./flexpart_data/bahamas/bias"
```

#### Installation
To install **HAMSTER**, do the following:

1. Clone the repository
    ```
    git clone https://github.com/h-cel/hamster
    cd hamster
    ```
2. Make an anaconda environment with the necessary python packages
    ```
    conda create -n _newenvironment_ --file requirements.txt
    ```
or install the packages listed in requirements.txt in your local environment. Note, however, that due to different versions, some errors might occur. It is thus recommended to load preinstalled environments, such as the one above. 


- - - -
## How do I run HAMSTER?

To run **HAMSTER**, change into the `src` directory
```
cd src
```
and run
```python
python main.py
```

Note that — without any flags — main.py is run with default values. Use 
```python
python main.py -h
```
for more details on setting dates, thresholds and other options. All user-specific paths are set in `paths.txt`. 

### Quick start.
To run all steps sequentially with the default settings, please proceed as follows. First, extract global 2-step trajectories (`steps--0  --maskval -999 --ctraj_len 0`) and perform the global diagnosis of fluxes (`--steps 1`), using
```
python main.py --steps 0 --maskval -999 --ctraj_len 0
python main.py --steps 1
```
Then, extract 10-day trajectories for a specific region (using a default maskvalue of 1 for the given maskfile; `--steps 0`), which are required to diagnose source regions of heat and moisture (`--steps 2`), and finally bias-correct these source regions using the global diagnosis data from above (`--steps 3`):
```
python main.py --steps 0
python main.py --steps 2
python main.py --steps 3
```

**The most important settings are**:

- `--steps` to select the part of hamster that is being executed (e.g., `--steps 0` runs flex2traj, `--steps 1` runs the diagnosis, `--steps 2` performs the attribution, ...)
- `--ayyyy` and `--am` to select the analysis year and month (e.g., `--ayyyy 2002 --am 1`)
- `--expid` to name a setting (e.g., `--expid "ALL-ABL"`)
- `--ctraj_len` to determine the maximum length of a trajectory for evaluation (e.g., `--ctraj_len 15` to select 15 days; 10 is the default)
- `--maskval` to filter for a value other than 1 using the maskfile from `paths.txt`(e.g., `--maskval 5001`)

- - - -
## What flags can I set? 

Here, we provide a short description of all flags that be used when running **HAMSTER**. 

#### The detection of fluxes (P, E and H) is set via a set of flags.

The detection of **precipitation** is set via `-–cprec_dqv` and `–-cprec_rh`:
- `-–cprec_dqv` to set the minimum loss of specific humidity (unit: kg kg-1; default: 0 kg kg-1)
- `–-cprec_rh` to the minimum average relative humidity (unit: %; default: 80% following the convection scheme from Emanuel, 1991)

The detection of **evaporation** is set via `--cevap_dqv`, `--fevap_drh`, `--cevap_drh`, `--cevap_hgt`:
- `--cevap_dqv` to set a minimum increase in specific humidity (unit: kg kg-1; default: 0.2 kg kg-1)
- `--cevap_hgt` to filter for specific heights (unit: m; default: 0)
- `--fevap_drh` to filter for relative humidity changes (False/True; default: False) employing a maximum change of `--cevap_drh` (unit: %, default: 15%)
- `--fallingdry` to shorten the trajectory length if a parcel *falls dry* (q<0.00005 kg kg-1; default: False)

The detection of **sensible heat fluxes** is set via `--cheat_dtemp`, `--cheat_hgt`, `--fheat_drh`, `--cheat_drh`, `--fheat_rdq` and `--cheat_rdq`:
- `--cheat_dtemp` to set a minimum increase in potential temperature (unit: K; default: 1 K)
- `--cheat_hgt` to filter for specific heights (unit: m; default: 0)
- `--fheat_drh` to filter for relative humidity changes (False/True; default: False) employing a maximum change of `--cheat_drh` (unit: %, default: 15%)
- `--fheat_rdq` to filter for relative specific humidity changes (False/True; default: False) employing a maximum change of `--cheat_rdq` (expressed as delta(qv)/qv; unit: %; default: 10%)

For E and H, the detection of fluxes can be limited to the atmospheric boundary layer (ABL):
- `--cpbl_method` describes the method used to determine the ABL height between two instances. Options are `actual`, `mean` and `max` (default: `max`).
- `--cpbl_factor` sets a factor that is used to increase (>1) or decrease (<1) the ABL heights (default: 1). 
- `--cpbl_strict` determines the 'strictness' of the ABL criteria (`--cpbl_strict 2` requires both instances to be within the actual/maximum/mean ABL, `--cpbl_strict 1` requires only one instance to be within the actual/mean/maximum ABL; `--cpbl_strict 0` does not filter for the ABL at all).

Note that the ABL criteria are set consistently for E and H.

Using these flags, a lot of the settings used in FLEXPART publications can be mimicked. 
- **Sodemann et al., 2008** for the detection of E (minimum threshold for dqv/dt; parcel has be reside in the vicinity of the ABL; note, however, that minor differences exists, e.g. through the application of the ABL factor 1.5 everywhere as in Keune and Miralles (2019), etc.): 
    ```
    --cevap_dqv 0.0002 --fallingdry True --fevap_drh False --cpbl_method "mean" --cpbl_factor 1.5
    ``` 
- **Fremme and Sodemann, 2019** and **Sodemann, 2020** (minimum threshold for dqv/dt; parcel does not have to reside in the ABL):
    ```
    --cevap_dqv 0.0001 --fevap_drh False --cpbl_strict 0
    ``` 
- **Schumacher et al., 2019** for the detection of H (minimum potential temperature increase; limitation by change in specific humidity content; parcel has to be within the maximum ABL at both time steps):
    ```
    --cheat_dtemp 1 --cheat_rdq 10 --fheat_rdq True --fheat_drh False --cpbl_strict 2 --cpbl_method "max" 
    ```


#### A few more notes on flags...
- Short flags available! See `python main.py -h` for details (e.g., `-–ayyyy`can be replaced with `-ay` etc.)
- Bias correction is, per default, performed using ERA-Interim data (i.e., the base for FLEXPART simulations). However, alternative data sets for P, E and/or H can be employed using `--pref_data others`, `--eref_data others` and `--href_data others`, respectively. The data has to be placed in the respective paths of `paths.txt` (`path_ref_p`, `path_ref_e`, `path_ref_h`) and **have to be in the correct format (netcdf), in the correct resolution on the correct (grid set via `--resolution`) with daily (or subdaily) time steps, and in the correct units (mm d-1 and W m-2)**. If `others` is employed, all files with the year (`--ayyyy`) from the respective directory are read in and used for bias correction.
- Analysis is performed on a monthly basis: for an independent analysis of months, the flag `--memento` is incorporated (default: True). This flag requires additional data of the previous month, that is extracted from the trajectories or, if not available from the trajectories, read in from the binary FLEXPART files in a preloop. Hence, the smoothest workflow is obtained if flex2traj dumps trajectories that are one day longer than what is needed in 02_attribution (e.g., run `python main.py --steps 0 --ctraj_len 16` but `python main.py --steps 0 --ctraj_len 15` -- the preloop is not needed in this case). 
- The `expid` has to be used consistently for the settings between steps 1-2-3. Otherwise, source-sink relationships may be bias-corrected with other criteria (DANGER!). There is no proper check for this – the user has to make sure they are using everything correctly. Various regions or attribution methods can be run using separate directories. 
- There are quite a few flags for 02_attribution (e.g., refering to settings concerning the random attribution) and 03_biascorrection (e.g., refering to the applied time scale and the aggregation of the output) available. Please use the help option for details for now. 
- While the output of flex2traj could be adjusted through modifications in 00_flex2traj.py, currently, all other steps require the following 9 variables (and in that specific order): `parcel id`, `lon`, `lat`, `ztra1`, `topo`, `qvi`, `rhoi`, `hmixi`, `tti`.
- If `--writestats True` is set for `--steps 2`, then the attribution statistics are written to a file `*_stats.csv` (absolute fraction of attributed precipitation, etc.). If `--writestats True` is set for `--steps 3`, then the validation statistics are written to a file `*_stats.csv` (bias in the sink region, the probability of detection etc.).  
- Use `--maskval -999` (or set maskfile=None in paths.txt) in combination with `--ctraj_len 0` to extract global 2-step trajectories for a global 'diagnosis' with flex2traj.

- - - -
## An example. 

Here, we provide a very basic example. 

1. Create a (global) netcdf file with a mask (value=1) for a specific region of interest, e.g., the Bahamas. 
2. Adjust the maskfile in `paths.txt`.
3. Construct trajectories for parcels whose arrival+midpoints are over the Bahamas (don't forget to untar the binary FLEXPART simulations for this and the previous month and the first day of the following month):
  ```python
  python main.py --steps 0 --ayyyy 2000 --am 6 --ctraj_len 11 --maskval 1 
  ```
4. Perform a global analysis of fluxes (and the previous month), and evaluate the bias and the reliability of detection for your region of interest and its (potential) source region, possibly selecting various diagnosis methods and fine tuning detection criteria (using the already available global data set on the VO), e.g.,  
  ```python
  python main.py --steps 1 --ayyyy 2000 --am 6 --expid "DEFAULT"
  ...
  python main.py --steps 1 --ayyyy 2000 --am 6 --cpbl_strict 2 --cpbl_method "max" --cevap_dqv 0 --cheat_dtemp 0 --expid "ALL-ABL"
  ```
5. Once you have fine-tuned your detection criteria, perform a first backward analysis considering a trajectory length of 10 days, e.g.
  ```python
  python main.py --steps 2 --ayyyy 2000 --am 6 --ctraj_len 10 --cpbl_strict 2 --cpbl_method "max" --cevap_dqv 0 --cheat_dtemp 0 --expid "ALL-ABL"
  ```
6. Bias-correct the established source and aggregate the results over the backward time dimension
  ```python
  python main.py --steps 3 --ayyyy 2000 --am 6 --expid "ALL-ABL" --bc_aggbwtime True
  ```
  The final netcdf file, `ALL-ABL_biascor-attr_r02_2002-06.nc` then contains all the source regions of heat and precipitation, both the raw and bias-corrected version (i.e., Had and Had_Hs, and E2P, E2P_Es, E2P_Ps, and E2P_EPs).  

- - - -
## Miscellaneous notes
- Everything is coded for a **backward** analysis (Where does the heat come from? What is the source region of precipitation?). Adjustments for a forward analysis can be easily made, but require code changes.
- Note that, however, flex2traj writes out data in a forward format (startdate --> enddate; but still filtering for the last step, i.e. in a backward manner), but that the time axis is swapped when reading this data in (enddate <-- startdate, see function `readtraj`) for all the remaining analysis steps.
- Everything is more or less hard-coded for (global) FLEXPART–ERA-Interim simulations with a 6-hourly time step and a maximum of ~2 million parcels. Any changes in resolution or input data require code adjustments!
- Note that regardless of the sink region size, 'flex2traj' reads in and temporarily stores data from all parcels during the backward analysis time period; in case of 15-day trajectories and 9 variables of interest, this translates to a numpy array with a size of ~ 7.2 GB (62 x 2e6 x 9 x 64 bit). For a small sink region with ~13'000 parcels (trajectory array: 62 x 13'000 x 9 x 64 bit ~ 0.5 GB), a total of 10 GB RAM is recommended to safely run flex2traj with a trajectory length of 15 days.
- flex2traj-related directories are currently assumed to have an annual structure (e.g., path_f2t_diag + "/2002") - these are created automatically.
- The minimum time scale for steps 1-2-3 is daily, which we assumed to be a reasonable limit for the FLEXPART–ERA-Interim simulations with 6-hourly time steps. 
- An additional file `*_warning.txt` is written, if a monthly bias-correction was required and daily data cannot be trusted (this is the case if, e.g., the reference data set contains precipitation for a specific day, but precipitation was not detected using FLEXPART and the selected detection criteria; and hence no trajectories were evaluated and no attribution for that specific day was performed, but the contribution of other precipitation days was upscaled to match the monthly precipitation amount). 

- - - -
## Known issues. 
- **Remote I/O error** when reading/writing a netcdf or h5 file: this only seems to happen with specific h5py / netCDF4 module versions (used in anaconda) and only on some clusters (victini). If you're on the UGent HPC cluster, please try to use the HPC modules listed above - then the error should disappear... 
- **Resolutions other than 1°**: the current version only supports grids of 1 degree - for ALL data (mask, reference data, hamster output) - the code should abort if that is not the case, stating that the grids are not identical, but it may not necessarily do so... so check your outputs carefully, if its runs through anyhow (at least step 3 - bias correction - should be erroneous!)

- - - -
## Epilogue
Keep in mind that... 
- **This code is not bug-free.** Please report any bugs through 'Issues' on https://github.com/h-cel/hamster/issues. 
- **This code is not intended to cover specific research-related tasks.** This code is intended to serve as a common base for the analysis of (FLEXPART) trajectories. Every user may create their own branch and adjust the code accordingly. Features of broad interest may be merged in future releases. 

### Contact and support
Dominik Schumacher (dominik.schumacher@ugent.be) and Jessica Keune (jessica.keune@ugent.be)

### Referencing
If you use HAMSTER, please cite:
Keune, J., Schumacher, D. L., & Miralles, D. G. (2022). A unified framework to estimate the origins of atmospheric moisture and heat using Lagrangian models. Geoscientific Model Development, 15, 1875–1898. https://doi.org/10.5194/gmd-15-1875-2022

### License
Copyright 2021 Dominik Schumacher, Jessica Keune, Diego G. Miralles. 

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

## References
- Fremme, A. and Sodemann, H.: The role of land and ocean evaporation on the variability of precipitation in the Yangtze River valley, Hydrol. Earth Syst. Sci., 23, 2525–2540, https://doi.org/10.5194/hess-23-2525-2019, 2019.
- Keune, J., and Miralles, D. G.: A precipitation recycling network to assess freshwater vulnerability: Challenging the watershed convention, Water Resour. Res., 55, 9947– 9961, https://doi.org/10.1029/2019WR025310, 2019.
- Schumacher, D.L., Keune, J., van Heerwaarden, C.C. et al.: Amplification of mega-heatwaves through heat torrents fuelled by upwind drought, Nat. Geosci. 12, 712–717, https://doi.org/10.1038/s41561-019-0431-6, 2019.
- Schumacher, D.L., Keune, J. and Miralles, D.G.: Atmospheric heat and moisture transport to energy- and water-limited ecosystems. Ann. N.Y. Acad. Sci., 1472, 123-138, https://doi.org/10.1111/nyas.14357, 2020.
- Sodemann, H., Schwierz, C., and Wernli, H.: Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res. Atmos., 113(D3), https://doi.org/10.1029/2007JD008503, 2008.
- Sodemann, H.: Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, J. Atmos. Sci., 77(2), 413-433, https://doi.org/10.1175/JAS-D-18-0336.1, 2020.

