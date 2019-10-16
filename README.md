# hamster

**HAMSTER - A Heat And MoiSturE Tracking fRamework**

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
## Run HAMSTER.
**HAMSTER** consists of 3 modules, 
1. Diagnosis
2. Attribution
3. Bias-correction

which build up on each other. It is suggested to run them sequentially to obtain the most efficient and informative workflow. 

### 1. Diagnosis
The diagnosis part of **HAMSTER** identifies atmospheric fluxes of humidity (precipitation and evaporation) or heat (sensible heat flux) using output from a Lagrangian model. There are several thresholds and criteria that can be set (see docs) to reduce the bias, increase the probability of detection and reduce the probability of false detection. The output from this part can be used to bias correct source-sink relationships. 

There are two options to perform the diagnosis for a single month, 
1. Run the python script interactively, use
```
python main.py
```
to run the script. 
2. Run the script as a job on a cluster. use
```
qsub job.pbs 
```
to submit the job, which executes main.py
Note that - without any flags - main.py is run with default values. Use 
```
python main.py -h
```
for more details on setting dates, thresholds and other options. 

The modular setup of main.py enables to run multiple runs in parallel using the same script (prerequisite: worker environment)
- [ ] adjust settings in work/dates.txt
- [ ] adjust the work environment settings in work/workerjob.pbs
- [ ] and submit the jobs. 
```
vi work/dates.py            # -- job settings (dates, thresholds, options, ...)
vi work/workerjob.pbs       # -- job settings
cd work 
wsub workerjob.pbs
``` 

### 2. Attribution
The attribution part of **HAMSTER** constructs mass- and energy-conserving trajectories of heat and moisture (e.g. using a linear discounting of changes en route), and establishes a first (biased) source-sink relationship. Multiple options to ensure mass- and energy conservation along trajectories are available (see docs). Various time and space-scales for attribution are possible (see docs). 

### 3. Bias-correction
The last module of **HAMSTER** uses information from the former two modules to bias-correct source-receptor relationships. Multiple options for bias-correction are available (see docs). 
