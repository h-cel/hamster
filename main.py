#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN SCRIPT MOISTURE AND HEAT DIAGNOSIS
@author: dominik and jessica

To execute interactively: 
> exec(open("./main.py").read())

"""

###########################################################################
##--- MODULES
###########################################################################

import os
import argparse
import imp
import time
from hamsterfunctions import *
from diagnosis import main_diagnosis
from attribution import main_attribution
from biascorrection import main_biascorrection

###########################################################################
##--- FUNCTIONS + COMMAND LINE ARGUMENTS
###########################################################################

## (1) LOADING FUNCTIONS
#exec(open("hamsterfunctions.py").read())
#exec(open("flex2traj.py").read())
#exec(open("diagnosis.py").read())
#exec(open("attribution.py").read())
#exec(open("biascorrection.py").read())

## (2) COMMAND LINE ARGUMENTS
# read command line arguments (dates, thresholds and other flags)
args    = read_cmdargs()
verbose = args.verbose
print(printsettings(args))

# just waiting a random number of seconds (max. 30s)
# to avoid overlap of path.exist and makedirs between parallel jobs (any better solution?)
if args.waiter:
    waiter  = random.randint(0,30)
    time.sleep(waiter)

###########################################################################
##--- PATHS
###########################################################################

## determine working directory
wpath = os.getcwd()
os.chdir(wpath)

## load input and output paths & input file name base(s)
print("Using paths from: "+ wpath+"/"+args.pathfile)
content = imp.load_source('',wpath+"/"+args.pathfile) # load like a python module
path_refp = content.path_ref_p
path_refe = content.path_ref_e
path_reft = content.path_ref_t
path_refh = content.path_ref_h
path_orig = content.path_orig
path_diag = content.path_diag
path_attr = content.path_attr
path_bias = content.path_bias
maskfile  = content.maskfile
path_f2t_diag = content.path_f2t_diag
base_f2t_diag = content.base_f2t_diag
path_f2t_traj = content.path_f2t_traj
base_f2t_traj = content.base_f2t_traj

# create output directories if they do not exist (in dependency of step)
if args.steps==0 and args.ctraj_len==0 and not os.path.exists(path_f2t_diag):
        os.makedirs(path_f2t_diag)
        os.makedirs(path_f2t_diag+"/"+str(args.ryyyy))
if args.steps==0 and args.ctraj_len>0 and not os.path.exists(path_f2t_traj):
        os.makedirs(path_f2t_traj)
        os.makedirs(path_f2t_traj+"/"+str(args.ryyyy))
if args.steps==1 and not os.path.exists(path_diag):
        os.makedirs(path_diag)
if args.steps==2 and not os.path.exists(path_attr):
        os.makedirs(path_attr)
if args.steps==3 and not os.path.exists(path_bias):
        os.makedirs(path_bias)

###########################################################################
##--- MAIN
###########################################################################
## (3) RUN main scripts with arguments
if args.steps ==0:
    if args.ctraj_len==0:
        path_f2t=path_f2t_diag
        base_f2t=base_f2t_diag
    elif args.ctraj_len>0:
        path_f2t=path_f2t_traj
        base_f2t=base_f2t_traj
    main_flex2traj(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
                   tml=args.ctraj_len,
                   maskfile=maskfile,
                   maskval=args.maskval,
                   idir=path_orig,
                   odir=path_f2t,
                   fout=base_f2t)

if args.steps == 1:
    main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
              ipath=path_f2t_diag,
              ifile_base=base_f2t_diag,
              opath=path_diag,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              verbose=args.verbose,
              veryverbose=args.veryverbose,
              fproc_npart=args.fproc_npart,
              # E criteria
              fevap=args.fevap,
              cevap_dqv=args.cevap_dqv,
              fevap_drh=args.fevap_drh,
              cevap_drh=args.cevap_drh,
              cevap_hgt=args.cevap_hgt, 
              # P criteria
              fprec=args.fprec,
              cprec_dqv=args.cprec_dqv, 
              cprec_rh=args.cprec_rh,
              # H criteria
              fheat=args.fheat,
              cheat_dtemp=args.cheat_dtemp,
              fheat_drh=args.fheat_drh,
              cheat_drh=args.cheat_drh,
              cheat_hgt=args.cheat_hgt,
              fheat_rdq=args.fheat_rdq,
              cheat_rdq=args.cheat_rdq,
              # pbl and height criteria
              cpbl_method=args.cpbl_method,
              cpbl_strict=args.cpbl_strict,
              cpbl_factor=args.cpbl_factor,
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis,
              fvariable_mass=args.variable_mass,
              strargs=printsettings(args))

if args.steps == 2:
    main_attribution(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad, 
              ipath=path_f2t_traj,
              ifile_base=base_f2t_traj,
              ipath_f2t=path_orig,
              opath=path_attr,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              maskfile=maskfile,
              maskval=args.maskval,
              verbose=args.verbose,
              veryverbose=args.veryverbose,
              ctraj_len=args.ctraj_len,
              # E criteria
              cevap_dqv=args.cevap_dqv,
              fevap_drh=args.fevap_drh,
              cevap_drh=args.cevap_drh,
              cevap_hgt=args.cevap_hgt, 
              # P criteria
              cprec_dqv=args.cprec_dqv, 
              cprec_rh=args.cprec_rh,
              # H criteria
              cheat_dtemp=args.cheat_dtemp,
              fheat_drh=args.fheat_drh,
              cheat_drh=args.cheat_drh,
              cheat_hgt=args.cheat_hgt,
              fheat_rdq=args.fheat_rdq,
              cheat_rdq=args.cheat_rdq,
              # pbl and height criteria
              cpbl_method=args.cpbl_method,
              cpbl_strict=args.cpbl_strict,
              cpbl_factor=args.cpbl_factor,
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis, 
              fdry=args.fallingdry,
              fmemento=args.memento,
              mattribution=args.mattribution,
              crandomnit=args.ratt_nit,
              randatt_forcall=args.ratt_forcall,
              explainp=args.explainp,
              fdupscale=args.dupscale,
              fmupscale=args.mupscale,
              fvariable_mass=args.variable_mass,
              fwritestats=args.writestats,
              strargs=printsettings(args))

if args.steps == 3:
    main_biascorrection(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am,
               opath_attr=path_attr, 
               opath_diag=path_diag, 
               ipath_refp=path_refp,
               ipath_refe=path_refe,
               ipath_reft=path_reft,
               ipath_refh=path_refh,
               opath=path_bias, 
               ofile_base=args.expid, # output
               mode=args.mode,
               maskfile=maskfile,
               maskval=args.maskval,
               verbose=args.verbose,
               veryverbose=args.veryverbose,
               fuseattp=args.bc_useattp,
               bcscale=args.bc_time,
               pref_data=args.pref_data,
               eref_data=args.eref_data,
               href_data=args.href_data,
               faggbwtime=args.bc_aggbwtime,
               fbc_e2p=args.bc_e2p,
               fbc_e2p_p=args.bc_e2p_p,
               fbc_e2p_e=args.bc_e2p_e,
               fbc_e2p_ep=args.bc_e2p_ep,
               fbc_t2p_ep=args.bc_t2p_ep,
               fbc_had=args.bc_had,
               fbc_had_h=args.bc_had_h,
               fdebug=args.debug,
               fwrite_netcdf=args.write_netcdf,
               fwrite_month=args.write_month,
               fwritestats=args.writestats,
               precision=args.precision,
               strargs=printsettings(args))
