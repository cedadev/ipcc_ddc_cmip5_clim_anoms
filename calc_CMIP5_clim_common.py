#!/usr/bin/env python
"""
Program     : Common functions for calculating climatologies of CMIP5 models for upload to the IPCC Data Distribution Centre
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 07/03/2017
Requires    : CDO, NCO, Python2.7, drslib (for Python)

IPCC DDC climatologies for AR5
==============================

Atmospheric variables
---------------------

- precipitation flux (kg m^{-2} s^{-1})
- air pressure at sea level (Pa)
- air temperature (K)
- air temperature daily max (K)
- air temperature daily min (K)
- eastward wind (m s^{-1}) (on pressure levels)
- northward wind (m s^{-1}) (on pressure levels)
- eastward wind (m s^{-1}) 
- northward wind (m s^{-1})
- surface downwelling shortwave flux in air (W m^{-2})
- specific humidity (1) (on pressure levels)
- specific humidity at the surface (1)

Scenarios
---------

- Historical (20thC)
- RCP2.6
- RCP4.5
- RCP6.0
- RCP8.5
- piControl  (pre-industrial control)
- 1pctCO2 (1 percent per year CO2 increase)

Anomaly specification
---------------------

All values are taken as anomalies from the 1986->2005 mean.
For simulations these base line values are taken from the historical run (20thC), for the historical, and RCP experiments.
For the piControl and 1pctCO2 experiments the reference period is the last 20 years of the piControl run.

Periods
-------

- periods (20x and 30a) are defined as being:
    - Relative to 2000 for future scenarios (RCPy.x)
    - Relative to 1900 for historical run (20thC)
    - Relative to the beginning of the experiment for the piControl and 1pctCO2 runs 

Periods covered are:

- 20x: Twenty year averages +20-30, +46-65, +80-99, +180-199
- 30a: Thirty year averages +01-30, +31-60, +61-90
- AR5: 2016-2035, 2046-2065, 2081-2100

Output
------

- Period mean of the monthly mean of the variables
- 12 monthly values per climatology
- Climatology is a 2D field
- One climatology per model per scenario. See AR5, WG1 Annex 1, Table AI.1
    - The information from this table is replicated in CMIP5_model_list
"""

import subprocess
import os
from drslib.drs import DRSFileSystem, CmipDRS
from drslib import cmip5
import sys
from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list

DEBUG=True

def get_output_base_path():
    return "/group_workspaces/jasmin/cedaproc/nrmassey/CMIP5"


def get_cmip5_base_path():
    return "/badc/cmip5/data/cmip5"


def cdo(cmd_params):
    """Run the cdo command using subprocess.check_output"""
    cmd = ["/usr/bin/cdo", "-s", "--no_warnings"]
    cmd.extend(cmd_params)
    if DEBUG:
        print " ".join(cmd)
    return subprocess.check_output(cmd) 


def get_drs_object(model_desc, experiment, variable):
    """Get a drs object from a model description, experiment and variable.  This can then be used to build paths to the files.
       :param list[] model_desc: description of CMIP5 model to use, a single member of CMIP5_model_list
       :param string experiment: which experiment (historical|rcp26|rcp45|rcp60|rcp85|piControl|1pctCO2)
       :param string variable: name of the variable, standard name or CMOR name
    """
   # get the drs object
    drs_obj = CmipDRS(activity="cmip5", product="output1",
                      institute=model_desc[0], model=model_desc[1],
                      experiment=experiment, ensemble=model_desc[2],
                      frequency="mon", table="Amon",
                      variable=variable)
    return drs_obj


def get_drs_obj_latest_version(drs_obj):
    """Get the latest version of a drs object in the CMIP5 archive.
       This function searches for the "latest" symbolic link in the drs filepath, resolves this link and sets the version to be that of the link
       :param CmipDRS drs_obj: drs_object as created by get_drs_object above
       :returns int version: integer version number
    """

    # get the path from the drs object
    if drs_obj.product == "derived":
        base_path = get_output_base_path()
    else:
        base_path = get_cmip5_base_path()

    # return a different version for the GFDL-CM3 model for rcp45 experiment
    if drs_obj.model == "CanESM2" and drs_obj.experiment == "rcp45":
        drs_obj.version = 20110829
        return 20110829
    # fake the version
    drs_obj.version = 1
    # get the fake path
    cmip5_trans = cmip5.make_translator(base_path)
    logical_path = cmip5_trans.drs_to_path(drs_obj)
    path = os.path.realpath(logical_path)
    # split the path
    split_path = path.split("/")
    # get the path up to the version number [-2]
    version_base_path = "/".join(split_path[:-2])
    # resolve the symbolic link
    if not os.path.exists(os.path.realpath(version_base_path)):
        return 0
    file_list = os.listdir(os.path.realpath(version_base_path))
    # if symbolic link called latest is present then resolve it as the path
    if "latest" in file_list:
        true_version_path = os.path.realpath(version_base_path + "/latest")
    else:
        # find the latest
        latest = file_list[0]
        for f in file_list[1:]:
            if int(f[1:]) > int(latest[1:]):
                latest = f
        true_version_path = os.path.realpath(version_base_path + "/" + latest)

    # get the version as an integer
    drs_obj.version = int(true_version_path.split("/")[-1][1:])
    return drs_obj.version


def get_drs_obj_nearest_version(drs_obj, ref_version):
    """Get the latest version of a drs object in the CMIP5 archive.
       This function searches for the "latest" symbolic link in the drs filepath, resolves this link and sets the version to be that of the li$
       :param CmipDRS drs_obj: drs_object as created by get_drs_object above
       :param integer ref_version: version of reference drs_obj, we want to get as close as possible to this version
       :returns int version: integer version number
    """

    # get the path from the drs object
    if drs_obj.product == "derived":
        base_path = get_output_base_path()
    else:
        base_path = get_cmip5_base_path()

    # fake the version
    drs_obj.version = 1
    # get the fake path
    cmip5_trans = cmip5.make_translator(base_path)
    path = os.path.realpath(cmip5_trans.drs_to_path(drs_obj))
    # split the path
    split_path = path.split("/")
    # get the path up to the version number [-2]
    version_base_path = "/".join(split_path[:-2])
    # resolve the symbolic link
    if not os.path.exists(os.path.realpath(version_base_path)):
        return -1
    file_list = os.listdir(os.path.realpath(version_base_path))
    # find the nearest version in the file list
    nv = 1e20
    for f in file_list:
        try:
            # current version looking at
            cv = int(f[1:])
            # calculate the difference
            if abs(cv - ref_version) < abs(nv - ref_version):
                nv = cv
        except:
            pass

    # set the version
    drs_obj.version = nv
    return nv
    

def get_files_between_year_period(drs_obj, start_year, end_year):
    """Get the path to a file(s) which encompasses start_year and end_year for a model, scenario, run and variable.
       This function is necessary as some models have a number of output files, each with different years in and some models have one output file with all years in
       :param CmipDRS drs_obj: drs object as created by get_drs_obj above
       :param integer start_year: start of the period
       :param integer end_year: end of the period
       :return: list of paths to the file(s) containing the start_year and / or end_year  
    """

    # could be in more than one file
    out_files = []
 
    # get the path from the drs object
    if drs_obj.product == "derived":
        base_path = get_output_base_path()
    else:
        base_path = get_cmip5_base_path()
    cmip5_trans = cmip5.make_translator(base_path)
    path = os.path.realpath(cmip5_trans.drs_to_path(drs_obj))
    
    # check the path exists
    if not os.path.exists(path):
        return ""

    # get a list of the files
    files = os.listdir(path)
    # loop over the files
    for f in files:
        # reject those with .nc extension
        if f[-3:] != ".nc":
            continue
        # resolve all symbolic links on the file
        real_f = os.path.realpath(path + "/" + f)
        if False:
            # This is the slow version which uses cdo to get the years from the files
            # It can fail on an empty file
            # run cdo showyear on the file
            try:
                out_yrs = cdo(["showyear", real_f])
                # convert the output into a list of integers
                out_yrs = map(lambda x: int(x), out_yrs.split())
            except:
                pass
        else:
            # This is the faster version which uses the drslib to get the start and end year from the filename
            drs_out_obj = cmip5_trans.filename_to_drs(f)
            out_yrs = range(drs_out_obj.subset[0][0], drs_out_obj.subset[1][0]+1)
        # check if start year or end year in out_yrs and add file to output list of files if it is
        if (len(out_yrs) != 0 and \
            ((out_yrs[0] >= start_year and out_yrs[-1] <= end_year) or \
             start_year in out_yrs or end_year in out_yrs)):
            out_files.append(real_f)
            
    return out_files


def get_min_max_years(drs_obj):
    """Get the minimum and maximum years in a model run described by model_desc, experiment, variable
       :param CmipDRS drs_obj: drs object as created by get_drs_obj above
       :return: min and max year in experiment
       :rtype: Tuple(integer, integer)
    """

    # get the path from the drs object
    if drs_obj.product == "derived":
        base_path = get_output_base_path()
    else:
        base_path = get_cmip5_base_path()

    # get the path from the drs object
    cmip5_trans = cmip5.make_translator(base_path)
    logical_path = cmip5_trans.drs_to_path(drs_obj)
    path = os.path.realpath(logical_path)

    # get a list of the files
    if not os.path.exists(path):
        return -1, -1
    files = os.listdir(path)
    # loop over the files to find the maximum and minimum year
    max_year = -1
    min_year = 1e20

    for f in files:
        # reject those with .nc extension
        if f[-3:] != ".nc":
            continue
        if False:
            # This is version one, which uses cdo to get the years from the files - i.e. we don't trust the filename.  It's thorough but slow
            # run cdo showyear on the file: need to resolve all symbolic links first!
            real_f = os.path.realpath(path + "/" + f)
            out_yrs = cdo(["showyear", real_f])
            # convert the output into a list of integers
            out_yrs = map(lambda x: int(x), out_yrs.split())
        else:
            # This is version two, which relies on the drs filename being correct - it's much quicker!
            drs_out_obj = cmip5_trans.filename_to_drs(f)
            out_yrs = [drs_out_obj.subset[0][0], drs_out_obj.subset[1][0]]
        # get the min and max year
        if max(out_yrs) > max_year:
            max_year = max(out_yrs)
        if min(out_yrs) < min_year:
            min_year = min(out_yrs)
    return min_year, max_year


def get_first_year(drs_obj):
    """Get the first year in a model_desc, experiment, variable experiment description
       This function is necessary to find the first year of the 1pct CO2 runs and piControl runs
       :param CmipDRS drs_obj: drs object as created by get_drs_obj above
       :return: First year of the experiment
       :rtype: integer
    """

    # get the min and max years in the experiment
    min_year, max_year = get_min_max_years(drs_obj)
    return min_year


def get_start_end_years(drs_obj, clim_period):
    """Get the start and end year for an experiment, variable and climate_meaning_period
       :param CmipDRS drs_obj: drs object as created by get_drs_obj above
       :param Tuple(int, int) clim_period: climatological period in a tuple: (start_year, end_year)
       :returns: the start and end year of the climatological period
       :rtype: Tuple(int, int)
    """
    # get the start year and end year
    if clim_period[0] > 1000:    # actual real years, rather than offsets
        start_year = clim_period[0]
        end_year = clim_period[1]
    else:
        if drs_obj.experiment == "historical":
            # ref year is 1900 for historical run 
            ref_year = 1900
        elif drs_obj.experiment in ["rcp26", "rcp45", "rcp60", "rcp85"]:
            # ref year is 2000 for future (rcp) runs
            ref_year = 2000
        elif drs_obj.experiment in ["1pctCO2", "piControl"]:
            ref_year = get_first_year(drs_obj)

        if ref_year > 0:
            start_year = ref_year + clim_period[0]
            end_year = ref_year + clim_period[1]
        else:
            start_year = -1
            end_year = -1

    return start_year, end_year
