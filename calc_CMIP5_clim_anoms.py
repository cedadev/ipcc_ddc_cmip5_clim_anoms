#!/usr/bin/env python
"""
Program     : Calculate mean climatology anomaly fields of CMIP5 models for upload to the IPCC Data Distribution Centre
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 14/03/2017
Modified    : 23/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import subprocess
import os, sys
import argparse
from drslib.drs import CmipDRS
from drslib import cmip5
from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list
from CMIP5_climatological_periods import CMIP5_climatological_periods

from calc_CMIP5_clim_common import *


def calculate_reference_period(input_drs_obj):
    """Calculate the reference period.  Anomalies are expressed with respect to these reference values.
       The reference period depends on the experiment.
       For 'historical' experiment it is 1986->2005.
       For 'piControl' experiment it is the last 20 years of the run
       :param CmipDRS input_drs_obj: A drs object describing the input to generate the reference file from
    """
    # reference start and end years
    if input_drs_obj.experiment == "historical":
        ref_sy = 1986
        ref_ey = 2005
    else:
        min_yr, max_yr = get_min_max_years(input_drs_obj)
        ref_sy = max_yr - 20  
        ref_ey = max_yr

    # get the files and drs_version
    ref_file_list = get_files_between_year_period(input_drs_obj, ref_sy, ref_ey)

    # create the output drs object
    # copy the input drs object but replace some information
    output_drs_obj = CmipDRS(input_drs_obj)
    output_drs_obj.product = "reference"
    output_drs_obj.subset = ((ref_sy, None, None, None, None, None),
                             (ref_ey, None, None, None, None, None), None)
    output_drs_obj.extended = "clim_mean"

    # get the path from the drs object
    cmip5_trans = cmip5.make_translator(get_output_base_path())

    # get the directory and file for the output
    ref_output_path = os.path.realpath(cmip5_trans.drs_to_path(output_drs_obj))
    ref_output_filepath = os.path.realpath(cmip5_trans.drs_to_filepath(output_drs_obj))

    # create the output file if it doesn't exist
    if not os.path.exists(ref_output_filepath):
        # create the output path if it doesn't exist
        if not os.path.exists(ref_output_path):
            os.makedirs(ref_output_path)
        # concatenate all of the files together in time and output to a temporary file in the output directory
        tmp_file = ref_output_path + "/tmp.nc"
        if os.path.exists(tmp_file):
            os.unlink(tmp_file)
        cdo(["mergetime", " ".join(ref_file_list), tmp_file])
        # calculate the climatologies, first selecting the ref_sy->ref_ey period
        cdo(["ymonmean", "-selyear,"+str(ref_sy)+"/"+str(ref_ey), tmp_file, ref_output_filepath])
        # delete the temporary file
        os.unlink(tmp_file)
    return ref_output_filepath


def calculate_climatological_mean_anomaly(model_desc, experiment, variable, clim_period):
    """Calculate the climatological mean anomaly for a climatological period for a particular CMIP5 model, scenario and variable
       :param List[string] model_desc: description of model, see 
       :param string experiment: name of the experiment, historical, rcp45, etc
       :param string variable: name of the variable, standard name or CMOR name
       :param Tuple(int, int) clim_period: climatological period in a tuple: (start_year, end_year)
    """

    # create the input drs object
    input_drs_obj = get_drs_object(model_desc, experiment, variable)
    get_drs_obj_latest_version(input_drs_obj)

    # get the files between the start_year and end_year for the model description, experiment and variable
    start_year, end_year = get_start_end_years(input_drs_obj, clim_period)
    input_file_list = get_files_between_year_period(input_drs_obj, start_year, end_year)

    # check if any files returned
    if len(input_file_list) == 0:
        return ""

    # the 30 year means (30a) are defined as being between 2001->2030 for the first period.
    # so, we need to add the historical file to the input_file_list
    if "rcp" in experiment and clim_period[0] == 1:
        # first copy the drs obj, and change the experiment to "historical"
        historical_drs_obj = CmipDRS(input_drs_obj)
        historical_drs_obj.experiment = "historical"
        historical_file_list = get_files_between_year_period(historical_drs_obj, 2001, 2005)
        input_file_list.extend(historical_file_list)

    # create the input reference files drs object
    # which reference to use?
    if experiment in ["1pctCO2", "piControl"]:
        ref_exp = "piControl"
    else:
        ref_exp = "historical"

    # create the input reference drs object and get the version of it which is the nearest to the input_drs_obj version
    input_ref_drs_obj = get_drs_object(model_desc, ref_exp, variable)
    get_drs_obj_nearest_version(input_ref_drs_obj, input_drs_obj.version)

    # calculate the reference period (if necessary, the function will check) and return its name
    ref_filepath = calculate_reference_period(input_ref_drs_obj)

    # get the drs output object - copy the input object and change some info
    output_drs_obj = CmipDRS(input_drs_obj)
    output_drs_obj.product = "derived"
    output_drs_obj.subset = ((start_year, None, None, None, None, None),
                             (end_year, None, None, None, None, None), None)
    output_drs_obj.extended = "clim_mean_anom"

    # get the path from the drs object
    anom_cmip5_trans = cmip5.make_translator(get_output_base_path())
    anom_output_path = os.path.realpath(anom_cmip5_trans.drs_to_path(output_drs_obj))
    anom_output_filepath = os.path.realpath(anom_cmip5_trans.drs_to_filepath(output_drs_obj))

    # delete the file if it exists - then we can do checking later, and after multiple passes
    if os.path.exists(anom_output_filepath):
        os.unlink(anom_output_filepath)

    # create the output path if it doesn't exist
    if not os.path.exists(anom_output_path):
        os.makedirs(anom_output_path)

    # create the output name of the temporary file
    tmp_file = anom_output_path + "/tmp.nc"
    # delete temporary file if it already exists
    if os.path.exists(tmp_file):
        os.unlink(tmp_file)

    # concatenate all of the files together in time and output to a temporary file in the output directory
    selyear_str = "-selyear,{:d}/{:d}".format(start_year, end_year)
    cdo_cmd = ["-mergetime"]
    for f in input_file_list:
        cdo_cmd.append(selyear_str)
        cdo_cmd.append(f)
    cdo_cmd.append(tmp_file)
    cdo(cdo_cmd)

    # calculate the climatologies and subtract the reference climatology
    cdo(["sub", "-ymonmean", selyear_str, tmp_file, ref_filepath, anom_output_filepath])

    # delete the temporary file
    os.unlink(tmp_file)
    return anom_output_filepath


if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", nargs=1, help="Model number 1-" + str(len(CMIP5_model_list)))
    parser.add_argument("-e", nargs=1, help="Experiment: " + ("|").join(CMIP5_scenario_list))
    parser.add_argument("-v", nargs=1, help="Variable:" + ("|").join(CMIP5_var_list.values()), default=None)
    parser.add_argument("-c", nargs=1, help="Climatological period: " + ("|").join(CMIP5_climatological_periods), default=None)
    parser.add_argument("-s", action="store_true", help="Display model numbers")
    args = parser.parse_args()

    if args.s:
        # print the model list
        for m in range(0, len(CMIP5_model_list)):
            print m+1, CMIP5_model_list[m][1]
        sys.exit()
    # ensure that the arguments are valid
    assert(args.m and int(args.m[0]) <= len(CMIP5_model_list))
    assert(args.e and (args.e[0] in CMIP5_scenario_list))
    assert(args.v == None or args.v[0] in CMIP5_var_list.keys() or args.v[0] in CMIP5_var_list.values())
    assert(args.c == None or (args.c[0] in CMIP5_climatological_periods))

    # get the model description and scenario
    model_desc = CMIP5_model_list[int(args.m[0])-1]
    experiment = args.e[0]
    # get the variable or use all
    if args.v == None:
        var_list = CMIP5_var_list.values()
    else:
        var_list = [args.v[0]]

    if args.c == None:
        clim_periods = CMIP5_climatological_periods.keys()
    else:
        clim_periods = [args.c[0]]

    print model_desc, experiment, clim_periods, var_list

    # subset CMIP5_climatological_periods to remove AR5 periods if experiment is "historical"
    for variable in var_list:
        for p in clim_periods:
            # skip AR5 period if historical or piControl or 1pctCO2
            if p == "AR5" and experiment in ["historical", "piControl", "1pctCO2"]:
                continue
            # now get the years:
            for clim_period in CMIP5_climatological_periods[p]:
                try:
                    clim_mean_anom_file = calculate_climatological_mean_anomaly(model_desc, experiment, variable, clim_period)
                except:
                    print "Failed on ", model_desc, experiment, variable
