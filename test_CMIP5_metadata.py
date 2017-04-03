#!/usr/bin/env python
"""
Program     : Test the metadata of the derived CMIP5 data that is to be uploaded to the IPCC DDC
              Data quality checks are:
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 03/04/2017
Requires    : CDO, NCO, Python, drslib
"""
import os, sys
import argparse
from drslib.drs import CmipDRS
from drslib import cmip5
from netCDF4 import Dataset

from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list
from CMIP5_climatological_periods import CMIP5_climatological_periods

from calc_CMIP5_clim_common import get_drs_obj_latest_version, get_output_base_path
from calc_CMIP5_clim_common import get_cmip5_base_path, get_drs_object, get_start_end_years, get_files_between_year_period


def produce_metadata_report(model, variable, experiment, clim_period, years):
    """Check the metadata of the derived data file matches the metadata of the original source file 
    """
    # get the model description, variable, experiment and climatological period
    model_desc = CMIP5_model_list[model-1]

    # initial empty string for the errors
    run_string = "r" + str(model_desc[2][0]) + "i" + str(model_desc[2][1]) + "p" + str(model_desc[2][2])
#    ERROR_STRING = "[" + model_desc[0] + " " + model_desc[1] + " " + run_string + ", " + variable + ", " + experiment + ", " + clim_period +"]\n"
    ERROR_STRING = ""

    # get the original drs object
    orig_drs = get_drs_object(model_desc, experiment, variable)
    # get the version
    orig_drs.version = get_drs_obj_latest_version(orig_drs)
    # get the years
    start_year, end_year = get_start_end_years(orig_drs, years)
    orig_drs.subset = ((start_year,None,None,None,None,None), (end_year,None,None,None,None,None), None)

    # CMIP5 files do not have a direct mapping from the dates to the files - so we just use the first of 
    # the files that are within the desired date range
    orig_files = get_files_between_year_period(orig_drs, start_year, end_year)

    if len(orig_files) == 0:
#        ERROR_STRING += "  [FILE ERROR] [missing file]: "+model_desc[0]+" "+model_desc[1]+" "+experiment+" ("+clim_period+") "+str(start_year)+"->"+str(end_year)+" "+variable+"\n"
        return ERROR_STRING
    orig_path = orig_files[0]
    
    # create the derived product drs
    derived_drs = CmipDRS(orig_drs)
    derived_drs.extended = "clim_mean_anom"
    derived_drs.product = "derived"
    derived_trans = cmip5.make_translator(get_output_base_path())
    derived_logical_path = derived_trans.drs_to_filepath(derived_drs)
    derived_path = os.path.realpath(derived_logical_path)

    # File error
    if not os.path.exists(derived_path):
#        ERROR_STRING += "[FILE ERROR] [missing file] " + derived_path + " [source file] " + orig_path +"\n"
        return ERROR_STRING

    # open the file at the original path and the derived path
    orig_fh = Dataset(orig_path, "r")
    derived_fh = Dataset(derived_path, "r")
    
    # get the variables
    orig_var = orig_fh.variables[variable]
    derived_var = orig_fh.variables[variable]

    # get the attributes
    orig_attrs = orig_var.ncattrs()
    derived_attrs = derived_var.ncattrs()

    # check attribute presence
    for a in orig_attrs:
        if not a in derived_attrs:
            ERROR_STRING += "[ATTRIBUTE ERROR] [missing attribute] " + a + " for variable " + variable + " in file " + derived_path + "\n"
    else:
        if derived_var.getncattr(a) != orig_var.getncattr(a):
            ERROR_STRING += "[ATTRIBUTE ERROR] [attribute value] " + a + " for variable " + variable + " in file " + derived_path + "\n"

    # close the files
    orig_fh.close()
    derived_fh.close()
    return ERROR_STRING

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
    assert(args.m == None or int(args.m[0]) <= len(CMIP5_model_list))
    assert(args.e == None or (args.e[0] in CMIP5_scenario_list))
    assert(args.v == None or args.v[0] in CMIP5_var_list.keys() or args.v[0] in CMIP5_var_list.values())
    assert(args.c == None or (args.c[0] in CMIP5_climatological_periods))

    # build the lists of the model, experiment, variable and climatological period
    if args.m == None:
        model_list = range(1, len(CMIP5_model_list)+1)
    else:
        model_list = [int(args.m[0])]

    if args.e == None:
        exp_list = CMIP5_scenario_list
    else:
        exp_list = [args.e[0]]

    if args.v == None:
        var_list = CMIP5_var_list.values()
    else:
        var_list = [args.v[0]]

    if args.c == None:
        clim_periods = CMIP5_climatological_periods.keys()
    else:
        clim_periods = [args.c[0]]

    fh = open("attribute_check.txt", "w")
    # subset CMIP5_climatological_periods to remove AR5 periods if experiment is "historical"
    for m in model_list:
        for v in var_list:
            for e in exp_list:
                for cp in clim_periods:
                    # skip AR5 period if historical or piControl or 1pctCO2
                    if cp == "AR5" and e in ["historical", "piControl", "1pctCO2"]:
                        continue
                    # now get the years:
                    for yrs in CMIP5_climatological_periods[cp]:
                        # produce the metadata report
                        md_report = produce_metadata_report(m, v, e, cp, yrs)
                        if md_report != "":
                            fh.write(md_report)
    fh.close()
