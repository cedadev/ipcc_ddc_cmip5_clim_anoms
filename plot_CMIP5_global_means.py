#!/usr/bin/env python
"""
Program     : Plot the derived CMIP5 data that is to be uploaded to the IPCC DDC
              This is part of the data quality check, as the plots can be compared directly to those in IPCC AR5 WG1 Annex 1
              The program will produce three plots for each variable: 20x (twenty year averages), 30a (thirty year averages) and AR5 periods
              Each plot is a timeseries of the anomalies from the reference period for each time period, for each experiment.
              Time periods are relative to: historical: 1900, rcpx.y: 2000, piControl and 1pctCO2: start of the run

Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 22/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import os
import argparse
from drslib.drs import CmipDRS
from drslib import cmip5
from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list
from CMIP5_climatological_periods import CMIP5_climatological_periods

from calc_CMIP5_clim_common import get_start_end_years, get_output_base_path, get_drs_version


def plot_CMIP5_global_mean_anoms(variable):
    """Plot the timeseries of global mean anomalies for each model from CMIP5.
       The historical, piControl and 1pctCO2 will be plotted on one figure and
       The historical, rcp26, rcp45, rcp60 and rcp85 plotted on another figure
       :param string variable: the variable to plot
    """
    pass


def calc_CMIP5_global_mean_anoms(experiment, clim_period, variable):
    """Calculate the timeseries of global mean anomalies for each model from CMIP5.
       :param string variable: the variable to plot
    """

    for m in range(0, len(CMIP5_model_list)):
        model_desc = CMIP5_model_list[m]
        print CMIP5_climatological_periods[clim_period]
        for c in CMIP5_climatological_periods[clim_period]:
            # get the start year and end year  
            start_year, end_year = get_start_end_years(model_desc, experiment, c, variable)
            print start_year, end_year
            if start_year < 0:
                continue
            # get the drs output object
            anom_drs_obj = CmipDRS(activity="cmip5", product="derived",
                                   institute=model_desc[0], model=model_desc[1],
                                   experiment=experiment, ensemble=model_desc[2],
                                   frequency="mon", table="Amon",
                                   variable=variable,
                                   subset=((start_year, None, None, None, None, None),
                                           (end_year, None, None, None, None, None), None),
                                   extended="clim_mean_anom")
            # get the version of the drs object
            anom_drs_obj.version = get_drs_version(anom_drs_obj, get_output_base_path())

            # get the path from the drs object
            anom_cmip5_trans = cmip5.make_translator(get_output_base_path())
            anom_output_path = os.path.realpath(anom_cmip5_trans.drs_to_path(anom_drs_obj))                
            anom_output_filepath = os.path.realpath(anom_cmip5_trans.drs_to_filepath(anom_drs_obj))
            
            # calculate the global (field) mean


if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", nargs=1, help="Variable:" + ("|").join(CMIP5_var_list), default=None)
    parser.add_argument("-e", nargs=1, help="Experiment: " + ("|").join(CMIP5_scenario_list))
    parser.add_argument("-c", nargs=1, help="Climatological period:" + ("|").join(CMIP5_climatological_periods), default=None)
    args = parser.parse_args()
    # ensure that the arguments are valid
    assert(args.v != None and (args.v[0] in CMIP5_var_list.keys() or args.v[0] in CMIP5_var_list.values()))
    assert(args.c != None and (args.c[0] in CMIP5_climatological_periods))
    assert(args.e and (args.e[0] in CMIP5_scenario_list))

    # calculate the GMTS first
    calc_CMIP5_global_mean_anoms(args.e[0], args.c[0], args.v[0])
