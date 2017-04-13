#!/usr/bin/env python
"""
Program     : Rename the CMIP5 files created by calc_CMIP5_clim_anoms to conform to the drs conventions
              and to add a reference to the reference period (ref5)
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 14/03/2017
Modified    : 23/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import os, sys

from drslib.drs import CmipDRS
from drslib import cmip5

from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list
from CMIP5_climatological_periods import CMIP5_climatological_periods
from calc_CMIP5_clim_common import get_drs_object, get_output_base_path, get_start_end_years, get_drs_obj_latest_version


def get_CMIP5_clim_anom_filename(model_desc, experiment, variable, clim_period, years):
    """Get the path to a file calculated by calc_CMIP5_clim_anom.py
       :param Tuple model_desc: model description from CMIP5_model_list
       :param string experiment: experiment description from CMIP5_scenario_list
       :param string variable: variable name, CMOR style, from CMIP5_var_list
       :param strin clim_period: climatological period
       :param Tuple(int, int) years: start year and end year of climatological period
    """
    # get the drs object, the start and end years of the climatological period and the version
    drs_obj = get_drs_object(model_desc, experiment, variable)
    drs_obj.product = "derived"
    drs_obj.version = get_drs_obj_latest_version(drs_obj)

    sy, ey = get_start_end_years(drs_obj, years)
    drs_obj.subset = ((sy, None, None, None, None, None),
                      (ey, None, None, None, None, None), None)
    drs_obj.extended = "clim_mean_anom"

    # get the path from the drs object
    cmip5_trans = cmip5.make_translator(get_output_base_path())

    # get the directory and file for the output
    output_filepath = os.path.realpath(cmip5_trans.drs_to_filepath(drs_obj))
    return output_filepath


if __name__ == "__main__":
    for model_desc in CMIP5_model_list:
        for experiment in CMIP5_scenario_list:
            for variable in CMIP5_var_list.values():
                for clim_period in CMIP5_climatological_periods:
                    # skip AR5 period if historical or piControl or 1pctCO2
                    if clim_period == "AR5" and experiment in ["historical", "piControl", "1pctCO2"]:
                        continue
                    # now get the years:
                    for years in CMIP5_climatological_periods[clim_period]:
                        orig_fname = get_CMIP5_clim_anom_filename(model_desc, experiment, variable, clim_period, years)
                        if os.path.exists(orig_fname):
                            new_fname = orig_fname[:-(len("clim_mean_anom.nc"))] + "clim-anom-ref5.nc"
                            os.rename(orig_fname, new_fname)
