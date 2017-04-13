#!/usr/bin/env python
"""
Program     : Fix the metadata in the files created by calc_CMIP5_clim_anoms
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

from netCDF4 import Dataset

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
    drs_obj.extended = "clim-anom-ref5"

    # get the path from the drs object
    cmip5_trans = cmip5.make_translator(get_output_base_path())

    # get the directory and file for the output
    output_filepath = os.path.realpath(cmip5_trans.drs_to_filepath(drs_obj))
    return output_filepath


def fix_metadata(fname, vname, years):
    """Fix the metadata for the CMIP5 anomaly file
       :param string fname: name of file to fix metadata for
    """
    # open the dataset
    nc_fh = Dataset(fname, 'a')
    # get the attributes
    nc_attrs = nc_fh.ncattrs()
    # fix the history attribute first
    hist_attr = nc_fh.getncattr("history").split("\n")
    for h in range(0, len(hist_attr)):
        # make the history string shorter
        if "cdo -s --no_warnings" in hist_attr[h]:
            s = hist_attr[h].split(" ")
            new_s = []
            for t in s:
                if "datacentre" in t or "group_workspaces" in t:
                    t = t.split("/")[-1]
                new_s.append(t)
                
            s = new_s[:]
            # fix an introduced mistake
            if "2005" in new_s:
                idx_2005 = new_s.index("2005")
                new_s[idx_2005] = "-selyear,1986/2005"
            if str(1900+years[1]) in new_s:
                idx_y1 = new_s.index(str(1900+years[1]))
                new_s[idx_y1] = "-selyear,"+str(1900+years[0])+"/"+str(1900+years[1])
            if str(2000+years[1]) in new_s:
                idx_y1 = new_s.index(str(2000+years[1]))
                new_s[idx_y1] = "-selyear,"+str(2000+years[0])+"/"+str(2000+years[1])

            cdo_pos = s.index("cdo")
            cdo_string = " ".join(new_s[cdo_pos:])
            hist_attr[h] = cdo_string
    new_hist_attr = "\n".join(hist_attr)

    # fix the cdo history
    nc_fh.setncattr("history", new_hist_attr)
    # add some other descriptive data about the content / processing
    nc_fh.setncattr("processing_description", "anomalies of monthly multi-year averages")
    nc_fh.setncattr("processing_type", "mean average")
    nc_fh.setncattr("average_dimension", "time [nested]")
    n_yrs = years[1] - years[0] + 1
    if (n_yrs == 10):
        avg_type = "decadal"
    elif (n_yrs == 20):
        avg_type = "20 years"
    elif (n_yrs == 30):
        avg_type = "30 years"
    else:
        avg_type = "yearly"
    nc_fh.setncattr("average_type_outer", avg_type)
    nc_fh.setncattr("average_type_inner", "monthly")
    nc_fh.setncattr("anomaly_reference_period", "1986-2005")
    nc_fh.setncattr("contact", "support@ceda.ac.uk")
    nc_fh.setncattr("reference", "IPCC DDC AR5 climatologies")
    nc_fh.close()


def fix_metadata_email(fname, vname, years):
    nc_fh = Dataset(fname, 'a')
    nc_attrs = nc_fh.ncattrs()
    nc_fh.setncattr("contact", "support@ceda.ac.uk")    
    nc_fh.close()

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
                        fname = get_CMIP5_clim_anom_filename(model_desc, experiment, variable, clim_period, years)
                        if os.path.exists(fname):
                            print fname
                            fix_metadata(fname, variable, years)
