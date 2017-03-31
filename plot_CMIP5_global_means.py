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

from calc_CMIP5_clim_common import get_start_end_years, get_output_base_path, get_drs_obj_latest_version
from netCDF4 import Dataset

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import scipy.stats as stats

def get_gmts_filename(experiment, clim_period, variable):
    out_fname = "gmts/gmts_" + experiment + "_" + clim_period + "_" + variable + ".npy"
    return out_fname


def plot_boxplots(sp, xpos, data, idx, color):
    # get the sub data from the data
    sub_data = data[:, idx]
    # get the valid values
    valid_data = sub_data[~sub_data.mask]
    # calculate the 5th, 95th, 25th, 75th and 50th percentiles
    p5  = numpy.percentile(valid_data, 5)
    p25 = numpy.percentile(valid_data, 25)
    p50 = numpy.percentile(valid_data, 50)
    p75 = numpy.percentile(valid_data, 75)
    p95 = numpy.percentile(valid_data, 95)
    # bar between 25th and 75th
    sp.bar(xpos+0.1, p75-p25, width=0.8, bottom=p25, ec=color, fc=color, alpha=0.5)
    # 5th / 95th lines
    sp.plot([xpos+0.1,xpos+0.9], [p5,p5], color=color, lw=1.0)
    sp.plot([xpos+0.5,xpos+0.5], [p5,p25], color=color, lw=1.0)
    sp.plot([xpos+0.1,xpos+0.9], [p95,p95], color=color, lw=1.0)
    sp.plot([xpos+0.5,xpos+0.5], [p75,p95], color=color, lw=1.0)
    # median
    sp.plot([xpos,xpos+1], [p50,p50], color=color, lw=2.0)


def plot_CMIP5_global_mean_anoms(clim_period, variable):
    """Plot the timeseries of global mean anomalies for each model from CMIP5.
       The historical, piControl and 1pctCO2 will be plotted on one figure and
       The historical, rcp26, rcp45, rcp60 and rcp85 plotted on one figure
       :param string clim_period: the climatalogical period to calculate (20x|30a|AR5)
       :param string variable: the variable to plot
    """
    # load the historical, rcp26, rcp45, rcp60 and rcp80 scenarios
    mv = 2e20
    if clim_period == "AR5":
        histo_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("historical", "20x", variable)), mv)
    else:
        histo_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("historical", clim_period, variable)), mv)
    rcp26_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("rcp26", clim_period, variable)), mv)
    rcp45_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("rcp45", clim_period, variable)), mv)
    rcp60_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("rcp60", clim_period, variable)), mv)
    rcp85_data = numpy.ma.masked_equal(numpy.load(get_gmts_filename("rcp85", clim_period, variable)), mv)

    # create the December->February means
    histo_DJF = 1.0/3 * (histo_data[11] + histo_data[0] + histo_data[1])
    rcp26_DJF = 1.0/3 * (rcp26_data[11] + rcp26_data[0] + rcp26_data[1])
    rcp45_DJF = 1.0/3 * (rcp45_data[11] + rcp45_data[0] + rcp45_data[1])
    rcp60_DJF = 1.0/3 * (rcp60_data[11] + rcp60_data[0] + rcp60_data[1])
    rcp85_DJF = 1.0/3 * (rcp85_data[11] + rcp85_data[0] + rcp85_data[1])

    # create the plot
    gs = gridspec.GridSpec(1,5)
    sp0 = plt.subplot(gs[0,:-1])
    sp1 = plt.subplot(gs[0,-1])
    gs.update(wspace=0.25, hspace=0.1) 

    # get the ranges
    histo_max = numpy.ma.max(histo_DJF, axis=0)
    histo_min = numpy.ma.min(histo_DJF, axis=0)
    histo_mean = numpy.ma.mean(histo_DJF, axis=0)

    rcp26_max = numpy.ma.max(rcp26_DJF, axis=0)
    rcp26_min = numpy.ma.min(rcp26_DJF, axis=0)
    rcp26_mean = numpy.ma.mean(rcp26_DJF, axis=0)

    rcp45_max = numpy.ma.max(rcp45_DJF, axis=0)
    rcp45_min = numpy.ma.min(rcp45_DJF, axis=0)
    rcp45_mean = numpy.ma.mean(rcp45_DJF, axis=0)

    rcp60_max = numpy.ma.max(rcp60_DJF, axis=0)
    rcp60_min = numpy.ma.min(rcp60_DJF, axis=0)
    rcp60_mean = numpy.ma.mean(rcp60_DJF, axis=0)

    rcp85_max = numpy.ma.max(rcp85_DJF, axis=0)
    rcp85_min = numpy.ma.min(rcp85_DJF, axis=0)
    rcp85_mean = numpy.ma.mean(rcp85_DJF, axis=0)

    # colors
    histo_c = "#AAAAAA"
    rcp26_c = "#191970"
    rcp45_c = "#6495ED"
    rcp60_c = "#FF8C00"
    rcp85_c = "#8B0000"

    for c in range(0, len(CMIP5_climatological_periods[clim_period])):
        cp = CMIP5_climatological_periods[clim_period][c]

        if clim_period == "AR5":
            # plot the histo first
            cph = CMIP5_climatological_periods["20x"][c]
            histo_sy = cph[0] + 1900
            histo_ey = cph[1] + 1900
            # rcps next
            rcp_sy = cp[0]
            rcp_ey = cp[1]
        else:        
            # plot the histo first
            histo_sy = cp[0] + 1900
            histo_ey = cp[1] + 1900
            # rcps next
            rcp_sy = cp[0] + 2000
            rcp_ey = cp[1] + 2000

        # do a square covering the period
        h_l = sp0.fill_between([histo_sy,histo_ey],[histo_min[c],histo_min[c]],[histo_max[c],histo_max[c]], color=histo_c, alpha=0.75)
        sp0.plot([histo_sy, histo_ey], [histo_mean[c], histo_mean[c]], color=histo_c, lw=2.0)
        
        # width of points between sy and ey
        rcp_width = (float(rcp_ey) - rcp_sy)/4.0

        rcp26_l = sp0.fill_between([rcp_sy,rcp_sy+rcp_width],[rcp26_min[c],rcp26_min[c]],[rcp26_max[c],rcp26_max[c]], color=rcp26_c, alpha=0.5)
        rcp45_l = sp0.fill_between([rcp_sy+rcp_width,rcp_sy+2*rcp_width],[rcp45_min[c],rcp45_min[c]],[rcp45_max[c],rcp45_max[c]], color=rcp45_c, alpha=0.5)
        rcp60_l = sp0.fill_between([rcp_sy+2*rcp_width,rcp_sy+3*rcp_width],[rcp60_min[c],rcp60_min[c]],[rcp60_max[c],rcp60_max[c]], color=rcp60_c, alpha=0.5)
        rcp85_l = sp0.fill_between([rcp_sy+3*rcp_width,rcp_sy+4*rcp_width],[rcp85_min[c],rcp85_min[c]],[rcp85_max[c],rcp85_max[c]], color=rcp85_c, alpha=0.5)

        # plot the min, max lines
        sp0.plot([rcp_sy, rcp_ey], [rcp26_min[c], rcp26_min[c]], color=rcp26_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp45_min[c], rcp45_min[c]], color=rcp45_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp60_min[c], rcp60_min[c]], color=rcp60_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp85_min[c], rcp85_min[c]], color=rcp85_c)

        sp0.plot([rcp_sy, rcp_ey], [rcp26_max[c], rcp26_max[c]], color=rcp26_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp45_max[c], rcp45_max[c]], color=rcp45_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp60_max[c], rcp60_max[c]], color=rcp60_c)
        sp0.plot([rcp_sy, rcp_ey], [rcp85_max[c], rcp85_max[c]], color=rcp85_c)

        # plot the means- thicker lines
        sp0.plot([rcp_sy, rcp_ey], [rcp26_mean[c], rcp26_mean[c]], color=rcp26_c, lw=2.0)
        sp0.plot([rcp_sy, rcp_ey], [rcp45_mean[c], rcp45_mean[c]], color=rcp45_c, lw=2.0)
        sp0.plot([rcp_sy, rcp_ey], [rcp60_mean[c], rcp60_mean[c]], color=rcp60_c, lw=2.0)
        sp0.plot([rcp_sy, rcp_ey], [rcp85_mean[c], rcp85_mean[c]], color=rcp85_c, lw=2.0)

    # remove the axes from sp1 plot and set the ticks to show on the right of sp0
    sp1.axis('off')
    sp0.tick_params(labelright=True)

    # do the box and whiskers plot on sp1
    # index for 2081->2100 in 20a is -2, 30x is -1 and AR5 is also -1
    if clim_period == "20x":
        last_idx = -2
    else:
        last_idx = -1

    # plot the boxplots on sp1 - custom function to take into account missing values
    plot_boxplots(sp1, 0, rcp26_DJF, last_idx, rcp26_c)
    plot_boxplots(sp1, 1, rcp45_DJF, last_idx, rcp45_c)
    plot_boxplots(sp1, 2, rcp60_DJF, last_idx, rcp60_c)
    plot_boxplots(sp1, 3, rcp85_DJF, last_idx, rcp85_c)

    # set the sp1 range to be the same as sp0
    sp1.set_ylim(sp0.get_ylim())
    sp1.set_xlim([0,4])
    # dotted line at 0
    sp1.plot([0,4],[0,0], '--k')

    # legend on sp0
    if variable == "rsds" or variable == "psl":
        leg_loc = 3
    else:
        leg_loc = 2
    sp0.legend([h_l, rcp26_l, rcp45_l, rcp60_l, rcp85_l], ["Hist.", "RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"], loc=leg_loc, ncol=2)

    # output the plot
    if not os.path.isdir("./images"):
        os.makedirs("./images")
    image_fname = "images/gmts_" + clim_period + "_" + variable + ".png"

    plt.gcf().set_size_inches(8,4)
#    plt.tight_layout()
    plt.savefig(image_fname)


def calc_CMIP5_global_mean_anoms(experiment, clim_period, variable):
    """Calculate the timeseries of global mean anomalies for each model from CMIP5.
       :param stirng experiment: the experiment to calculate the values for (piControl|historical|rcp26|rcp46|rcp60|rcp85|1pctCO2)
       :param string variable: the variable to calculate (rsds|va|tasmin|huss|hus|tas|ua|vas|pr|psl|uas|tasmax)
       :param string clim_period: the climatalogical period to calculate (20x|30a|AR5)
    """

    # don't do historical for AR5
    if experiment == "historical" and clim_period == "AR5":
        return

    # create the output directory
    if not os.path.isdir("./gmts"):
        os.makedirs("./gmts")

    # create the output filename
    out_fname = get_gmts_filename(experiment, clim_period, variable)

    # do not calculate if already calculated
    if os.path.exists(out_fname):
        return

    # create the output array of the shape [12, len(CMIP5_model_list), len(clim_period)]
    mv = 2e20
    out_array = numpy.ones([12, len(CMIP5_model_list), len(CMIP5_climatological_periods[clim_period])], 'f') * mv

    for m in range(0, len(CMIP5_model_list)):
        model_desc = CMIP5_model_list[m]
        c_count = 0
        for c in CMIP5_climatological_periods[clim_period]:
            # get the drs output object
            anom_drs_obj = CmipDRS(activity="cmip5", product="derived",
                                   institute=model_desc[0], model=model_desc[1],
                                   experiment=experiment, ensemble=model_desc[2],
                                   frequency="mon", table="Amon",
                                   variable=variable,
                                   extended="clim_mean_anom")
            # get the version of the drs object
            anom_drs_obj.version = get_drs_obj_latest_version(anom_drs_obj)

            # get the start year and end year
            start_year, end_year = get_start_end_years(anom_drs_obj, c)
            if start_year < 0:
                continue

            anom_drs_obj.subset=((start_year, None, None, None, None, None),
                                 (end_year, None, None, None, None, None), None)

            # get the path from the drs object
            anom_cmip5_trans = cmip5.make_translator(get_output_base_path())
            anom_output_filepath = os.path.realpath(anom_cmip5_trans.drs_to_filepath(anom_drs_obj))
            # file may not exist!
            if not os.path.exists(anom_output_filepath):
                continue

            # open the netCDF input file / dataset
            nc_fh = Dataset(anom_output_filepath, "r")

            # get the variable and data
            var_data = nc_fh.variables[variable][:]

            # get the variable missing value and reassign to global mv
            var_data = numpy.ma.masked_outside(var_data, -1e6, 1e6)

            # ua and va have height data as well as lat, lon, time - i.e. 4D rather than 3D
            # just get the 2nd (index=1) field
            if var_data.ndim == 4:
                var_data = var_data[:,1,:,:]

            # get the latitude variable data
            lats = nc_fh.variables["lat"][:]
            lat_weights = numpy.cos(numpy.radians(lats))

            # calculate the global (field) mean 
            gmts = numpy.ma.mean(numpy.ma.average(var_data, axis=1, weights=lat_weights), axis=1)
            out_array[:, m, c_count] = gmts
            c_count += 1

            # close the netCDF input file
            nc_fh.close()
    # save the output file
    numpy.save(out_fname, out_array)

if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", nargs=1, help="Variable:" + ("|").join(CMIP5_var_list), default=None)
    parser.add_argument("-c", nargs=1, help="Climatological period:" + ("|").join(CMIP5_climatological_periods), default=None)
    args = parser.parse_args()
    # ensure that the arguments are valid
    assert(args.v != None and (args.v[0] in CMIP5_var_list.keys() or args.v[0] in CMIP5_var_list.values()))
    assert(args.c != None and (args.c[0] in CMIP5_climatological_periods))

    # calculate the GMTS first
    for e in CMIP5_scenario_list:
        calc_CMIP5_global_mean_anoms(e, args.c[0], args.v[0])

    plot_CMIP5_global_mean_anoms(args.c[0], args.v[0])
