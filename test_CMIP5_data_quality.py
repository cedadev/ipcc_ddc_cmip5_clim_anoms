#!/usr/bin/env python
"""
Program     : Test the quality of the derived CMIP5 data that is to be uploaded to the IPCC DDC
              Data quality checks are:
                1. File presence.  The files calculated match those in Table AI.1 in IPCC AR5 Annex 1
                2. Metadata integrity.  The metadata in the derived file matches that in the original file. (The derived file will have some extra attributes, of course)
                
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 09/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import os
import argparse
from drslib.drs import CmipDRS
from drslib import cmip5
from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list

from calc_CMIP5_clim_common import get_output_base_path, get_drs_version

def load_CMIP5_AR5_table():
    """ load the CMIP5 WG1 AR5 AI.1 table
        note that the table HAS to be in the order Model Name, piControl, historical, rcp26, rcp45, rcp60, rcp85"""
    fh = open("CMIP5_AR5_table.csv")
    lines = fh.readlines()
    # output is a dictionary
    table = {}
    for l in lines[1:]:
        spl = l.split(",")
        ele = spl[1:]
        ele[-1] = ele[-1].strip()
        table[spl[0]] = ele
    fh.close()
    return table

def produce_model_run_scenario_table():
    """Produce a table cataloging which runs are used for each model and scenario, and output it in a manner similar to
       table AI.1 in IPCC AR5 Annex 1"""
    
    # output table
    table = {}
    # loop over all model (descriptions)
    for m in CMIP5_model_list:
        # loop over all scenarios
        sc_list = []
        for s in CMIP5_scenario_list:
            var_count = 0
            for v in CMIP5_var_list:
                # create a drs object containing the model, scenario, run, etc. info
                drs_obj = CmipDRS(activity="cmip5", product="derived",
                                  institute=m[0], model=m[1],
                                  experiment=s, ensemble=m[2],
                                  frequency="mon", table="Amon",
                                  variable=CMIP5_var_list[v],
                                  extended="clim_mean")
                drs_obj.version = get_drs_version(drs_obj, get_output_base_path())
                # get the path from the drs object
                cmip5_trans = cmip5.make_translator(get_output_base_path())
                output_path = os.path.realpath(cmip5_trans.drs_to_path(drs_obj))
                # get the number of files
                if os.path.exists(output_path) and len(os.listdir(output_path)) != 0:
                    var_count += 1
            if var_count > 0:
                sc_list.append(str(m[2][0]))
            else:
                sc_list.append('')
        # append pertubed physics name to GISS models
        model_name = m[1]
        if model_name == "GISS-E2-H" or model_name == "GISS-E2-R":
            model_name += " p"+str(m[2][2])
        table[model_name] = sc_list
    return table


def produce_variables_table():
    """Print the variables that exist for each model, scenario pair as a table"""
    table = {}
    for m in CMIP5_model_list:
        sc_list = []
        for s in CMIP5_scenario_list:
            out_vars = []
            for v in CMIP5_var_list:
                # create a drs object containing the model, scenario, run, etc. info
                drs_obj = CmipDRS(activity="cmip5", product="derived",
                                  institute=m[0], model=m[1],
                                  experiment=s, ensemble=m[2],
                                  frequency="mon", table="Amon",
                                  variable=CMIP5_var_list[v],
                                  extended="clim_mean")
                drs_obj.version = get_drs_version(drs_obj, get_output_base_path())
                # get the path from the drs object
                cmip5_trans = cmip5.make_translator(get_output_base_path())
                output_path = os.path.realpath(cmip5_trans.drs_to_path(drs_obj))
                # get the number of files that are non_zero and compare to the number of files listed
                if os.path.exists(output_path):
                    n_files = len(os.listdir(output_path))
                    c_files = 0
                    for f in os.listdir(output_path):
                        if os.path.getsize(output_path + "/" + f) > 0:
                            c_files += 1
                    # if all files are non-zero then append the variable to the output variable list
                    if n_files == c_files:
                        out_vars.append(CMIP5_var_list[v])
            sc_list.append(out_vars)
        # append pertubed physics name to GISS models
        model_name = m[1]
        if model_name == "GISS-E2-H" or model_name == "GISS-E2-R":
            model_name += " p"+str(m[2][2])
        table[model_name] = sc_list
    return table


def print_variables_table(table):
    # load the CMIP5 table so as to replicate the order of the models
    cmip5_table = load_CMIP5_AR5_table()
    # sort the keys alphabetically but not by case
    k = sorted(cmip5_table.keys(), key=str.lower)
    for m in k:
        print m + "," ,
        for sc in table[m]:
            var_string = " ".join(sc)
            var_string += " ["+str(len(sc))+"]"
            print var_string + "," ,
        print

def compare_tables(cmip5_table, model_table):
    """Compare two tables produced by the above functions
     We are interested in two cases:
     1. An entry is in the cmip5_table but not in the model_table (missing run)
     2. An entry is not in the cmip5_table but is in the model_table (extra run)"""
    missing_runs = []
    extra_runs = []
    for m in cmip5_table:
        cmip5_runs = cmip5_table[m]
        model_runs = model_table[m]
        for s in range(0, len(cmip5_runs)):
            if cmip5_runs[s] == "" and model_runs[s] != "":
                extra_runs.append((m, CMIP5_scenario_list[s]))
            if cmip5_runs[s] != "" and model_runs[s] == "":
                missing_runs.append((m, CMIP5_scenario_list[s]))
    # print out the missing runs
    print "Missing runs"
    print "============"
    for m in missing_runs:
        print m[0], m[1]
    print
    print "Extra runs"
    print "=========="
    for e in extra_runs:
        print e[0], e[1]
    print


def print_missing_variables_table(cmip5_table, model_table):
    """Compare the cmip5 table with the variables model_table"""
    missing_runs = []
    missing_vars = {}
    cmip5_set = set(['hus', 'tasmin', 'ua', 'rsds', 'va', 'tasmax', 'tas', 'pr', 'psl', 'uas', 'vas', 'huss'])
    for m in cmip5_table:
        cmip5_runs = cmip5_table[m]
        model_runs = model_table[m]
        for s in range(0, len(cmip5_runs)):
            if cmip5_runs[s] != "" and len(model_runs[s]) == 0:
                missing_runs.append((m, CMIP5_scenario_list[s]))
            if cmip5_runs[s] != "" and len(model_runs[s]) != len(cmip5_set):
                var_set = set(model_runs[s])
                item = (m, CMIP5_scenario_list[s], cmip5_set - var_set)
            else:
                item = ()
            if not m in missing_vars.keys():
                missing_vars[m] = [item]
            else:
                missing_vars[m].append(item)
    # print the table heading
    print "CMIP5 model name, piControl, historical, rcp26, rcp45, rcp60, rcp85, 1pctCO2"
    # print out the missing runs
    keys = sorted(missing_vars.keys(), key=str.lower)
    for m in keys:
        item = missing_vars[m]
        print m + ",", 
        for x in item:
            if len(x) != 0:
                print " ".join(x[2]) + ", ",
            else:
                print ", ",
        print
    

if __name__ == "__main__":
#    cmip5_table = load_CMIP5_AR5_table()
#    model_table = produce_model_run_scenario_table()
#    compare_tables(cmip5_table, model_table)
    model_table = produce_variables_table()
    cmip5_table = load_CMIP5_AR5_table()
#    print_missing_variables_table(cmip5_table, model_table)
    print_variables_table(model_table)
