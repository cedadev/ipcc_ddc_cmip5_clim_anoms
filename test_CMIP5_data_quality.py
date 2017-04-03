#!/usr/bin/env python
"""
Program     : Test the quality of the derived CMIP5 data that is to be uploaded to the IPCC DDC
              Data quality checks are:
                1. File presence.  The files calculated match those in Table AI.1 in IPCC AR5 Annex 1
                2. File size >= than a minimum field size (HadCM3 field size)
                
Author      : Neil Massey
Organisation: CEDA (STFC)
Email       : neil.massey@stfc.ac.uk
Date        : 09/03/2017
Requires    : CDO, NCO, Python, drslib
"""
import os, sys
import argparse
from drslib.drs import CmipDRS
from drslib import cmip5
from CMIP5_model_list import CMIP5_model_list
from CMIP5_var_list import CMIP5_var_list
from CMIP5_scenario_list import CMIP5_scenario_list
from CMIP5_climatological_periods import CMIP5_climatological_periods

from calc_CMIP5_clim_common import get_drs_obj_latest_version, get_start_end_years, get_files_between_year_period, get_output_base_path

# minimum field size for every month plus latitude and longitude to check file size - equal to HadCM3 resolution
MIN_FIELD_SIZE=(96*73*4*12)+96*4+73*4

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

def produce_model_run_scenario_table(model):
    """Produce a table cataloging which runs are used for each model and scenario, and output it in a manner similar to
       table AI.1 in IPCC AR5 Annex 1"""

    if model == -1:
        model_list = CMIP5_model_list
    else:
        model_list = [CMIP5_model_list[model]]

    # output table
    table = {}
    # loop over all model (descriptions)
    for m in model_list:
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
                drs_obj.version = get_drs_obj_latest_version(drs_obj)
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


def produce_variables_table(model):
    """Print the variables that exist for each model, scenario pair as a table"""
    table = {}
    if model == -1:
        model_list = CMIP5_model_list
    else:
        model_list = [CMIP5_model_list[model]]
    for m in model_list:
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
                drs_obj.version = get_drs_obj_latest_version(drs_obj)
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
    """Print the variables that exist for each model, scenario pair as a table"""
    # sort the keys alphabetically but not by case
    k = sorted(table.keys(), key=str.lower)
    for m in k:
        print m + "," ,
        for sc in table[m]:
            var_string = " ".join(sc)
            var_string += " ["+str(len(sc))+"]"
            print var_string + "," ,
        print

def print_missing_runs_table(cmip5_table, model_table):
    """Compare two tables produced by the above functions
     We are interested in two cases:
     1. An entry is in the cmip5_table but not in the model_table (missing run)
     2. An entry is not in the cmip5_table but is in the model_table (extra run)"""
    missing_runs = []
    extra_runs = []
    for m in model_table:
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
    cmip5_var_set = set(CMIP5_var_list.values())
    for m in model_table:
        cmip5_runs = cmip5_table[m]
        model_runs = model_table[m]
        for s in range(0, len(cmip5_runs)):
            if cmip5_runs[s] != "" and len(model_runs[s]) == 0:
                missing_runs.append((m, CMIP5_scenario_list[s]))
            if cmip5_runs[s] != "" and len(model_runs[s]) != len(cmip5_var_set):
                var_set = set(model_runs[s])
                item = (m, CMIP5_scenario_list[s], cmip5_var_set - var_set)
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
    mn = 1
    for m in keys:
        item = missing_vars[m]
        print m + ",", 
        for x in item:
            if len(x) != 0:
                print " ".join(x[2]) + ", ",
            else:
                print ", ",
        print
        mn += 1


def print_missing_climatological_periods_table(model, ignore_vars):
    """Print a table of model, variable, experiment and climatological periods covered.
       This is quite a large table!"""

    if model == -1:
        model_list = CMIP5_model_list
    else:
        model_list = [CMIP5_model_list[model]]

    # path translator for use later
    cmip5_trans = cmip5.make_translator(get_output_base_path())

    # get the model list
    cmip5_table = load_CMIP5_AR5_table()

    for m in model_list:
        # append pertubed physics name to GISS models
        model_name = m[1]
        if model_name == "GISS-E2-H" or model_name == "GISS-E2-R":
            model_name += " p"+str(m[2][2])

        # get the cmip5 model - order is piControl, historical, rcp26, rcp45, rcp60, rcp85
        cmip5_model = cmip5_table[model_name]
        cmip5_experiment_order = ["piControl", "historical", "rcp26", "rcp45", "rcp60", "rcp85", "1pctCO2"]

        # print the model name and scenarios as a header
        print model_name +", ",
        for e in CMIP5_scenario_list:
            print e +", ",
        print

        # loop over the variables
        for v in CMIP5_var_list:
            print CMIP5_var_list[v] + ",",
            if ignore_vars and v in ["vas", "uas", "huss"]:
                print ",,,,,,,"
            # loop over the experiments
            for e in CMIP5_scenario_list:
                # loop over the climatological periods
                cp_str = ""
                # create a drs object containing the model, scenario, run, etc. info
                drs_obj = CmipDRS(activity="cmip5", product="derived",
                                  institute=m[0], model=m[1],
                                  experiment=e, ensemble=m[2],
                                  frequency="mon", table="Amon",
                                  variable=CMIP5_var_list[v],
                                  extended="clim_mean")
                drs_obj.version = get_drs_obj_latest_version(drs_obj)

                # get the experiment details from the cmip5 model table
                cmip5_experiment_index = cmip5_experiment_order.index(e)
                cmip5_experiment = cmip5_model[cmip5_experiment_index]
                # don't do the checks if this experiment is not in the CMIP5 table
                if cmip5_experiment == "":
                    print ",",
                    continue

                for cp in CMIP5_climatological_periods:
                    # only perform the AR5 test for the rcp scenarios and not for historical, piControl and 1pctCO2
                    if cp == "AR5" and not "rcp" in e:
                        continue
                    # get the start/end year offsets tuple
                    years = CMIP5_climatological_periods[cp]
                    # keep a count of the good files
                    good_files = 0
                    bad_file_string = ""
                    for y in years:
                        # don't check the +180 offset for all the models as only some models have that offset
                        if cp == "20x" and y[0] == 180:
                            continue
                        # need a matching input drs to get the actual start and end years from the offsets
                        input_drs_obj = CmipDRS(drs_obj)
                        input_drs_obj.extended = None
                        input_drs_obj.product = "output1"
                        start_year, end_year = get_start_end_years(input_drs_obj, y)
                        # if the start / end year could not be found
                        if start_year == -1 or end_year == -1:
                            bad_file_string += "xx"+str(y[0])+" "                            
                            continue
                        # check to see if the exact file has been created
                        drs_file_obj = CmipDRS(drs_obj)
                        drs_file_obj.subset = ((start_year,None,None,None,None,None), (end_year,None,None,None,None,None), None)
                        drs_file_obj.extended = "clim_mean_anom"
                        file = os.path.realpath(cmip5_trans.drs_to_filepath(drs_file_obj))
                        if os.path.exists(file) and os.path.getsize(file) >= MIN_FIELD_SIZE:
                            good_files += 1
                        else:
                            # if a bad file is encountered, set the counter to something very large, so that good_files is never > 0 after this point
                            good_files = -1e10
                            bad_file_string += str(start_year)+" "
                    # if the number of good files is not zero then add the climatic period to the output string
                    if good_files <= 0:
                        cp_str +=  cp + "[" + bad_file_string[0:-1] + "] "
                print cp_str + ",",
            print
        print


tests = {"print_variables_table" : print_variables_table.__doc__,
         "print_missing_variables_table" : print_missing_variables_table.__doc__,
         "print_missing_runs_table" : print_missing_runs_table.__doc__,
         "print_missing_climatological_periods_table" : print_missing_climatological_periods_table.__doc__,
        }


if __name__ == "__main__":
    # create the argument parser and the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", nargs=1, help="Model number 1-" + str(len(CMIP5_model_list)), default=None)
    parser.add_argument("-t", nargs=1, help="Test number 1-" + str(len(tests)))
    parser.add_argument("-T", action="store_true", help="Show test numbers")
    parser.add_argument("-s", action="store_true", help="Show model numbers")
    parser.add_argument("-i", action="store_true", help="Ignore huss, vas, uas variables")

    args = parser.parse_args()
    if args.T:
        # print the test list
        tkeys = tests.keys()
        for t in range(0, len(tkeys)):
            print str(t+1)+".", tkeys[t], tests[tkeys[t]]
        sys.exit()

    if args.s:
        # print the model list
        for m in range(0, len(CMIP5_model_list)):
            print m+1, CMIP5_model_list[m][1]
        sys.exit()

    if args.i:
        ignore_vars = True
    else:
        ignore_vars = False

    assert(args.t and int(args.t[0]) <= len(tests))
    assert(args.m == None or int(args.m[0]) <= len(CMIP5_model_list))

    t = int(args.t[0])

    if args.m:
        model = int(args.m[0])-1
    else:
        model = -1

    if t == 1 + tests.keys().index("print_variables_table"):     # print_variables_table
        model_table = produce_variables_table(model)
        print_variables_table(model_table)
    elif t == 1 + tests.keys().index("print_missing_variables_table"):   # print_missing_variables_table
        cmip5_table = load_CMIP5_AR5_table()
        model_table = produce_variables_table(model)
        print_missing_variables_table(cmip5_table, model_table)
    elif t == 1 + tests.keys().index("print_missing_runs_table"):   # print_missing_runs_table
        cmip5_table = load_CMIP5_AR5_table()
        model_table = produce_model_run_scenario_table(model)
        print_missing_runs_table(cmip5_table, model_table)
    elif t == 1 + tests.keys().index("print_missing_climatological_periods_table"):   # print_climatological_periods_table
        print_missing_climatological_periods_table(model, ignore_vars)
