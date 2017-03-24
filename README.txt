CEDA IPCC DDC climatologies for CMIP5 / AR5
===========================================

This repository contains the files used to calculate the climatologies 
of the models in CMIP5 / used in AR5 for the IPCC DDC (data distribution 
centre).

The docstring at the top of calc_CMIP5_clim_common.py contains a summation of 
the task performed by this repository.

Requires:
    - CDO (climate data operators)
    - NCO (netCDF operators)
    - Python2.7
    - drslib (for Python)

MIP tables
----------

The CMIP5 / CMOR MIP tables need to be downloaded via git from:
    https://github.com/PCMDI/cmip5-cmor-tables

and placed in a directory called cmip5-cmor-tables. 
