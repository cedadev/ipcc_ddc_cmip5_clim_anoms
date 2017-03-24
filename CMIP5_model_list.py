# List of CMIP5 models that are used in the code to generate the IPCC DDC climatologies
# These match table AI1 in IPCC AR5 WG1 Annex I
# The model details are encoded as (in drs language):
# [institute, model, ensemble member]
# ensemble member is a tuple (L,M,N) which converts to the string rLiMpN

CMIP5_model_list = [\
["CSIRO-BOM",    "ACCESS1-0",      (1,1,1)],
["CSIRO-BOM",    "ACCESS1-3",      (1,1,1)],
["BCC",          "bcc-csm1-1",     (1,1,1)],
["BCC",          "bcc-csm1-1-m",   (1,1,1)],
["BNU",          "BNU-ESM",        (1,1,1)],
["CCCma",        "CanESM2",        (1,1,1)],
["NCAR",         "CCSM4",          (1,1,1)],
["NSF-DOE-NCAR", "CESM1-BGC",      (1,1,1)],
["NSF-DOE-NCAR", "CESM1-CAM5",     (1,1,1)],
["CMCC",         "CMCC-CM",        (1,1,1)],
["CMCC",         "CMCC-CMS",       (1,1,1)],
["CNRM-CERFACS", "CNRM-CM5",       (1,1,1)],
["CSIRO-QCCCE",  "CSIRO-Mk3-6-0",  (1,1,1)],
["ICHEC",        "EC-EARTH",       (8,1,1)],
["LASG-CESS",    "FGOALS-g2",      (1,1,1)],
["FIO",          "FIO-ESM",        (1,1,1)],
["NOAA-GFDL",    "GFDL-CM3",       (1,1,1)],
["NOAA-GFDL",    "GFDL-ESM2G",     (1,1,1)],
["NOAA-GFDL",    "GFDL-ESM2M",     (1,1,1)],
["NASA-GISS",    "GISS-E2-H",      (1,1,1)],
["NASA-GISS",    "GISS-E2-H",      (1,1,2)],
["NASA-GISS",    "GISS-E2-H",      (1,1,3)],
["NASA-GISS",    "GISS-E2-H-CC",   (1,1,1)],
["NASA-GISS",    "GISS-E2-R",      (1,1,1)],
["NASA-GISS",    "GISS-E2-R",      (1,1,2)],
["NASA-GISS",    "GISS-E2-R",      (1,1,3)],
["NASA-GISS",    "GISS-E2-R-CC",   (1,1,1)],
["NIMR-KMA",     "HadGEM2-AO",     (1,1,1)],
["MOHC",         "HadGEM2-CC",     (1,1,1)],
["MOHC",         "HadGEM2-ES",     (2,1,1)],
["INM",          "inmcm4",         (1,1,1)],
["IPSL",         "IPSL-CM5A-LR",   (1,1,1)],  
["IPSL",         "IPSL-CM5A-MR",   (1,1,1)],
["IPSL",         "IPSL-CM5B-LR",   (1,1,1)],
["MIROC",        "MIROC5",         (1,1,1)],
["MIROC",        "MIROC-ESM",      (1,1,1)],
["MIROC",        "MIROC-ESM-CHEM", (1,1,1)],
["MPI-M",        "MPI-ESM-LR",     (1,1,1)],
["MPI-M",        "MPI-ESM-MR",     (1,1,1)],
["MPI-M",        "MPI-ESM-P",      (1,1,1)],
["MRI",          "MRI-CGCM3",      (1,1,1)],
["NCC",          "NorESM1-M",      (1,1,1)],
["NCC",          "NorESM1-ME",     (1,1,1)],
]
