#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
import matplotlib.pyplot as plt
import cdo as cdo_module
from baseline import global_sum_baselines
from cmip_files import get_filename

cdo = cdo_module.Cdo(tempdir='.')
#cdo.debug = True

# Constants.
MIPS = ['CMIP', 'LUMIP']
EXPERIMENTS = {
        'CMIP': ['esm-hist'],
        'LUMIP': ['esm-ssp585-ssp126Lu']}
TABLES = ['Lmon', 'Emon']
VARIABLES = {
        'Lmon': ['gpp', 'npp', 'ra', 'rh', 'nbp', 'cVeg', 'cLitter'],
        # ACCESS does not have nep
        'Emon': ['cSoil']}
ENSEMBLES = [str(e) for e in range(1,11)]
LAND_FRAC_FILE = '/g/data/fs38/publications/CMIP6/LUMIP/CSIRO/ACCESS-ESM1-5/esm-ssp585-ssp126Lu/'\
        'r10i1p1f1/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_r10i1p1f1_gn.nc'
KG_IN_PG = 1000000000000
SEC_IN_YEAR = 60*60*24*365
SEC_IN_DAY = 60*60*24
NEW_UNITS_FACTOR = SEC_IN_YEAR/KG_IN_PG


# Calc global sum and load variable
gpp_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', '1', 'Lmon', 'gpp')
time_units = SEC_IN_DAY
gpp = cdo.divc(str(KG_IN_PG),
        input=f'-yearsum -muldpm -mulc,{time_units} -fldsum -mul -mul '+gpp_files[0]+\
                f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
        options='-L',
        returnCdf=True).variables['gpp'][:].squeeze()

# Plot the graph
years = [y for y in range(2015,2101)]
plt.plot(years, gpp-global_sum_baselines['gpp'], label="ACCESS-ESM1-5")
plt.hlines(0, years[0], years[-1], colors='black', linestyles='dotted')
plt.xlabel('Year')
plt.ylabel('GPP anomaly (PgC/year)')
plt.legend()
plt.show()
