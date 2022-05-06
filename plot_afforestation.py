#!/usr/bin/env python3
# plot_afforestation.py generates plots for the afforestation simulation.
import glob

import cdo as cdo_module
import matplotlib.pyplot as plt
import numpy as np
from baseline import global_sum_baselines
from cmip_files import get_filename
from constants import TABLES, VARIABLES, ENSEMBLES, SEC_IN_YEAR, KG_IN_PG

cdo = cdo_module.Cdo(tempdir='.')
#cdo.debug = True

# Local constants
LAND_FRAC_FILE = '/g/data/fs38/publications/CMIP6/LUMIP/CSIRO/ACCESS-ESM1-5/esm-ssp585-ssp126Lu/'\
        'r10i1p1f1/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_r10i1p1f1_gn.nc'
SEC_IN_DAY = 60*60*24
COLORS = {'gpp':'green',
        'npp':'olive',
        'ra':'orange',
        'rh':'saddlebrown',
        'nbp':'purple',
        'cVeg':'darkgreen',
        'cLitter':'chocolate',
        'cSoil':'black'}

for table in TABLES:
    for var in VARIABLES[table]:
        # Calc global sum and load variable
        data = np.ones((len(ENSEMBLES), 86))*np.nan
        for ens in ENSEMBLES:
            data_files = get_filename('LUMIP', 'esm-ssp585-ssp126Lu', ens, table, var)
            if var in ['gpp', 'npp', 'ra', 'rh', 'nbp']:
                data[int(ens)-1,:] = cdo.divc(str(KG_IN_PG),
                        input=f'-yearsum -muldpm -mulc,{SEC_IN_DAY} -fldsum -mul -mul '\
                                +data_files[0]+\
                                f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
                        options='-L',
                        returnCdf=True).variables[var][:].squeeze()
            else:
                data[int(ens)-1,:] = cdo.divc(str(KG_IN_PG),
                        input=f'-yearmonmean -fldsum -mul -mul '\
                                +data_files[0]+\
                                f' -divc,100 {LAND_FRAC_FILE} -gridarea {LAND_FRAC_FILE}',
                        options='-L',
                        returnCdf=True).variables[var][:].squeeze()

        # Anomaly, mean and std
        data_anomaly = data - global_sum_baselines[var]
        data_ens_mean = np.mean(data_anomaly, axis=0)
        data_ens_std = np.std(data_anomaly, axis=0, ddof=1)
        
        # Plot the graph
        years = [y for y in range(2015,2101)]
        plt.figure()
        for ens in ENSEMBLES:
            plt.plot(years, data_anomaly[int(ens)-1], color='lightgray', linewidth=0.6)
        plt.plot(years, data_ens_mean, color=COLORS[var], label="Ensemble mean")
        plt.plot(years, data_ens_mean+data_ens_std, color=COLORS[var], linewidth=0.8,
                label="+-1$\sigma$")
        plt.plot(years, data_ens_mean-data_ens_std, color=COLORS[var], linewidth=0.8)
        plt.hlines(0, years[0], years[-1], colors='black', linestyles='dotted')
        plt.xlabel('Year')
        plt.ylabel(f'{var.upper()} anomaly (PgC/year)')
        plt.title(f"ACCESS-ESM1-5 {var.upper()}")
        plt.savefig(f'{var}_ACCESS-ESM1-5_esm-ssp585-ssp126Lu_ensembles_anomalies.svg')

# Clean up
#temp_files = glob.glob('./cdoPy*')
#os.remove(temp_files)
#plt.show()

